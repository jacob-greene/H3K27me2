#!/bin/bash
set -euo pipefail

PROJECT="."
OUT="$PROJECT/data/results/replitag_sphase_nucleation"
LOG_DIR="$OUT/logs"
TASKS="$OUT/config/matrix_tasks.tsv"
GENERATOR="$PROJECT/processing_scripts/replitag_nucleation_generate_matrices.sh"
TILE_SIZES="${TILE_SIZES:-1000}"
export TILE_SIZES

cd "$PROJECT"
mkdir -p "$LOG_DIR"

source /app/lmod/lmod/init/profile
module load BEDTools/2.31.0-GCC-12.3.0
module load deepTools/3.5.4.post1-gfbf-2022b

bash "$GENERATOR" prepare

ntasks="$(( $(wc -l < "$TASKS") - 1 ))"
if (( ntasks < 1 )); then
  echo "No matrix tasks found in $TASKS" >&2
  exit 1
fi

array_script="$OUT/config/replitag_nucleation_matrix_array.sbatch"
collate_script="$OUT/config/replitag_nucleation_matrix_collate.sbatch"

cat > "$array_script" <<SBATCH
#!/bin/bash
#SBATCH --job-name=replitag_nuc_matrix
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --time-min=01:00:00
#SBATCH --partition=campus-new
#SBATCH --output=$LOG_DIR/matrix_array_%A_%a.out
#SBATCH --error=$LOG_DIR/matrix_array_%A_%a.err

set -euo pipefail
export TILE_SIZES="$TILE_SIZES"
source /app/lmod/lmod/init/profile
module load BEDTools/2.31.0-GCC-12.3.0
module load deepTools/3.5.4.post1-gfbf-2022b
cd "$PROJECT"
bash "$GENERATOR" coverage "\${SLURM_ARRAY_TASK_ID}"
SBATCH

cat > "$collate_script" <<SBATCH
#!/bin/bash
#SBATCH --job-name=replitag_nuc_collate
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --partition=campus-new
#SBATCH --output=$LOG_DIR/collate_%j.out
#SBATCH --error=$LOG_DIR/collate_%j.err

set -euo pipefail
export TILE_SIZES="$TILE_SIZES"
source /app/lmod/lmod/init/profile
module load BEDTools/2.31.0-GCC-12.3.0
module load deepTools/3.5.4.post1-gfbf-2022b
cd "$PROJECT"
bash "$GENERATOR" collate
SBATCH

array_job="$(sbatch --parsable --array=1-"$ntasks"%80 "$array_script")"
collate_job="$(sbatch --parsable --dependency=afterok:"$array_job" "$collate_script")"

echo "Submitted matrix array job: $array_job ($ntasks tasks, max 80 concurrent; TILE_SIZES=$TILE_SIZES)"
echo "Submitted dependent collation job: $collate_job"
echo "Monitor with: squeue -j $array_job,$collate_job"
