#!/usr/bin/env bash
#SBATCH --job-name=feature_counts
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=logs/feature_counts_%A_%a.log  # %a for array index

set -euo pipefail

command -v featureCounts >/dev/null || module load Subread/2.0.3-GCC-11.2.0

# --- Config ---
export outdir="data/2024/cCREs/intersected_annos/FC_out"
export sections_dir="data/2024/SH_all_data/sections4"

# Directory containing all SAFs you want to run
export saf_dir="data/2024/LAD_data"

# SAF selection: all *.saf in saf_dir
mapfile -t SAFS < <(find "$saf_dir" -maxdepth 1 -type f -name "*.saf" | sort)
if [[ ${#SAFS[@]} -eq 0 ]]; then
  echo "No SAF files found in $saf_dir" >&2
  exit 1
fi

# BAM list (expand your sets here); one path per line in the *_K.txt lists
mapfile -t BAMS < <(cat "$sections_dir"/{K27ac,K27me1,K27me2,K27me3,P2}_K.txt)
if [[ ${#BAMS[@]} -eq 0 ]]; then
  echo "No BAMs found from lists in $sections_dir" >&2
  exit 1
fi

mkdir -p "$outdir"

run_one_saf () {
  local saf="$1"
  local base
  base="$(basename "$saf" .saf)"

  local out_txt="$outdir/${base}_featureCounts_bulk_counts.txt"
  local out_tsv="$outdir/${base}_featureCounts_bulk_counts.tsv"
  local out_sum_txt="$outdir/${base}_featureCounts_bulk_counts.txt.summary"
  local out_sum_tsv="$outdir/${base}_featureCounts_bulk_counts.summary.tsv"

  echo "[$(date)] Running featureCounts for SAF: $saf"
  featureCounts \
    -T "${SLURM_CPUS_PER_TASK:-16}" \
    -p --countReadPairs \
    -B -d 1 -D 1000 \
    --primary -O \
    -F SAF -a "$saf" \
    -o "$out_txt" \
    "${BAMS[@]}"

  # Clean to TSV and rename summary
  awk 'BEGIN{FS=OFS="\t"} !/^#/{print}' "$out_txt" > "$out_tsv"
  rm -f "$out_txt"
  mv "$out_sum_txt" "$out_sum_tsv"

  echo "[$(date)] Done: $out_tsv ; summary: $out_sum_tsv"
}

# --- Array mode: process a single SAF per task ---
if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  idx="${SLURM_ARRAY_TASK_ID}"
  if (( idx < 0 || idx >= ${#SAFS[@]} )); then
    echo "SLURM_ARRAY_TASK_ID $idx out of range (0..$(( ${#SAFS[@]} - 1 )))" >&2
    exit 1
  fi
  run_one_saf "${SAFS[$idx]}"
else
  # Fallback: process all SAFs sequentially in one job
  echo "SLURM_ARRAY_TASK_ID not set; processing all ${#SAFS[@]} SAFs sequentially."
  for saf in "${SAFS[@]}"; do
    run_one_saf "$saf"
  done
fi
