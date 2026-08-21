#!/bin/bash
#SBATCH --job-name=Reads_in_feat   # Job name
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16            # Number of CPU cores per task
#SBATCH --mem=32G                     # Job memory request
#SBATCH --time=48:00:00               # Time limit hrs:min:sec
#SBATCH --output=logs/Reads_in_feat_%j.log   # Standard output and error log

# Specify directories and files
export data_dir="data/2024/K562_annotations/processed_beds/sorted"
# export bw_dir="data/2024/RepliTag/bws_nodup"
# export bw_dir_spike="data/2024/RepliTag/bws_nodup_spike"
export bw_dir_bulk="data/2024/RepliTag/bws_bulk"
export bw_dir_drug="data/2024/RepliTag/trackHub/bw_drug"
export outdir="data/2024/RepliTag/Matrices"
# export peakdir="data/2024/SH_all_data/peakPRINT_50"
peakdir="data/2024/SH_all_data/transitions"
# Create output directories
mkdir -p logs $outdir


command -v computeMatrix >/dev/null || module load deepTools/3.5.4.post1-gfbf-2022b

# for bed in "$peakdir"/*K_regions.bed; do
#     base_name=$(basename "$bed" .txt)
#     echo $bed
    
#     multiBigwigSummary BED-file --BED $bed \
#     -b $(ls $bw_dir_drug/*.bw) \
#     --smartLabels \
#     -out $outdir/coverage_scores_${base_name}_no_spike.npz \
#     --outRawCounts $outdir/coverage_scores_${base_name}_no_spike.tab \
#     -p $SLURM_CPUS_PER_TASK
    
    # multiBigwigSummary BED-file --BED $bed \
    # -b $(ls $bw_dir_bulk/*.bw) \
    # --smartLabels \
    # -out $outdir/coverage_scores_${base_name}_bulk.npz \
    # --outRawCounts $outdir/coverage_scores_${base_name}_bulk.tab \
    # -p $SLURM_CPUS_PER_TASK
# done


# for bed in "$data_dir"/gencode.v44.basic.TSS_1000-300_sorted.bed; do
#     base_name=$(basename "$bed" .txt)
#     echo $bed
#     # multiBigwigSummary BED-file --BED $bed \
#     # -b $(ls $bw_dir/*.bw) \
#     # --smartLabels \
#     # -out $outdir/coverage_scores_${base_name}_no_spike.npz \
#     # --outRawCounts $outdir/coverage_scores_${base_name}_no_spike.tab \
#     # -p $SLURM_CPUS_PER_TASK
    
#     # multiBigwigSummary BED-file --BED $bed \
#     # -b $(ls $bw_dir_spike/*.bw) \
#     # --smartLabels \
#     # -out $outdir/coverage_scores_${base_name}.npz \
#     # --outRawCounts $outdir/coverage_scores_${base_name}.tab \
#     # -p $SLURM_CPUS_PER_TASK

#     multiBigwigSummary BED-file --BED $bed \
#     -b $(ls $bw_dir_bulk/*.bw) \
#     --smartLabels \
#     -out $outdir/coverage_scores_${base_name}_bulk.npz \
#     --outRawCounts $outdir/coverage_scores_${base_name}_bulk.tab \
#     -p $SLURM_CPUS_PER_TASK
# done

for bed in "$data_dir"/gencode.v44.basic.TSS-1000_TES_sorted.bed; do
    base_name=$(basename "$bed" .txt)
#     echo $bed
#     multiBigwigSummary BED-file --BED $bed \
#     -b $(ls $bw_dir/*.bw) \
#     --smartLabels \
#     -out $outdir/coverage_scores_${base_name}_no_spike.npz \
#     --outRawCounts $outdir/coverage_scores_${base_name}_no_spike.tab \
#     -p $SLURM_CPUS_PER_TASK
    
#     multiBigwigSummary BED-file --BED $bed \
#     -b $(ls $bw_dir_spike/*.bw) \
#     --smartLabels \
#     -out $outdir/coverage_scores_${base_name}.npz \
#     --outRawCounts $outdir/coverage_scores_${base_name}.tab \
#     -p $SLURM_CPUS_PER_TASK

    multiBigwigSummary BED-file --BED $bed \
    -b $(ls $bw_dir_bulk/*.bw) \
    --smartLabels \
    -out $outdir/coverage_scores_${base_name}_bulk.npz \
    --outRawCounts $outdir/coverage_scores_${base_name}_bulk.tab \
    -p $SLURM_CPUS_PER_TASK
done

#mod to use no spike bws

# multiBigwigSummary bins --binSize 10000 \
# -b $(ls $bw_dir_drug/*.bw) \
# --smartLabels \
# -out $outdir/coverage_scores_no_spike_10kb.npz \
# --outRawCounts $outdir/coverage_scores_no_spike_10kb.tab \
# -p $SLURM_CPUS_PER_TASK

# multiBigwigSummary bins --binSize 10000 \
# -b $(ls $bw_dir_spike/*.bw) \
# --smartLabels \
# -out $outdir/coverage_scores_10kb.npz \
# --outRawCounts $outdir/coverage_scores_10kb.tab \
# -p $SLURM_CPUS_PER_TASK

# multiBigwigSummary bins --binSize 10000 \
# -b $(ls $bw_dir_bulk/*.bw) \
# --smartLabels \
# -out $outdir/coverage_scores_bulk_10kb.npz \
# --outRawCounts $outdir/coverage_scores_bulk_10kb.tab \
# -p $SLURM_CPUS_PER_TASK
