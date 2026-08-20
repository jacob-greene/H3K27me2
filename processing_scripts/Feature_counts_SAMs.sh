#!/bin/bash
#SBATCH --job-name=feature_counts   # Job name
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16            # Number of CPU cores per task
#SBATCH --mem=32G                     # Job memory request
#SBATCH --time=48:00:00               # Time limit hrs:min:sec
#SBATCH --output=logs/feature_counts_%j.log   # Standard output and error log

# Load required module
module load Subread/2.0.3-GCC-11.2.0

# Set output directory and GTF file path
export outdir="data/2024/cCREs/intersected_annos/FC_out"
mkdir -p "$outdir"
export sections_dir="data/2024/SH_all_data/sections4"
export saf_file="data/2024/K562_annotations/processed_beds/sorted/gencode.v44.basic.TSS-1000_TES_sorted.saf"

# Count reads in genes using featureCounts, using the output of the cat command as input
featureCounts -T $SLURM_CPUS_PER_TASK -p --countReadPairs -B -d 1 -D 1000 --primary -O -F SAF -a "$saf_file" \
-o "$outdir/hg38.gencodev44_HKP_genes-1000_featureCounts_bulk_counts.txt" \
$(cat $sections_dir/{K27ac,K27me1,K27me2,K27me3,P2}_{H,K,P}.txt)

# Clean up featureCounts output
awk 'BEGIN{FS=OFS="\t"} !/^#/{print}' "$outdir/hg38.gencodev44_HKP_genes-1000_featureCounts_bulk_counts.txt" > "$outdir/hg38.gencodev44_HKP_genes-1000_featureCounts_bulk_counts.tsv"
rm "$outdir/hg38.gencodev44_HKP_genes-1000_featureCounts_bulk_counts.txt"
mv "$outdir/hg38.gencodev44_HKP_genes-1000_featureCounts_bulk_counts.txt.summary" "$outdir/hg38.gencodev44_HKP_genes-1000_featureCounts_bulk_counts.summary.tsv"
