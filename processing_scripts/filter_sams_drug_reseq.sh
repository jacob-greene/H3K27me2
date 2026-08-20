#!/bin/bash
#SBATCH --job-name=filter_blacklist
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --output=logs/filter_blacklist_%j.log

# Load required modules
module load SAMtools/1.19.2-GCC-13.2.0
module load BEDTools/2.30.0-GCC-12.2.0
module load parallel

# Set paths
export BLACKLIST="data/Kc_merged_blacklist_hg38.bed"
export OUT_DIR="data/2024/RepliTag/filtered_sams_drug_1.0/reseq"
export INPUT_DIR1="/shared/ngs/illumina/henikoff/250915_bowtie2/JG_HsDm/marked/sams/"
export INPUT_DIR2="/shared/ngs/illumina/henikoff/250910_bowtie2/JG_HsDm/marked/sams/"

export saf_file="data/2024/K562_annotations/processed_beds/sorted/gencode.v44.basic.TSS-1000_TES_sorted.saf"

mkdir -p "$OUT_DIR" logs

# Function for filtering SAM files
filter_sam() {
    sam="$1"
    base=$(basename "$sam" .sam)

    # Extract blacklisted read names
    samtools sort -n "$sam" -O SAM | \
      samtools view -h - | \
      bedtools bamtobed -bedpe -i stdin | \
      awk -v OFS="\t" '$1 == $4 {print $1, $2, $6, $7}' | \
      bedtools intersect -a stdin -b "$BLACKLIST" -f 1.0 -wa | \
      cut -f4 | sort -u > "$OUT_DIR/${base}_blacklisted.txt"

    # Remove blacklisted reads from SAM
    grep -Fv -f "$OUT_DIR/${base}_blacklisted.txt" "$sam" | \
    samtools sort -n -O SAM -o "$OUT_DIR/${base}_filtered.sam" -

}
export -f filter_sam

# Run in parallel
find "$INPUT_DIR1" -name "*.sam" | \
  parallel --bar -j $SLURM_CPUS_PER_TASK filter_sam {}

find "$INPUT_DIR2" -name "*.sam" | \
  parallel --bar -j $SLURM_CPUS_PER_TASK filter_sam {}

# Count reads in genes using featureCounts
module load Subread/2.0.3-GCC-11.2.0

# Generate list of SAM files
sam_files=$(ls "$OUT_DIR"/*.sam 2>/dev/null)
if [[ -z "$sam_files" ]]; then
    echo "Error: No SAM files found in $OUT_DIR."
    exit 1
fi

featureCounts -T $SLURM_CPUS_PER_TASK -p --countReadPairs --ignoreDup -B -d 1 -D 1000 --primary -O -F SAF -a "$saf_file" \
-o "$OUT_DIR/hg38.gencodev44_RepliTag_genes-1000_featureCounts_counts_noDUP.txt" \
$sam_files

featureCounts -T $SLURM_CPUS_PER_TASK -p --countReadPairs -B -d 1 -D 1000 --primary -O -F SAF -a "$saf_file" \
-o "$OUT_DIR/hg38.gencodev44_RepliTag_genes-1000_featureCounts_counts_plusDUP.txt" \
$sam_files