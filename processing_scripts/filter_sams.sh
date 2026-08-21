#!/bin/bash
#SBATCH --job-name=filter_blacklist
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --output=logs/filter_blacklist_%j.log

# Load required modules
command -v samtools >/dev/null || module load SAMtools/1.19.2-GCC-13.2.0
command -v bedtools >/dev/null || module load BEDTools/2.30.0-GCC-12.2.0
command -v parallel >/dev/null || module load parallel

# Set paths
export BLACKLIST="data/Kc_merged_blacklist_hg38.bed"
export OUT_DIR="data/2024/RepliTag/filtered_sams_WT_1.0"
# Duplicate-marked bowtie2 SAMs. geo_to_sams.sh writes them here by default; override
# INPUT_DIR in the environment to point at your own.
export INPUT_DIR="${INPUT_DIR:-data/2024/RepliTag/geo_sams}"
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
find "$INPUT_DIR" -name "*.sam" | \
  parallel --bar -j $SLURM_CPUS_PER_TASK filter_sam {}

# # Count reads in genes using featureCounts
command -v featureCounts >/dev/null || module load Subread/2.0.3-GCC-11.2.0

# # Generate list of SAM files
sam_files=$(ls "$OUT_DIR"/*.sam 2>/dev/null | grep -v 'Kc')
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

# Reshape featureCounts' output into what 251113_DGE_RUVseq.R reads.
#
# featureCounts writes a leading '#' comment line, names each count column with the full
# path of its SAM, and calls its own summary '<counts>.txt.summary'. The R script reads
# one table per mark, with the six featureCounts meta columns (it needs Geneid and Length)
# and count columns named for the bare sample, plus the summary as '.txt.summary.tsv'.
counts_txt="$OUT_DIR/hg38.gencodev44_RepliTag_genes-1000_featureCounts_counts_plusDUP.txt"

mv "$counts_txt.summary" "$counts_txt.summary.tsv"

for mark in P2S25p K27me3 K27me2 K27me1 K27ac; do
  awk -v OFS='\t' -v mark="$mark" '
    /^#/ { next }
    !seen_header {
      seen_header = 1
      for (i = 1; i <= NF; i++) {
        name = $i
        sub(/^.*\//, "", name)                        # drop the directory
        sub(/\.DTmarked_filtered\.sam$/, "", name)    # drop the SAM suffix
        if (i <= 6 || name ~ ("_" mark "_")) { keep[++n] = i; hdr[i] = name }
      }
      line = hdr[keep[1]]
      for (k = 2; k <= n; k++) line = line OFS hdr[keep[k]]
      print line
      next
    }
    {
      line = $(keep[1])
      for (k = 2; k <= n; k++) line = line OFS $(keep[k])
      print line
    }
  ' "$counts_txt" > "${counts_txt%.txt}_${mark}.tsv"
done