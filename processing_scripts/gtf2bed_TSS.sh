#!/bin/bash

# Extracts TSS-1000..TES coordinates from the GENCODE v44 basic GTF, adds the "chr" prefix,
# and emits the sorted SAF that filter_sams*.sh and Feature_counts_SAMs.sh consume and the
# sorted BED that RIFs_RepliTag.sh and two notebooks read.
#
# Run from the repository root. Requires gtftools on PATH.
#
# Input:  data/2024/K562_annotations/gencode.v44.basic.annotation.coding.gtf   (see README)
# Output: data/2024/K562_annotations/processed_beds/sorted/
#             gencode.v44.basic.TSS-1000_TES_sorted.bed
#             gencode.v44.basic.TSS-1000_TES_sorted.saf

set -euo pipefail

anno_dir="data/2024/K562_annotations"
gtf="$anno_dir/gencode.v44.basic.annotation.coding.gtf"
bed_dir="$anno_dir/processed_beds"
sorted_dir="$bed_dir/sorted"

[[ -s "$gtf" ]] || { echo "ERROR: missing $gtf" >&2; exit 1; }
mkdir -p "$bed_dir" "$sorted_dir"

# #TSS-1000-TES
# gtftools writes one row per gene: chr, gene start, gene end, strand, gene_id, gene_name,
# gene_type, then the TSS-1000 start and TES end.
gtftools -f 1000-0 -g "$bed_dir/gencode.v44.basic.TSS-1000_TES.bed" "$gtf"
gtftools_sorted="$bed_dir/gencode.v44.basic.TSS-1000_TES_gtftools_sorted.bed"
awk 'BEGIN {OFS="\t"} {$1="chr"$1; print}' "$bed_dir/gencode.v44.basic.TSS-1000_TES.bed" \
  | LC_COLLATE=C sort -k1,1 -k2,2n > "$gtftools_sorted"

output_bed="$sorted_dir/gencode.v44.basic.TSS-1000_TES_sorted.bed"
output_saf="$sorted_dir/gencode.v44.basic.TSS-1000_TES_sorted.saf"

# Truncate first: the loop below appends, so re-running would otherwise duplicate every row.
: > "$output_saf"
while IFS=$'\t' read -r chr start end strand gene_id gene_name gene_type start1 end1; do
    echo -e "$gene_name\t$chr\t$start1\t$end1\t$strand" >> "$output_saf"
done < "$gtftools_sorted"

# The BED the notebooks read (251203_*, 260803_*) has the TSS-1000..TES interval in
# BED6 order, then gene_id, gene_name, gene_type:
#     chr, TSS-1000 start, TES end, (empty name), (empty score), strand, gene_id, gene_name, gene_type
# This matches the published file byte for byte.
awk -F'\t' 'BEGIN {OFS="\t"} {print $1, $8, $9, "", "", $4, $5, $6, $7}' \
  "$gtftools_sorted" > "$output_bed"
