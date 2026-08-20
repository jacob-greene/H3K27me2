#!/bin/bash

#extracts coordinates from gtf and adds chr prefix

# #TSS-1000-TES
gtftools -f 1000-0 -g processed_beds/gencode.v44.basic.TSS-1000_TES.bed gencode.v44.basic.annotation.coding.gtf
awk 'BEGIN {OFS="\t"} {$1="chr"$1; print}' processed_beds/gencode.v44.basic.TSS-1000_TES.bed | LC_COLLATE=C sort -k1,1 -k2,2n > processed_beds/sorted/gencode.v44.basic.TSS-1000_TES_sorted.bed

output_bed=processed_beds/sorted/gencode.v44.basic.TSS-1000_TES_sorted.bed
output_saf=processed_beds/sorted/gencode.v44.basic.TSS-1000_TES_sorted.saf
while IFS=$'\t' read -r chr start end strand gene_id gene_name gene_type start1 end1; do
    echo -e "$gene_name\t$chr\t$start1\t$end1\t$strand" >> "$output_saf"
done < "$output_bed"