Kc1=/shared/ngs/illumina/henikoff/241106_bowtie2/JG_HsDm/JG_HsDm_K562_K27ac_241101_Kc_d2pcr.bed
Kc2=/shared/ngs/illumina/henikoff/250605_bowtie2/BT_Dm/hg38/BT_Dm_3836.bed
Kc3=/shared/ngs/illumina/henikoff/211210_bowtie2/BT_Dm/hg38/BT_Dm_BT419_H29.bed
Kc4=/shared/ngs/illumina/henikoff/211210_bowtie2/BT_Dm/hg38/BT_Dm_BT420_H29.bed
OUT=data/2024/RepliTag/refs
module load BEDTools/2.31.0-GCC-12.3.0

# Clean, combine, sort, and merge
cat "$Kc1" "$Kc2" "$Kc3" "$Kc4" | \
  awk 'NF>=3 {print $1, $2, $3}' OFS='\t' | \
  sort -k1,1 -k2,2n | \
  bedtools merge -i stdin | \
  awk 'BEGIN{OFS="\t"} {print $1, $2, $3, $3 - $2}' > "${OUT}/Kc_merged_blacklist_hg38.bed"