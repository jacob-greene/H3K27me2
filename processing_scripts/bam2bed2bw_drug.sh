#!/bin/bash
set -euo pipefail

# Directories
sam_files_dir='data/2024/RepliTag/filtered_sams_drug_1.0'
sam_files_dir1='data/2024/RepliTag/filtered_sams_drug_1.0/reseq'

bw_files_dir='data/2024/RepliTag/filtered_sams_drug_1.0_merged'
mkdir -p "$bw_files_dir" "$bw_files_dir/scripts" logs

# Genome parameters
export scale=3210038034  # hg38 genome size
export flens="data/chr_lens.txt"
export min_len=1
export max_len=1000

# Extract sample base names (e.g., JG_HsDm_K562_K27ac_241101_G1_d2pcr)
samples=$(find "$sam_files_dir" -name "*_[ab].DTmarked_filtered.sam" | \
    sed -E 's|.*/(.+)_([ab])\.DTmarked_filtered\.sam|\1|' | sort -u)

# Process each sample
for sample in $samples; do
    export sam_a="$sam_files_dir/${sample}_a.DTmarked_filtered.sam"
    export sam_b="$sam_files_dir/${sample}_b.DTmarked_filtered.sam"
    export sam_a1="$sam_files_dir1/${sample}_a.DTmarked_filtered.sam"   # resequenced (optional)
    export sam_b1="$sam_files_dir1/${sample}_b.DTmarked_filtered.sam"   # resequenced (optional)

    export bed_a="$bw_files_dir/${sample}_a.bed"
    export bed_b="$bw_files_dir/${sample}_b.bed"
    export bed_a1="$bw_files_dir/${sample}_a_reseq.bed"
    export bed_b1="$bw_files_dir/${sample}_b_reseq.bed"

    export merged_bed="$bw_files_dir/${sample}_merged_no_chrM.bed"
    export slurm_script="$bw_files_dir/scripts/${sample}_merged_bam2bw.sh"

    # Create SLURM script
    cat > "$slurm_script" << 'EOF'
#!/bin/bash
#SBATCH --job-name=bed2bam2bw
#SBATCH --output=logs/bed2bam2bw_%j.out
#SBATCH --time=01:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G

set -euo pipefail

ml BEDTools/2.30.0-GCC-11.2.0
ml SAMtools/1.17-GCC-12.2.0

# helper: SAM -> fragment BED (chr, start, end, strand, len), drop chrM
sam2bed () {
  local in_sam="$1"
  local out_bed="$2"
  samtools view -B -h "$in_sam" \
    | bedtools bamtobed -bedpe -i stdin \
    | awk -v OFS='\t' '$1!="chrM"{
        s=$2; e=$6; if ($5 < s) s=$5; if ($3 > e) e=$3;
        len=e-s; print $1,s,e,$9,len
      }' \
    | sort -k1,1 -k2,2n -k3,3n > "$out_bed"
}

# Convert original A/B
sam2bed "$sam_a" "$bed_a"
sam2bed "$sam_b" "$bed_b"

# If resequenced SAMs exist, convert them too (otherwise skip)
if [[ -f "$sam_a1" ]]; then sam2bed "$sam_a1" "$bed_a1"; fi
if [[ -f "$sam_b1" ]]; then sam2bed "$sam_b1" "$bed_b1"; fi

# Concatenate whatever exists and sort once
tmp_cat="$(mktemp)"
: > "$tmp_cat"
cat "$bed_a" >> "$tmp_cat"
cat "$bed_b" >> "$tmp_cat"
[[ -f "$bed_a1" ]] && cat "$bed_a1" >> "$tmp_cat"
[[ -f "$bed_b1" ]] && cat "$bed_b1" >> "$tmp_cat"

sort -k1,1 -k2,2n -k3,3n "$tmp_cat" > "$merged_bed"
rm -f "$tmp_cat"

# Make normalized BigWig
processing_scripts/fraction_norm.csh "$merged_bed" $scale bw $flens $min_len $max_len
# processing_scripts/fraction_norm.csh "$merged_bed" $scale bg $flens $min_len $max_len
EOF

    # Submit SLURM job
    sbatch "$slurm_script"
done
