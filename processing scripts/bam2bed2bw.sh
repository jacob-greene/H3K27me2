#!/bin/bash

# Directories
sam_files_dir='data/2024/RepliTag/filtered_sams_WT_1.0'

bw_files_dir='data/2024/RepliTag/filtered_sams_WT_1.0_merged'
mkdir -p "$sam_files_dir" "$bw_files_dir" "$bw_files_dir/scripts" logs

# Genome parameters
export scale=3210038034  # hg38 genome size
export flens="/shared/ngs/illumina/henikoff/databases/human/hg38/chr_lens.txt"
export min_len=1
export max_len=1000

# Extract sample base names (e.g., JG_HsDm_K562_K27ac_241101_G1_d2pcr)
samples=$(find "$sam_files_dir" -name "*_[ab].DTmarked_filtered.sam" | \
    sed -E 's|.*/(.+)_([ab])\.DTmarked_filtered\.sam|\1|' | sort -u)

# Process each sample
for sample in $samples; do
    export sam_a="$sam_files_dir/${sample}_a.DTmarked_filtered.sam"
    export sam_b="$sam_files_dir/${sample}_b.DTmarked_filtered.sam"
    export bed_a="$sam_files_dir/${sample}_a.bed"
    export bed_b="$sam_files_dir/${sample}_b.bed"
    export merged_bed="$bw_files_dir/${sample}_merged_no_chrM.bed"
    export slurm_script="$bw_files_dir/scripts/${sample}_merged_bam2bw.sh"

    # Create SLURM script
    cat > "$slurm_script" << 'EOF'
#!/bin/bash
#SBATCH --job-name=bed2bam2bw
#SBATCH --output=logs/bed2bam2bw_%j.out
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G

ml BEDTools/2.30.0-GCC-11.2.0
ml SAMtools/1.17-GCC-12.2.0

# Convert SAM to BED (BEDPE -> fragment)
samtools view -B -h "$sam_a" | bedtools bamtobed -bedpe -i stdin | cut -f1,2,6,7 | grep -v chrM | sort -k1,1 -k2n,2n -k3n,3n | awk -v OFS='\t' '{len = $3 - $2 ; print $0, len }' > $bed_a
samtools view -B -h "$sam_b" | bedtools bamtobed -bedpe -i stdin | cut -f1,2,6,7 | grep -v chrM | sort -k1,1 -k2n,2n -k3n,3n | awk -v OFS='\t' '{len = $3 - $2 ; print $0, len }' > $bed_b

# Merge and sort A and B
cat "$bed_a" "$bed_b" | sort -k1,1 -k2,2n -k3,3n > "$merged_bed"

# Generate normalized BigWig
/shared/ngs/illumina/henikoff/bin/aligners/fraction_norm.csh "$merged_bed" $scale bw $flens $min_len $max_len
# /shared/ngs/illumina/henikoff/bin/aligners/fraction_norm.csh "$merged_bed" $scale bg $flens $min_len $max_len

EOF

    # Submit SLURM job
    sbatch "$slurm_script"
done
