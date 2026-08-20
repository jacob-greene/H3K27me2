#!/bin/bash

# Define directories
export data_dir="data/2024/SH_all_data"
export input_file="$data_dir/All_sams4.txt"
export sections_dir="$data_dir/sections4"  # Directory to store section files
export out_dir="$data_dir/merged_bams4"  # Output directory for merged BAM files

# Ensure sections and output directories exist
rm -r "$sections_dir"
mkdir -p "$sections_dir" "$out_dir" "$out_dir/logs"|| { echo "Failed to create necessary directories"; exit 1; }

# Read the input file line by line and create section files
while IFS= read -r line; do
    # Check if the line is a section header
    if [[ "$line" =~ ^[A-Za-z0-9_]+$ ]]; then
        # Start a new section file
        current_section_file="${sections_dir}/${line}.txt"
        touch "$current_section_file"
    elif [ -n "$current_section_file" ]; then
        # Add the file path to the current section file
        echo "$line" >> "$current_section_file"
    fi
done < "$input_file"

# Loop over each section file and submit a job to Slurm
for section_file in "$sections_dir"/*.txt; do
    base_name=$(basename "$section_file" .txt)

    # Create a Slurm batch script for each section
    cat <<EOT >"${sections_dir}/${base_name}_job.sh"
#!/bin/bash
#SBATCH --job-name=process_${base_name}
#SBATCH --output=${out_dir}/logs/${base_name}_%j.out
#SBATCH --error=${out_dir}/logs/${base_name}_%j.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G

# Load necessary modules
ml SAMtools/1.17-GCC-12.2.0

# Define directories
export data_dir="${data_dir}"
export sections_dir="${sections_dir}"
export out_dir="${out_dir}"
export section_file="${section_file}"

# Read SAM file paths from the section file
readarray -t sam_files < "\$section_file"

# Define a temporary merged BAM file path
merged_bam="\$out_dir/${base_name}_merged.bam"

# First, use samtools to merge SAM/BAM files
samtools merge -@ 8 -o "\$merged_bam" "\${sam_files[@]}"

# Then, use samtools to sort the merged BAM file and output the final sorted BAM
samtools sort -n -@ 8 -O bam -o "\$out_dir/${base_name}.bam" -T "\$out_dir/${base_name}_tmp" "\$merged_bam"

# Optionally, remove the intermediate merged BAM file to save space
rm "\$merged_bam"

echo "Processed \$section_file into \$out_dir/${base_name}.bam"
EOT

    # Submit the job to Slurm
    sbatch "${sections_dir}/${base_name}_job.sh"
    rm "${sections_dir}/${base_name}_job.sh"
done

echo "All section files submitted to Slurm."
