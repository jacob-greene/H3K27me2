#!/bin/bash

# Script to lift over BigWig files from hg19 to hg38 using CrossMap

# Set input and output directories
input_dir='data/2024/RepliSeq'
output_dir='data/2024/RepliTag/RepliSeq/hg38_bws_crossmap'
chain_file="hg19ToHg38.over.chain.gz"

# Create output directory if it doesn't exist
mkdir -p "$output_dir"

# Download chain file if not already present
if [ ! -f "$chain_file" ]; then
  echo "Downloading chain file for hg19 to hg38..."
  wget -q http://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
fi


# Loop through all BigWig files in the input directory
for bw_file in "$input_dir"/*.bigWig; do
  # Get the base filename without extension
  base_name=$(basename "$bw_file" .bigWig)

  # Define output file
  lifted_bw="$output_dir/${base_name}_hg38.bigWig"

  # Perform liftover directly on BigWig files using CrossMap
  echo "Running CrossMap for $base_name..."
  CrossMap bigwig "$chain_file" "$bw_file" "$lifted_bw"

  if [ -s "$lifted_bw" ]; then
    echo "Liftover completed for $base_name. Output saved to $lifted_bw."
  else
    echo "No regions lifted for $base_name."
  fi

done

echo "All files processed. Output saved to $output_dir."
