#!/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -ne 5 ]; then
    echo "Usage: $0 <bedgraph_file> <GENOME> <CHRM_SIZES> <BINWIDTH> <MIN_SIZE>"
    exit 1
fi

# Input files
BEDGRAPH_FILE=$1
GENOME=$2
CHRM_SIZES=$3
BINWIDTH=$4
MIN_SIZE=$5
BEDGRAPH_BASENAME=$(basename "$BEDGRAPH_FILE" .bedgraph)

# Output files
DOWNSAMPLED_FILE="${BEDGRAPH_BASENAME}_downsampled_bins.bed"
MERGED_INTERVALS_FILE="${BEDGRAPH_BASENAME}_merged_intervals.bed"
NORM_CUMSUM_TEMP="${BEDGRAPH_BASENAME}_normcumsum_temp.bed"
MAPPED_BED="${BEDGRAPH_BASENAME}_mapped.bed"

command -v bedtools >/dev/null || module load BEDTools/2.30.0-GCC-11.2.0

# Set LC_COLLATE for sorting
export LC_COLLATE=C

# Bin width argument
GENOMIC_BINS_FILE="${GENOME}_${BINWIDTH}bp_bins.bed"

# Ensure chromosome sizes file exists
if [[ ! -s "$CHRM_SIZES" ]]; then
    echo "Error: Chromosome sizes file ($CHRM_SIZES) is missing or empty!"
    exit 1
fi

# Check if genomic bins file exists, if not then create it
if [[ ! -s "$GENOMIC_BINS_FILE" ]]; then
    echo "Creating genomic bins file: $GENOMIC_BINS_FILE"
    bedtools makewindows -g "$CHRM_SIZES" -w "$BINWIDTH" | sort -k1,1 -k2,2n > "$GENOMIC_BINS_FILE"

    # Verify the bins file was created successfully
    if [[ ! -s "$GENOMIC_BINS_FILE" ]]; then
        echo "Error: Failed to create $GENOMIC_BINS_FILE"
        exit 1
    fi
    echo "File $GENOMIC_BINS_FILE successfully created."
else
    echo "File $GENOMIC_BINS_FILE already exists. No action taken."
fi

# Step 1: Map bedgraph values to bins
bedtools map -a "$GENOMIC_BINS_FILE" -b "$BEDGRAPH_FILE" -c 4 -o mean -null 0 > "$MAPPED_BED"
head $MAPPED_BED

# Check if MAPPED_BED was created
if [[ ! -s "$MAPPED_BED" ]]; then
    echo "Error: bedtools map did not create $MAPPED_BED"
    exit 1
fi

# Step 2: Run Python script
# Only fall back to the Lmod Python if the active interpreter cannot run the rank script.
# Loading it unconditionally would prepend it ahead of an activated conda environment
# (envs/processing.yml), and the Lmod Python/3.9.6 has no numpy.
python3 -c 'import numpy' 2>/dev/null || module load Python/3.9.6-GCCcore-11.2.0
rank_script="processing_scripts/compute_rank_threshold.py"
python3 "$rank_script" "$MAPPED_BED" "$NORM_CUMSUM_TEMP"

# Check if Python script produced output
if [[ ! -s "$NORM_CUMSUM_TEMP" ]]; then
    echo "Error: Python script did not generate $NORM_CUMSUM_TEMP"
    exit 1
fi
head $NORM_CUMSUM_TEMP
# Step 3: Find Rank Threshold
RANK_THRESHOLD=$(awk 'BEGIN {OFS="\t"} {if ($7 != "NA" && $7 > 1) {print $5; exit}}' "$NORM_CUMSUM_TEMP")

# Check if RANK_THRESHOLD is empty
if [[ -z "$RANK_THRESHOLD" ]]; then
    echo "Warning: No rank threshold found in $NORM_CUMSUM_TEMP."
    head "$NORM_CUMSUM_TEMP"  # Show first few lines for debugging
    exit 1
fi

# Print results
MAX_RANK=$(awk 'END {print NR}' "$NORM_CUMSUM_TEMP")
echo "Rank threshold is $(echo "$RANK_THRESHOLD / $MAX_RANK" | bc -l)"

# Step 5: Output all lines above rank, merge, then remove bins less than min_size
awk -v rank_threshold="$RANK_THRESHOLD" '$5 > rank_threshold' "$NORM_CUMSUM_TEMP" | \
sort -k1,1 -k2,2n | \
bedtools merge -i stdin | \
awk -v min_size="$MIN_SIZE" '$3 - $2 > min_size' > "$MERGED_INTERVALS_FILE"

# Step 6: Downsample output to yield 10,000 evenly spaced bins
total_lines=$(wc -l < "$NORM_CUMSUM_TEMP")
step=$(( total_lines / 10000 ))
if [[ "$step" -lt 1 ]]; then step=1; fi

awk -v step="$step" 'NR % step == 0' "$NORM_CUMSUM_TEMP" > "$DOWNSAMPLED_FILE"

# Cleanup temporary files unless debug flag is enabled
# rm "$MAPPED_BED" "$NORM_CUMSUM_TEMP"

echo "Pipeline completed. Downsampled output written to $DOWNSAMPLED_FILE"
echo "Merged intervals written to $MERGED_INTERVALS_FILE"
