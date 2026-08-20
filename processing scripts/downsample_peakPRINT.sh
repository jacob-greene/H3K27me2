#!/bin/bash
#set error debugging, print each line
# set -e
# set -x

# Script arguments
inputBed=$1
currentPct=$2
tmpdir=$3
outputDir=$4

# Load required modules or software environments
module load BEDTools/2.30.0-GCC-11.2.0

# Define additional variables
GENOME=hg38
BINWIDTH=50
MIN_SIZE=300
scale=3298912062
min_len=1
max_len=1000
PEAKscript="processing scripts/bedtools_fingerprint_slope3.sh"
CHRM_SIZES="/shared/ngs/illumina/henikoff/databases/human/hg38/chr_lens.txt"
basename=$(basename "$inputBed" .bed)

# Calculate number of lines to sample based on percentage
totalLines=$(wc -l < $inputBed)
linesToSample=$(echo "($totalLines * $currentPct) / 100" | bc)

# 1. Randomly downsample and remove chrM
downsampledBed="$tmpdir/$basename/${basename}_${currentPct}pct.bed"
grep -v chrM $inputBed | shuf -n $linesToSample | LC_COLLATE=C sort -k1,1 -k2,2n > $downsampledBed
totalReads=$(cat $downsampledBed | wc -l)
echo "downsampling ${currentPct} complete."

# 2. Compute a BEDGraph file scaled to read count over the hg38 genome size for all fragments between 1-1000 bp
bedGraphFile="${downsampledBed%.bed}.${min_len}-${max_len}.bedGraph"

# Filter based on fragment length, assuming column 5 has the length (adjust if necessary)
filteredBed=$tmpdir/$basename/filtered_${currentPct}pct.bed
awk -v min=$min_len -v max=$max_len '{if ($5 >= min && $5 <= max) print}' $downsampledBed > $filteredBed

# Calculate scale factor and generate BEDGraph
totCnt=$(awk '{sum += $5} END {print sum}' $filteredBed)
scaleFactor=$(echo "$scale / $totCnt" | bc -l)
bedtools genomecov -bg -scale $scaleFactor -i $filteredBed -g $CHRM_SIZES > $bedGraphFile
echo "bedgraph ${currentPct} complete."

# 3. Run PEAKprint on the BEDGraph file
BEDGRAPH_BASENAME=$(basename "$bedGraphFile" .bedgraph)
PEAKoutput="${BEDGRAPH_BASENAME}_merged_intervals.bed"

bash $PEAKscript $bedGraphFile $GENOME $CHRM_SIZES $BINWIDTH $MIN_SIZE

echo "peakPRINT ${currentPct} complete."
echo "$PEAKoutput"

# Check if SEACR output file is empty (i.e., 0 peaks)
if [ ! -s "$PEAKoutput" ]; then
    # If SEACR output file is empty, set all dependent values to 0
    genomeCovered=0
    fractionReadsInPeaks=0
    nPeaks=0
    averagePeakWidth=0
    peakWidths=""
    echo "SEACR found 0 peaks; setting dependent values to 0."
else
    # 4. Compute the percent of the genome covered by peaks
    totalPeaksCoverage=$(bedtools merge -i "${PEAKoutput}" | awk '{sum+=$3-$2} END {print sum}')
    genomeCovered=$(echo "scale=4; $totalPeaksCoverage / $scale * 100" | bc)
    echo "coverage ${currentPct} complete."

    # 5. Compute fractions of reads in peaks
    readsInPeaks=$(bedtools intersect -u -a $downsampledBed -b "${PEAKoutput}" | wc -l)
    fractionReadsInPeaks=$(echo "scale=4; $readsInPeaks / $totalReads" | bc)
    echo "FRiPs ${currentPct} complete."

    # 6. Compute total peaks
    nPeaks=$(cat $PEAKoutput | wc -l)
    echo "nPeaks ${currentPct} complete."

    # Calculate average peak width
    totalPeakWidth=$(awk '{sum+=$3-$2} END {print sum}' $PEAKoutput)
    averagePeakWidth=$(echo "scale=2; $totalPeakWidth / $nPeaks" | bc)
    echo "Average Peak Width: $averagePeakWidth"
    
    # Calculate all peak widths and saving as a comma-separated string
    # peakWidths=$(awk '{print $3-$2}' $PEAKoutput | paste -sd, -)
    # echo "Peak Widths: $peakWidths"
fi

# Append collected data to TSV
echo -e "${inputBed}\t${currentPct}\t${totalReads}\t${genomeCovered}\t${fractionReadsInPeaks}\t${nPeaks}\t${averagePeakWidth}\t${peakWidths}" >> $outputDir/${basename}_${currentPct}pct.tsv

echo "Analysis complete. Results saved to $outputDir/${basename}_${currentPct}pct.tsv"
exit