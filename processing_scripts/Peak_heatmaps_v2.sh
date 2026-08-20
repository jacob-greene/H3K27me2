#!/bin/bash
#SBATCH --job-name=Peak_heatmap   # Job name
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16            # Number of CPU cores per task
#SBATCH --mem=32G                     # Job memory request
#SBATCH --time=48:00:00               # Time limit hrs:min:sec
#SBATCH --output=logs/TSS_1.5-2.5kb_heatmap_%j.log   # Standard output and error log

export bw_dir="data/2024/SH_all_data/merged_bws4"
export peak_dir='data/2024/SH_all_data/downsample_peakPRINT_slope1_min300/min350/top95'
export fig_dir="figures/"
export mat_dir="data/2024/SH_all_data/deeptools_mats"
mkdir -p $mat_dir

module load deepTools/3.5.4.post1-gfbf-2022b

# run compute matrix to collect the data needed for plotting
computeMatrix reference-point -S $bw_dir/K27me3_K.1-1000.bw $bw_dir/K27me2_K.1-1000.bw $bw_dir/K27me1_K.1-1000.bw $bw_dir/K27ac_K.1-1000.bw $bw_dir/P2_K.1-1000.bw \
                              -R $peak_dir/K27me3_K_top95regions.bed $peak_dir/K27me2_K_top95regions.bed \
                              $peak_dir/K27me1_K_top95regions.bed $peak_dir/K27ac_K_top95regions.bed $peak_dir/P2_K_top95regions.bed \
                              --referencePoint center -b 5000 -a 5000 \
                              --outFileSortedRegions $mat_dir/K_top95regions_10kb.bed \
                              --missingDataAsZero -o $mat_dir/K_top95regions_10kb.mat.gz -p 16  

plotHeatmap -m $mat_dir/K_top95regions_10kb.mat.gz \
            --plotFileFormat pdf --linesAtTickMarks \
            --perGroup \
            --colorList '#ffffff, #C73131' '#ffffff, #000000' '#ffffff, #FF00FF' '#ffffff, #22D0FB' '#ffffff, #FFC60B' \
            --missingDataColor "#ffffff" \
            --averageTypeSummaryPlot median \
            -o $fig_dir/20250418_ALL_K_top95regions_10kb.pdf

plotProfile -m $mat_dir/K_top95regions_10kb.mat.gz  \
    -out $fig_dir/20250418_ALL_K_top95regions_10kb_profile.pdf \
    --perGroup --numPlotsPerRow 5 \
    --plotFileFormat pdf --yMax 11 \
    --averageType median --plotWidth 9 --plotHeight 10 \
    --colors '#C73131' '#000000' '#FF00FF' '#22D0FB' '#FFC60B'
  