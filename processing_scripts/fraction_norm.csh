#!/bin/csh 

#2/19/2022 output bw format instead of bg
#Normalize as a fraction of total counts
#                   1           2            3            4              5       6
#fraction_norm.csh genome.bed  scale output(bg|bga|d) genome_chr_lens  min_len max_len

#	min_len and max_len refer to the fragment lengths in genome.bed

#Format for genome_chr_lens file is (information is in the sam file headers):
#chr1	249250621
#chr2	243199373
#...

#One way to get bed files from sam alignment files is to convert sam to bam format and then run
#sam to bam:  picard SamFormatConverter
#bedtools bam2bed -bedpe

unalias rm
set sbin = /shared/ngs/illumina/henikoff/bin

#echo $#argv
if ($#argv < 6) then
	echo "USAGE fraction_norm.csh genome.bed scale output(bg|bga|d) genome_chr_lens min_len max_len"
	echo "Fraction calibration using bedtools genomecov"
	echo "Possible scale=10000 or genome size"
	echo "Fragment length must be in 5th column of genome.bed file"
	echo "min_len and max_len refer to the lengths of fragments in genome.bed to normalize"
	echo "To normalize all fragments use min_len = 1 and max_len = 1000"
	echo "Output will be placed in the current directory"
	echo "This script should be run on the gizmo system using sbatch"
	exit(-1)
endif

set genome_bed = $1
set scale = $2
set report = $3
set genome_len = $4
set min_len = $5
set max_len = $6

echo $genome_bed $scale $report $min_len $max_len
if (!(-e $genome_bed) || (-z $genome_bed)) then
	echo "$genome_bed not found or is empty"
	exit(-1)
endif
if (!(-e $genome_len) || (-z $genome_len)) then
	echo "$genome_len not found or is empty"
	exit(-1)
endif

set temp = $genome_bed:t
set name = $temp:r
set output = $name.${min_len}-${max_len}.$report

#Select rows within the fragment length range, assumes fragment length
#is in column 5 of the bed file
#$sbin/select_col.pl $genome_bed 5 $min_len $max_len > $$.temp.bed
cat $genome_bed | awk -v min=$min_len -v max=$max_len '{if ($5 >= min && $5 <= max) print}' > $$.temp.bed

#Get total counts for genomecov, col 5 has the lengths, sum them 
set totcnt = `cat $$.temp.bed | awk '{sum += $5} END {print sum}'`
set scale_factor = `echo "$scale / $totcnt" | bc -l`
echo "totcnt=$totcnt scale_factor=$scale_factor"

#ml bedtools
ml BEDTools
which bedtools
ml Kent_tools
#genomecov
#see http://bedtools.readthedocs.io/en/latest/content/tools/genomecov.html
#The first position is start-1 and the last is end for bed files, -bg and -bga
#-bga prints zero intervals -bg doesn't, -d prints each bp starting from 1 (not 0)
#bedtools doesn't recognize type bw
if ($report == "bw") then
	set bedrep = "bg"
else
	set bedrep = $report
endif
bedtools genomecov -$bedrep -scale $scale_factor -i $$.temp.bed -g $genome_len > $$.output
if ($report == "bg" || $report == "bga") then
	#This is for IGV which doesn't recognize .bg or .bga files
	set output = $name.${min_len}-${max_len}.bedgraph
	echo track type=bedGraph name=$name > $$.track
 	cat $$.track $$.output > $output
else
	if ($report == "bw") then
		set output = $name.${min_len}-${max_len}.bw
 		bedGraphToBigWig $$.output $genome_len $output
	else
		mv $$.output $output
	endif
endif
echo "Output is in file $output in the current directory"

rm $$.*
exit
track type=bedGraph name=SB_ScDm_Sth_H2A_0713 description="SB_ScDm_Sth_H2A_0713 yeast fly_spiked"
