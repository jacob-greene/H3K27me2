•	Next-generation sequencing data (CUT&Tag) are available under GEO accession number GSE327802 as of the date of publication.
•	Processed data is available at Zenodo under https://doi.org/10.5281/zenodo.19489947.
•	Data processing and plotting scripts are publicly available on GitHub at https://github.com/jacob-greene/H3K27me2.
•	An interactive UCSC browser session with the manuscript data is available at https://genome.ucsc.edu/s/jegreene/H3K27me2_manuscript.

DNA sequencing data processing
The size distributions and molar concentration of libraries were determined using an Agilent 4200 TapeStation and libraries were mixed to achieve equal representation as desired aiming for a final concentration as recommended by the manufacturer. Paired-end 50x50 bp sequencing on the Illumina NextSeq 2000 platform was performed by the Fred Hutchinson Cancer Center Genomics Shared Resources and data were analyzed as described (https://www.protocols.io/view/cut-amp-tag-data-processing-and-analysis-tutorial-e6nvw93x7gmk/v1).
For processing sequencing data, we used cutadapt 4.472 to trim adapters from 50 bp paired-end reads with parameters “-j 8 --nextseq-trim 20 -m 20 -a
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -Z”
We used Bowtie2 2.5.173 to map the paired-end 50 bp reads to the hg38 human genome reference sequence from UCSC with parameters "--very-sensitive-local --soft-clipped-unmapped-tlen --dovetail --no-mixed --no-discordant -q --phred33 -I 10 -X 1000". 
To remove spike-in reads, hg38 positions mapped to by reads from S2 profiles were blacklisted and reads that completely overlapped with these intervals were removed using bedtools intersect. The spike-in reads were not used for data analysis.



