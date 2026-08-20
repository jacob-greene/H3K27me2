#!/bin/bash
#SBATCH --job-name=geo_to_sams
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=logs/geo_to_sams_%A_%a.log
#
# ===========================================================================================
#  GEO -> duplicate-marked SAMs.  PROVIDED FOR CONVENIENCE — NOT RUN, NOT TESTED.
# ===========================================================================================
#
#  This script was written from the methods described in `paper.md` so that someone starting
#  from GEO `GSE327802` can reach the input that `filter_sams.sh` expects. It is NOT the
#  script the authors ran: the manuscript data were trimmed and aligned by a lab-internal
#  pipeline that is not part of this repository. **It has never been executed.** Treat it
#  as a transcription of the published methods into runnable form, and check its output
#  before trusting anything downstream of it.
#
#  Every place where `paper.md` does not determine a parameter is marked
#  `# AMBIGUOUS (paper):` with the choice made and why. Grep for that string to find them all:
#
#      grep -n 'AMBIGUOUS (paper)' processing_scripts/geo_to_sams.sh
#
#  ---------------------------------------------------------------------------------------
#  PIPELINE
#      fastq-dump          SRA run -> paired FASTQ
#      cutadapt 4.4        adapter + NextSeq quality trim   (parameters quoted in paper.md)
#      bowtie2 2.5.1       -> hg38                          (parameters quoted in paper.md)
#      samtools            coordinate-sort + mark duplicates
#
#  OUTPUT
#      $OUT_DIR/<sample>.DTmarked.sam   — duplicate-MARKED, not deduplicated
#
#      That is the input to `filter_sams.sh`, which removes the Drosophila spike-in
#      blacklist. Point its `INPUT_DIR` at `$OUT_DIR` and it will pick these up. Downstream
#      of that, `bam2bed2bw.sh` globs `*_[ab].DTmarked_filtered.sam`, so <sample> must end
#      in `_a` or `_b` — see the sample sheet below.
#
#  REQUIRED INPUTS
#      1. $SAMPLE_SHEET  — TSV, no header, two columns:
#
#             <SRA run accession>   <sample>
#             SRRxxxxxxx            JG_HsDm_K562_K27me3_241101_G1_d2pcr_a
#             SRRxxxxxxx            JG_HsDm_K562_K27me3_241101_G1_d2pcr_b
#
#         Build it from the SRA Run Selector for GSE327802:
#             https://www.ncbi.nlm.nih.gov/Traces/study/?acc=GSE327802
#         or, with NCBI Entrez Direct installed:
#             esearch -db sra -query GSE327802 | efetch -format runinfo > runinfo.csv
#
#         AMBIGUOUS (paper): the GEO sample titles are not reproduced in paper.md, so the
#         mapping from run accession to the <sample> names the rest of this repository uses
#         (mark, date, cell-cycle fraction, replicate a/b) cannot be derived here. You must
#         supply it. Getting it wrong silently mislabels every downstream figure.
#
#      2. $BT2_INDEX — a bowtie2 index of hg38.
#         paper.md: "the hg38 human genome reference sequence from UCSC".
#         AMBIGUOUS (paper): it does not say whether alt/random/chrUn contigs were included,
#         nor which UCSC snapshot. The blacklist shipped in this repository
#         (data/Kc_merged_blacklist_hg38.bed) is on primary `chr*` contigs, and every
#         downstream script works in `chr`-prefixed hg38 coordinates. Build the index from
#         the UCSC `hg38.fa` and keep the UCSC contig names.
#
#  DEPENDENCIES
#      sra-tools (fastq-dump), cutadapt 4.4, bowtie2 2.5.1, samtools
#      The `module load` lines below are Fred Hutch Lmod names; substitute your own.
#
#  USAGE
#      # one Slurm task per sample sheet line
#      sbatch --array=1-$(wc -l < samples.tsv) processing_scripts/geo_to_sams.sh samples.tsv
#
#      # or, without a scheduler, every line in sequence
#      bash processing_scripts/geo_to_sams.sh samples.tsv
# ===========================================================================================

set -euo pipefail

SAMPLE_SHEET="${1:-samples.tsv}"
OUT_DIR="${OUT_DIR:-data/2024/RepliTag/geo_sams}"
FASTQ_DIR="${FASTQ_DIR:-data/2024/RepliTag/geo_fastq}"
BT2_INDEX="${BT2_INDEX:?set BT2_INDEX to the hg38 bowtie2 index prefix}"
THREADS="${THREADS:-${SLURM_CPUS_PER_TASK:-8}}"   # paper.md quotes cutadapt -j 8
KEEP_FASTQ="${KEEP_FASTQ:-0}"                     # 1 = keep the trimmed/raw FASTQs

module load SRA-Toolkit 2>/dev/null || true
module load cutadapt 2>/dev/null || true
module load Bowtie2 2>/dev/null || true
module load SAMtools 2>/dev/null || true

for tool in fastq-dump cutadapt bowtie2 samtools; do
    command -v "$tool" >/dev/null || { echo "ERROR: $tool not on PATH" >&2; exit 1; }
done
[[ -s "$SAMPLE_SHEET" ]] || { echo "ERROR: no sample sheet at $SAMPLE_SHEET" >&2; exit 1; }

mkdir -p "$OUT_DIR" "$FASTQ_DIR" logs

# A sample name must be unique: it becomes an output filename, and a collision would have
# one run silently overwrite another.
if [[ $(cut -f2 "$SAMPLE_SHEET" | sort | uniq -d | wc -l) -gt 0 ]]; then
    echo "ERROR: duplicate sample names in $SAMPLE_SHEET:" >&2
    cut -f2 "$SAMPLE_SHEET" | sort | uniq -d >&2
    exit 1
fi
# AMBIGUOUS (paper): where a GEO sample is split across several SRA runs, the methods do not
# say whether the runs were concatenated before trimming or aligned separately and merged.
# This script assumes ONE run per sample, which the uniqueness check above enforces. If your
# runinfo has several runs for a sample, concatenate their FASTQs first and give the merged
# pair a single sample name.

process_one() {
    local srr="$1" sample="$2"
    local sam="$OUT_DIR/${sample}.DTmarked.sam"

    if [[ -s "$sam" ]]; then
        echo "[$sample] $sam exists, skipping"
        return 0
    fi
    echo "[$sample] === $srr ==="

    # -- 1. GEO/SRA -> paired FASTQ ---------------------------------------------------------
    # --split-3 puts a truly-paired run in _1/_2 and any orphan reads in a third file; a run
    # that is not paired therefore produces no _1/_2 and is caught below rather than being
    # aligned as if it were single-end. paper.md: "Paired-end 50x50 bp".
    local r1="$FASTQ_DIR/${srr}_1.fastq.gz" r2="$FASTQ_DIR/${srr}_2.fastq.gz"
    if [[ ! -s "$r1" || ! -s "$r2" ]]; then
        fastq-dump --split-3 --gzip --outdir "$FASTQ_DIR" "$srr"
    fi
    [[ -s "$r1" && -s "$r2" ]] || {
        echo "ERROR: $srr did not yield a paired FASTQ ($r1 / $r2)" >&2; return 1; }

    # -- 2. cutadapt ------------------------------------------------------------------------
    # Verbatim from paper.md: cutadapt 4.4, "-j 8 --nextseq-trim 20 -m 20
    #   -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -Z"
    # -j is taken from $THREADS here so the script respects its Slurm allocation; the paper's
    # value (8) is the default. -Z is cutadapt's fast-compression shorthand and affects only
    # the speed/size of the intermediate FASTQ, not the reads.
    local t1="$FASTQ_DIR/${sample}_trim_1.fastq.gz" t2="$FASTQ_DIR/${sample}_trim_2.fastq.gz"
    cutadapt -j "$THREADS" --nextseq-trim 20 -m 20 \
        -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
        -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
        -Z -o "$t1" -p "$t2" "$r1" "$r2" \
        > "logs/${sample}.cutadapt.log"

    # -- 3. bowtie2 -------------------------------------------------------------------------
    # Verbatim from paper.md: Bowtie2 2.5.1, "--very-sensitive-local
    #   --soft-clipped-unmapped-tlen --dovetail --no-mixed --no-discordant -q --phred33
    #   -I 10 -X 1000".
    # No spike-in genome here on purpose: paper.md removes Drosophila reads AFTER alignment,
    # by blacklisting the hg38 positions that Drosophila profiles map to
    # (data/Kc_merged_blacklist_hg38.bed, applied by filter_sams.sh), not by aligning to a
    # combined index.
    bowtie2 -p "$THREADS" -x "$BT2_INDEX" -1 "$t1" -2 "$t2" \
        --very-sensitive-local --soft-clipped-unmapped-tlen --dovetail \
        --no-mixed --no-discordant -q --phred33 -I 10 -X 1000 \
        -S "$OUT_DIR/${sample}.raw.sam" \
        2> "logs/${sample}.bowtie2.log"

    # -- 4. mark duplicates -----------------------------------------------------------------
    # AMBIGUOUS (paper): paper.md does not mention duplicate handling at all. That duplicates
    # are MARKED and not removed is inferred from this repository, not from the paper:
    # `filter_sams.sh` consumes SAMs from a directory called `marked/`, its filenames carry a
    # `DTmarked` token, and it runs featureCounts twice — once with `--ignoreDup` and once
    # without — which is only meaningful if the 0x400 flag is set and the reads are still
    # present. The tool and its parameters are NOT recoverable; `samtools markdup` is used
    # here as the plain equivalent. If the original used Picard MarkDuplicates its optical-
    # duplicate handling will differ slightly, which shifts the `noDUP` counts but not the
    # `plusDUP` ones.
    samtools collate -@ "$THREADS" -O -u "$OUT_DIR/${sample}.raw.sam" \
      | samtools fixmate -@ "$THREADS" -m -u - - \
      | samtools sort -@ "$THREADS" -u - \
      | samtools markdup -@ "$THREADS" - "$sam.tmp" -O SAM \
            -f "logs/${sample}.markdup.log"
    mv "$sam.tmp" "$sam"

    rm -f "$OUT_DIR/${sample}.raw.sam"
    [[ "$KEEP_FASTQ" == "1" ]] || rm -f "$t1" "$t2" "$r1" "$r2"
    echo "[$sample] wrote $sam"
}

if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    line=$(awk -v n="$SLURM_ARRAY_TASK_ID" 'NR == n' "$SAMPLE_SHEET")
    [[ -n "$line" ]] || { echo "ERROR: no line $SLURM_ARRAY_TASK_ID in $SAMPLE_SHEET" >&2; exit 1; }
    IFS=$'\t' read -r srr sample <<< "$line"
    process_one "$srr" "$sample"
else
    while IFS=$'\t' read -r srr sample; do
        [[ -z "${srr:-}" || "$srr" == \#* ]] && continue
        process_one "$srr" "$sample"
    done < "$SAMPLE_SHEET"
fi

echo "done. Next: point filter_sams.sh's INPUT_DIR at $OUT_DIR"
