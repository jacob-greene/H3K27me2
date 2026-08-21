#!/usr/bin/env bash
#
# Thin driver for processing_scripts/. Run it from the repository root:
#
#     bash processing_scripts/run.sh list           # show the stages
#     bash processing_scripts/run.sh <stage>        # run one
#     DRY_RUN=1 bash processing_scripts/run.sh <stage>   # print, don't execute
#
# Each stage just invokes the underlying script(s) with the repo-root-relative paths they
# already expect. This driver adds no analysis of its own and reorders nothing.
#
# Stages are listed in pipeline order, but most need the Zenodo archive unpacked into
# data/ first (see README). `selftest` is the exception: it needs no external data.

set -uo pipefail

usage() {
    cat <<'EOF'
usage: bash processing_scripts/run.sh <stage>

  selftest    peakPRINT chain on synthetic data — no external data needed
  align       geo_to_sams.sh          GEO/SRA -> duplicate-marked SAMs   [NEVER RUN BY THE AUTHORS]
  filter      filter_sams*.sh         apply the Drosophila spike-in blacklist
  merge       merge_sams2bam_slurm.sh split All_sams4.txt into sections, merge BAMs
  tracks      bam2bed2bw*.sh          filtered SAMs -> normalised bigWigs
  annotation  gtf2bed_TSS.sh          GENCODE v44 GTF -> TSS-1000..TES BED + SAF   [needs gtftools]
  counts      Feature_counts_SAMs*.sh featureCounts over the SAF(s)
  dge         *_DGE_RUVseq*.R         RUVSeq/edgeR differential expression
  peakprint   downsample_peakPRINT.sh <inputBed> <pct> <tmpdir> <outdir>
  heatmaps    Peak_heatmaps_v2.sh     deepTools matrices + heatmaps
  rifs        RIFs_RepliTag.sh        reads-in-features matrices
  repliseq    crossmaphg19tohg38.sh   UW Repli-seq hg19 bigWigs -> hg38
  nucleation  run_replitag_nucleation_matrices_slurm.sh, then build_per_bin_covcorr.py
  blacklist   blacklist_dm6.sh        rebuild the spike-in blacklist   [LAB-INTERNAL INPUTS]

DRY_RUN=1 prints the commands instead of running them.
EOF
}

[[ $# -ge 1 ]] || { usage; exit 2; }

# Every script below uses repo-root-relative paths, so refuse to run from anywhere else.
[[ -d processing_scripts && -f README.md ]] || {
    echo "ERROR: run this from the repository root (see README)." >&2; exit 2; }

# Fail fast and propagate the real exit code. A stage that reports success while the
# script it wrapped has failed is worse than no driver at all.
run() {
    echo "+ $*" >&2
    [[ "${DRY_RUN:-0}" == "1" ]] && return 0
    "$@" || { rc=$?; echo "FAILED (exit $rc): $*" >&2; exit "$rc"; }
}

stage="$1"; shift

case "$stage" in
  list|-h|--help|help) usage; exit 0 ;;

  selftest)
    # Exercises downsample_peakPRINT.sh -> bedtools_fingerprint_slope3.sh ->
    # compute_rank_threshold.py, i.e. every intra-repo script-to-script reference, using
    # synthetic fragments and a two-chromosome sizes file so it finishes in ~1 minute.
    # It runs in a scratch directory so the repository stays clean.
    if [[ "${DRY_RUN:-0}" == "1" ]]; then
        echo "+ bash processing_scripts/downsample_peakPRINT.sh selftest.bed 100 tmpd outd" >&2
        exit 0
    fi
    for t in bedtools bc; do
        command -v "$t" >/dev/null || { echo "ERROR: $t not on PATH (envs/processing.yml)" >&2; exit 1; }
    done
    # Check the interpreter the PIPELINE will use, not this one: bedtools_fingerprint_slope3.sh
    # falls back to an Lmod Python when the active one has no numpy.
    python3 -c 'import numpy' 2>/dev/null \
        || { echo "ERROR: python3 has no numpy (envs/processing.yml); the pipeline would" >&2
             echo "       fall back to an Lmod Python that also lacks it." >&2; exit 1; }
    work="$(mktemp -d)"; trap 'rm -rf "$work"' EXIT
    ln -s "$PWD/processing_scripts" "$work/processing_scripts"
    mkdir -p "$work/data" "$work/tmpd/selftest" "$work/outd"
    printf 'chr21\t2000000\nchr22\t2000000\n' > "$work/data/chr_lens.txt"
    awk 'BEGIN{srand(7); for(i=0;i<20000;i++){c=(i%2)?"chr21":"chr22";
         s=int(rand()*1990000); l=100+int(rand()*300); print c"\t"s"\t"s+l"\tfrag"i"\t"l}}' \
      | LC_COLLATE=C sort -k1,1 -k2,2n > "$work/selftest.bed"
    ( cd "$work" && run bash processing_scripts/downsample_peakPRINT.sh selftest.bed 100 tmpd outd ) \
        || { echo "selftest FAILED" >&2; exit 1; }
    [[ "${DRY_RUN:-0}" == "1" ]] && exit 0
    # downsample_peakPRINT.sh exits 0 even when its inner pipeline fails, so check the
    # result rather than the exit code: a real run calls peaks, a broken one writes zeros.
    row="$(cat "$work/outd/selftest_100pct.tsv" 2>/dev/null)"
    [[ -n "$row" ]] || { echo "selftest FAILED: no result row written" >&2; exit 1; }
    echo "--- result row (inputBed pct reads genomeCovered FRiP nPeaks avgWidth) ---"
    printf '%s\n' "$row"
    npeaks="$(printf '%s' "$row" | cut -f6)"
    [[ "${npeaks:-0}" -gt 0 ]] || { echo "selftest FAILED: 0 peaks called" >&2; exit 1; }
    echo "selftest OK: $npeaks peaks called through processing_scripts/ -> processing_scripts/" >&2
    exit 0 ;;

  align)
    echo "NOTE: geo_to_sams.sh was transcribed from the manuscript Methods and has never" >&2
    echo "      been executed by the authors. Set BT2_INDEX and pass a sample sheet." >&2
    run bash processing_scripts/geo_to_sams.sh "$@" ;;

  filter)
    run bash processing_scripts/filter_sams.sh
    run bash processing_scripts/filter_sams_drug.sh
    run bash processing_scripts/filter_sams_drug_reseq.sh ;;

  merge)      run bash processing_scripts/merge_sams2bam_slurm.sh ;;

  tracks)
    run bash processing_scripts/bam2bed2bw.sh
    run bash processing_scripts/bam2bed2bw_drug.sh ;;

  annotation) run bash processing_scripts/gtf2bed_TSS.sh ;;

  counts)
    run bash processing_scripts/Feature_counts_SAMs.sh
    run bash processing_scripts/Feature_counts_SAMs-allFeatures_k562.sh ;;

  dge)
    run Rscript processing_scripts/251113_DGE_RUVseq.R
    run Rscript processing_scripts/251202_DGE_RUVseq_drug.R ;;

  peakprint)
    [[ $# -eq 4 ]] || {
        echo "usage: run.sh peakprint <inputBed> <percent> <tmpdir> <outdir>" >&2
        echo "note: \$tmpdir/\$(basename inputBed .bed) must already exist" >&2; exit 2; }
    run bash processing_scripts/downsample_peakPRINT.sh "$@" ;;

  heatmaps)   run bash processing_scripts/Peak_heatmaps_v2.sh ;;
  rifs)       run bash processing_scripts/RIFs_RepliTag.sh ;;
  repliseq)   run bash processing_scripts/crossmaphg19tohg38.sh ;;

  nucleation)
    # The first script submits a Slurm array and its dependent collation job, then returns.
    # build_per_bin_covcorr.py needs those matrices, so run it only after they finish.
    run bash processing_scripts/run_replitag_nucleation_matrices_slurm.sh
    echo "NOTE: matrices are built by the Slurm jobs above. When they complete, run:" >&2
    echo "      python3 processing_scripts/build_per_bin_covcorr.py" >&2 ;;

  blacklist)
    echo "NOTE: blacklist_dm6.sh reads four lab-internal files that are not redistributed." >&2
    echo "      Its output, data/Kc_merged_blacklist_hg38.bed, already ships in this repo." >&2
    run bash processing_scripts/blacklist_dm6.sh ;;

  *) echo "unknown stage: $stage" >&2; usage; exit 2 ;;
esac
