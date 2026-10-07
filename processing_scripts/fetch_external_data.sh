#!/usr/bin/env bash
#
# Fetches the public third-party inputs the pipeline needs. Run it from the repository
# root, before any processing stage:
#
#     bash processing_scripts/fetch_external_data.sh list      # show the targets
#     bash processing_scripts/fetch_external_data.sh all       # fetch everything
#     bash processing_scripts/fetch_external_data.sh <target>  # fetch one
#
#     DRY_RUN=1   print the URLs instead of downloading
#     FORCE=1     re-download targets whose output file already exists
#
# Everything lands under data/, which .gitignore excludes, so a fetch never dirties the
# working tree. It does NOT fetch the raw reads, GEO GSE327802; `run.sh align` does.
#
# Total download is roughly 2.8 GB, dominated by the eight Repli-seq bigWigs.

set -euo pipefail

[[ -d processing_scripts && -f README.md ]] || {
    echo "ERROR: run this from the repository root (see README)." >&2; exit 2; }

GENCODE=https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release
UCSC_RS=https://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwRepliSeq
MSIGDB=https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2025.1.Hs
ENCODE=https://www.encodeproject.org/files
GEO=https://ftp.ncbi.nlm.nih.gov/geo/samples
ZENODO=https://zenodo.org/records/19489948/files   # this study's processed data, v1

# Zenodo record 19489948: "md5  file-name  destination". The record is flat; the
# destinations are where the notebooks and the DGE scripts read and write these files.
ZENODO_FILES=(
    "b227e2638cb072363bc57b48add9233f  K27ac_drug_RUVr_series.tsv  data/2024/RepliTag/Matrices/TAZEED/K27ac_drug_RUVr_series.tsv"
    "ad3bba96c76b027a139dcad043203c82  K27me1_drug_RUVr_series.tsv  data/2024/RepliTag/Matrices/TAZEED/K27me1_drug_RUVr_series.tsv"
    "632c377a53050a5f1dda1f43c99d8110  K27me2_drug_RUVr_series.tsv  data/2024/RepliTag/Matrices/TAZEED/K27me2_drug_RUVr_series.tsv"
    "f7a1267399e97fea2b797cffebc87991  K27me3_drug_RUVr_series.tsv  data/2024/RepliTag/Matrices/TAZEED/K27me3_drug_RUVr_series.tsv"
    "0d2f8c00e8709a5e3a2adf589d14fd07  P2S25p_drug_RUVr_series.tsv  data/2024/RepliTag/Matrices/TAZEED/P2S25p_drug_RUVr_series.tsv"
    "baa3b20d3c515365ad675027b0288771  K27ac_K_top95regions.bed  data/2024/SH_all_data/downsample_peakPRINT_slope1_min300/min350/top95/K27ac_K_top95regions.bed"
    "d1795d29df95cb37e7c201a6adc6e768  K27me1_K_top95regions.bed  data/2024/SH_all_data/downsample_peakPRINT_slope1_min300/min350/top95/K27me1_K_top95regions.bed"
    "5baf9dfad8cfa506ba80e893a7a613c3  K27me2_K_top95regions.bed  data/2024/SH_all_data/downsample_peakPRINT_slope1_min300/min350/top95/K27me2_K_top95regions.bed"
    "1019c35c89b8471449852fc8b706c4b8  K27me3_K_top95regions.bed  data/2024/SH_all_data/downsample_peakPRINT_slope1_min300/min350/top95/K27me3_K_top95regions.bed"
    "d3b5fee60bf9803af91fb82e7781cbd8  P2_K_top95regions.bed  data/2024/SH_all_data/downsample_peakPRINT_slope1_min300/min350/top95/P2_K_top95regions.bed"
    "4ee2c8e0839b538516329c74b70dd582  K27ac_RUVr_series.tsv  data/2024/RepliTag/filtered_sams_WT_1.0/251113_DGE/K27ac/output.tsv"
    "598cdfeaa72cd9a6dcda089ac2f6cfd8  K27me1_RUVr_series.tsv  data/2024/RepliTag/filtered_sams_WT_1.0/251113_DGE/K27me1/output.tsv"
    "4f138280dcf01d86f9a94677e61389e4  K27me2_RUVr_series.tsv  data/2024/RepliTag/filtered_sams_WT_1.0/251113_DGE/K27me2/output.tsv"
    "62695239937eb304c1d0e844bdbb5f07  K27me3_RUVr_series.tsv  data/2024/RepliTag/filtered_sams_WT_1.0/251113_DGE/K27me3/output.tsv"
    "f9d8e73b7b6b1fadb658d7b88cefd48c  P2S25p_RUVr_series.tsv  data/2024/RepliTag/filtered_sams_WT_1.0/251113_DGE/P2S25p/output.tsv"
)

REPLISEQ_FILES=(G1PctSignalRep1 G2PctSignalRep1 S1PctSignalRep1 S2PctSignalRep1
                S3PctSignalRep1 S4PctSignalRep1 SumSignalRep1 WaveSignalRep1)

usage() {
    cat <<'EOF'
usage: bash processing_scripts/fetch_external_data.sh <target> [<target> ...]

  all         every target below
  zenodo      this study's processed data (15 files, ~280 MB): the differential
              H3K27 time series (WT and drug) and the top-95% peak regions
                                     -> data/2024/RepliTag/, data/2024/SH_all_data/
  gencode27   GENCODE v27 GTF        -> data/2024/LAD_data/gencode.v27.annotation.gtf.gz
  gencode44   GENCODE v44 basic GTF, protein-coding genes only, for `run.sh annotation`
                                     -> data/2024/K562_annotations/gencode.v44.basic.annotation.coding.gtf
  msigdb      MSigDB GOBP_CELL_CYCLE -> data/2024/RepliTag/refs/GOBP_CELL_CYCLE.v2025.1.Hs.json
  repliseq    UW Repli-seq K562 hg19 bigWigs (8 files, ~2 GB), the input to
              `run.sh repliseq`      -> data/2024/RepliSeq/
  encode      ENCODE EZH2 + SUZ12 K562 fold-change bigWigs
                                     -> data/2024/MTF2_GSE164804/
  mtf2        GEO GSE164804 MTF2 shCT bedGraph (hg19)
                                     -> data/2024/MTF2_GSE164804/

DRY_RUN=1 prints the URLs instead of downloading. FORCE=1 re-downloads existing files.
EOF
}

# curl and wget are both fine; environment.yml ships wget, most systems have curl.
if command -v curl >/dev/null; then
    get() { curl -fsSL --retry 3 -o "$1" "$2"; }
elif command -v wget >/dev/null; then
    get() { wget -q -O "$1" "$2"; }
else
    echo "ERROR: neither curl nor wget on PATH (environment.yml ships wget)." >&2; exit 1
fi

# Download $2 to $1 unless it is already there. Downloads to a .part file first, so an
# interrupted fetch never leaves a truncated file that looks complete on the next run.
fetch() {
    local dest="$1" url="$2"
    if [[ -s "$dest" && "${FORCE:-0}" != "1" ]]; then
        echo "  have $dest" >&2; return 0
    fi
    if [[ "${DRY_RUN:-0}" == "1" ]]; then
        echo "  $url -> $dest" >&2; return 0
    fi
    mkdir -p "$(dirname "$dest")"
    echo "  fetching $dest" >&2
    get "$dest.part" "$url"
    mv "$dest.part" "$dest"
}

do_gencode27() {
    fetch data/2024/LAD_data/gencode.v27.annotation.gtf.gz \
          "$GENCODE"_27/gencode.v27.annotation.gtf.gz
}

do_gencode44() {
    local out=data/2024/K562_annotations/gencode.v44.basic.annotation.coding.gtf
    if [[ -s "$out" && "${FORCE:-0}" != "1" ]]; then echo "  have $out" >&2; return 0; fi
    fetch data/2024/K562_annotations/gencode.v44.basic.annotation.gtf.gz \
          "$GENCODE"_44/gencode.v44.basic.annotation.gtf.gz
    if [[ "${DRY_RUN:-0}" == "1" ]]; then return 0; fi
    # gtf2bed_TSS.sh wants the protein-coding subset. Keep the headers and every line
    # belonging to a protein-coding gene; gtftools strips the "chr" prefix itself and
    # gtf2bed_TSS.sh adds it back, so leave the GENCODE naming alone.
    echo "  filtering to protein-coding genes -> $out" >&2
    gzip -dc data/2024/K562_annotations/gencode.v44.basic.annotation.gtf.gz \
      | awk '/^#/ || /gene_type "protein_coding"/' > "$out.part"
    mv "$out.part" "$out"
}

do_msigdb() {
    local out=data/2024/RepliTag/refs/GOBP_CELL_CYCLE.v2025.1.Hs.json
    if [[ -s "$out" && "${FORCE:-0}" != "1" ]]; then echo "  have $out" >&2; return 0; fi
    # The single-gene-set endpoint on gsea-msigdb.org serves the CURRENT release, which is
    # not the v2025.1.Hs the notebook pins: 1705 gene symbols in v2025.1.Hs against 1754
    # in the release that endpoint served when this script was written.
    # Take the release-pinned collection file and cut the one gene set out of it instead.
    fetch data/2024/RepliTag/refs/c5.go.bp.v2025.1.Hs.json "$MSIGDB"/c5.go.bp.v2025.1.Hs.json
    if [[ "${DRY_RUN:-0}" == "1" ]]; then return 0; fi
    echo "  extracting GOBP_CELL_CYCLE -> $out" >&2
    python3 -c '
import json, sys
src, dst = sys.argv[1], sys.argv[2]
gs = json.load(open(src))["GOBP_CELL_CYCLE"]
json.dump({"GOBP_CELL_CYCLE": gs}, open(dst, "w"), indent=2)
print("  %d gene symbols" % len(gs["geneSymbols"]))
' data/2024/RepliTag/refs/c5.go.bp.v2025.1.Hs.json "$out"
}

do_zenodo() {
    local row md5 key dest
    for row in "${ZENODO_FILES[@]}"; do
        read -r md5 key dest <<< "$row"
        fetch "$dest" "$ZENODO/$key?download=1"
        [[ "${DRY_RUN:-0}" == "1" ]] && continue
        echo "$md5  $dest" | md5sum -c --quiet - \
          || { echo "ERROR: checksum mismatch for $dest" >&2; exit 1; }
    done
}

do_repliseq() {
    local f
    for f in "${REPLISEQ_FILES[@]}"; do
        fetch "data/2024/RepliSeq/wgEncodeUwRepliSeqK562${f}.bigWig" \
              "$UCSC_RS/wgEncodeUwRepliSeqK562${f}.bigWig"
    done
}

do_encode() {
    # Both are already GRCh38 fold-change-over-control bigWigs; only the name changes.
    fetch data/2024/MTF2_GSE164804/ENCODE_EZH2_K562_ENCFF587SWK_fc_hg38.bw \
          "$ENCODE"/ENCFF587SWK/@@download/ENCFF587SWK.bigWig
    fetch data/2024/MTF2_GSE164804/ENCODE_SUZ12_K562_ENCFF065KGU_fc_hg38.bw \
          "$ENCODE"/ENCFF065KGU/@@download/ENCFF065KGU.bigWig
}

do_mtf2() {
    fetch data/2024/MTF2_GSE164804/GSM5019814_ChIPseq_MTF2_shCT.bedgraph.gz \
          "$GEO"/GSM5019nnn/GSM5019814/suppl/GSM5019814_ChIPseq_MTF2_shCT.bedgraph.gz
    # GSE164804 is mapped to hg19. build_per_bin_covcorr.py and 260801_nucleation_pub.ipynb
    # read data/2024/MTF2_GSE164804/MTF2_shCT_hg38_coverage.bw, which is this bedGraph
    # lifted to hg38 and converted to bigWig. This script does not produce that lifted
    # copy, and it is not in the Zenodo record.
    echo "  NOTE: MTF2_shCT_hg38_coverage.bw (hg19->hg38 liftOver of the file above) is" >&2
    echo "        not produced here." >&2
}

[[ $# -ge 1 ]] || { usage; exit 2; }

for target in "$@"; do
    case "$target" in
      list|-h|--help|help) usage; exit 0 ;;
      all)
        for t in zenodo gencode27 gencode44 msigdb repliseq encode mtf2; do
            echo "== $t" >&2; "do_$t"
        done ;;
      zenodo|gencode27|gencode44|msigdb|repliseq|encode|mtf2)
        echo "== $target" >&2; "do_$target" ;;
      *) echo "unknown target: $target" >&2; usage; exit 2 ;;
    esac
done
