#!/bin/bash
set -euo pipefail

PROJECT="."
BED_DIR="$PROJECT/data/2024/RepliTag/filtered_sams_WT_1.0"
RT_BIGWIG="$PROJECT/data/2024/RepliTag/RepliSeq/hg38_bws_crossmap/wgEncodeUwRepliSeqK562WaveSignalRep1_hg38.bigWig.bw"
OUT="$PROJECT/data/results/replitag_sphase_nucleation"
CONFIG_DIR="$OUT/config"
TILE_DIR="$OUT/tiles"
MERGED_DIR="$OUT/merged_sample_beds"
COV_DIR="$OUT/coverage_columns"
MATRIX_DIR="$OUT/matrices"
TMP_DIR="$OUT/tmp_collate"
DOMAIN_META_DIR="$OUT/domain_metadata"
TILE_SIZES="${TILE_SIZES:-1000}"
DOMAIN_MERGE_DISTANCE="${DOMAIN_MERGE_DISTANCE:-10000}"

mkdir -p "$CONFIG_DIR" "$TILE_DIR" "$MERGED_DIR" "$COV_DIR" "$MATRIX_DIR" "$TMP_DIR" "$DOMAIN_META_DIR"

sample_manifest="$CONFIG_DIR/samples.tsv"
task_manifest="$CONFIG_DIR/matrix_tasks.tsv"

usage() {
  cat <<'USAGE'
Usage:
  TILE_SIZES="2000 5000" replitag_nucleation_generate_matrices.sh prepare
  replitag_nucleation_generate_matrices.sh coverage TASK_ID
  replitag_nucleation_generate_matrices.sh collate

Modes:
  prepare       Build tiled domain BEDs, RT metadata, and sample/task manifests.
  coverage      Run one Slurm-array task: one domain set x one replicate BED.
  collate       Join per-sample read-count columns into wide count and RPKM matrices.

The script expects BEDTools and deepTools to be available in PATH. The Slurm wrapper loads them.
USAGE
}

make_tiles() {
  local domain_name="$1"
  local domains_bed="$2"
  local tile_size="$3"
  local out_bed="$TILE_DIR/${domain_name}.${tile_size}bp.tiles.bed"

  awk -v OFS='\t' -v size="$tile_size" '
    $1 !~ /^#/ && NF >= 3 {
      chr=$1; start=$2; end=$3;
      if (end <= start) next;
      domain_id=chr "_" start "_" end;
      idx=0;
      for (s=start; s<end; s+=size) {
        e=s+size;
        if (e>end) e=end;
        printf "%s\t%d\t%d\t%s\t%s__bin%06d\t%d\t%d\n", chr, s, e, domain_id, domain_id, idx, idx, size;
        idx++;
      }
    }
  ' "$domains_bed" > "$out_bed"
}

write_domains_bed() {
  local domain_name="$1"
  local domain_bed="$2"
  local out_bed="$DOMAIN_META_DIR/${domain_name}.domains.bed"
  local tmp_bed="$DOMAIN_META_DIR/${domain_name}.domains.unmerged.tmp.bed"

  awk -v OFS='\t' '
    $1 !~ /^#/ && NF >= 3 {
      chr=$1; start=$2; end=$3;
      if (end <= start) next;
      print chr, start, end;
    }
  ' "$domain_bed" | sort -u -k1,1 -k2,2n > "$tmp_bed"

  bedtools merge -i "$tmp_bed" -d "$DOMAIN_MERGE_DISTANCE" \
    | awk -v OFS='\t' '
      {
        domain_id=$1 "_" $2 "_" $3;
        print $1, $2, $3, domain_id;
      }
    ' > "$out_bed"

  rm -f "$tmp_bed"
}

assign_rt() {
  local domain_name="$1"
  local domains_bed="$DOMAIN_META_DIR/${domain_name}.domains.bed"
  local rt_npz="$DOMAIN_META_DIR/${domain_name}.rt.npz"
  local rt_raw="$DOMAIN_META_DIR/${domain_name}.rt.raw.tsv"
  local rt_meta="$DOMAIN_META_DIR/${domain_name}.domains.rt.tsv"

  multiBigwigSummary BED-file \
    --BED "$domains_bed" \
    --bwfiles "$RT_BIGWIG" \
    --labels RT \
    --outFileName "$rt_npz" \
    --outRawCounts "$rt_raw"

  awk -v OFS='\t' '
    BEGIN {
      split("0.16605208 0.27145242 0.37685275 0.48225309 0.58765342 0.69305376 0.79845409", cut, " ");
      split("Early RT|Mid-early RT|Mid RT|Late-mid RT|Mid-late RT|Late RT", lab, "|");
      print "domain_id", "chr", "start", "end", "RT_raw", "RT", "RT_inverted", "RT_bin", "RT_label";
    }
    NR == FNR {
      id[$1 FS $2 FS $3]=$4;
      next;
    }
    FNR == 1 { next; }
    {
      key=$1 FS $2 FS $3;
      raw=$4;
      rt=raw/100;
      if (rt > 1) rt=1;
      if (rt < 0) rt=0;
      inv=1-rt;
      bin=1;
      for (i=2; i<=6; i++) {
        if (inv > cut[i]) bin=i;
      }
      print id[key], $1, $2, $3, raw, rt, inv, bin, lab[bin];
    }
  ' "$domains_bed" "$rt_raw" > "$rt_meta"
}

remove_old_tile_sizes() {
  find "$TILE_DIR" -maxdepth 1 -type f \( -name '*.25bp.*' -o -name '*.250bp.*' -o -name '*.1500bp.*' \) -delete
  find "$MATRIX_DIR" -maxdepth 1 -type f \( -name '*.25bp.*' -o -name '*.250bp.*' -o -name '*.1500bp.*' -o -name '*.1000bp.bp_covered.*' -o -name '*.1000bp.frac_covered.*' \) -delete
  rm -rf "$COV_DIR"/encode_coremarks/25bp "$COV_DIR"/encode_coremarks/250bp "$COV_DIR"/encode_coremarks/1500bp
  rm -rf "$COV_DIR"/k27me3_top95/25bp "$COV_DIR"/k27me3_top95/250bp "$COV_DIR"/k27me3_top95/1500bp
  rm -rf "$TMP_DIR"/*.25bp "$TMP_DIR"/*.250bp "$TMP_DIR"/*.1500bp
}

remove_requested_tile_size_outputs() {
  local tile_size
  for tile_size in $TILE_SIZES; do
    rm -f "$TILE_DIR"/*.${tile_size}bp.tiles.bed
    rm -f "$MATRIX_DIR"/*.${tile_size}bp.read_count.matrix.tsv.gz
    rm -f "$MATRIX_DIR"/*.${tile_size}bp.rpkm.matrix.tsv.gz
    rm -rf "$COV_DIR"/encode_coremarks/${tile_size}bp
    rm -rf "$COV_DIR"/k27me3_top95/${tile_size}bp
    rm -rf "$TMP_DIR"/*.${tile_size}bp
  done
}

prepare() {
  remove_old_tile_sizes
  remove_requested_tile_size_outputs

  {
    printf "sample\tmark\tphase\treplicate\tn_reads\tbed_path\n"
    find "$BED_DIR" -maxdepth 1 -type f -name '*.bed' | LC_ALL=C sort | while read -r bed; do
      name="$(basename "$bed")"
      if [[ "$name" =~ ^JG_HsDm_K562_(K27me3|K27me2|K27me1|K27ac|P2S25p)_241101_(G1|S1|S2|S3|S4|G2)_d2pcr_([ab])\.bed$ ]]; then
        mark="${BASH_REMATCH[1]}"
        phase="${BASH_REMATCH[2]}"
        rep="${BASH_REMATCH[3]}"
        n_reads="$(wc -l < "$bed")"
        printf "%s_%s_%s\t%s\t%s\t%s\t%s\t%s\n" "$mark" "$phase" "$rep" "$mark" "$phase" "$rep" "$n_reads" "$bed"
      fi
    done
  } > "$sample_manifest"

  write_domains_bed "encode_coremarks" "$PROJECT/data/2024/ENCODE_states/K562_E123_15_coreMarks_domains.bed"
  write_domains_bed "k27me3_top95" "$PROJECT/data/2024/SH_all_data/downsample_peakPRINT_slope1_min300/min350/top95/K27me3_K_top95regions.bed"

  assign_rt "encode_coremarks"
  assign_rt "k27me3_top95"

  for tile_size in $TILE_SIZES; do
    make_tiles \
      "encode_coremarks" \
      "$DOMAIN_META_DIR/encode_coremarks.domains.bed" \
      "$tile_size"
    make_tiles \
      "k27me3_top95" \
      "$DOMAIN_META_DIR/k27me3_top95.domains.bed" \
      "$tile_size"
  done

  {
    printf "domain_set\ttile_size\tsample\tbed_path\n"
    while IFS=$'\t' read -r sample mark phase replicate n_reads bed_path; do
      [[ "$sample" == "sample" ]] && continue
      for domain_set in encode_coremarks k27me3_top95; do
        for tile_size in $TILE_SIZES; do
          printf "%s\t%s\t%s\t%s\n" "$domain_set" "$tile_size" "$sample" "$bed_path"
        done
      done
    done < "$sample_manifest"
  } > "$task_manifest"

  printf "Wrote %s (%d samples)\n" "$sample_manifest" "$(( $(wc -l < "$sample_manifest") - 1 ))"
  printf "Wrote %s (%d tasks)\n" "$task_manifest" "$(( $(wc -l < "$task_manifest") - 1 ))"
}

coverage() {
  local task_id="${1:?coverage mode requires TASK_ID}"
  local line
  line="$(awk -v n="$task_id" 'NR == n + 1 {print; exit}' "$task_manifest")"
  if [[ -z "$line" ]]; then
    echo "No task line for task id $task_id in $task_manifest" >&2
    exit 2
  fi

  local domain_set tile_size sample bed_path
  IFS=$'\t' read -r domain_set tile_size sample bed_path <<< "$line"

  local tile_bed="$TILE_DIR/${domain_set}.${tile_size}bp.tiles.bed"
  local coverage_dir="$COV_DIR/$domain_set/${tile_size}bp"
  local coverage_tsv_gz="$coverage_dir/${sample}.reads.tsv.gz"
  local done_file="$coverage_tsv_gz.done"

  mkdir -p "$coverage_dir"

  if [[ -s "$coverage_tsv_gz" && -e "$done_file" ]]; then
    gzip -t "$coverage_tsv_gz"
    echo "Output exists, skipping: $coverage_tsv_gz"
    return
  fi

  local tmp_cov="$coverage_tsv_gz.tmp.$$"
  bedtools intersect -a "$tile_bed" -b "$bed_path" -c \
    | awk -v OFS='\t' '{print $5, $NF}' \
    | LC_ALL=C sort -k1,1 \
    | gzip -c > "$tmp_cov"

  gzip -t "$tmp_cov"
  mv "$tmp_cov" "$coverage_tsv_gz"
  touch "$done_file"
  echo "Wrote $coverage_tsv_gz"
}

collate_one() {
  local domain_set="$1"
  local tile_size="$2"
  local tile_bed="$TILE_DIR/${domain_set}.${tile_size}bp.tiles.bed"
  local work_dir="$TMP_DIR/${domain_set}.${tile_size}bp"
  local matrix_gz="$MATRIX_DIR/${domain_set}.${tile_size}bp.read_count.matrix.tsv.gz"
  local rpkm_gz="$MATRIX_DIR/${domain_set}.${tile_size}bp.rpkm.matrix.tsv.gz"

  if [[ -s "$matrix_gz" && -s "$rpkm_gz" ]]; then
    gzip -t "$matrix_gz"
    gzip -t "$rpkm_gz"
    echo "Matrices exist, skipping: $domain_set ${tile_size}bp"
    return
  fi

  rm -rf "$work_dir"
  mkdir -p "$work_dir"
  awk -v OFS='\t' '{print $5,$4,$1,$2,$3,$6,$7}' "$tile_bed" | LC_ALL=C sort -k1,1 > "$work_dir/current.tsv"

  mapfile -t samples < <(awk 'NR > 1 {print $1}' "$sample_manifest")
  for sample in "${samples[@]}"; do
    local cov_gz="$COV_DIR/$domain_set/${tile_size}bp/${sample}.reads.tsv.gz"
    if [[ ! -s "$cov_gz" ]]; then
      echo "Missing coverage column: $cov_gz" >&2
      exit 3
    fi
    gzip -t "$cov_gz"
    zcat "$cov_gz" > "$work_dir/${sample}.tsv"
    LC_ALL=C join -t $'\t' -a 1 -e 0 -o auto "$work_dir/current.tsv" "$work_dir/${sample}.tsv" > "$work_dir/next.tsv"
    mv "$work_dir/next.tsv" "$work_dir/current.tsv"
    rm -f "$work_dir/${sample}.tsv"
  done

  {
    printf "domain_id\tchr\tstart\tend\tbin_id\tbin_index\ttile_size"
    for sample in "${samples[@]}"; do
      printf "\t%s" "$sample"
    done
    printf "\n"
    awk -v OFS='\t' '{
      printf "%s\t%s\t%s\t%s\t%s\t%s\t%s", $2,$3,$4,$5,$1,$6,$7;
      for (i=8; i<=NF; i++) printf "\t%s", $i;
      printf "\n";
    }' "$work_dir/current.tsv"
  } | gzip -c > "$matrix_gz.tmp"
  gzip -t "$matrix_gz.tmp"
  mv "$matrix_gz.tmp" "$matrix_gz"

  zcat "$matrix_gz" \
    | awk -v OFS='\t' -v samples_file="$sample_manifest" '
      BEGIN {
        while ((getline line < samples_file) > 0) {
          if (line ~ /^sample\t/) continue;
          n=split(line, a, "\t");
          lib[a[1]]=a[5];
        }
      }
      NR == 1 {
        for (i=8; i<=NF; i++) header[i]=$i;
        print;
        next;
      }
      {
        kb=($4-$3)/1000;
        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s", $1,$2,$3,$4,$5,$6,$7;
        for (i=8; i<=NF; i++) {
          sample_name=header[i];
          v=0;
          if (kb > 0 && lib[sample_name] > 0) v=$i / (kb * lib[sample_name] / 1000000);
          printf "\t%.7g", v;
        }
        printf "\n";
      }
    ' | gzip -c > "$rpkm_gz.tmp"
  gzip -t "$rpkm_gz.tmp"
  mv "$rpkm_gz.tmp" "$rpkm_gz"

  rm -rf "$work_dir"
  echo "Wrote $matrix_gz and $rpkm_gz"
}

collate() {
  for domain_set in encode_coremarks k27me3_top95; do
    for tile_size in $TILE_SIZES; do
      collate_one "$domain_set" "$tile_size"
    done
  done
}

case "${1:-}" in
  prepare)
    prepare
    ;;
  coverage)
    coverage "${2:-}"
    ;;
  collate)
    collate
    ;;
  *)
    usage >&2
    exit 1
    ;;
esac
