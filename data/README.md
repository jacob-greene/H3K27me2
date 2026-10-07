# data/

Inputs for the notebooks and the processing scripts. A few small files ship with the
repository. `bash processing_scripts/run.sh fetch` downloads the public ones, and the
processing stages write the rest. Git ignores everything here except the files below.

## Files that ship with the repository

| File | Size | Provenance |
|---|---|---|
| `250214_K27me123ac_P2S25p_TAZEEDUTX_cell_counts.xlsx` | 13 KB | this study: cell counts and viability |
| `20260206_drug_8only_aurora/*.csv` (15 files) | 6.3 MB | this study: flow-cytometry exports |
| `Kc_merged_blacklist_hg38.bed` | 7.9 MB | this study: *Drosophila* spike-in blacklist, 286,386 intervals |
| `2024/LAD_data/JG_coding_genes_maxRPKM.tsv` | 630 KB | this study: per-gene maximum RPKM and replication timing |
| `chr_lens.txt` | 12 KB | UCSC hg38 chromosome sizes |
| `2024/LAD_data/Dataset_S1_Promoter_SuRE_Classification.tsv` | 3.2 MB | **published supplement, not this study** |
| `2024/LAD_data/Dataset_S2_TRIP.tsv` | 1.1 MB | **published supplement, not this study** |

`Dataset_S1` and `Dataset_S2` are from Leemans et al., *Cell* 177:852–864 (2019). They are
redistributed unmodified under **CC BY-NC-ND 4.0**, so that
`260803_WT_timeseries_clean_chromHMM_pub.ipynb` runs without a manual download. **Cite that
paper, not this one, for those two files.** See
[`THIRD_PARTY_NOTICES.md`](../THIRD_PARTY_NOTICES.md) for the terms.

## Full layout after `fetch` and the processing stages

Each note names the producer: `fetch` (downloaded by `run.sh fetch`), a `run.sh` stage, or
a notebook. Entries marked `no script` were built in the authors' working copy, and no
script in this repository writes them.

```
data/
├── 2024/
│   ├── cCREs/intersected_annos/FC_out/*_featureCounts_bulk_*.tsv   (counts)
│   ├── ENCODE_states/K562_E123_15_coreMarks_domains.bed   (K27me3_states_111epigenomes_pub)
│   ├── K562_annotations/
│   │   ├── gencode.v44.basic.annotation.coding.gtf          (fetch)
│   │   └── processed_beds/sorted/gencode.v44.basic.TSS-1000_TES_sorted.{bed,saf}   (annotation)
│   ├── LAD_data/
│   │   ├── gencode.v27.annotation.gtf.gz                    (fetch)
│   │   └── *.saf, for the counts stage                      (no script)
│   ├── MTF2_GSE164804/
│   │   ├── ENCODE_{EZH2,SUZ12}_K562_*_fc_hg38.bw            (fetch)
│   │   └── MTF2_shCT_hg38_coverage.bw                       (no script: hg38 liftover of the fetched GSE164804 bedGraph)
│   ├── RepliSeq/*.bigWig                                    (fetch: UW Repli-seq, hg19)
│   ├── RepliTag/
│   │   ├── bws_bulk/*.bw                                    (no script)
│   │   ├── filtered_sams_WT_1.0/, filtered_sams_drug_1.0/   (filter)
│   │   ├── filtered_sams_WT_1.0/251113_DGE/*/output.tsv     (fetch from Zenodo, or dge)
│   │   ├── Matrices/coverage_scores_*.tab                   (rifs)
│   │   ├── Matrices/TAZEED/*_drug_RUVr_series.tsv           (fetch from Zenodo, or dge)
│   │   ├── refs/GOBP_CELL_CYCLE.v2025.1.Hs.json             (fetch)
│   │   └── RepliSeq/hg38_bws_crossmap/                      (repliseq)
│   └── SH_all_data/
│       ├── sections4/, merged_bams4/                        (merge)
│       ├── merged_bws4/                                     (no script)
│       ├── peakPRINT_50/MERGED_downsampled_bins.tsv         (no script)
│       └── downsample_peakPRINT_slope1_min300/
│           └── min350/{*.bed, *_meta_v*.tsv}                (peakprint; the meta tables: no script)
│               └── top95/*_top95regions.bed                 (fetch from Zenodo; see the note below)
└── results/
    ├── mtf2_nucleation/per_bin_covcorr_full.tsv.gz          (nucleation)
    ├── mtf2_nucleation/nuc_cpgonly_spread_bins.tsv.gz       (no script)
    └── replitag_sphase_nucleation/matrices/*.matrix.tsv.gz  (nucleation)
```

### Note on `top95/*_top95regions.bed`

`250411_peakPrint_clean_pub.ipynb` writes all fifteen `{mark}_{cell}_top95regions.bed`
files, but as published it reads only `K27me3_K_peak_coverage.bed`. It therefore
regenerates `K27me3_K_top95regions.bed` (identical to the Zenodo copy) and writes the other
fourteen files empty. `Peak_heatmaps_v2.sh` needs the five `*_K_*` files, so after running
that notebook, run `bash processing_scripts/run.sh fetch zenodo` again: it re-downloads any
empty file and checks every checksum.
