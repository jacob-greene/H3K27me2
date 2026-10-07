# Processing scripts

These scripts turn the raw CUT&Tag reads (GEO
[GSE327802](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE327802)) into the
processed tables the notebooks read. [`run.sh`](run.sh) is the single entry point. Run it
from the repository root, inside the `h3k27me2` environment:

```bash
conda activate h3k27me2

bash processing_scripts/run.sh list        # the stages, in dependency order
bash processing_scripts/run.sh selftest    # ~1 min check on synthetic data; needs no downloads
bash processing_scripts/run.sh fetch       # Zenodo processed data + public reference files
bash processing_scripts/run.sh <stage>     # DRY_RUN=1 prints the commands without running them
```

Every path in every script is relative to the repository root, and every output lands
under `data/`. `run.sh` refuses to run from anywhere else.

## Stages

| Stage | Script(s) | What it does |
|---|---|---|
| `selftest` | `downsample_peakPRINT.sh` and helpers | Calls peaks on synthetic fragments, to check a fresh clone is wired up |
| `fetch` | `fetch_external_data.sh` | Downloads the Zenodo processed data and the public references (GENCODE, MSigDB, ENCODE, UW Repli-seq, GEO GSE164804) into `data/`, ~2.8 GB; `fetch_external_data.sh list` names each target |
| `align` | `geo_to_sams.sh` | GEO/SRA reads → trimmed, aligned, duplicate-marked SAMs (see the caveat below) |
| `annotation` | `gtf2bed_TSS.sh` | GENCODE v44 → TSS−1 kb..TES gene BED and SAF |
| `filter` | `filter_sams*.sh` | Removes reads in the *Drosophila* spike-in blacklist; counts reads per gene |
| `merge` | `merge_sams2bam_slurm.sh` | Merges replicate SAMs into BAMs (Slurm) |
| `tracks` | `bam2bed2bw*.sh`, `fraction_norm.csh` | Filtered SAMs → normalised bigWig tracks (Slurm) |
| `counts` | `Feature_counts_SAMs*.sh` | featureCounts over gene and genomic-feature SAFs |
| `dge` | `251113_DGE_RUVseq.R`, `251202_DGE_RUVseq_drug.R` | RUVSeq/edgeR time series, wild type and PRC2 inhibitors |
| `peakprint` | `downsample_peakPRINT.sh` | Downsampled peak calling for one fragment BED |
| `heatmaps` | `Peak_heatmaps_v2.sh` | deepTools matrices and heatmaps over peaks |
| `rifs` | `RIFs_RepliTag.sh` | Reads-in-features matrices (the `coverage_scores_*.tab` tables) |
| `repliseq` | `crossmaphg19tohg38.sh` | Lifts the UW Repli-seq bigWigs from hg19 to hg38 |
| `nucleation` | `run_replitag_nucleation_matrices_slurm.sh`, then `build_per_bin_covcorr.py` | Nucleation-site matrices (Slurm array) and the per-bin table `260801_nucleation_pub` reads |
| `blacklist` | `blacklist_dm6.sh` | Rebuilds the spike-in blacklist from four lab-internal BEDs; its output already ships, so skip it |

Running the stages in the order `list` prints satisfies every dependency between them. In
particular `annotation` writes the SAF that `filter` and `counts` read, and `repliseq`
writes the lifted bigWig that `nucleation` reads.

## Caveats

- **`align` has never been run by the authors.** The authors aligned reads with a
  lab-internal pipeline. `geo_to_sams.sh` reimplements it from the manuscript Methods
  (cutadapt 4.4, bowtie2 2.5.1, samtools duplicate marking) so that the GEO reads can be
  processed outside the lab.
- **`nucleation` returns before its jobs finish.** It submits a Slurm array and a dependent
  collation job. Run `python3 processing_scripts/build_per_bin_covcorr.py` once they finish.
- **`merge`, `tracks` and `nucleation` submit Slurm jobs.** The other stages run in the
  foreground.
- On an Lmod cluster the scripts load a module only when the tool is not already on `PATH`,
  so the activated conda environment always takes precedence.
