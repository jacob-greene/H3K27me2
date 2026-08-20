#!/usr/bin/env Rscript
# edgeR_paired_then_series_hg38_total_RUVr_oneTSV.R
#
# Drug time series:
#   - hg38 only, main + reseq counts summed per gene
#   - lib.size = total hg38 reads from featureCounts summary (main + reseq),
#                with reseq runs collapsed by basename
#   - RUVr (k=1) using ~ Rep + Time*Drug residuals, all genes as controls
#   - Statistical tests (only for genes with DMSO-based max RPKM >= 0.5):
#       * Per-timepoint, paired:   EED vs DMSO, TAZ vs DMSO  (design: ~ Rep + Drug + W)
#       * Series-wide, factorial:  ~ Rep + Time*Drug + W with
#           - OMNIBUS F-test (all time-specific effects)
#           - AVG (mean effect across time)
#           - INTERACTION F-test (Time*Drug terms)
#
# Output: one TSV per mark with:
#   - logFC / PValue / FDR for each per-timepoint contrast
#   - logFC / PValue / FDR for AVG contrasts
#   - F / PValue / FDR for OMNIBUS and INTERACTION
#   - AveExpr (logCPM from series model)
#   - DMSO-based RPKM summary
#   - Per-sample backend CPM (RUVr-normalized)
#   - meanCPM_<Time>_<Drug> from backend CPM
#   - Length and raw counts per sample

suppressPackageStartupMessages({
  library(edgeR)
  library(RUVSeq)
  library(stringr)
  library(dplyr)
  library(readr)
  library(tibble)
})

`%||%` <- function(a,b) if(!is.null(a)) a else b

# ====== 0) Parameters ======
hg38_data_dir <- "data/2024/RepliTag/filtered_sams_drug_1.0/"

# counts (per-gene) – main + reseq
hg38_counts_main  <- file.path(hg38_data_dir, "hg38.gencodev44_RepliTag_genes-1000_featureCounts_counts_plusDUP.tsv")
hg38_counts_reseq <- file.path(hg38_data_dir, "reseq/hg38.gencodev44_RepliTag_genes-1000_featureCounts_counts_plusDUP.tsv")

# total hg38 reads (from featureCounts summaries) – main + reseq
hg38_summary_main  <- file.path(hg38_data_dir, "hg38.gencodev44_RepliTag_genes-1000_featureCounts_counts_plusDUP.txt.summary.tsv")
hg38_summary_reseq <- file.path(hg38_data_dir, "reseq/hg38.gencodev44_RepliTag_genes-1000_featureCounts_counts_plusDUP.txt.summary.tsv")
#hg38_summary_reseq <- "data/2024/RepliTag/featureCounts_drug_dm6/genes_v5_KABK_230607_featureCounts_counts_plusDUP.txt.summary"

output_dir    <- "data/2024/RepliTag/Matrices/TAZEED"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

time_levels <- c("G1","S1","S2","S3","S4","G2")
drug_levels <- c("DMSO","EED","TAZ")
marks_keep  <- c("K27ac","K27me1","K27me2","K27me3","P2S25p")
k_ruv       <- 1   # RUVr dimension

# ====== 1) Helpers ======
strip_hg38_suffix <- function(x) sub("\\.DTmarked_filtered\\.sam$", "", basename(x))

parse_coldata <- function(snames) {
  # JG_HsDm_K562_<MARK>_<DATE>_(G1|S1-4|G2)_(DMSO|EED|TAZ)[_CUTAC]?_d2pcr_(a|b)
  pat <- "^JG_HsDm_K562_([^_]+)_[0-9]+_(G1|S[1-4]|G2)_(DMSO|EED|TAZ)(?:_CUTAC)?_d2pcr_([ab])$"
  m <- stringr::str_match(snames, pat)
  out <- data.frame(
    Sample = snames,
    Mark   = m[,2],
    Time   = factor(m[,3], levels = time_levels),
    Drug   = factor(m[,4], levels = drug_levels),
    Rep    = factor(m[,5]),
    stringsAsFactors = FALSE
  )
  rownames(out) <- out$Sample
  out
}

safe_read_counts <- function(fn) {
  if (is.null(fn) || !file.exists(fn)) return(NULL)
  x <- read.delim(fn, check.names = FALSE, stringsAsFactors = FALSE)
  rownames(x) <- x[[1]]
  x[[1]] <- NULL
  as.data.frame(x, check.names = FALSE)
}

safe_read_summary <- function(fn) {
  if (is.null(fn) || !file.exists(fn)) return(NULL)
  x <- read.delim(fn, check.names = FALSE, stringsAsFactors = FALSE)
  rownames(x) <- x[[1]]
  x[[1]] <- NULL
  as.data.frame(x, check.names = FALSE)
}

strip_meta <- function(x) {
  meta_drop <- intersect(colnames(x), c("Chr","Start","End","Strand","Length"))
  if (length(meta_drop)) x <- x[ , setdiff(colnames(x), meta_drop), drop = FALSE]
  x
}

sum_two_mats <- function(m1, m2) {
  if (is.null(m1) && is.null(m2)) return(NULL)
  if (is.null(m1)) return(m2)
  if (is.null(m2)) return(m1)
  genes <- union(rownames(m1), rownames(m2))
  samps <- union(colnames(m1), colnames(m2))
  out <- matrix(0, nrow = length(genes), ncol = length(samps),
                dimnames = list(genes, samps))
  out[rownames(m1), colnames(m1)] <- out[rownames(m1), colnames(m1)] + as.matrix(m1)
  out[rownames(m2), colnames(m2)] <- out[rownames(m2), colnames(m2)] + as.matrix(m2)
  as.data.frame(out, check.names = FALSE)
}

# collapse duplicate sample columns in summary files *after* stripping dirs
collapse_summary_cols <- function(x) {
  if (is.null(x)) return(NULL)
  cn <- colnames(x)
  if (!any(duplicated(cn))) return(x)
  groups <- split(seq_along(cn), cn)
  collapsed_mat <- sapply(groups, function(idx) {
    rowSums(as.matrix(x[, idx, drop = FALSE]))
  })
  collapsed <- as.data.frame(collapsed_mat, check.names = FALSE)
  rownames(collapsed) <- rownames(x)
  collapsed
}

get_total_from_summary <- function(x) {
  if (is.null(x)) return(NULL)
  # Prefer "Assigned" if present, else "Total", else sum of all rows
  if ("Assigned" %in% rownames(x)) {
    tot <- x["Assigned", , drop = TRUE]
  } else if ("Total" %in% rownames(x)) {
    tot <- x["Total", , drop = TRUE]
  } else {
    tot <- colSums(x)
  }
  tot
}

coef_names_for_time_drug <- function(coln, drug, times, baseTime) {
  # Effect at baseline time = main Drug coefficient
  # At other times = main Drug + interaction(Time_t:Drug)
  out <- lapply(times, function(ti) {
    pieces <- c()
    main <- paste0("Drug", drug)
    if (main %in% coln) pieces <- c(pieces, main)
    if (ti != baseTime) {
      intr     <- paste0("Time", ti, ":Drug", drug)
      intr_alt <- paste0("Drug", drug, ":Time", ti)
      if (intr %in% coln) pieces <- c(pieces, intr)
      else if (intr_alt %in% coln) pieces <- c(pieces, intr_alt)
    }
    pieces
  })
  names(out) <- times
  out
}

avg_contrast_vector <- function(coln, pieces_list) {
  K <- length(pieces_list)
  v <- setNames(numeric(length(coln)), coln)
  for (set in pieces_list) for (p in set) v[p] <- v[p] + 1/K
  matrix(v, ncol = 1)
}

omnibus_contrast_matrix <- function(coln, pieces_list) {
  mats <- lapply(names(pieces_list), function(ti) {
    v <- setNames(numeric(length(coln)), coln)
    for (p in pieces_list[[ti]]) v[p] <- v[p] + 1
    v
  })
  do.call(cbind, mats)
}

interaction_contrast_matrix <- function(coln, drug, times, baseTime) {
  inter_cols <- c()
  for (ti in setdiff(times, baseTime)) {
    inter_cols <- c(inter_cols,
                    paste0("Time", ti, ":Drug", drug),
                    paste0("Drug", drug, ":Time", ti))
  }
  inter_cols <- unique(inter_cols[inter_cols %in% coln])
  if (!length(inter_cols)) return(NULL)
  M <- matrix(0, nrow = length(coln), ncol = length(inter_cols),
              dimnames = list(coln, paste0(drug, "_interaction_", seq_along(inter_cols))))
  for (j in seq_along(inter_cols)) M[inter_cols[j], j] <- 1
  M
}

# ====== 2) Load & aggregate hg38 counts and total reads ======
hg38_main_raw  <- safe_read_counts(hg38_counts_main)
hg38_reseq_raw <- safe_read_counts(hg38_counts_reseq)

if (is.null(hg38_main_raw) && is.null(hg38_reseq_raw)) {
  stop("No hg38 counts loaded from either main or reseq.")
}

# extract gene lengths from main counts file (WT-style)
len_tab <- read.delim(hg38_counts_main, check.names = FALSE, stringsAsFactors = FALSE)
gene_ids_len <- len_tab[[1]]
gene_length  <- len_tab$Length
names(gene_length) <- gene_ids_len

# strip meta from counts
hg38_main  <- if (!is.null(hg38_main_raw))  strip_meta(hg38_main_raw)  else NULL
hg38_reseq <- if (!is.null(hg38_reseq_raw)) strip_meta(hg38_reseq_raw) else NULL

if (!is.null(hg38_main))  colnames(hg38_main)  <- strip_hg38_suffix(colnames(hg38_main))
if (!is.null(hg38_reseq)) colnames(hg38_reseq) <- strip_hg38_suffix(colnames(hg38_reseq))

hg38 <- sum_two_mats(hg38_main, hg38_reseq)
if (is.null(hg38)) stop("No hg38 counts after summing main and reseq.")

# total hg38 reads from summaries (main + reseq), with duplicate basenames collapsed
sum_main  <- safe_read_summary(hg38_summary_main)
sum_reseq <- safe_read_summary(hg38_summary_reseq)

if (!is.null(sum_main)) {
  colnames(sum_main) <- strip_hg38_suffix(colnames(sum_main))
  sum_main <- collapse_summary_cols(sum_main)
}
if (!is.null(sum_reseq)) {
  colnames(sum_reseq) <- strip_hg38_suffix(colnames(sum_reseq))
  sum_reseq <- collapse_summary_cols(sum_reseq)
}

tot_main  <- get_total_from_summary(sum_main)
tot_reseq <- get_total_from_summary(sum_reseq)

if (!is.null(tot_main))  names(tot_main)  <- colnames(sum_main)
if (!is.null(tot_reseq)) names(tot_reseq) <- colnames(sum_reseq)

samples_all <- colnames(hg38)
total_reads <- setNames(rep(0, length(samples_all)), samples_all)

if (!is.null(tot_main)) {
  for (s in intersect(samples_all, names(tot_main))) {
    total_reads[s] <- total_reads[s] + as.numeric(tot_main[s])
  }
}
if (!is.null(tot_reseq)) {
  for (s in intersect(samples_all, names(tot_reseq))) {
    total_reads[s] <- total_reads[s] + as.numeric(tot_reseq[s])
  }
}

# fallback: if any samples missing totals, use column sum of counts
zero <- is.na(total_reads) | total_reads <= 0
if (any(zero)) {
  warn_samps <- names(total_reads)[zero]
  warning("Using column sums as total hg38 reads for samples: ",
          paste(warn_samps, collapse = ", "))
  total_reads[zero] <- colSums(hg38)[warn_samps]
}

# align gene_length to hg38
gene_length <- gene_length[rownames(hg38)]

# ====== 3) Parse coldata and filter samples ======
coldata_all <- parse_coldata(colnames(hg38))
bad <- which(is.na(coldata_all$Mark) | is.na(coldata_all$Time) | is.na(coldata_all$Drug))
if (length(bad)) {
  warning("Dropping unparsable samples: ", paste(coldata_all$Sample[bad], collapse=", "))
  keep <- setdiff(seq_len(nrow(coldata_all)), bad)
  good_samples <- coldata_all$Sample[keep]
  coldata_all  <- coldata_all[keep, , drop = FALSE]
  hg38         <- hg38[, good_samples, drop = FALSE]
  total_reads  <- total_reads[good_samples]
}

marks <- intersect(sort(unique(coldata_all$Mark)), marks_keep)
message("Marks detected: ", paste(marks, collapse = ", "))

# ====== 4) Per-mark loop (with RUVr, one TSV per mark) ======
for (mark in marks) {
  message("=== Processing mark: ", mark, " ===")
  cds <- dplyr::filter(coldata_all, Mark == mark)
  if (nrow(cds) < 6) {
    warning("Too few samples for mark ", mark, "; skipping.")
    next
  }
  cols <- cds$Sample
  hg38_m <- hg38[, cols, drop = FALSE]
  
  hg38_m[is.na(hg38_m)] <- 0
  hg38_m <- pmax(round(hg38_m), 0)
  
  rn_all <- rownames(hg38_m)
  
  # library sizes for this mark (from summaries, with reseq collapsed)
  lib_sizes_mark <- total_reads[cols]
  if (any(is.na(lib_sizes_mark) | lib_sizes_mark <= 0)) {
    bad_lib <- names(lib_sizes_mark)[is.na(lib_sizes_mark) | lib_sizes_mark <= 0]
    warning("Library sizes <=0 or NA for ", mark, " samples: ",
            paste(bad_lib, collapse = ", "),
            ". Using column sums instead.")
    lib_sizes_mark[bad_lib] <- colSums(hg38_m)[bad_lib]
  }
  
  # ===============================
  # RPKM summary based on DMSO only
  # ===============================
  gene_len_sub <- as.numeric(gene_length[rn_all])
  if (any(is.na(gene_len_sub))) {
    warning("NA gene lengths for some genes in mark ", mark,
            "; leaving them as NA. RPKM will be NA for those genes.")
  }
  
  dge_rpkm <- DGEList(counts = hg38_m, lib.size = lib_sizes_mark)
  dge_rpkm <- calcNormFactors(dge_rpkm, method = "RLE")
  
  rpkm_all <- rpkm(
    dge_rpkm,
    gene.length          = gene_len_sub,
    normalized.lib.sizes = TRUE
  )
  
  # mean DMSO RPKM per timepoint
  cond_means <- sapply(time_levels, function(tp) {
    idx <- cds$Time == tp & cds$Drug == "DMSO"
    if (!any(idx)) {
      rep(NA_real_, nrow(rpkm_all))
    } else {
      rowMeans(rpkm_all[, cds$Sample[idx], drop = FALSE])
    }
  })
  colnames(cond_means) <- time_levels
  rownames(cond_means) <- rn_all
  
  max_mean_rpkm <- apply(cond_means, 1, function(z) {
    if (all(is.na(z))) NA_real_ else max(z, na.rm = TRUE)
  })
  min_mean_rpkm <- apply(cond_means, 1, function(z) {
    if (all(is.na(z))) NA_real_ else min(z, na.rm = TRUE)
  })
  
  high_expr <- !is.na(max_mean_rpkm) & (max_mean_rpkm >= 0.5)
  
  if (!any(high_expr)) {
    warning("No genes meet high_expr RPKM threshold for mark ", mark,
            "; statistical tests will be skipped.")
  }
  
  # ===============================
  # DGEList for RUVr and modeling
  # ===============================
  y_all <- DGEList(counts = hg38_m, lib.size = lib_sizes_mark)
  y_all$samples$norm.factors <- 1
  
  # -----------------------------
  # RUVr
  # -----------------------------
  design0_df <- cds
  design0_df$Time <- factor(as.character(design0_df$Time), levels = time_levels)
  design0_df$Drug <- factor(as.character(design0_df$Drug), levels = drug_levels)
  design0 <- model.matrix(~ Rep + Time * Drug, data = design0_df)
  
  y0   <- estimateDisp(y_all, design0)
  fit0 <- glmQLFit(y0, design0)
  
  resids <- residuals(fit0)
  ctrls  <- rn_all
  ruv_r  <- RUVr(as.matrix(hg38_m), ctrls, k = k_ruv, res = resids)
  
  W_ruv <- ruv_r$W
  if (is.null(colnames(W_ruv))) {
    colnames(W_ruv) <- paste0("W", seq_len(ncol(W_ruv)))
  }
  if (is.null(rownames(W_ruv))) {
    rownames(W_ruv) <- cols
  }
  W_ruv <- W_ruv[cols, , drop = FALSE]
  
  # ==========================================================
  # STEP 1 — Per-timepoint, PAIRED (design: ~ Rep + Drug + W)
  # ==========================================================
  res_tp_EED <- setNames(vector("list", length(time_levels)), time_levels)
  res_tp_TAZ <- setNames(vector("list", length(time_levels)), time_levels)
  names(res_tp_EED) <- time_levels
  names(res_tp_TAZ) <- time_levels
  
  if (any(high_expr)) {
    for (tp in time_levels) {
      cds_tp <- dplyr::filter(cds, Time == tp)
      if (nrow(cds_tp) < 4) next  # need at least 2 drugs × 2 reps
      cols_tp <- cds_tp$Sample
      
      y_tp_full <- y_all[, cols_tp, keep.lib.sizes = TRUE]
      y_tp <- y_tp_full[high_expr, , keep.lib.sizes = TRUE]
      if (nrow(y_tp) == 0) next
      
      W_tp <- W_ruv[cols_tp, , drop = FALSE]
      
      cds_tp$Drug <- factor(as.character(cds_tp$Drug), levels = drug_levels)
      
      design_df <- cds_tp
      if (ncol(W_tp)) design_df <- cbind(design_df, W_tp)
      
      rhs <- paste0("Rep + Drug",
                    if (ncol(W_tp)) paste0(" + ", paste(colnames(W_tp), collapse = " + ")) else "")
      design <- model.matrix(as.formula(paste("~", rhs)), data = design_df)
      
      y_tp <- estimateDisp(y_tp, design)
      fit_tp <- glmQLFit(y_tp, design)
      
      if ("DrugEED" %in% colnames(design)) {
        tt <- topTags(glmQLFTest(fit_tp, coef = "DrugEED"), n = Inf, sort.by = "none")$table
        res_tp_EED[[tp]] <- tt
      }
      if ("DrugTAZ" %in% colnames(design)) {
        tt <- topTags(glmQLFTest(fit_tp, coef = "DrugTAZ"), n = Inf, sort.by = "none")$table
        res_tp_TAZ[[tp]] <- tt
      }
    }
  }
  
  # ==========================================================
  # STEP 2 — Series-wide model: ~ Rep + Time*Drug (+ W)
  # ==========================================================
  res_EED_omni <- res_TAZ_omni <- NULL
  res_EED_avg  <- res_TAZ_avg  <- NULL
  res_EED_int  <- res_TAZ_int  <- NULL
  
  if (any(high_expr)) {
    design_df <- cds
    design_df$Time <- factor(as.character(design_df$Time), levels = time_levels)
    design_df$Drug <- factor(as.character(design_df$Drug), levels = drug_levels)
    
    W_full <- W_ruv[design_df$Sample, , drop = FALSE]
    if (ncol(W_full)) design_df <- cbind(design_df, W_full)
    
    rhs_series <- paste0("Rep + Time * Drug",
                         if (ncol(W_full)) paste0(" + ", paste(colnames(W_full), collapse = " + ")) else "")
    design_series <- model.matrix(as.formula(paste("~", rhs_series)), data = design_df)
    
    y_series_full <- y_all[, design_df$Sample, keep.lib.sizes = TRUE]
    y_series <- y_series_full[high_expr, , keep.lib.sizes = TRUE]
    
    if (nrow(y_series) > 0) {
      y_series <- estimateDisp(y_series, design_series)
      fit <- glmQLFit(y_series, design_series)
      
      coln <- colnames(design_series)
      baseTime <- time_levels[1]
      
      pieces_EED <- coef_names_for_time_drug(coln, "EED", time_levels, baseTime)
      pieces_TAZ <- coef_names_for_time_drug(coln, "TAZ", time_levels, baseTime)
      
      C_EED_omni <- omnibus_contrast_matrix(coln, pieces_EED)
      C_TAZ_omni <- omnibus_contrast_matrix(coln, pieces_TAZ)
      
      if (!is.null(C_EED_omni)) {
        res_EED_omni <- topTags(glmQLFTest(fit, contrast = C_EED_omni),
                                n = Inf, sort.by = "none")$table
      }
      if (!is.null(C_TAZ_omni)) {
        res_TAZ_omni <- topTags(glmQLFTest(fit, contrast = C_TAZ_omni),
                                n = Inf, sort.by = "none")$table
      }
      
      C_EED_avg <- avg_contrast_vector(coln, pieces_EED)
      C_TAZ_avg <- avg_contrast_vector(coln, pieces_TAZ)
      
      res_EED_avg <- topTags(glmQLFTest(fit, contrast = C_EED_avg),
                             n = Inf, sort.by = "none")$table
      res_TAZ_avg <- topTags(glmQLFTest(fit, contrast = C_TAZ_avg),
                             n = Inf, sort.by = "none")$table
      
      C_EED_int <- interaction_contrast_matrix(coln, "EED", time_levels, baseTime)
      C_TAZ_int <- interaction_contrast_matrix(coln, "TAZ", time_levels, baseTime)
      
      if (!is.null(C_EED_int)) {
        res_EED_int <- topTags(glmQLFTest(fit, contrast = C_EED_int),
                               n = Inf, sort.by = "none")$table
      }
      if (!is.null(C_TAZ_int)) {
        res_TAZ_int <- topTags(glmQLFTest(fit, contrast = C_TAZ_int),
                               n = Inf, sort.by = "none")$table
      }
    }
  }
  
  # ==============================
  # Assemble wide output per mark
  # ==============================
  df_out <- data.frame(row.names = rn_all)
  
  for (drug in c("EED","TAZ")) {
    for (tp in time_levels) {
      base <- paste0(drug, "_vs_DMSO_", tp)
      df_out[[paste0("logFC_", base)]]  <- NA_real_
      df_out[[paste0("PValue_", base)]] <- NA_real_
      df_out[[paste0("FDR_", base)]]    <- NA_real_
    }
    df_out[[paste0("logFC_",  drug, "_AVG")]]   <- NA_real_
    df_out[[paste0("PValue_", drug, "_AVG")]]   <- NA_real_
    df_out[[paste0("FDR_",    drug, "_AVG")]]   <- NA_real_
    df_out[[paste0("F_",      drug, "_OMNIBUS")]]      <- NA_real_
    df_out[[paste0("PValue_", drug, "_OMNIBUS")]]      <- NA_real_
    df_out[[paste0("FDR_",    drug, "_OMNIBUS")]]      <- NA_real_
    df_out[[paste0("F_",      drug, "_INTERACTION")]]  <- NA_real_
    df_out[[paste0("PValue_", drug, "_INTERACTION")]]  <- NA_real_
    df_out[[paste0("FDR_",    drug, "_INTERACTION")]]  <- NA_real_
  }
  df_out[["AveExpr"]] <- NA_real_
  
  # per-timepoint tests
  for (tp in time_levels) {
    tabE <- res_tp_EED[[tp]]
    if (!is.null(tabE)) {
      rn <- rownames(tabE)
      base <- paste0("EED_vs_DMSO_", tp)
      df_out[rn, paste0("logFC_",  base)] <- tabE$logFC
      df_out[rn, paste0("PValue_", base)] <- tabE$PValue
      df_out[rn, paste0("FDR_",    base)] <- tabE$FDR
      df_out[rn, "AveExpr"] <- tabE$logCPM
    }
    tabT <- res_tp_TAZ[[tp]]
    if (!is.null(tabT)) {
      rn <- rownames(tabT)
      base <- paste0("TAZ_vs_DMSO_", tp)
      df_out[rn, paste0("logFC_",  base)] <- tabT$logFC
      df_out[rn, paste0("PValue_", base)] <- tabT$PValue
      df_out[rn, paste0("FDR_",    base)] <- tabT$FDR
      df_out[rn, "AveExpr"] <- tabT$logCPM
    }
  }
  
  # AVG tests
  if (!is.null(res_EED_avg)) {
    rn <- rownames(res_EED_avg)
    df_out[rn, "logFC_EED_AVG"]   <- res_EED_avg$logFC
    df_out[rn, "PValue_EED_AVG"]  <- res_EED_avg$PValue
    df_out[rn, "FDR_EED_AVG"]     <- res_EED_avg$FDR
    df_out[rn, "AveExpr"]         <- res_EED_avg$logCPM
  }
  if (!is.null(res_TAZ_avg)) {
    rn <- rownames(res_TAZ_avg)
    df_out[rn, "logFC_TAZ_AVG"]   <- res_TAZ_avg$logFC
    df_out[rn, "PValue_TAZ_AVG"]  <- res_TAZ_avg$PValue
    df_out[rn, "FDR_TAZ_AVG"]     <- res_TAZ_avg$FDR
    df_out[rn, "AveExpr"]         <- res_TAZ_avg$logCPM
  }
  
  # OMNIBUS
  if (!is.null(res_EED_omni)) {
    rn <- rownames(res_EED_omni)
    df_out[rn, "F_EED_OMNIBUS"]      <- res_EED_omni$F
    df_out[rn, "PValue_EED_OMNIBUS"] <- res_EED_omni$PValue
    df_out[rn, "FDR_EED_OMNIBUS"]    <- res_EED_omni$FDR
    df_out[rn, "AveExpr"]            <- res_EED_omni$logCPM
  }
  if (!is.null(res_TAZ_omni)) {
    rn <- rownames(res_TAZ_omni)
    df_out[rn, "F_TAZ_OMNIBUS"]      <- res_TAZ_omni$F
    df_out[rn, "PValue_TAZ_OMNIBUS"] <- res_TAZ_omni$PValue
    df_out[rn, "FDR_TAZ_OMNIBUS"]    <- res_TAZ_omni$FDR
    df_out[rn, "AveExpr"]            <- res_TAZ_omni$logCPM
  }
  
  # INTERACTION
  if (!is.null(res_EED_int)) {
    rn <- rownames(res_EED_int)
    df_out[rn, "F_EED_INTERACTION"]      <- res_EED_int$F
    df_out[rn, "PValue_EED_INTERACTION"] <- res_EED_int$PValue
    df_out[rn, "FDR_EED_INTERACTION"]    <- res_EED_int$FDR
    df_out[rn, "AveExpr"]                <- res_EED_int$logCPM
  }
  if (!is.null(res_TAZ_int)) {
    rn <- rownames(res_TAZ_int)
    df_out[rn, "F_TAZ_INTERACTION"]      <- res_TAZ_int$F
    df_out[rn, "PValue_TAZ_INTERACTION"] <- res_TAZ_int$PValue
    df_out[rn, "FDR_TAZ_INTERACTION"]    <- res_TAZ_int$FDR
    df_out[rn, "AveExpr"]                <- res_TAZ_int$logCPM
  }
  
  # ==========================
  # RPKM summaries for all genes (DMSO-based)
  # ==========================
  df_out[["max_RPKM"]]      <- max_mean_rpkm[rownames(df_out)]
  df_out[["min_RPKM"]]      <- min_mean_rpkm[rownames(df_out)]
  df_out[["max_RPKM_DMSO"]] <- max_mean_rpkm[rownames(df_out)]
  df_out[["min_RPKM_DMSO"]] <- min_mean_rpkm[rownames(df_out)]
  
  # per-timepoint DMSO mean RPKM
  for (tp in time_levels) {
    col_name <- paste0("meanRPKM_", tp, "_DMSO")
    df_out[[col_name]] <- cond_means[rownames(df_out), tp]
  }
  
  # ==========================
  # Backend CPM from RUVr
  # ==========================
  backend_cpm <- cpm(
    ruv_r$normalizedCounts,
    log         = FALSE,
    prior.count = 0
  )
  backend_cpm <- backend_cpm[rownames(df_out), cols, drop = FALSE]
  
  for (s in colnames(backend_cpm)) {
    df_out[[paste0(s, "_CPM")]] <- backend_cpm[rownames(df_out), s]
  }
  
  group_labels <- with(cds, paste0(Time, "_", Drug))
  unique_groups <- unique(group_labels)
  for (g in unique_groups) {
    idx <- group_labels == g
    samp_ids <- cds$Sample[idx]
    mat <- backend_cpm[, samp_ids, drop = FALSE]
    mean_vec <- if (ncol(mat) == 1) mat[,1] else rowMeans(mat)
    df_out[[paste0("meanCPM_", g)]] <- mean_vec[rownames(df_out)]
  }
  
  # Length and raw counts
  df_out[["Length"]] <- gene_length[rownames(df_out)]
  for (s in cols) {
    df_out[[s]] <- hg38_m[rownames(df_out), s]
  }
  
  df_out <- cbind(Geneid = rownames(df_out), df_out)
  
  out_file <- file.path(output_dir, paste0(mark, "_drug_RUVr_series.tsv"))
  write_tsv(df_out, out_file)
  message("Wrote: ", out_file)
}

message("All done. Results in: ", normalizePath(output_dir))
