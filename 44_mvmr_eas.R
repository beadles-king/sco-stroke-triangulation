# =============================================================================
# Reproducibility repository - Wang et al., transportability-aware triangulation
# of metabolic stroke evidence across the Shanghai Cooperation Organisation.
# Paths are resolved by config.R; see README.md.
# =============================================================================
if (!exists("RAWDAT")) {
  if (file.exists("config.R")) source("config.R") else
    stop("Run this script from the code/ directory (the folder containing config.R).")
}

# =============================================================================
# Multivariable Mendelian Randomization: East Asian Ancestry Cohort  (v2)
# Exposures (k = 5): FG, LDL, SBP, BMI (TPMI) + eGFR (BBJ proxy)
# Outcomes  (3-4)  : Stroke (any), IS, ICH, SAH (GIGASTROKE EAS)
# Methods :
#   - Primary  : MVMR-IVW
#   - Sensitivity : MVMR-Egger (intercept tests pleiotropy)
#   - Diagnostics : Conditional F per exposure; conditional Q (Sanderson-Spiller)
# Weak-instrument threshold : cF >= 5 (Sanderson 2021 Stat Med)
# Multiple testing : BH-FDR within each outcome on IVW p-values
# Clumping :
#   - Within each exposure : r2 < 0.001, kb = 10000 (EAS panel)
#   - Across pooled exposures (2nd clump): r2 < 0.005 (looser, removes only
#     between-exposure LD; matches the original EAS analysis)
#
# Data sources & quirks (vs. EUR):
#   - All exposures and outcomes are LOCAL files (no OpenGWAS download).
#   - TPMI exposures (FG/LDL/SBP/BMI) lack rsIDs natively; you have already
#     pre-mapped chr:pos:A1/A2 via SNPnexus into *_snp_AS.txt files.
#   - BBJ eGFR already has rsIDs and uses different column names.
#   - Outcomes are GWAS Catalog harmonised TSVs (no rsID needed; matched by
#     SNP after exposure rsIDs are resolved).
#
# Bug-prevention measures:
#   1. Each exposure has a `name` field; never index columns by [,1]..[,5].
#   2. After mv_harmonise_data, columns of exposure_beta are reordered to our
#      canonical name list, with a SNP-level numeric cross-check vs. the
#      pre-harmonisation expo_dat_AS.
#   3. id.outcome is set to a clear, unambiguous label per outcome.
# =============================================================================

# ------------------------------------------------------------------
# 0. Setup
# ------------------------------------------------------------------
library(TwoSampleMR)
library(MVMR)
library(MendelianRandomization)
library(data.table)
library(dplyr)
library(ggplot2)
library(tidyr)

# OpenGWAS authentication: set your own token before running, e.g.
#   Sys.setenv(OPENGWAS_JWT = "<your token>")
# Tokens are not distributed with this repository.

work_dir <- file.path(GWAS_DIR, "AS")
setwd(work_dir)

out_dir  <- file.path(RAWDAT, "MR/EAS/MVMR_res")
plot_dir <- file.path(out_dir, "plots")
proc_dir <- file.path(out_dir, "processed")
invisible(lapply(c(out_dir, plot_dir, proc_dir),
                 dir.create, showWarnings = FALSE, recursive = TRUE))

# --- Constants ------------------------------------------------------
P_THRESH        <- 5e-8
R2_WITHIN       <- 0.001    # clumping within each exposure
R2_BETWEEN      <- 0.005    # 2nd clump on pooled exposures (looser)
KB_THRESH       <- 10000
F_THRESH        <- 5
GENCOV          <- 0
POP             <- "EAS"

# Path to the package-cached working folder for clump_data writes
PKG_TMP <- file.path(find.package("TwoSampleMR"), "1.gbd_str")
dir.create(PKG_TMP, showWarnings = FALSE, recursive = TRUE)

# --- Exposure metadata ----------------------------------------------
# Each entry encodes the source format (TPMI vs BBJ) and column mappings.
# `snpnexus_file` is the SNPnexus output for that exposure (already prepared
# offline). Set to NA for sources that already carry rsIDs (e.g. BBJ eGFR).
exposure_meta <- list(
  list(name = "FG",
       source = "TPMI",
       gwas_file = "FG_AS.txt",
       snpnexus_file = "FG_snp_AS.txt",
       cols = list(beta = "BETA", se = "SE", pval = "p.value",
                   ea = "Allele2", oa = "Allele1", eaf = "AF_Allele2",
                   chr = "CHR", pos = "POS"),
       qc_filter = quote(QC == "PASS")),
  list(name = "LDL",
       source = "TPMI",
       gwas_file = "LDL_AS.txt",
       snpnexus_file = "LDL_snp_AS.txt",
       cols = list(beta = "BETA", se = "SE", pval = "p.value",
                   ea = "Allele2", oa = "Allele1", eaf = "AF_Allele2",
                   chr = "CHR", pos = "POS"),
       qc_filter = quote(QC == "PASS")),
  list(name = "SBP",
       source = "TPMI",
       gwas_file = "SBP_AS.txt",
       snpnexus_file = "SBP_snp_AS.txt",
       cols = list(beta = "BETA", se = "SE", pval = "p.value",
                   ea = "Allele2", oa = "Allele1", eaf = "AF_Allele2",
                   chr = "CHR", pos = "POS"),
       qc_filter = quote(QC == "PASS")),
  list(name = "BMI",
       source = "TPMI",
       gwas_file = "BMI_AS.txt",
       snpnexus_file = "BMI_snp_AS.txt",
       cols = list(beta = "BETA", se = "SE", pval = "p.value",
                   ea = "Allele2", oa = "Allele1", eaf = "AF_Allele2",
                   chr = "CHR", pos = "POS"),
       qc_filter = quote(QC == "PASS")),
  list(name = "eGFR",
       source = "BBJ",
       gwas_file = "eGFR_AS.txt.gz",
       snpnexus_file = NA,                    # BBJ already has rsID
       cols = list(beta = "BETA", se = "SE", pval = "P_BOLT_LMM_INF",
                   ea = "ALLELE1", oa = "ALLELE0", eaf = "A1FREQ",
                   snp = "SNP"),
       qc_filter = NULL)
)
exposure_names <- sapply(exposure_meta, `[[`, "name")

# --- Outcome metadata -----------------------------------------------
# All outcomes are GIGASTROKE EAS GWAS Catalog harmonised TSVs.
# Columns are uniform across these files.
OUT_COLS <- list(
  snp  = "variant_id",        beta = "beta",       se   = "standard_error",
  ea   = "effect_allele",     oa   = "other_allele",
  eaf  = "effect_allele_frequency", pval = "p_value"
)

outcome_meta <- list(
  list(name = "AllStroke", file = "Stroke.tsv", label = "Any stroke (GIGASTROKE-EAS)"),
  list(name = "IS",        file = "isc.tsv",    label = "Ischemic stroke (GIGASTROKE-EAS)"),
  list(name = "ICH",       file = "ich.tsv",    label = "Intracerebral hemorrhage (GIGASTROKE-EAS)"),
  list(name = "SAH",       file = "sah.tsv",    label = "Subarachnoid hemorrhage (GIGASTROKE-EAS)")
)

# ------------------------------------------------------------------
# 1. Helper: load TPMI-format exposure, map rsID via SNPnexus, clump
# ------------------------------------------------------------------
prepare_exposure_TPMI <- function(meta) {
  message(sprintf("\n[Exposure-TPMI] %s  (file = %s)", meta$name, meta$gwas_file))
  cache <- file.path(proc_dir, sprintf("%s_clumped.rds", tolower(meta$name)))
  if (file.exists(cache)) {
    message("  Loading cached clumped data."); return(readRDS(cache))
  }
  cc <- meta$cols

  # Read full summary stats
  raw <- fread(meta$gwas_file, sep = "auto", header = TRUE,
               nThread = getDTthreads(), na.strings = c("", "NA"),
               showProgress = TRUE)

  # Filter by p-value + QC
  sub <- raw[get(cc$pval) < P_THRESH]
  if (!is.null(meta$qc_filter)) sub <- sub[eval(meta$qc_filter, sub)]
  message(sprintf("  After p<%.0e + QC: %d SNPs", P_THRESH, nrow(sub)))

  # Build MarkerID for SNPnexus join
  sub[, MarkerID := paste0(get(cc$chr), ":", get(cc$pos), ":",
                           get(cc$oa), "/", get(cc$ea), ":1")]

  # Read SNPnexus output (already prepared offline)
  snpx <- read.table(meta$snpnexus_file, header = FALSE)
  joined <- inner_join(sub, snpx, by = c("MarkerID" = "V1"))
  message(sprintf("  After SNPnexus rsID mapping: %d SNPs", nrow(joined)))

  # Write to package tmp (clump_data needs file path)
  csv_path <- file.path(PKG_TMP, sprintf("%s_cor_AS.csv", tolower(meta$name)))
  fwrite(joined, csv_path)

  # Read as exposure for clumping
  exp_dat <- read_exposure_data(
    filename = csv_path, sep = ",", snp_col = "V2",
    beta_col = cc$beta, se_col = cc$se, pval_col = cc$pval,
    effect_allele_col = cc$ea, other_allele_col = cc$oa, eaf_col = cc$eaf
  )
  exp_dat <- as.data.frame(exp_dat)

  # Clump in EAS panel
  clumped_keys <- clump_data(exp_dat, clump_kb = KB_THRESH, clump_r2 = R2_WITHIN,
                             clump_p1 = 1, clump_p2 = 1, pop = POP)
  exp_dat <- exp_dat[exp_dat$SNP %in% clumped_keys$SNP, ]
  exp_dat$id.exposure <- meta$name
  exp_dat$exposure    <- meta$name

  # F-statistic
  exp_dat$F_value <- (exp_dat$beta.exposure / exp_dat$se.exposure)^2
  message(sprintf("  Final clumped: %d SNPs (mean F = %.1f)",
                  nrow(exp_dat), mean(exp_dat$F_value)))

  saveRDS(exp_dat, cache)
  exp_dat
}

# ------------------------------------------------------------------
# 2. Helper: load BBJ-format exposure (eGFR), already has rsID
# ------------------------------------------------------------------
prepare_exposure_BBJ <- function(meta) {
  message(sprintf("\n[Exposure-BBJ] %s  (file = %s)", meta$name, meta$gwas_file))
  cache <- file.path(proc_dir, sprintf("%s_clumped.rds", tolower(meta$name)))
  if (file.exists(cache)) {
    message("  Loading cached clumped data."); return(readRDS(cache))
  }
  cc <- meta$cols

  raw <- fread(meta$gwas_file, sep = "auto", header = TRUE,
               nThread = getDTthreads(), na.strings = c("", "NA"),
               showProgress = TRUE)
  sub <- raw[get(cc$pval) < P_THRESH]
  sub <- sub[grepl("^rs", get(cc$snp))]
  message(sprintf("  After p<%.0e + rsID filter: %d SNPs", P_THRESH, nrow(sub)))

  csv_path <- file.path(PKG_TMP, sprintf("%s_cor_AS.csv", tolower(meta$name)))
  fwrite(sub, csv_path)

  exp_dat <- read_exposure_data(
    filename = csv_path, sep = ",", snp_col = cc$snp,
    beta_col = cc$beta, se_col = cc$se, pval_col = cc$pval,
    effect_allele_col = cc$ea, other_allele_col = cc$oa, eaf_col = cc$eaf
  )
  exp_dat <- as.data.frame(exp_dat)

  clumped_keys <- clump_data(exp_dat, clump_kb = KB_THRESH, clump_r2 = R2_WITHIN,
                             clump_p1 = 1, clump_p2 = 1, pop = POP)
  exp_dat <- exp_dat[exp_dat$SNP %in% clumped_keys$SNP, ]
  exp_dat$id.exposure <- meta$name
  exp_dat$exposure    <- meta$name
  exp_dat$F_value     <- (exp_dat$beta.exposure / exp_dat$se.exposure)^2
  message(sprintf("  Final clumped: %d SNPs (mean F = %.1f)",
                  nrow(exp_dat), mean(exp_dat$F_value)))

  saveRDS(exp_dat, cache)
  exp_dat
}

# Dispatcher
prepare_exposure <- function(meta) {
  switch(meta$source,
         "TPMI" = prepare_exposure_TPMI(meta),
         "BBJ"  = prepare_exposure_BBJ(meta),
         stop("Unknown exposure source: ", meta$source))
}

# ------------------------------------------------------------------
# 3. Pool & cross-clump exposures, build SNP_uni with chr:pos
# ------------------------------------------------------------------
pool_and_clump <- function(exposure_list) {
  pooled <- bind_rows(lapply(exposure_list, function(d)
    d[, c("SNP", "pval.exposure")]))
  pooled <- pooled[!duplicated(pooled$SNP), ]
  message(sprintf("\n[Pool] %d unique SNPs across all exposures (pre-clump)",
                  nrow(pooled)))

  # Second clump (looser r2 to remove between-exposure LD only)
  cache <- file.path(proc_dir, "SNP_uni_AS.rds")
  if (file.exists(cache)) {
    SNP_uni <- readRDS(cache)
    message(sprintf("  Loaded cached SNP_uni: %d SNPs", nrow(SNP_uni)))
  } else {
    pooled_renamed <- pooled
    colnames(pooled_renamed) <- c("SNP", "pval")  # clump_data expects 'pval'
    SNP_uni <- clump_data(pooled_renamed, clump_kb = KB_THRESH,
                          clump_r2 = R2_BETWEEN, clump_p1 = 1, clump_p2 = 1,
                          pop = POP)
    saveRDS(SNP_uni, cache)
    message(sprintf("  After 2nd clump (r2<%.3f): %d SNPs", R2_BETWEEN,
                    nrow(SNP_uni)))
  }
  SNP_uni
}

# Attach chr/pos to SNP_uni via the user-prepared Marker_ID.txt
attach_chr_pos <- function(SNP_uni) {
  marker_file <- "Marker_ID.txt"
  if (!file.exists(marker_file)) {
    stop("Marker_ID.txt not found. Run SNPnexus on SNP_uni and save as Marker_ID.txt before continuing.")
  }
  marker_ID <- read.table(marker_file, header = FALSE)
  marker_ID <- data.frame(SNP = marker_ID$V1,
                          CHR = paste0("chr", marker_ID$V3),
                          POS = marker_ID$V4)
  merge(SNP_uni, marker_ID, by = "SNP")
}

# ------------------------------------------------------------------
# 4. Build wide exposure_dat by re-querying each raw GWAS at SNP_uni
# ------------------------------------------------------------------
build_exposure_dat <- function(SNP_uni, exposure_meta) {
  out <- vector("list", length(exposure_meta))
  for (i in seq_along(exposure_meta)) {
    meta <- exposure_meta[[i]]; cc <- meta$cols
    message(sprintf("  Re-querying %s for SNP_uni...", meta$name))
    raw <- fread(meta$gwas_file, sep = "auto", header = TRUE,
                 nThread = getDTthreads(), na.strings = c("", "NA"))

    if (meta$source == "TPMI") {
      hit <- merge(SNP_uni, raw, by.x = c("CHR", "POS"),
                   by.y = c(cc$chr, cc$pos))
    } else {
      hit <- merge(SNP_uni, raw, by.x = "SNP", by.y = cc$snp)
      # add SNP back as column so downstream is uniform
      if (!"SNP" %in% colnames(hit)) hit$SNP <- hit[[cc$snp]]
    }
    if (nrow(hit) == 0) {
      warning(sprintf("0 SNPs for %s after merge with SNP_uni", meta$name))
      out[[i]] <- NULL; next
    }

    df <- data.frame(
      SNP                    = as.character(hit$SNP),
      exposure               = meta$name,
      id.exposure            = meta$name,
      pval.exposure          = as.numeric(hit[[cc$pval]]),
      beta.exposure          = as.numeric(hit[[cc$beta]]),
      se.exposure            = as.numeric(hit[[cc$se]]),
      effect_allele.exposure = as.character(hit[[cc$ea]]),
      other_allele.exposure  = as.character(hit[[cc$oa]]),
      eaf.exposure           = as.numeric(hit[[cc$eaf]]),
      stringsAsFactors = FALSE
    )
    # Sanity check: pval underflow detection
    n_zero_pval <- sum(df$pval.exposure == 0, na.rm = TRUE)
    if (n_zero_pval > 0) {
      message(sprintf("    Warning: %d SNPs have pval=0 after numeric conversion (underflow). Setting to 1e-300.",
                      n_zero_pval))
      df$pval.exposure[df$pval.exposure == 0] <- 1e-300
    }
    out[[i]] <- df
    message(sprintf("    %s: %d SNPs", meta$name, nrow(df)))
  }
  bind_rows(out)
}

# ------------------------------------------------------------------
# 5. Helper: align mvmr_dat columns to canonical exposure name order
#    + SNP-level numeric verification (same logic as EUR v3)
# ------------------------------------------------------------------
align_and_verify_mvmr <- function(mvmr_dat, exposure_names,
                                  expo_dat_AS, test_snp_check = TRUE) {
  # CRITICAL: column order in exposure_beta is given by colnames, NOT by 
  # expname row order. mv_harmonise_data sorts colnames alphabetically 
  # while expname may keep input order. We use colnames as the truth source.
  cur_cols <- colnames(mvmr_dat$exposure_beta)
  if (is.null(cur_cols)) {
    stop("colnames(exposure_beta) is NULL - cannot determine column order.")
  }
  
  miss <- setdiff(exposure_names, cur_cols)
  if (length(miss) > 0) {
    stop(sprintf("Exposure names missing in exposure_beta colnames: %s\n  Found: %s",
                 paste(miss, collapse = ", "),
                 paste(cur_cols, collapse = ", ")))
  }
  
  # Reorder by matching against colnames (the true column order)
  col_ord <- match(exposure_names, cur_cols)
  mvmr_dat$exposure_beta <- mvmr_dat$exposure_beta[, col_ord, drop = FALSE]
  mvmr_dat$exposure_se   <- mvmr_dat$exposure_se[,   col_ord, drop = FALSE]
  mvmr_dat$exposure_pval <- mvmr_dat$exposure_pval[, col_ord, drop = FALSE]
  
  # Reorder expname independently — it has its own row order
  exp_ord <- match(exposure_names, mvmr_dat$expname$id.exposure)
  if (anyNA(exp_ord)) {
    stop(sprintf("Exposure names missing in expname$id.exposure: %s",
                 paste(exposure_names[is.na(exp_ord)], collapse = ", ")))
  }
  mvmr_dat$expname <- mvmr_dat$expname[exp_ord, , drop = FALSE]
  
  colnames(mvmr_dat$exposure_beta) <-
    colnames(mvmr_dat$exposure_se)   <-
    colnames(mvmr_dat$exposure_pval) <- exposure_names
  
  # SNP-level numeric check
  if (test_snp_check) {
    test_snps <- rownames(mvmr_dat$exposure_beta)[1:min(3, nrow(mvmr_dat$exposure_beta))]
    for (snp in test_snps) {
      for (i in seq_along(exposure_names)) {
        ex <- expo_dat_AS$beta.exposure[
          expo_dat_AS$SNP == snp &
          expo_dat_AS$id.exposure == exposure_names[i]]
        if (length(ex) == 0 || is.na(ex[1])) next
        ob <- mvmr_dat$exposure_beta[snp, i]
        # mv_harmonise_data may flip the allele coding -> beta sign flip;
        # accept either exact match or sign-flipped match
        if (!(abs(ex[1] - ob) < 1e-6 || abs(ex[1] + ob) < 1e-6)) {
          stop(sprintf(
            "Numeric mismatch: SNP=%s, exposure=%s\n  expo_dat beta=%.6f, mvmr_dat beta=%.6f",
            snp, exposure_names[i], ex[1], ob))
        }
      }
    }
  }

  message("\n  Column alignment check (PASSED):")
  chk <- data.frame(
    pos        = seq_along(exposure_names),
    name       = exposure_names,
    expname_id = mvmr_dat$expname$id.exposure,
    matrix_col = colnames(mvmr_dat$exposure_beta)
  )
  print(chk, row.names = FALSE)
  mvmr_dat
}

# ------------------------------------------------------------------
# 6. Helper: full MVMR pipeline for one outcome (local file)
# ------------------------------------------------------------------
run_mvmr_pair_local <- function(expo_dat_AS, SNP_uni, exposure_names,
                                outcome_name, outcome_file, outcome_label) {
  message(sprintf("\n=== Outcome: %s (%s) ===", outcome_name, outcome_file))

  # 6a. Read outcome, restrict to SNP_uni
  raw_out <- fread(outcome_file, sep = "auto", header = TRUE,
                   nThread = getDTthreads(), na.strings = c("", "NA"))
  out_sub <- merge(SNP_uni, raw_out, by.x = "SNP", by.y = OUT_COLS$snp)
  if (nrow(out_sub) == 0) {
    warning(sprintf("0 SNPs for outcome %s", outcome_name)); return(NULL)
  }
  csv_path <- file.path(proc_dir, sprintf("%s_outcome.csv", outcome_name))
  fwrite(out_sub, csv_path)

  out_dat <- read_outcome_data(
    snps = SNP_uni$SNP, filename = csv_path, sep = ",",
    snp_col = "SNP",
    beta_col = OUT_COLS$beta, se_col = OUT_COLS$se,
    effect_allele_col = OUT_COLS$ea, other_allele_col = OUT_COLS$oa,
    eaf_col = OUT_COLS$eaf, pval_col = OUT_COLS$pval
  )
  out_dat$id.outcome <- outcome_name
  out_dat$outcome    <- outcome_label

  # 6b. Harmonise + align
  mvmr_dat <- mv_harmonise_data(expo_dat_AS, out_dat)
  mvmr_dat <- align_and_verify_mvmr(mvmr_dat, exposure_names, expo_dat_AS,
                                    test_snp_check = TRUE)
  k_snp <- nrow(mvmr_dat$exposure_beta)
  message(sprintf("  SNPs after harmonisation: %d", k_snp))

  # 6c. MVMR-IVW (primary)
  mr_input <- mr_mvinput(
    bx = mvmr_dat$exposure_beta, bxse = mvmr_dat$exposure_se,
    by = mvmr_dat$outcome_beta,  byse = mvmr_dat$outcome_se,
    exposure = exposure_names, outcome = outcome_name
  )
  ivw <- mr_mvivw(mr_input, model = "default", correl = FALSE,
                  distribution = "normal", alpha = 0.05)
  ivw_tbl <- data.frame(
    outcome = outcome_name, exposure = exposure_names, method = "MVMR-IVW",
    nsnp = k_snp,
    b = ivw@Estimate, se = ivw@StdError, pval = ivw@Pvalue,
    OR = exp(ivw@Estimate),
    OR_lci = exp(ivw@CILower), OR_uci = exp(ivw@CIUpper),
    stringsAsFactors = FALSE
  )

  # 6d. MVMR-Egger (sensitivity)
  egger <- mr_mvegger(mr_input, orientate = 1, correl = FALSE,
                      distribution = "normal", alpha = 0.05)
  egger_tbl <- data.frame(
    outcome = outcome_name, exposure = exposure_names, method = "MVMR-Egger",
    nsnp = k_snp,
    b = egger@Estimate, se = egger@StdError.Est, pval = egger@Pvalue.Est,
    OR = exp(egger@Estimate),
    OR_lci = exp(egger@CILower.Est), OR_uci = exp(egger@CIUpper.Est),
    stringsAsFactors = FALSE
  )
  egger_intercept <- data.frame(
    outcome = outcome_name,
    intercept = egger@Intercept, se = egger@StdError.Int,
    pval = egger@Pvalue.Int,
    ci_lower = egger@CILower.Int, ci_upper = egger@CIUpper.Int,
    interpretation = ifelse(egger@Pvalue.Int < 0.05,
                            "Evidence of directional pleiotropy",
                            "No evidence of directional pleiotropy"),
    stringsAsFactors = FALSE
  )

  # 6e. cF + Q
  F_data <- format_mvmr(BXGs = mvmr_dat$exposure_beta,
                        BYG  = mvmr_dat$outcome_beta,
                        seBXGs = mvmr_dat$exposure_se,
                        seBYG  = mvmr_dat$outcome_se,
                        RSID = rownames(mvmr_dat$exposure_beta))
  cF_raw <- strength_mvmr(r_input = F_data, gencov = GENCOV)
  cF_tbl <- data.frame(
    outcome = outcome_name, exposure = exposure_names,
    cond_F = as.numeric(cF_raw),
    weak_F5 = as.numeric(cF_raw) < 5,
    weak_F10 = as.numeric(cF_raw) < 10,
    stringsAsFactors = FALSE
  )
  Q_raw <- pleiotropy_mvmr(r_input = F_data, gencov = GENCOV)
  Q_tbl <- data.frame(
    outcome = outcome_name,
    Q = Q_raw$Qstat, df = k_snp - length(exposure_names),
    Q_pval = Q_raw$Qpval,
    interpretation = ifelse(Q_raw$Qpval < 0.05,
                            "Heterogeneity / horizontal pleiotropy detected",
                            "No detectable heterogeneity / pleiotropy"),
    stringsAsFactors = FALSE
  )

  fwrite(F_data, file.path(proc_dir, sprintf("%s_F_data.csv", outcome_name)))
  saveRDS(mvmr_dat, file.path(proc_dir, sprintf("%s_mvmr_dat.rds", outcome_name)))

  list(ivw = ivw_tbl, egger = egger_tbl, intercept = egger_intercept,
       cF = cF_tbl, Q = Q_tbl, nsnp = k_snp)
}

# ==================================================================
#  EXECUTE PIPELINE
# ==================================================================

# ------------------------------------------------------------------
# 7. Prepare all exposures
# ------------------------------------------------------------------
exposures <- lapply(exposure_meta, prepare_exposure)
names(exposures) <- exposure_names

instrument_summary <- bind_rows(lapply(seq_along(exposures), function(i) {
  d <- exposures[[i]]
  data.frame(
    exposure = exposure_names[i],
    source   = exposure_meta[[i]]$source,
    n_snp    = nrow(d),
    mean_F   = round(mean(d$F_value), 1),
    min_F    = round(min(d$F_value),  1),
    stringsAsFactors = FALSE
  )
}))
fwrite(instrument_summary, file.path(out_dir, "00_instrument_summary_EAS.csv"))
print(instrument_summary)

# ------------------------------------------------------------------
# 8. Pool, cross-clump, attach chr/pos, build wide expo_dat_AS
# ------------------------------------------------------------------
SNP_uni    <- pool_and_clump(exposures)
SNP_uni    <- attach_chr_pos(SNP_uni)
expo_dat_AS <- build_exposure_dat(SNP_uni, exposure_meta)

message(sprintf("\n[Final exposure_dat]  total rows = %d  (across %d exposures)",
                nrow(expo_dat_AS), length(exposures)))
saveRDS(expo_dat_AS, file.path(proc_dir, "expo_dat_AS.rds"))
saveRDS(SNP_uni,     file.path(proc_dir, "SNP_uni_with_chrpos.rds"))

# ------------------------------------------------------------------
# 9. Run MVMR for each outcome
# ------------------------------------------------------------------
all_ivw <- list(); all_egger <- list(); all_intercept <- list()
all_cF  <- list(); all_Q     <- list()

for (om in outcome_meta) {
  res <- run_mvmr_pair_local(expo_dat_AS, SNP_uni, exposure_names,
                             om$name, om$file, om$label)
  if (is.null(res)) next
  all_ivw[[length(all_ivw)+1]]             <- res$ivw
  all_egger[[length(all_egger)+1]]         <- res$egger
  all_intercept[[length(all_intercept)+1]] <- res$intercept
  all_cF[[length(all_cF)+1]]               <- res$cF
  all_Q[[length(all_Q)+1]]                 <- res$Q
}

ivw_tbl       <- do.call(rbind, all_ivw)
egger_tbl     <- do.call(rbind, all_egger)
intercept_tbl <- do.call(rbind, all_intercept)
cF_tbl        <- do.call(rbind, all_cF)
Q_tbl         <- do.call(rbind, all_Q)

# ------------------------------------------------------------------
# 10. FDR + flag weak instruments
# ------------------------------------------------------------------
ivw_tbl <- ivw_tbl %>%
  mutate(CI_text = sprintf("%.3f (%.3f, %.3f)", OR, OR_lci, OR_uci)) %>%
  group_by(outcome) %>%
  mutate(p_fdr = p.adjust(pval, method = "BH")) %>%
  ungroup() %>%
  left_join(cF_tbl[, c("outcome", "exposure", "cond_F", "weak_F5", "weak_F10")],
            by = c("outcome", "exposure")) %>%
  mutate(flag = case_when(
    weak_F5  ~ "WEAK (cF<5)",
    weak_F10 ~ "Borderline (cF 5-10)",
    TRUE     ~ "OK (cF>=10)"
  ))

egger_tbl <- egger_tbl %>%
  mutate(CI_text = sprintf("%.3f (%.3f, %.3f)", OR, OR_lci, OR_uci))

# ------------------------------------------------------------------
# 11. Save outputs
# ------------------------------------------------------------------
fwrite(ivw_tbl,       file.path(out_dir, "01_MVMR_IVW_EAS.csv"))
fwrite(egger_tbl,     file.path(out_dir, "02_MVMR_Egger_EAS.csv"))
fwrite(intercept_tbl, file.path(out_dir, "03_Egger_intercept_pleiotropy_EAS.csv"))
fwrite(cF_tbl,        file.path(out_dir, "04_conditional_F_EAS.csv"))
fwrite(Q_tbl,         file.path(out_dir, "05_Q_heterogeneity_EAS.csv"))

pub_tbl <- ivw_tbl %>%
  left_join(Q_tbl[, c("outcome", "Q", "Q_pval")], by = "outcome") %>%
  left_join(intercept_tbl[, c("outcome", "intercept", "pval")] %>%
              rename(egger_intercept = intercept, egger_intercept_pval = pval),
            by = "outcome") %>%
  select(outcome, exposure, nsnp, b, se, OR, OR_lci, OR_uci, CI_text,
         pval, p_fdr, cond_F, flag,
         Q, Q_pval, egger_intercept, egger_intercept_pval) %>%
  arrange(outcome, exposure)
fwrite(pub_tbl, file.path(out_dir, "00_MVMR_publication_table_EAS.csv"))

# ------------------------------------------------------------------
# 12. Forest plot per outcome
# ------------------------------------------------------------------
plot_mvmr_forest <- function(outcome_name, ivw_sub, egger_sub) {
  combined <- bind_rows(
    ivw_sub   %>% mutate(method = "MVMR-IVW"),
    egger_sub %>% mutate(method = "MVMR-Egger")
  ) %>% mutate(exposure = factor(exposure, levels = rev(exposure_names)))

  p <- ggplot(combined, aes(x = OR, y = exposure, color = method)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
    geom_pointrange(aes(xmin = OR_lci, xmax = OR_uci),
                    position = position_dodge(width = 0.5), size = 0.4) +
    scale_x_log10() +
    scale_color_manual(values = c("MVMR-IVW" = "#0072B2", "MVMR-Egger" = "#D55E00")) +
    labs(title = sprintf("MVMR results for %s (EAS)", outcome_name),
         x = "OR (95% CI), log scale", y = NULL, color = NULL) +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank())

  ggsave(file.path(plot_dir, sprintf("%s_MVMR_forest.pdf", outcome_name)),
         p, width = 6.5, height = 4)
}

for (om in outcome_meta) {
  plot_mvmr_forest(
    om$name,
    ivw_tbl   %>% filter(outcome == om$name),
    egger_tbl %>% filter(outcome == om$name)
  )
}

# ------------------------------------------------------------------
# 13. Console summary
# ------------------------------------------------------------------
message("\n=== EAS MVMR analysis complete ===")
message(sprintf("Results saved to: %s", out_dir))
cat("\n--- Conditional F (weak instrument flags) ---\n");  print(cF_tbl)
cat("\n--- Q heterogeneity / pleiotropy test ---\n");      print(Q_tbl)
cat("\n--- Egger intercept pleiotropy test ---\n");        print(intercept_tbl)
cat("\n--- IVW publication table ---\n");                  print(pub_tbl)
