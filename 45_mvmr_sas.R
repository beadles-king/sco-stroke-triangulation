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
# Multivariable Mendelian Randomization: South Asian Ancestry Cohort  (v2)
# Exposures (k = 5): FG, SBP, LDL, BMI, eGFR  (pan-UKBB CSA via OpenGWAS)
# Outcomes (2)    : IS and any-stroke (GIGASTROKE-SAS)
#                   NB: SAS any-stroke has 178 additional cases vs IS (mostly
#                   ICH/SAH); the two outcomes are highly correlated but
#                   provide complementary signal especially for FG (β=0.09 vs
#                   0.11 in univariate IVW). Both reported for completeness.
# Methods :
#   - Primary  : MVMR-IVW (MendelianRandomization::mr_mvivw)
#   - Sensitivity : MVMR-Egger (intercept tests pleiotropy)
#   - Diagnostics : Conditional F-statistic per exposure (MVMR::strength_mvmr)
#                   Conditional Q-statistic (MVMR::pleiotropy_mvmr)
# Weak-instrument threshold : cF >= 5 (Sanderson 2021 Stat Med, lenient cutoff)
# Multiple testing : BH-FDR within outcome on IVW p-values
# Clumping : r2 < 0.001 (MVMR-recommended strict threshold; consistent with
#            EUR/EAS MVMR pipelines; default of mv_extract_exposures)
# Instrument p threshold : 1e-5 (SAS-specific relaxation due to small CSA
#                           sample size; matches univariate SAS MR pipeline)
#
# IMPORTANT — SAS-SPECIFIC LIMITATIONS:
#   pan-UKBB CSA exposure GWAS (n ~ 8,000) yield underpowered MVMR.
#   In typical runs, all five conditional F-statistics fall below 10, with
#   most below the lenient threshold of 5 — indicating substantial weak
#   instrument bias. SAS MVMR results should therefore be presented as
#   EXPLORATORY / SUPPLEMENTARY, not as primary causal evidence. The script
#   automatically flags this in 04_conditional_F.csv (column `flag`).
#
# Outcome quirk: GIGASTROKE-SAS files have NO rsID column. We use a
# pre-built chr:pos -> rsID map derived from the per-exposure clumped
# univariate MR objects (FG_GWAS, BMI_GWAS, ...) which DO carry chr/pos.
# This map is then applied to the outcome to attach rsIDs before harmonise.
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

work_dir <- file.path(GWAS_DIR, "SAS")
setwd(work_dir)

out_dir  <- file.path(RAWDAT, "MR/SAS/MVMR_res")
plot_dir <- file.path(out_dir, "plots")
proc_dir <- file.path(out_dir, "processed")
invisible(lapply(c(out_dir, plot_dir, proc_dir),
                 dir.create, showWarnings = FALSE, recursive = TRUE))

# --- Constants ------------------------------------------------------
P_THRESH  <- 1e-5     # SAS-specific (small sample; matches univariate SAS MR)
R2_THRESH <- 0.001    # MVMR-recommended strict clumping
KB_THRESH <- 10000
F_THRESH  <- 5        # Sanderson 2021 lenient cutoff
GENCOV    <- 0
POP       <- "SAS"

# --- Exposure metadata ----------------------------------------------
# Same five pan-UKBB CSA IDs used in the univariate SAS pipeline.
# The verify_label_contains regex is checked against OpenGWAS server-returned
# trait labels — guards against any future ID/name swap.
exposure_meta <- list(
  list(name = "FG",   opengwas_id = "ebi-a-GCST90002234",
       verify_label_contains = "[Ff]asting"),
  list(name = "LDL",  opengwas_id = "ukb-e-30780_CSA",
       verify_label_contains = "LDL|Low.density"),
  list(name = "SBP",  opengwas_id = "ukb-e-4080_CSA",
       verify_label_contains = "[Ss]ystolic|blood pressure"),
  list(name = "BMI",  opengwas_id = "ukb-e-23104_CSA",
       verify_label_contains = "[Bb]ody mass|BMI"),
  list(name = "eGFR", opengwas_id = "ukb-e-recode6_CSA",
       verify_label_contains = "glomerular|creatinine|eGFR")
)
exposure_ids   <- sapply(exposure_meta, `[[`, "opengwas_id")
exposure_names <- sapply(exposure_meta, `[[`, "name")
exposure_regex <- sapply(exposure_meta, `[[`, "verify_label_contains")

# --- Outcome metadata -----------------------------------------------
# GIGASTROKE-SAS files lack rsID; matched to exposures via chr:pos.
# All-stroke outcome (GCST90104559) is omitted in the main MVMR because
# it is near-identical to IS for SAS (~95% of stroke cases are ischemic).
OUT_COLS <- list(
  chr  = "chromosome",       pos  = "base_pair_location",
  beta = "beta",             se   = "standard_error",
  ea   = "effect_allele",    oa   = "other_allele",
  eaf  = "effect_allele_frequency", pval = "p_value"
)

outcome_meta <- list(
  list(name = "AllStroke",
       file = file.path(RAWDAT, "GWAS_C/GIGA_stroke_SAS.tsv.gz"),
       label = "Any stroke (GIGASTROKE-SAS)"),
  list(name = "IS",
       file = file.path(RAWDAT, "GWAS_C/GIGA_IS_SAS.tsv.gz"),
       label = "Ischemic stroke (GIGASTROKE-SAS)")
)

# ------------------------------------------------------------------
# 1. Build chr:pos -> rsID map from existing univariate exposure objects
#    (FG_GWAS, LDL_GWAS, SBP_GWAS, BMI_GWAS, eGFR_GWAS already in env)
# ------------------------------------------------------------------
# ------------------------------------------------------------------
# Helper: build a single univariate exposure (mirrors univariate SAS MR)
# Uses the SAME parameters: p<5e-8, clump r2=0.001, kb=10000, F>10, SAS panel
# Note: these parameters match the published-default SAS univariate MR script;
# they are independent of the MVMR p threshold (1e-5) — the snp_map only
# needs chr/pos for *some* SNPs that overlap with the MVMR exposure pool.
# ------------------------------------------------------------------
build_univariate_one <- function(name, opengwas_id) {
  cache <- file.path(proc_dir, sprintf("univariate_%s.rds", tolower(name)))
  if (file.exists(cache)) {
    d <- readRDS(cache)
    message(sprintf("  [%s] loaded cached univariate exposure: %d SNPs",
                    name, nrow(d)))
    return(d)
  }
  message(sprintf("  [%s] downloading instruments from %s ...", name, opengwas_id))
  d <- tryCatch(
    extract_instruments(outcomes = opengwas_id, p1 = P_THRESH, clump = FALSE),
    error = function(e) { message("    failed: ", e$message); NULL })
  if (is.null(d) || nrow(d) == 0) return(NULL)
  
  d <- clump_data(d, clump_r2 = R2_THRESH, clump_kb = KB_THRESH, pop = POP)
  d$F_value <- (d$beta.exposure / d$se.exposure)^2
  d <- subset(d, F_value > 10)
  saveRDS(d, cache)
  message(sprintf("  [%s] final: %d SNPs after clumping + F>10", name, nrow(d)))
  d
}

# ------------------------------------------------------------------
# Build chr:pos -> rsID map by running univariate exposure extraction
# This is self-contained: does NOT depend on the univariate SAS MR script
# being sourced first. Each exposure is cached in proc_dir.
# ------------------------------------------------------------------
build_snp_map <- function() {
  message("[snp_map] Building univariate exposures for chr:pos lookup...")
  univ_list <- list()
  for (m in exposure_meta) {
    d <- build_univariate_one(m$name, m$opengwas_id)
    if (!is.null(d)) univ_list[[m$name]] <- d
  }
  if (length(univ_list) == 0) {
    stop("Failed to retrieve any univariate exposure for snp_map construction.")
  }
  
  snp_map <- bind_rows(lapply(univ_list, function(d) {
    needed <- c("SNP", "chr.exposure", "pos.exposure")
    if (!all(needed %in% colnames(d))) {
      stop("Univariate exposure missing chr.exposure/pos.exposure columns")
    }
    data.frame(SNP = d$SNP, chr = d$chr.exposure, pos = d$pos.exposure,
               stringsAsFactors = FALSE)
  }))
  snp_map <- unique(snp_map)
  snp_map$id <- paste0(snp_map$chr, ":", snp_map$pos)
  message(sprintf("[snp_map] %d unique SNP-chr-pos triples built", nrow(snp_map)))
  
  # Save to disk for inspection
  saveRDS(snp_map, file.path(proc_dir, "snp_map.rds"))
  fwrite(snp_map, file.path(proc_dir, "snp_map.csv"))
  
  snp_map
}

# ------------------------------------------------------------------
# 2. Extract & harmonise MV exposures (cached; reused across outcomes)
# ------------------------------------------------------------------
extract_mv_exposures <- function() {
  exp_cache <- file.path(proc_dir, "mv_exposure_dat.rds")
  if (file.exists(exp_cache)) {
    d <- readRDS(exp_cache)
    message(sprintf("[exposures] Loaded cached: %d unique SNPs",
                    length(unique(d$SNP))))
    return(d)
  }
  message("[exposures] Extracting MV exposures from OpenGWAS (SAS panel)...")
  d <- mv_extract_exposures(
    id_exposure          = exposure_ids,
    pval_threshold       = P_THRESH,
    clump_r2             = R2_THRESH,
    clump_kb             = KB_THRESH,
    harmonise_strictness = 2,
    find_proxies         = TRUE,
    pop                  = POP
  )
  saveRDS(d, exp_cache)
  message(sprintf("  Done. %d unique SNPs across %d exposures.",
                  length(unique(d$SNP)), length(exposure_ids)))
  d
}

# ------------------------------------------------------------------
# 3. Helper: align mvmr_dat columns to canonical name order
#    + label-regex verification + SNP-level numeric core check
#    (Same logic as EUR v3 / EAS v2 with bug fix.)
# ------------------------------------------------------------------
align_and_verify_mvmr <- function(mvmr_dat,
                                  exposure_ids,
                                  exposure_names,
                                  exposure_regex,
                                  exposure_dat,
                                  test_snp_check = TRUE) {

  # Truth source: column order is given by colnames(exposure_beta), NOT by
  # expname row order. Always rely on colnames for matrix reordering.
  cur_cols <- colnames(mvmr_dat$exposure_beta)
  if (is.null(cur_cols)) {
    stop("colnames(exposure_beta) is NULL - cannot determine column order.")
  }

  # ---- 3a. Verify all requested IDs are present in matrix columns ----
  miss <- setdiff(exposure_ids, cur_cols)
  if (length(miss) > 0) {
    stop(sprintf("Exposure IDs missing in exposure_beta colnames: %s\n  Found: %s",
                 paste(miss, collapse = ", "),
                 paste(cur_cols, collapse = ", ")))
  }

  # ---- 3b. Verify each ID's expname label matches expected phenotype ----
  cur_ids_in_expname    <- mvmr_dat$expname$id.exposure
  cur_labels_in_expname <- mvmr_dat$expname$exposure
  mismatch_msgs <- character(0)
  for (i in seq_along(exposure_ids)) {
    idx <- match(exposure_ids[i], cur_ids_in_expname)
    if (is.na(idx)) {
      mismatch_msgs <- c(mismatch_msgs, sprintf(
        "  - '%s' (%s) not found in expname$id.exposure",
        exposure_names[i], exposure_ids[i]))
      next
    }
    server_lbl <- cur_labels_in_expname[idx]
    if (!grepl(exposure_regex[i], server_lbl)) {
      mismatch_msgs <- c(mismatch_msgs, sprintf(
        "  - Expected '%s' (%s) to match /%s/ but server returned: '%s'",
        exposure_names[i], exposure_ids[i],
        exposure_regex[i], server_lbl))
    }
  }
  if (length(mismatch_msgs) > 0) {
    stop("\nID->NAME MISMATCH DETECTED in exposure_meta:\n",
         paste(mismatch_msgs, collapse = "\n"),
         "\n  Fix the opengwas_id/name pairing in exposure_meta and re-run.")
  }

  # ---- 3c. Reorder matrix columns based on colnames (truth) ----
  col_ord <- match(exposure_ids, cur_cols)
  mvmr_dat$exposure_beta <- mvmr_dat$exposure_beta[, col_ord, drop = FALSE]
  mvmr_dat$exposure_se   <- mvmr_dat$exposure_se[,   col_ord, drop = FALSE]
  mvmr_dat$exposure_pval <- mvmr_dat$exposure_pval[, col_ord, drop = FALSE]

  # ---- 3d. Reorder expname independently (its own row order) ----
  exp_ord <- match(exposure_ids, cur_ids_in_expname)
  mvmr_dat$expname <- mvmr_dat$expname[exp_ord, , drop = FALSE]

  # ---- 3e. Apply human-readable column names ----
  colnames(mvmr_dat$exposure_beta) <-
    colnames(mvmr_dat$exposure_se)   <-
    colnames(mvmr_dat$exposure_pval) <- exposure_names

  # ---- 3f. SNP-level numeric core check ----
  if (test_snp_check) {
    test_snps <- rownames(mvmr_dat$exposure_beta)[1:min(3, nrow(mvmr_dat$exposure_beta))]
    for (snp in test_snps) {
      for (i in seq_along(exposure_ids)) {
        ex <- exposure_dat$beta.exposure[
          exposure_dat$SNP == snp &
          exposure_dat$id.exposure == exposure_ids[i]]
        if (length(ex) == 0 || is.na(ex[1])) next
        ob <- mvmr_dat$exposure_beta[snp, i]
        # Allow exact match or sign-flipped match (allele alignment)
        if (!(abs(ex[1] - ob) < 1e-6 || abs(ex[1] + ob) < 1e-6)) {
          stop(sprintf(
            "Numeric mismatch: SNP=%s, exposure=%s [%s]\n  exposure_dat=%.6f, mvmr_dat=%.6f",
            snp, exposure_names[i], exposure_ids[i], ex[1], ob))
        }
      }
    }
  }

  # ---- 3g. Print confirmation ----
  message("\n  Column alignment check (PASSED):")
  chk <- data.frame(
    pos        = seq_along(exposure_names),
    name       = exposure_names,
    id         = exposure_ids,
    server_lbl = mvmr_dat$expname$exposure,
    matrix_col = colnames(mvmr_dat$exposure_beta)
  )
  print(chk, row.names = FALSE)
  mvmr_dat
}

# ------------------------------------------------------------------
# 4. Helper: read local outcome, attach rsID via chr:pos map, format
# ------------------------------------------------------------------
prepare_outcome_local <- function(outcome_meta_one, snp_map, exposure_snps) {
  message(sprintf("\n[outcome] Loading %s from %s",
                  outcome_meta_one$name, outcome_meta_one$file))

  raw <- fread(outcome_meta_one$file, sep = "auto", header = TRUE,
               nThread = getDTthreads(), na.strings = c("", "NA"),
               showProgress = TRUE)

  # Build chr:pos id, attach rsID via snp_map (left join on id)
  raw[, id := paste0(get(OUT_COLS$chr), ":", get(OUT_COLS$pos))]
  hit <- merge(raw, snp_map[, c("SNP", "id")], by = "id", all.x = FALSE)
  if (nrow(hit) == 0) {
    warning(sprintf("0 SNPs after chr:pos match for outcome %s",
                    outcome_meta_one$name))
    return(NULL)
  }
  message(sprintf("  Outcome SNPs after chr:pos match: %d", nrow(hit)))

  # Deduplicate: same SNP could be matched twice if multiple chr:pos rows
  # share an rsID upstream
  hit <- hit[!duplicated(hit$SNP), ]

  out_dat <- format_data(
    dat   = as.data.frame(hit),
    type  = "outcome",
    snps  = exposure_snps,
    snp_col = "SNP",
    beta_col = OUT_COLS$beta, se_col = OUT_COLS$se,
    effect_allele_col = OUT_COLS$ea, other_allele_col = OUT_COLS$oa,
    eaf_col = OUT_COLS$eaf, pval_col = OUT_COLS$pval
  )
  out_dat$id.outcome <- outcome_meta_one$name
  out_dat$outcome    <- outcome_meta_one$label

  message(sprintf("  Outcome SNPs after format_data: %d", nrow(out_dat)))
  out_dat
}

# ------------------------------------------------------------------
# 5. Helper: full MVMR pipeline for one outcome
# ------------------------------------------------------------------
run_mvmr_pair <- function(exposure_dat, exposure_ids, exposure_names,
                          exposure_regex, snp_map,
                          outcome_meta_one) {
  outcome_name  <- outcome_meta_one$name
  outcome_label <- outcome_meta_one$label
  message(sprintf("\n=== Outcome: %s ===", outcome_name))

  # 5a. Outcome
  out_dat <- prepare_outcome_local(outcome_meta_one, snp_map, exposure_dat$SNP)
  if (is.null(out_dat) || nrow(out_dat) == 0) return(NULL)

  # 5b. Harmonise + align
  mvmr_dat <- mv_harmonise_data(exposure_dat, out_dat, harmonise_strictness = 2)
  mvmr_dat <- align_and_verify_mvmr(mvmr_dat, exposure_ids, exposure_names,
                                    exposure_regex, exposure_dat,
                                    test_snp_check = TRUE)
  k_snp <- nrow(mvmr_dat$exposure_beta)
  message(sprintf("  SNPs after harmonisation: %d", k_snp))

  # 5c. MendelianRandomization input
  mr_input <- mr_mvinput(
    bx = mvmr_dat$exposure_beta, bxse = mvmr_dat$exposure_se,
    by = mvmr_dat$outcome_beta,  byse = mvmr_dat$outcome_se,
    exposure = exposure_names, outcome = outcome_name
  )

  # 5d. MVMR-IVW (primary)
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

  # 5e. MVMR-Egger (sensitivity + intercept pleiotropy test)
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

  # 5f. cF + Q
  F_data <- format_mvmr(BXGs = mvmr_dat$exposure_beta,
                        BYG  = mvmr_dat$outcome_beta,
                        seBXGs = mvmr_dat$exposure_se,
                        seBYG  = mvmr_dat$outcome_se,
                        RSID = rownames(mvmr_dat$exposure_beta))
  cF_raw <- strength_mvmr(r_input = F_data, gencov = GENCOV)
  cF_tbl <- data.frame(
    outcome = outcome_name, exposure = exposure_names,
    cond_F = as.numeric(cF_raw),
    weak_F5  = as.numeric(cF_raw) < 5,
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
# 6. Build chr:pos -> rsID map and extract MV exposures
# ------------------------------------------------------------------
snp_map      <- build_snp_map()
exposure_dat <- extract_mv_exposures()

matched_n <- sum(exposure_dat$SNP %in% snp_map$SNP)
message(sprintf("\n[chr:pos coverage] %d / %d MV exposure SNPs found in snp_map (%.1f%%)",
                matched_n, length(unique(exposure_dat$SNP)),
                100 * matched_n / length(unique(exposure_dat$SNP))))

# ------------------------------------------------------------------
# 7. Run MVMR for each outcome
# ------------------------------------------------------------------
all_ivw <- list(); all_egger <- list(); all_intercept <- list()
all_cF  <- list(); all_Q     <- list()

for (om in outcome_meta) {
  res <- run_mvmr_pair(exposure_dat, exposure_ids, exposure_names,
                       exposure_regex, snp_map, om)
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
# 8. FDR + flag weak instruments
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
# 9. SAS-specific weak-instrument advisory
# ------------------------------------------------------------------
n_weak  <- sum(cF_tbl$weak_F5)
n_total <- nrow(cF_tbl)
if (n_weak >= n_total / 2) {
  cat("\n", strrep("=", 70), "\n", sep = "")
  cat("SAS MVMR ADVISORY:\n",
      sprintf("  %d / %d exposure-outcome pairs have cF < 5 (severe weak instrument).\n",
              n_weak, n_total),
      "  Per Sanderson 2021, results should be treated as EXPLORATORY.\n",
      "  Recommended: present in Supplementary, not main text.\n",
      "  See Methods/Limitations templates in this script's header.\n", sep = "")
  cat(strrep("=", 70), "\n", sep = "")
}

# ------------------------------------------------------------------
# 10. Save outputs
# ------------------------------------------------------------------
fwrite(ivw_tbl,       file.path(out_dir, "01_MVMR_IVW_SAS.csv"))
fwrite(egger_tbl,     file.path(out_dir, "02_MVMR_Egger_SAS.csv"))
fwrite(intercept_tbl, file.path(out_dir, "03_Egger_intercept_pleiotropy_SAS.csv"))
fwrite(cF_tbl,        file.path(out_dir, "04_conditional_F_SAS.csv"))
fwrite(Q_tbl,         file.path(out_dir, "05_Q_heterogeneity_SAS.csv"))

pub_tbl <- ivw_tbl %>%
  left_join(Q_tbl[, c("outcome", "Q", "Q_pval")], by = "outcome") %>%
  left_join(intercept_tbl[, c("outcome", "intercept", "pval")] %>%
              rename(egger_intercept = intercept, egger_intercept_pval = pval),
            by = "outcome") %>%
  select(outcome, exposure, nsnp, b, se, OR, OR_lci, OR_uci, CI_text,
         pval, p_fdr, cond_F, flag,
         Q, Q_pval, egger_intercept, egger_intercept_pval) %>%
  arrange(outcome, exposure)
fwrite(pub_tbl, file.path(out_dir, "00_MVMR_publication_table_SAS.csv"))

# ------------------------------------------------------------------
# 11. Forest plot per outcome (IVW vs Egger)
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
    labs(title = sprintf("MVMR results for %s (SAS) — exploratory", outcome_name),
         subtitle = "Note: most exposures have cF < 10; interpret with caution",
         x = "OR (95% CI), log scale", y = NULL, color = NULL) +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank(),
          plot.subtitle = element_text(color = "grey40", size = 9))

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
# 12. Console summary
# ------------------------------------------------------------------
message("\n=== SAS MVMR analysis complete ===")
message(sprintf("Results saved to: %s", out_dir))
cat("\n--- Conditional F (weak instrument flags) ---\n");  print(cF_tbl)
cat("\n--- Q heterogeneity / pleiotropy test ---\n");      print(Q_tbl)
cat("\n--- Egger intercept pleiotropy test ---\n");        print(intercept_tbl)
cat("\n--- IVW publication table ---\n");                  print(pub_tbl)
