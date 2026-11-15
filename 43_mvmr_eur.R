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
# Multivariable Mendelian Randomization: European Ancestry Cohort  (v3)
# Exposures (k = 5): FG, SBP, LDL, BMI, eGFR  (jointly modeled)
# Outcomes  (4)    : All stroke, IS, ICH, SAH
# Methods :
#   - Primary  : MVMR-IVW (MendelianRandomization::mr_mvivw)
#   - Sensitivity : MVMR-Egger (intercept tests pleiotropy)
#   - Diagnostics : Conditional F-statistic per exposure (MVMR::strength_mvmr)
#                   Conditional Q-statistic (MVMR::pleiotropy_mvmr)
#                     -> Joint heterogeneity/pleiotropy test in MVMR framework
# Weak-instrument threshold : cF >= 5 (Sanderson 2021 Stat Med, lenient cutoff)
#   Note: gencov = 0 is assumed (no genetic covariance between exposures);
#         this is a conservative simplification when phenotypic covariance is
#         not available from summary data.
# Multiple testing : BH-FDR within each outcome, applied to IVW p-values
# Clumping : r2 < 0.001 (MVMR default, stricter than univariate r2 < 0.01;
#            consistent with Sanderson 2019, 2021 recommendations)
#
# CRITICAL FIX in v3:
#   v2 had wrong OpenGWAS ID->name mapping (SBP and LDL were swapped).
#   Verified IDs (from OpenGWAS server response):
#     ebi-a-GCST90002232 = Fasting glucose
#     ieu-b-4846 = LDL cholesterol     <-- was mislabeled SBP in v2
#     ieu-b-4818 = Systolic blood pressure  <-- was mislabeled LDL in v2
#     ieu-b-4816 = Body mass index
#     ieu-b-4832 = Estimated glomerular filtration rate
#   The validation block below now flags this kind of error automatically.
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

work_dir <- file.path(GWAS_DIR, "EU/MVMR")
setwd(work_dir)

out_dir  <- file.path(RAWDAT, "MR/EUR/MVMR_res")
plot_dir <- file.path(out_dir, "plots")
proc_dir <- file.path(out_dir, "processed")
invisible(lapply(c(out_dir, plot_dir, proc_dir),
                 dir.create, showWarnings = FALSE, recursive = TRUE))

# --- Constants ------------------------------------------------------
P_THRESH  <- 5e-8     # SNP-exposure p-value cutoff
R2_THRESH <- 0.001    # MVMR-recommended strict clumping (default)
KB_THRESH <- 10000    # clumping window
F_THRESH  <- 5        # conditional F threshold (Sanderson 2021 lenient)
GENCOV    <- 0        # genetic covariance assumed zero (limitation)

# --- Exposure metadata ----------------------------------------------
# IDs verified against OpenGWAS server labels.
# Each entry has a `verify_label_contains` regex used by the validation
# block below to confirm the server returns the expected phenotype.
exposure_meta <- list(
  list(name = "FG",   opengwas_id = "ebi-a-GCST90002232",
       verify_label_contains = "[Ff]asting"),
  list(name = "SBP",  opengwas_id = "ieu-b-4818",
       verify_label_contains = "[Ss]ystolic"),
  list(name = "LDL",  opengwas_id = "ieu-b-4846",
       verify_label_contains = "LDL"),
  list(name = "BMI",  opengwas_id = "ieu-b-4816",
       verify_label_contains = "[Bb]ody mass"),
  list(name = "eGFR", opengwas_id = "ieu-b-4832",
       verify_label_contains = "glomerular")
)
exposure_ids   <- sapply(exposure_meta, `[[`, "opengwas_id")
exposure_names <- sapply(exposure_meta, `[[`, "name")
exposure_regex <- sapply(exposure_meta, `[[`, "verify_label_contains")

# --- Outcome metadata -----------------------------------------------
outcome_meta <- list(
  list(name = "AllStroke", opengwas_id = "finn-b-C_STROKE",     label = "Any stroke (FinnGen)"),
  list(name = "IS",        opengwas_id = "ebi-a-GCST90018864",  label = "Ischemic stroke"),
  list(name = "ICH",       opengwas_id = "ebi-a-GCST90018870",  label = "Intracerebral hemorrhage"),
  list(name = "SAH",       opengwas_id = "ebi-a-GCST90018923",  label = "Subarachnoid hemorrhage")
)

# ------------------------------------------------------------------
# 1. Extract & harmonise exposures (cached; reused across outcomes)
# ------------------------------------------------------------------
exp_cache <- file.path(proc_dir, "mv_exposure_dat.rds")
if (file.exists(exp_cache)) {
  exposure_dat <- readRDS(exp_cache)
  message(sprintf("Loaded cached exposure_dat: %d SNPs",
                  length(unique(exposure_dat$SNP))))
} else {
  message("Extracting MV exposures from OpenGWAS ...")
  exposure_dat <- mv_extract_exposures(
    id_exposure          = exposure_ids,
    pval_threshold       = P_THRESH,
    clump_r2             = R2_THRESH,
    clump_kb             = KB_THRESH,
    harmonise_strictness = 2,
    find_proxies         = FALSE,
    pop                  = "EUR"
  )
  saveRDS(exposure_dat, exp_cache)
  message(sprintf("  Done. %d unique SNPs across %d exposures.",
                  length(unique(exposure_dat$SNP)), length(exposure_ids)))
}

# ------------------------------------------------------------------
# 2. Helper: align column order of mvmr_dat to our exposure_meta order
#    AND verify ID->name mapping against OpenGWAS server labels.
#    This catches the v2 bug class permanently.
# ------------------------------------------------------------------
align_and_verify_mvmr <- function(mvmr_dat,
                                  exposure_ids,
                                  exposure_names,
                                  exposure_regex,
                                  test_snp_check = TRUE) {
  
  # Truth source: column order is given by colnames(exposure_beta), NOT by 
  # expname row order. The two may differ (mv_harmonise_data sorts colnames 
  # alphabetically while expname keeps a different internal order).
  cur_cols <- colnames(mvmr_dat$exposure_beta)
  if (is.null(cur_cols)) {
    stop("colnames(exposure_beta) is NULL - cannot determine column order.")
  }
  
  # ---- 2a. Verify all requested IDs are present in the matrix columns ----
  miss <- setdiff(exposure_ids, cur_cols)
  if (length(miss) > 0) {
    stop(sprintf("Exposure IDs missing in exposure_beta colnames: %s\n  Found: %s",
                 paste(miss, collapse = ", "),
                 paste(cur_cols, collapse = ", ")))
  }
  
  # ---- 2b. Verify each ID's expname label matches expected phenotype ----
  # Match server-returned labels to our regex; uses expname (independent of matrix).
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
  
  # ---- 2c. Reorder MATRIX columns based on colnames (the truth) ----
  col_ord <- match(exposure_ids, cur_cols)
  mvmr_dat$exposure_beta <- mvmr_dat$exposure_beta[, col_ord, drop = FALSE]
  mvmr_dat$exposure_se   <- mvmr_dat$exposure_se[,   col_ord, drop = FALSE]
  mvmr_dat$exposure_pval <- mvmr_dat$exposure_pval[, col_ord, drop = FALSE]
  
  # ---- 2d. Reorder expname INDEPENDENTLY (it has its own row order) ----
  exp_ord <- match(exposure_ids, cur_ids_in_expname)
  mvmr_dat$expname <- mvmr_dat$expname[exp_ord, , drop = FALSE]
  
  # ---- 2e. Set human-readable column names ----
  colnames(mvmr_dat$exposure_beta) <-
    colnames(mvmr_dat$exposure_se)   <-
    colnames(mvmr_dat$exposure_pval) <- exposure_names
  
  # ---- 2f. SNP-level numeric check against original exposure_dat ----
  # Pull a beta value from mvmr_dat and compare to the same SNP-exposure
  # combination in the original exposure_dat. They MUST match exactly.
  if (test_snp_check && exists("exposure_dat", envir = .GlobalEnv)) {
    ed <- get("exposure_dat", envir = .GlobalEnv)
    test_snps <- rownames(mvmr_dat$exposure_beta)[1:min(3, nrow(mvmr_dat$exposure_beta))]
    for (snp in test_snps) {
      for (i in seq_along(exposure_ids)) {
        expected <- ed$beta.exposure[ed$SNP == snp &
                                       ed$id.exposure == exposure_ids[i]]
        observed <- mvmr_dat$exposure_beta[snp, i]
        if (length(expected) == 0 || is.na(expected[1])) next
        if (abs(expected[1] - observed) > 1e-6) {
          stop(sprintf(
            "Numeric mismatch: SNP=%s, exposure=%s [%s]\n  expected=%.6f, observed=%.6f",
            snp, exposure_names[i], exposure_ids[i],
            expected[1], observed))
        }
      }
    }
  }
  
  # ---- 2g. Print confirmation table ----
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
# 3. Helper: full MVMR pipeline for one outcome
# ------------------------------------------------------------------
run_mvmr_pair <- function(exposure_dat, exposure_ids, exposure_names,
                          exposure_regex,
                          outcome_id, outcome_name, outcome_label) {
  message(sprintf("\n=== Outcome: %s (%s) ===", outcome_name, outcome_id))
  
  out_dat <- extract_outcome_data(exposure_dat$SNP, outcome_id, proxies = FALSE)
  if (is.null(out_dat) || nrow(out_dat) == 0) {
    warning(sprintf("No outcome data for %s", outcome_name)); return(NULL)
  }
  mvmr_dat <- mv_harmonise_data(exposure_dat, out_dat, harmonise_strictness = 2)
  mvmr_dat <- align_and_verify_mvmr(mvmr_dat, exposure_ids, exposure_names,
                                    exposure_regex, test_snp_check = TRUE)
  k_snp <- nrow(mvmr_dat$exposure_beta)
  message(sprintf("  SNPs after harmonisation: %d", k_snp))
  
  # 3a. MendelianRandomization input
  mr_input <- mr_mvinput(
    bx       = mvmr_dat$exposure_beta,
    bxse     = mvmr_dat$exposure_se,
    by       = mvmr_dat$outcome_beta,
    byse     = mvmr_dat$outcome_se,
    exposure = exposure_names,
    outcome  = outcome_name
  )
  
  # 3b. MVMR-IVW (primary)
  ivw <- mr_mvivw(mr_input, model = "default", correl = FALSE,
                  distribution = "normal", alpha = 0.05)
  ivw_tbl <- data.frame(
    outcome  = outcome_name,
    exposure = exposure_names,
    method   = "MVMR-IVW",
    nsnp     = k_snp,
    b        = ivw@Estimate,
    se       = ivw@StdError,
    pval     = ivw@Pvalue,
    OR       = exp(ivw@Estimate),
    OR_lci   = exp(ivw@CILower),
    OR_uci   = exp(ivw@CIUpper),
    stringsAsFactors = FALSE
  )
  
  # 3c. MVMR-Egger (sensitivity; intercept = pleiotropy test)
  egger <- mr_mvegger(mr_input, orientate = 1, correl = FALSE,
                      distribution = "normal", alpha = 0.05)
  egger_tbl <- data.frame(
    outcome  = outcome_name,
    exposure = exposure_names,
    method   = "MVMR-Egger",
    nsnp     = k_snp,
    b        = egger@Estimate,
    se       = egger@StdError.Est,
    pval     = egger@Pvalue.Est,
    OR       = exp(egger@Estimate),
    OR_lci   = exp(egger@CILower.Est),
    OR_uci   = exp(egger@CIUpper.Est),
    stringsAsFactors = FALSE
  )
  egger_intercept <- data.frame(
    outcome      = outcome_name,
    intercept    = egger@Intercept,
    se           = egger@StdError.Int,
    pval         = egger@Pvalue.Int,
    ci_lower     = egger@CILower.Int,
    ci_upper     = egger@CIUpper.Int,
    interpretation = ifelse(egger@Pvalue.Int < 0.05,
                            "Evidence of directional pleiotropy",
                            "No evidence of directional pleiotropy"),
    stringsAsFactors = FALSE
  )
  
  # 3d. MVMR diagnostics: conditional F + Q statistic
  F_data <- format_mvmr(
    BXGs   = mvmr_dat$exposure_beta,
    BYG    = mvmr_dat$outcome_beta,
    seBXGs = mvmr_dat$exposure_se,
    seBYG  = mvmr_dat$outcome_se,
    RSID   = rownames(mvmr_dat$exposure_beta)
  )
  
  cF_raw <- strength_mvmr(r_input = F_data, gencov = GENCOV)
  cF_tbl <- data.frame(
    outcome  = outcome_name,
    exposure = exposure_names,
    cond_F   = as.numeric(cF_raw),
    weak_F5  = as.numeric(cF_raw) < 5,
    weak_F10 = as.numeric(cF_raw) < 10,
    stringsAsFactors = FALSE
  )
  
  Q_raw <- pleiotropy_mvmr(r_input = F_data, gencov = GENCOV)
  Q_tbl <- data.frame(
    outcome = outcome_name,
    Q       = Q_raw$Qstat,
    df      = k_snp - length(exposure_names),
    Q_pval  = Q_raw$Qpval,
    interpretation = ifelse(Q_raw$Qpval < 0.05,
                            "Heterogeneity / horizontal pleiotropy detected",
                            "No detectable heterogeneity / pleiotropy"),
    stringsAsFactors = FALSE
  )
  
  # 3e. Save SNP-level data and harmonised matrices
  fwrite(F_data, file.path(proc_dir, sprintf("%s_F_data.csv", outcome_name)))
  saveRDS(mvmr_dat, file.path(proc_dir, sprintf("%s_mvmr_dat.rds", outcome_name)))
  
  list(ivw = ivw_tbl, egger = egger_tbl, intercept = egger_intercept,
       cF = cF_tbl, Q = Q_tbl, nsnp = k_snp)
}

# ------------------------------------------------------------------
# 4. Run for all outcomes
# ------------------------------------------------------------------
all_ivw <- list(); all_egger <- list(); all_intercept <- list()
all_cF  <- list(); all_Q     <- list()

for (om in outcome_meta) {
  res <- run_mvmr_pair(exposure_dat, exposure_ids, exposure_names,
                       exposure_regex,
                       om$opengwas_id, om$name, om$label)
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
# 5. FDR (within outcome, on IVW p-values across 5 exposures)
# ------------------------------------------------------------------
ivw_tbl <- ivw_tbl %>%
  mutate(CI_text = sprintf("%.3f (%.3f, %.3f)", OR, OR_lci, OR_uci)) %>%
  group_by(outcome) %>%
  mutate(p_fdr = p.adjust(pval, method = "BH")) %>%
  ungroup()

egger_tbl <- egger_tbl %>%
  mutate(CI_text = sprintf("%.3f (%.3f, %.3f)", OR, OR_lci, OR_uci))

# ------------------------------------------------------------------
# 6. Annotate IVW with cF and weak-instrument flag
# ------------------------------------------------------------------
ivw_tbl <- ivw_tbl %>%
  left_join(cF_tbl[, c("outcome", "exposure", "cond_F", "weak_F5", "weak_F10")],
            by = c("outcome", "exposure")) %>%
  mutate(
    flag = case_when(
      weak_F5  ~ "WEAK (cF<5)",
      weak_F10 ~ "Borderline (cF 5-10)",
      TRUE     ~ "OK (cF>=10)"
    )
  )

# ------------------------------------------------------------------
# 7. Save outputs
# ------------------------------------------------------------------
fwrite(ivw_tbl,       file.path(out_dir, "01_MVMR_IVW_EUR.csv"))
fwrite(egger_tbl,     file.path(out_dir, "02_MVMR_Egger_EUR.csv"))
fwrite(intercept_tbl, file.path(out_dir, "03_Egger_intercept_pleiotropy_EUR.csv"))
fwrite(cF_tbl,        file.path(out_dir, "04_conditional_F_EUR.csv"))
fwrite(Q_tbl,         file.path(out_dir, "05_Q_heterogeneity_EUR.csv"))

# Combined publication table
pub_tbl <- ivw_tbl %>%
  left_join(Q_tbl[, c("outcome", "Q", "Q_pval")], by = "outcome") %>%
  left_join(intercept_tbl[, c("outcome", "intercept", "pval")] %>%
              rename(egger_intercept = intercept, egger_intercept_pval = pval),
            by = "outcome") %>%
  select(outcome, exposure, nsnp, b, se, OR, OR_lci, OR_uci, CI_text,
         pval, p_fdr, cond_F, flag,
         Q, Q_pval, egger_intercept, egger_intercept_pval) %>%
  arrange(outcome, exposure)

fwrite(pub_tbl, file.path(out_dir, "00_MVMR_publication_table_EUR.csv"))

# ------------------------------------------------------------------
# 8. Forest plot per outcome (IVW vs Egger side-by-side)
# ------------------------------------------------------------------
plot_mvmr_forest <- function(outcome_name, ivw_sub, egger_sub) {
  combined <- bind_rows(
    ivw_sub   %>% mutate(method = "MVMR-IVW"),
    egger_sub %>% mutate(method = "MVMR-Egger")
  ) %>%
    mutate(exposure = factor(exposure, levels = rev(exposure_names)))
  
  p <- ggplot(combined, aes(x = OR, y = exposure, color = method)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
    geom_pointrange(aes(xmin = OR_lci, xmax = OR_uci),
                    position = position_dodge(width = 0.5), size = 0.4) +
    scale_x_log10() +
    scale_color_manual(values = c("MVMR-IVW" = "#0072B2", "MVMR-Egger" = "#D55E00")) +
    labs(title = sprintf("MVMR results for %s (EUR)", outcome_name),
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
# 9. Console summary
# ------------------------------------------------------------------
message("\n=== MVMR analysis complete ===")
message(sprintf("Results saved to: %s", out_dir))
cat("\n--- Conditional F (weak instrument flags) ---\n");  print(cF_tbl)
cat("\n--- Q heterogeneity / pleiotropy test ---\n");      print(Q_tbl)
cat("\n--- Egger intercept pleiotropy test ---\n");        print(intercept_tbl)
cat("\n--- IVW publication table ---\n");                  print(pub_tbl)
