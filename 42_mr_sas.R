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
# Two-Sample MR: South Asian Ancestry Cohort  (v2)
# Exposures :
#   - FG, LDL, SBP, BMI, eGFR : pan-UKBB (CSA = Central/South Asian)
#                               Downloaded via extract_instruments(); already has rsID
# Outcomes  : GIGASTROKE South Asian (Mishra et al. 2022)
#             stroke (any), IS (ischemic)
#             NB: outcome files have NO rsID -> match by chr:pos
# Main method  : IVW
# Sensitivity  : MR-Egger, Weighted Median,
#                Heterogeneity (Cochran Q), Pleiotropy (Egger intercept),
#                Leave-one-out
# Multiple testing : BH-FDR stratified by outcome
# Instrument QC    : p < 1e-5 (small SAS samples), r2 < 0.001, kb = 10000, F >= 10
# Sample overlap   : pan-UKBB-CSA exposures vs GIGASTROKE SAS outcomes
#                    GIGASTROKE-SAS includes some UKBB-CSA participants -> 
#                    potential modest overlap; flagged in limitations.
# =============================================================================

# ------------------------------------------------------------------
# 0. Setup
# ------------------------------------------------------------------
# OpenGWAS authentication: set your own token before running, e.g.
#   Sys.setenv(OPENGWAS_JWT = "<your token>")
# Tokens are not distributed with this repository.

library(TwoSampleMR)
library(data.table)
library(dplyr)
library(ggplot2)

# --- Working directory & output structure ---------------------------
work_dir <- file.path(GWAS_DIR, "SAS")
setwd(work_dir)

out_dir   <- file.path(RAWDAT, "MR/SAS/res")
plot_dir  <- file.path(out_dir, "plots")
snp_dir   <- file.path(out_dir, "snp_tables")
sens_dir  <- file.path(out_dir, "sensitivity")
proc_dir  <- file.path(out_dir, "processed")     # cached merged exposure-outcome
invisible(lapply(c(out_dir, plot_dir, snp_dir, sens_dir, proc_dir),
                 dir.create, showWarnings = FALSE, recursive = TRUE))

# --- Constants ------------------------------------------------------
P_THRESH  <- 1e-5      # SAS samples are small, relax from 5e-8
R2_THRESH <- 0.01
KB_THRESH <- 10000
F_THRESH  <- 10
POP       <- "SAS"

# --- Exposure metadata (pan-UKBB CSA via OpenGWAS) ------------------
# All have rsID directly from extract_instruments(); no chr:pos mapping needed
exposure_meta <- list(
  FG   = list(name = "FG",   opengwas_id = "ebi-a-GCST90002234", source = "pan-UKBB-CSA"),
  LDL  = list(name = "LDL",  opengwas_id = "ukb-e-30780_CSA",     source = "pan-UKBB-CSA"),
  SBP  = list(name = "SBP",  opengwas_id = "ukb-e-4080_CSA",   source = "pan-UKBB-CSA"),
  BMI  = list(name = "BMI",  opengwas_id = "ukb-e-23104_CSA",    source = "pan-UKBB-CSA"),
  eGFR = list(name = "eGFR", opengwas_id = "ukb-e-recode6_CSA",  source = "pan-UKBB-CSA")
)

# --- Outcome metadata (GIGASTROKE SAS) ------------------------------
# IMPORTANT: these files have NO rsID column; matched by chr:pos
outcome_meta <- list(
  stroke = list(name = "stroke",
                file = file.path(RAWDAT, "GWAS_C/GIGA_stroke_SAS.tsv.gz"),
                source = "GIGASTROKE-SAS"),
  IS     = list(name = "IS",
                file = file.path(RAWDAT, "GWAS_C/GIGA_IS_SAS.tsv.gz"),
                source = "GIGASTROKE-SAS")
)

# Outcome columns (GIGASTROKE GWAS Catalog standard)
OUT_COLS <- list(
  chr  = "chromosome",       pos  = "base_pair_location",
  beta = "beta",             se   = "standard_error",
  ea   = "effect_allele",    oa   = "other_allele",
  eaf  = "effect_allele_frequency", pval = "p_value"
)

# ------------------------------------------------------------------
# 1. Helper: prepare exposure from OpenGWAS (pan-UKBB CSA)
# ------------------------------------------------------------------
prepare_exposure_opengwas <- function(meta) {
  message(sprintf("\n[Exposure] %s (%s, OpenGWAS id = %s)",
                  meta$name, meta$source, meta$opengwas_id))
  
  # Cache: avoid re-downloading each run
  cache_path <- file.path(proc_dir, sprintf("%s_exposure_raw.rds", tolower(meta$name)))
  if (file.exists(cache_path)) {
    exp_dat <- readRDS(cache_path)
    message(sprintf("  Loaded cached download: %d SNPs", nrow(exp_dat)))
  } else {
    exp_dat <- tryCatch(
      extract_instruments(outcomes = meta$opengwas_id,
                          p1 = P_THRESH, clump = FALSE),
      error = function(e) {
        message("  extract_instruments failed: ", e$message); NULL
      })
    if (is.null(exp_dat) || nrow(exp_dat) == 0) {
      warning(sprintf("No instruments returned for %s", meta$name)); return(NULL)
    }
    saveRDS(exp_dat, cache_path)
    message(sprintf("  Downloaded & cached: %d SNPs at p<%.0e", nrow(exp_dat), P_THRESH))
  }
  
  exp_dat$exposure <- meta$name
  
  # Clump in SAS reference panel
  clumped <- tryCatch(
    clump_data(exp_dat, clump_kb = KB_THRESH, clump_r2 = R2_THRESH,
               clump_p1 = 1, clump_p2 = 1, pop = POP),
    error = function(e) {
      message("  clump_data failed: ", e$message); NULL
    })
  if (is.null(clumped) || nrow(clumped) == 0) return(NULL)
  message(sprintf("  After clumping (r2<%.3f, kb=%d, %s): %d SNPs",
                  R2_THRESH, KB_THRESH, POP, nrow(clumped)))
  
  # F-statistic and F>=10 filter
  clumped$F_value <- (clumped$beta.exposure / clumped$se.exposure)^2
  n_pre_F <- nrow(clumped)
  clumped <- clumped[clumped$F_value >= F_THRESH, ]
  message(sprintf("  After F>=%d: %d SNPs (was %d; min F = %.2f)",
                  F_THRESH, nrow(clumped), n_pre_F, min(clumped$F_value)))
  
  # Build chr:pos key for outcome matching (outcomes lack rsID)
  clumped$id <- paste0(clumped$chr.exposure, ":", clumped$pos.exposure)
  
  # R^2 (proportion of variance explained) - useful for power calc
  if (!is.null(clumped$samplesize.exposure) &&
      !all(is.na(clumped$samplesize.exposure))) {
    n_eff <- clumped$samplesize.exposure[1]
    clumped$r2 <- clumped$F_value / (clumped$F_value + n_eff - 2)
    message(sprintf("  Total R^2 (sum): %.4f  (N = %d)", sum(clumped$r2), n_eff))
  }
  
  clumped$exposure_name <- meta$name
  clumped$source        <- meta$source
  clumped
}

# ------------------------------------------------------------------
# 2. Helper: extract outcome SNPs by chr:pos from in-memory data.table
#    (GIGASTROKE-SAS files have no rsID, so we match positionally and
#     then graft the exposure rsIDs back onto the outcome rows.)
# ------------------------------------------------------------------
extract_outcome_local_bypos <- function(exposure_dat, outcome_dt, tag) {
  # Build chr:pos in outcome (already done at load time, but defensive)
  if (!"id" %in% colnames(outcome_dt)) {
    outcome_dt[, id := paste0(get(OUT_COLS$chr), ":", get(OUT_COLS$pos))]
  }
  
  # Subset outcome to exposure positions
  out_sub <- outcome_dt[id %in% exposure_dat$id]
  if (nrow(out_sub) == 0) {
    warning(sprintf("No chr:pos overlap in outcome for %s", tag))
    return(NULL)
  }
  
  # Graft rsID from exposure (id <-> SNP map)
  rs_map <- unique(exposure_dat[, c("id", "SNP")])
  out_sub <- merge(out_sub, rs_map, by = "id", all.x = TRUE)
  
  # Some chr:pos may map to >1 rsID upstream; keep first occurrence
  out_sub <- out_sub[!duplicated(out_sub$SNP), ]
  
  csv_path <- file.path(proc_dir, sprintf("%s_outcome.csv", tag))
  fwrite(out_sub, csv_path)
  
  message(sprintf("    Outcome SNPs matched by chr:pos: %d (of %d exposure SNPs)",
                  nrow(out_sub), nrow(exposure_dat)))
  
  read_outcome_data(
    snps = exposure_dat$SNP, filename = csv_path, sep = ",",
    snp_col = "SNP",
    beta_col = OUT_COLS$beta, se_col = OUT_COLS$se,
    effect_allele_col = OUT_COLS$ea, other_allele_col = OUT_COLS$oa,
    eaf_col = OUT_COLS$eaf, pval_col = OUT_COLS$pval
  )
}

# ------------------------------------------------------------------
# 3. Helper: full MR pipeline for one exposure-outcome pair
# ------------------------------------------------------------------
run_mr_pair <- function(exposure_dat, exposure_name,
                        outcome_dt, outcome_name) {
  tag <- sprintf("%s_on_%s", exposure_name, outcome_name)
  message(sprintf("  -> %s", tag))
  
  out_dat <- extract_outcome_local_bypos(exposure_dat, outcome_dt, tag)
  if (is.null(out_dat) || nrow(out_dat) == 0) {
    warning(sprintf("No outcome data for %s", tag)); return(NULL)
  }
  out_dat$id.outcome <- sprintf("%s_SAS", outcome_name)
  out_dat$outcome    <- sprintf("%s (South Asian)", outcome_name)
  
  dat <- harmonise_data(exposure_dat = exposure_dat,
                        outcome_dat  = out_dat, action = 2)
  fwrite(dat, file.path(snp_dir, sprintf("%s_snps.csv", tag)))
  
  n_keep <- sum(dat$mr_keep)
  if (n_keep < 3) {
    warning(sprintf("Only %d SNPs pass mr_keep for %s; skipping.", n_keep, tag))
    return(NULL)
  }
  
  # Main + sensitivity MR methods
  mr_res <- mr(dat, method_list = c("mr_ivw",
                                    "mr_egger_regression",
                                    "mr_weighted_median"))
  mr_res$exposure_name <- exposure_name
  mr_res$outcome_name  <- outcome_name
  
  # Heterogeneity
  het <- tryCatch(mr_heterogeneity(dat), error = function(e) NULL)
  if (!is.null(het)) {
    het$exposure_name <- exposure_name; het$outcome_name <- outcome_name
  }
  
  # Egger pleiotropy
  pleio <- tryCatch(mr_pleiotropy_test(dat), error = function(e) NULL)
  if (!is.null(pleio)) {
    pleio$exposure_name <- exposure_name; pleio$outcome_name <- outcome_name
  }
  
  # Leave-one-out
  loo <- tryCatch(mr_leaveoneout(dat), error = function(e) NULL)
  if (!is.null(loo)) {
    loo$exposure_name <- exposure_name; loo$outcome_name <- outcome_name
    fwrite(loo, file.path(sens_dir, sprintf("%s_leaveoneout.csv", tag)))
  }
  
  # Plots
  tryCatch({
    p_scatter <- mr_scatter_plot(mr_res, dat)[[1]]
    ggsave(file.path(plot_dir, sprintf("%s_scatter.pdf", tag)),
           p_scatter, width = 6, height = 5)
  }, error = function(e) message("    scatter failed: ", e$message))
  
  tryCatch({
    ss <- mr_singlesnp(dat)
    p_forest <- mr_forest_plot(ss)[[1]]
    ggsave(file.path(plot_dir, sprintf("%s_forest.pdf", tag)),
           p_forest, width = 6, height = max(4, 0.18 * nrow(ss)))
    p_funnel <- mr_funnel_plot(ss)[[1]]
    ggsave(file.path(plot_dir, sprintf("%s_funnel.pdf", tag)),
           p_funnel, width = 6, height = 5)
  }, error = function(e) message("    forest/funnel failed: ", e$message))
  
  tryCatch({
    if (!is.null(loo)) {
      p_loo <- mr_leaveoneout_plot(loo)[[1]]
      ggsave(file.path(plot_dir, sprintf("%s_leaveoneout.pdf", tag)),
             p_loo, width = 6, height = max(4, 0.18 * nrow(loo)))
    }
  }, error = function(e) message("    LOO plot failed: ", e$message))
  
  list(mr = mr_res, het = het, pleio = pleio)
}

# ------------------------------------------------------------------
# 4. Load all outcomes once into memory
# ------------------------------------------------------------------
message("=== Loading outcomes (GIGASTROKE SAS) ===")
outcomes <- lapply(outcome_meta, function(m) {
  message(sprintf("  Loading %s from %s ...", m$name, m$file))
  dt <- fread(m$file, sep = "auto", header = TRUE,
              nThread = getDTthreads(), na.strings = c("", "NA"),
              showProgress = TRUE)
  # Pre-build chr:pos key (outcomes have no rsID)
  dt[, id := paste0(get(OUT_COLS$chr), ":", get(OUT_COLS$pos))]
  dt
})

# Validate outcome columns
required_cols <- unlist(OUT_COLS)
for (nm in names(outcomes)) {
  missing <- setdiff(required_cols, colnames(outcomes[[nm]]))
  if (length(missing) > 0) {
    stop(sprintf(
      "Outcome '%s' missing columns: %s\nFound: %s\nUpdate OUT_COLS to match.",
      nm, paste(missing, collapse = ", "),
      paste(colnames(outcomes[[nm]]), collapse = ", ")))
  }
  message(sprintf("  %s: %d variants loaded", nm, nrow(outcomes[[nm]])))
}

# ------------------------------------------------------------------
# 5. Process all exposures
# ------------------------------------------------------------------
message("\n=== Processing exposures ===")
exposures <- list()
for (nm in names(exposure_meta)) {
  exposures[[nm]] <- prepare_exposure_opengwas(exposure_meta[[nm]])
}

# Instrument summary
instrument_summary <- do.call(rbind, lapply(names(exposures), function(nm) {
  d <- exposures[[nm]]; m <- exposure_meta[[nm]]
  if (is.null(d)) return(NULL)
  data.frame(
    exposure  = nm,
    source    = m$source,
    opengwas  = m$opengwas_id,
    n_snp     = nrow(d),
    mean_F    = round(mean(d$F_value), 2),
    min_F     = round(min(d$F_value),  2),
    max_F     = round(max(d$F_value),  2),
    total_R2  = if (!is.null(d$r2)) round(sum(d$r2), 4) else NA_real_,
    stringsAsFactors = FALSE
  )
}))
fwrite(instrument_summary, file.path(out_dir, "00_instrument_summary_SAS.csv"))
print(instrument_summary)

# ------------------------------------------------------------------
# 6. Run MR for every exposure x outcome combination
# ------------------------------------------------------------------
all_mr <- list(); all_het <- list(); all_pleio <- list()

for (ex_name in names(exposures)) {
  ex_dat <- exposures[[ex_name]]
  if (is.null(ex_dat)) next
  message(sprintf("\n=== Exposure: %s ===", ex_name))
  for (out_name in names(outcomes)) {
    res <- run_mr_pair(ex_dat, ex_name, outcomes[[out_name]], out_name)
    if (is.null(res)) next
    all_mr[[length(all_mr)+1]] <- res$mr
    if (!is.null(res$het))   all_het[[length(all_het)+1]]     <- res$het
    if (!is.null(res$pleio)) all_pleio[[length(all_pleio)+1]] <- res$pleio
  }
}

mr_tbl    <- do.call(rbind, all_mr)
het_tbl   <- do.call(rbind, all_het)
pleio_tbl <- do.call(rbind, all_pleio)

# ------------------------------------------------------------------
# 7. Format + FDR (stratified by outcome)
# ------------------------------------------------------------------
mr_tbl <- mr_tbl %>%
  mutate(
    OR      = exp(b),
    OR_lci  = exp(b - 1.96 * se),
    OR_uci  = exp(b + 1.96 * se),
    CI_text = sprintf("%.3f (%.3f, %.3f)", OR, OR_lci, OR_uci)
  ) %>%
  group_by(outcome_name) %>%
  mutate(p_fdr = {
    q <- rep(NA_real_, n())
    i <- method == "Inverse variance weighted"
    q[i] <- p.adjust(pval[i], method = "BH")
    q
  }) %>%
  ungroup()

# ------------------------------------------------------------------
# 8. Publication-ready primary table (IVW only)
# ------------------------------------------------------------------
primary_tbl <- mr_tbl %>%
  filter(method == "Inverse variance weighted") %>%
  select(exposure_name, outcome_name, nsnp,
         b, se, OR, OR_lci, OR_uci, CI_text, pval, p_fdr) %>%
  arrange(outcome_name, exposure_name)

fwrite(primary_tbl, file.path(out_dir, "01_primary_IVW_results_SAS.csv"))
fwrite(mr_tbl,      file.path(out_dir, "02_all_MR_methods_SAS.csv"))
fwrite(het_tbl,     file.path(out_dir, "03_heterogeneity_SAS.csv"))
fwrite(pleio_tbl,   file.path(out_dir, "04_pleiotropy_egger_intercept_SAS.csv"))

message("\n=== Analysis complete ===")
message(sprintf("Results saved to: %s", out_dir))
print(primary_tbl)
