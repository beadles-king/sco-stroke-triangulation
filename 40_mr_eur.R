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
# Two-Sample MR: European Ancestry Cohort
# Exposures: FG, LDL, SBP, BMI, eGFR
# Outcomes : Total Stroke (FinnGen), IS, ICH, SAH (UKB+FinnGen meta)
# Main method : IVW
# Sensitivity : MR-Egger, Weighted Median, MR-PRESSO,
#               Heterogeneity (Cochran Q), Pleiotropy (Egger intercept),
#               Leave-one-out
# Multiple testing : BH-FDR stratified by outcome (5 exposures per outcome)
# Instrument QC : p < 5e-8, r2 < 0.01, kb = 10000, F >= 10
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
work_dir <- file.path(RAWDAT, "MR/EUR")
setwd(work_dir)

out_dir       <- file.path(work_dir, "res")
plot_dir      <- file.path(out_dir, "plots")
snp_dir       <- file.path(out_dir, "snp_tables")
sens_dir      <- file.path(out_dir, "sensitivity")
dir.create(out_dir,  showWarnings = FALSE, recursive = TRUE)
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(snp_dir,  showWarnings = FALSE, recursive = TRUE)
dir.create(sens_dir, showWarnings = FALSE, recursive = TRUE)

# --- Exposure and outcome GWAS IDs ----------------------------------
exposure_ids <- list(
  FG   = "ebi-a-GCST90002232",  # MAGIC
  LDL  = "ieu-b-4846",          # Within Family GWAS
  SBP  = "ieu-b-4818",          # Within Family GWAS
  BMI  = "ieu-b-4816",          # Within Family GWAS
  eGFR = "ieu-b-4832"           # Within Family GWAS
)

outcome_ids <- list(
  stroke = "finn-b-C_STROKE",         # Total stroke (FinnGen)
  IS     = "ebi-a-GCST90018864",     # Ischaemic stroke (UKB + FinnGen meta)
  ICH    = "ebi-a-GCST90018870",     # Intracerebral hemorrhage
  SAH    = "ebi-a-GCST90018923"      # Subarachnoid hemorrhage
)

# Analysis-wide parameters
P_THRESH  <- 5e-8
R2_THRESH <- 0.01
KB_THRESH <- 10000
F_THRESH  <- 10

# ------------------------------------------------------------------
# 1. Helper functions
# ------------------------------------------------------------------

get_exposure <- function(exposure_name, exposure_id) {
  message(sprintf("\n[Exposure] %s  (%s)", exposure_name, exposure_id))
  inst <- extract_instruments(outcomes = exposure_id,
                              p1 = P_THRESH, p2 = P_THRESH,
                              r2 = R2_THRESH, kb = KB_THRESH)
  if (is.null(inst) || nrow(inst) == 0) {
    warning(sprintf("No instruments retrieved for %s", exposure_name))
    return(NULL)
  }
  inst$F_value <- (inst$beta.exposure / inst$se.exposure)^2
  
  if (!"eaf.exposure" %in% names(inst)) inst$eaf.exposure <- NA_real_
  inst$R2 <- 2 * inst$eaf.exposure * (1 - inst$eaf.exposure) *
    inst$beta.exposure^2
  
  n_before <- nrow(inst)
  inst_f   <- inst[inst$F_value >= F_THRESH, ]
  n_after  <- nrow(inst_f)
  message(sprintf("  SNPs: %d -> %d after F>=%d filter (min F = %.2f)",
                  n_before, n_after, F_THRESH,
                  min(inst_f$F_value, na.rm = TRUE)))
  inst_f$exposure_name <- exposure_name
  inst_f
}

run_mr_pair <- function(exposure_dat, exposure_name,
                        outcome_id,   outcome_name) {
  tag <- sprintf("%s_on_%s", exposure_name, outcome_name)
  message(sprintf("  -> %s", tag))
  
  out_dat <- extract_outcome_data(snps = exposure_dat$SNP,
                                  outcomes = outcome_id,
                                  proxies = FALSE)
  if (is.null(out_dat) || nrow(out_dat) == 0) {
    warning(sprintf("No outcome data for %s", tag))
    return(NULL)
  }
  
  dat <- harmonise_data(exposure_dat = exposure_dat,
                        outcome_dat  = out_dat,
                        action       = 2)
  
  fwrite(dat, file.path(snp_dir, sprintf("%s_snps.csv", tag)))
  
  n_keep <- sum(dat$mr_keep)
  if (n_keep < 3) {
    warning(sprintf("Only %d SNPs pass mr_keep for %s; skipping.", n_keep, tag))
    return(NULL)
  }
  
  # --- Main + sensitivity MR methods ------------------------
  mr_res <- mr(dat, method_list = c("mr_ivw",
                                    "mr_egger_regression",
                                    "mr_weighted_median"))
  mr_res$exposure_name <- exposure_name
  mr_res$outcome_name  <- outcome_name
  
  # --- Heterogeneity (Cochran Q) ----------------------------
  het <- tryCatch(mr_heterogeneity(dat), error = function(e) NULL)
  if (!is.null(het)) {
    het$exposure_name <- exposure_name
    het$outcome_name  <- outcome_name
  }
  
  # --- Horizontal pleiotropy (Egger intercept) --------------
  pleio <- tryCatch(mr_pleiotropy_test(dat), error = function(e) NULL)
  if (!is.null(pleio)) {
    pleio$exposure_name <- exposure_name
    pleio$outcome_name  <- outcome_name
  }
  
  # --- Leave-one-out ---------------------------------------
  loo <- tryCatch(mr_leaveoneout(dat), error = function(e) NULL)
  if (!is.null(loo)) {
    loo$exposure_name <- exposure_name
    loo$outcome_name  <- outcome_name
    fwrite(loo, file.path(sens_dir, sprintf("%s_leaveoneout.csv", tag)))
  }
  
  # --- Save plots (scatter, forest, funnel, LOO) -----------
  tryCatch({
    p_scatter <- mr_scatter_plot(mr_res, dat)[[1]]
    ggsave(file.path(plot_dir, sprintf("%s_scatter.pdf", tag)),
           p_scatter, width = 6, height = 5)
  }, error = function(e) message("    scatter plot failed: ", e$message))
  
  tryCatch({
    ss        <- mr_singlesnp(dat)
    p_forest  <- mr_forest_plot(ss)[[1]]
    ggsave(file.path(plot_dir, sprintf("%s_forest.pdf", tag)),
           p_forest, width = 6, height = max(4, 0.18 * nrow(ss)))
    p_funnel  <- mr_funnel_plot(ss)[[1]]
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
# 2. Retrieve all exposure instruments
# ------------------------------------------------------------------
exposure_list <- lapply(names(exposure_ids),
                        function(x) get_exposure(x, exposure_ids[[x]]))
names(exposure_list) <- names(exposure_ids)

instrument_summary <- do.call(rbind, lapply(names(exposure_list), function(nm) {
  d <- exposure_list[[nm]]
  if (is.null(d)) return(NULL)
  data.frame(
    exposure   = nm,
    gwas_id    = exposure_ids[[nm]],
    n_snp      = nrow(d),
    mean_F     = round(mean(d$F_value, na.rm = TRUE), 2),
    min_F      = round(min(d$F_value,  na.rm = TRUE), 2),
    max_F      = round(max(d$F_value,  na.rm = TRUE), 2),
    total_R2   = round(sum(d$R2, na.rm = TRUE), 4),
    stringsAsFactors = FALSE
  )
}))
fwrite(instrument_summary, file.path(out_dir, "00_instrument_summary_EUR.csv"))
print(instrument_summary)

# ------------------------------------------------------------------
# 3. Run MR for every exposure x outcome combination
# ------------------------------------------------------------------
all_mr    <- list()
all_het   <- list()
all_pleio <- list()

for (ex_name in names(exposure_list)) {
  ex_dat <- exposure_list[[ex_name]]
  if (is.null(ex_dat)) next
  message(sprintf("\n=== Exposure: %s ===", ex_name))
  for (out_name in names(outcome_ids)) {
    res <- run_mr_pair(ex_dat, ex_name,
                       outcome_ids[[out_name]], out_name)
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
# 4. Format main results + FDR (stratified by outcome)
# ------------------------------------------------------------------
mr_tbl <- mr_tbl %>%
  mutate(
    OR       = exp(b),
    OR_lci   = exp(b - 1.96 * se),
    OR_uci   = exp(b + 1.96 * se),
    CI_text  = sprintf("%.3f (%.3f, %.3f)", OR, OR_lci, OR_uci)
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
# 5. Publication-ready primary table (IVW only)
# ------------------------------------------------------------------
primary_tbl <- mr_tbl %>%
  filter(method == "Inverse variance weighted") %>%
  select(exposure_name, outcome_name, nsnp,
         b, se, OR, OR_lci, OR_uci, CI_text, pval, p_fdr) %>%
  arrange(outcome_name, exposure_name)

fwrite(primary_tbl, file.path(out_dir, "01_primary_IVW_results_EUR.csv"))
fwrite(mr_tbl,      file.path(out_dir, "02_all_MR_methods_EUR.csv"))
fwrite(het_tbl,     file.path(out_dir, "03_heterogeneity_EUR.csv"))
fwrite(pleio_tbl,   file.path(out_dir, "04_pleiotropy_egger_intercept_EUR.csv"))

message("\n=== Analysis complete ===")
message(sprintf("Results saved to: %s", out_dir))
print(primary_tbl)
