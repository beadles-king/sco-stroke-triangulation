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
# Two-Sample MR: East Asian Ancestry Cohort  (v2)
# Exposures :
#   - FG, LDL, SBP, BMI : TPMI (Taiwan) - SAIGE format, needs SNPnexus rsID mapping
#   - eGFR              : BBJ (Japan)   - BOLT-LMM format, has rsID
# Outcomes  : CKB pheweb (https://pheweb.ckbiobank.org)
#             stroke (i61/total), IS (i63), ICH (i61), SAH (i60)
# Main method  : IVW
# Sensitivity  : MR-Egger, Weighted Median,
#                Heterogeneity (Cochran Q), Pleiotropy (Egger intercept),
#                Leave-one-out
# Multiple testing : BH-FDR stratified by outcome
# Instrument QC    : p < 5e-8, r2 < 0.01, kb = 10000, F >= 10
# Sample overlap   : NONE (TPMI Taiwanese / BBJ Japanese vs CKB Chinese mainland)
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
work_dir <- file.path(GWAS_DIR, "AS")
setwd(work_dir)

out_dir   <- file.path(RAWDAT, "MR/EAS/res")
plot_dir  <- file.path(out_dir, "plots")
snp_dir   <- file.path(out_dir, "snp_tables")
sens_dir  <- file.path(out_dir, "sensitivity")
proc_dir  <- file.path(out_dir, "processed")     # cached processed exposures/outcomes
invisible(lapply(c(out_dir, plot_dir, snp_dir, sens_dir, proc_dir),
                 dir.create, showWarnings = FALSE, recursive = TRUE))

# --- Constants ------------------------------------------------------
P_THRESH  <- 5e-8
R2_THRESH <- 0.01
KB_THRESH <- 10000
F_THRESH  <- 10
POP       <- "EAS"

# --- Exposure metadata ----------------------------------------------
# format = "saige" : CHR, POS, Allele1, Allele2, AF_Allele2, BETA, SE, p.value, QC
#                    Needs SNPnexus rsID conversion (chr:pos -> rs)
# format = "bolt"  : SNP (rsID), ALLELE1, ALLELE0, A1FREQ, BETA, SE, P_BOLT_LMM_INF
exposure_meta <- list(
  FG = list(
    name = "FG", source = "TPMI", format = "saige",
    file = "FG_AS.txt",
    id_file  = "FG_id_AS.txt",   # SNPnexus input (auto-generated if missing)
    snp_file = "FG_snp_AS.txt"   # SNPnexus output (rsID lookup)
  ),
  LDL = list(
    name = "LDL", source = "TPMI", format = "saige",
    file = "LDL_AS.txt",
    id_file = "LDL_id_AS.txt", snp_file = "LDL_snp_AS.txt"
  ),
  SBP = list(
    name = "SBP", source = "TPMI", format = "saige",
    file = "SBP_AS.txt",
    id_file = "SBP_id_AS.txt", snp_file = "SBP_snp_AS.txt"
  ),
  BMI = list(
    name = "BMI", source = "TPMI", format = "saige",
    file = "BMI_AS.txt",
    id_file = "BMI_id_AS.txt", snp_file = "BMI_snp_AS.txt"
  ),
  eGFR = list(
    name = "eGFR", source = "BBJ", format = "bolt",
    file = "eGFR_AS.txt.gz"
  )
)

# --- Outcome metadata (CKB pheweb) ----------------------------------
outcome_meta <- list(
  stroke = list(name = "stroke", file = "Stroke.tsv", source = "CKB"),
  IS     = list(name = "IS",     file = "isc.tsv",    source = "CKB"),
  ICH    = list(name = "ICH",    file = "ich.tsv",    source = "CKB"),
  SAH    = list(name = "SAH",    file = "sah.tsv",    source = "CKB-pheweb")
)

# Outcome columns (CKB pheweb / GWAS Catalog standard)
OUT_COLS <- list(
  snp  = "variant_id",  beta = "beta",  se = "standard_error",
  ea   = "effect_allele", oa = "other_allele",
  eaf  = "effect_allele_frequency", pval = "p_value"
)

# ------------------------------------------------------------------
# 1. Helper: prepare SAIGE-format exposure (TPMI: FG/LDL/SBP/BMI)
# ------------------------------------------------------------------
prepare_exposure_saige <- function(meta) {
  message(sprintf("\n[Exposure] %s (%s, SAIGE format)", meta$name, meta$source))
  
  # Read raw GWAS, filter by p and QC
  dat <- fread(meta$file, sep = "auto", header = TRUE,
               nThread = getDTthreads(), na.strings = c("", "NA"))
  dat$p.value <- as.numeric(dat$p.value)   # in case stored as string
  dat <- dat[p.value < P_THRESH & QC == "PASS"]
  message(sprintf("  After p<%.0e & QC==PASS: %d SNPs", P_THRESH, nrow(dat)))
  if (nrow(dat) == 0) return(NULL)
  
  # SNPnexus rsID mapping: write input if not exists
  if (!file.exists(meta$snp_file)) {
    id_input <- dat %>%
      mutate(Type = "chromosome",
             Name = gsub("chr", "", CHR),
             Position = POS,
             Strand = 1) %>%
      select(Type, Name, Position, Allele1, Allele2, Strand)
    write.table(id_input, meta$id_file, sep = " ",
                row.names = FALSE, col.names = FALSE, quote = FALSE)
    stop(sprintf(
      "SNPnexus input written to '%s'.\n  STEPS:\n  1. Submit it to https://www.snp-nexus.org\n  2. Save the rsID output as '%s'\n  3. Re-run this script.",
      meta$id_file, meta$snp_file))
  }
  
  # Read SNPnexus output (V1 = chr:pos:A/A:1, V2 = rsID)
  snp_map <- read.table(meta$snp_file, header = FALSE, stringsAsFactors = FALSE)
  message(sprintf("  SNPnexus rsID map: %d entries", nrow(snp_map)))
  
  # Build MarkerID and join
  dat[, MarkerID := paste0(CHR, ":", POS, ":", Allele1, "/", Allele2, ":1")]
  joined <- inner_join(dat, snp_map, by = c("MarkerID" = "V1"))
  
  # Deduplicate: one chr:pos can map to multiple rsIDs (or vice versa)
  n_before <- nrow(joined)
  joined <- joined %>% distinct(V2, .keep_all = TRUE)
  if (nrow(joined) < n_before) {
    message(sprintf("  Removed %d duplicate rsID mappings", n_before - nrow(joined)))
  }
  message(sprintf("  After rsID mapping & dedup: %d SNPs", nrow(joined)))
  
  # Save to disk and read in TwoSampleMR exposure format
  csv_path <- file.path(proc_dir, sprintf("%s_cor.csv", tolower(meta$name)))
  fwrite(joined, csv_path)
  
  exp_dat <- read_exposure_data(
    filename = csv_path, sep = ",",
    snp_col = "V2", beta_col = "BETA", se_col = "SE",
    effect_allele_col = "Allele2", other_allele_col = "Allele1",
    eaf_col = "AF_Allele2", pval_col = "p.value"
  )
  exp_dat$exposure <- meta$name
  
  # Clump (EAS reference panel)
  clumped <- clump_data(exp_dat,
                        clump_kb = KB_THRESH, clump_r2 = R2_THRESH,
                        clump_p1 = 1, clump_p2 = 1, pop = POP)
  message(sprintf("  After clumping: %d SNPs", nrow(clumped)))
  
  # F-value and F>=10 filter
  clumped$F_value <- (clumped$beta.exposure / clumped$se.exposure)^2
  n_pre_F <- nrow(clumped)
  clumped <- clumped[clumped$F_value >= F_THRESH, ]
  message(sprintf("  After F>=%d: %d SNPs (was %d; min F = %.2f)",
                  F_THRESH, nrow(clumped), n_pre_F, min(clumped$F_value)))
  
  clumped$exposure_name <- meta$name
  clumped$source <- meta$source
  clumped
}

# ------------------------------------------------------------------
# 2. Helper: prepare BOLT-LMM exposure (BBJ: eGFR)
# ------------------------------------------------------------------
prepare_exposure_bolt <- function(meta) {
  message(sprintf("\n[Exposure] %s (%s, BOLT-LMM format)", meta$name, meta$source))
  
  dat <- fread(meta$file, sep = "auto", header = TRUE,
               nThread = getDTthreads(), na.strings = c("", "NA"))
  dat[[ "P_BOLT_LMM_INF" ]] <- as.numeric(dat[[ "P_BOLT_LMM_INF" ]])
  dat <- dat[P_BOLT_LMM_INF < P_THRESH]
  dat <- dat[grepl("^rs", SNP)]
  message(sprintf("  After p<%.0e & rsID filter: %d SNPs", P_THRESH, nrow(dat)))
  
  csv_path <- file.path(proc_dir, sprintf("%s_cor.csv", tolower(meta$name)))
  fwrite(dat, csv_path)
  
  exp_dat <- read_exposure_data(
    filename = csv_path, sep = ",",
    snp_col = "SNP", beta_col = "BETA", se_col = "SE",
    effect_allele_col = "ALLELE1", other_allele_col = "ALLELE0",
    eaf_col = "A1FREQ", pval_col = "P_BOLT_LMM_INF"
  )
  exp_dat$exposure <- meta$name
  
  clumped <- clump_data(exp_dat, clump_kb = KB_THRESH, clump_r2 = R2_THRESH,
                        clump_p1 = 1, clump_p2 = 1, pop = POP)
  message(sprintf("  After clumping: %d SNPs", nrow(clumped)))
  
  clumped$F_value <- (clumped$beta.exposure / clumped$se.exposure)^2
  n_pre_F <- nrow(clumped)
  clumped <- clumped[clumped$F_value >= F_THRESH, ]
  message(sprintf("  After F>=%d: %d SNPs (was %d; min F = %.2f)",
                  F_THRESH, nrow(clumped), n_pre_F, min(clumped$F_value)))
  
  clumped$exposure_name <- meta$name
  clumped$source <- meta$source
  clumped
}

# ------------------------------------------------------------------
# 3. Helper: extract outcome SNPs from in-memory data.table
# ------------------------------------------------------------------
extract_outcome_local <- function(exposure_snps, outcome_dt, tag) {
  out_sub <- outcome_dt[get(OUT_COLS$snp) %in% exposure_snps]
  if (nrow(out_sub) == 0) return(NULL)
  csv_path <- file.path(proc_dir, sprintf("%s_outcome.csv", tag))
  fwrite(out_sub, csv_path)
  read_outcome_data(
    snps = exposure_snps, filename = csv_path, sep = ",",
    snp_col = OUT_COLS$snp,
    beta_col = OUT_COLS$beta, se_col = OUT_COLS$se,
    effect_allele_col = OUT_COLS$ea, other_allele_col = OUT_COLS$oa,
    eaf_col = OUT_COLS$eaf, pval_col = OUT_COLS$pval
  )
}

# ------------------------------------------------------------------
# 4. Helper: full MR pipeline for one exposure-outcome pair
# ------------------------------------------------------------------
run_mr_pair <- function(exposure_dat, exposure_name,
                        outcome_dt, outcome_name) {
  tag <- sprintf("%s_on_%s", exposure_name, outcome_name)
  message(sprintf("  -> %s", tag))
  
  out_dat <- extract_outcome_local(exposure_dat$SNP, outcome_dt, tag)
  if (is.null(out_dat) || nrow(out_dat) == 0) {
    warning(sprintf("No outcome data for %s", tag)); return(NULL)
  }
  out_dat$outcome <- outcome_name
  
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
# 5. Load all outcomes once into memory
# ------------------------------------------------------------------
message("=== Loading outcomes (CKB pheweb) ===")
outcomes <- lapply(outcome_meta, function(m) {
  message(sprintf("  Loading %s from %s ...", m$name, m$file))
  fread(m$file, sep = "auto", header = TRUE,
        nThread = getDTthreads(), na.strings = c("", "NA"))
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
}
message("  All outcome files have expected columns.")

# ------------------------------------------------------------------
# 6. Process all exposures
# ------------------------------------------------------------------
message("\n=== Processing exposures ===")
exposures <- list()
for (nm in names(exposure_meta)) {
  m <- exposure_meta[[nm]]
  exposures[[nm]] <- if (m$format == "saige") prepare_exposure_saige(m)
  else                      prepare_exposure_bolt(m)
}

# Instrument summary
instrument_summary <- do.call(rbind, lapply(names(exposures), function(nm) {
  d <- exposures[[nm]]; m <- exposure_meta[[nm]]
  if (is.null(d)) return(NULL)
  data.frame(
    exposure = nm,
    source   = m$source,
    n_snp    = nrow(d),
    mean_F   = round(mean(d$F_value), 2),
    min_F    = round(min(d$F_value),  2),
    max_F    = round(max(d$F_value),  2),
    stringsAsFactors = FALSE
  )
}))
fwrite(instrument_summary, file.path(out_dir, "00_instrument_summary_EAS.csv"))
print(instrument_summary)

# ------------------------------------------------------------------
# 7. Run MR for every exposure x outcome combination
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
# 8. Format + FDR (stratified by outcome)
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
# 9. Publication-ready primary table (IVW only)
# ------------------------------------------------------------------
primary_tbl <- mr_tbl %>%
  filter(method == "Inverse variance weighted") %>%
  select(exposure_name, outcome_name, nsnp,
         b, se, OR, OR_lci, OR_uci, CI_text, pval, p_fdr) %>%
  arrange(outcome_name, exposure_name)

fwrite(primary_tbl, file.path(out_dir, "01_primary_IVW_results_EAS.csv"))
fwrite(mr_tbl,      file.path(out_dir, "02_all_MR_methods_EAS.csv"))
fwrite(het_tbl,     file.path(out_dir, "03_heterogeneity_EAS.csv"))
fwrite(pleio_tbl,   file.path(out_dir, "04_pleiotropy_egger_intercept_EAS.csv"))

message("\n=== Analysis complete ===")
message(sprintf("Results saved to: %s", out_dir))
print(primary_tbl)
