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
# D5 v8 - Full pipeline for B5
# 依赖：05_paf_attribution.R
# CHANGE FROM v7
#   [1] The structurally explicable category is emitted as SED, not SD.
#       Was:  class = if (length(adm)) "SD" else "UD")
#       Now:  class = if (length(adm)) "SED" else "UD")
#       Why:  SD is already defined in the manuscript as standard deviation, and
#             both meanings appear in figure 1b. The sealed output is the single
#             source of truth for the classification, so the final label has to
#             originate here rather than be translated downstream.
#       Re-run after this change to regenerate B5_out_classification.csv.
#
# (v7 header follows)
# D5 v7 - Full pipeline for B5
#
# Supersedes D5_v6. Every number quoted in the B5 inserts is emitted here.
#
# -----------------------------------------------------------------------------
# CHANGES FROM v6
#
# [7] Header, directories.
#     Was:  EUR_DIR / EAS_DIR only.
#     Now:  MVMR_DIR added, pointing at rawdat/MR/EUR/MVMR_res.
#     Why:  the multivariable results sit in a directory parallel to res/, not
#           inside it.
#
# [8] Section 6, multivariable input.
#     Was:  mv_file <- file.path(EUR_DIR, "02_all_MR_methods_EUR.csv")
#           ... mv[grepl("MVMR|multivariable", method, ignore.case = TRUE) ...]
#           if (!"mvmr_persists" %in% names(gen)) stop(...)
#     Now:  reads MVMR_res/01_MVMR_IVW_EUR.csv, treats every row as a
#           multivariable estimate, normalises column names and back-fills pval
#           and p_fdr when absent, and stops when no row matches rather than
#           filling the column with FALSE.
#     Why:  the dedicated MVMR file carries no `method` column, so the pattern
#           matched nothing. merge(all.x = TRUE) then produced an all-NA column
#           that the old guard accepted and the next line filled with FALSE, so
#           criterion 2 of the genetic ordering and mechanism M3 were both
#           skipped in silence. In the v6 run LDL was ordered above BMI by
#           criterion 3 instead, in the same direction, which is why the failure
#           did not show up in the classification.
#
# [9] Section 6, definition of criterion 2.
#     Was:  mvmr_persists = any multivariable estimate for that exposure passing
#           FDR, over all subtypes.
#     Now:  restricted to the exposure-by-subtype cells counted under criterion
#           1 in European ancestry, and requiring the direction to be preserved.
#     Why:  the rule as written in the Methods breaks ties by whether the
#           association counted under the first criterion persists, not by
#           whether some other association exists. Unchanged for SBP, LDL and
#           BMI; the code now says what the text says.
#
# [10] Section 7, M3().
#     Was:  m <- mv[grepl("MVMR|multivariable", method, ...) & ...]
#     Now:  m <- mv[exposure_name == ex & outcome_name == "IS"]
#     Why:  same missing `method` column as [8].
#
# -----------------------------------------------------------------------------
# CHANGES FROM v5
#
# [1] Section 6, factor keys.  NEW.
#     Was:  orderings <- rbind(gen[...], pop[...], ind[...])
#           pos <- dcast(orderings, exposure_name ~ stream, value.var = "position")
#     Now:  a CANON map is applied to every stream and to the evidence register
#           before any merge, and a guard stops the run if a factor is missing
#           from a stream.
#     Why:  the population stream labels kidney dysfunction KD while the cohort
#           and MR streams label the same construct eGFR. The merge produced six
#           rows instead of five, two of them half NA, and detect_divergence()
#           failed on `if (any(nds) && !all(nds))` with a missing value. A name
#           mismatch should stop the run at the merge, not surface three steps
#           later as a type error.
#
# [2] Section 6, genetic ordering.
#     Was:  gen[, position := fifelse(nd, max(which(!nd)) + 1L, seq_len(.N))]
#           gen[nd == TRUE, position := sum(!gen$nd) + 1L]
#     Now:  n_nonnull computed first, then a plain seq_len and an override.
#     Why:  max(which(!nd)) returns -Inf when no factor is flagged, and the
#           fifelse evaluated both branches regardless. The old line also
#           assigned positions before sorting placed the flagged factors last.
#
# [3] Section 6, individual ordering. Same restructuring as [2], same reason.
#
# [4] Section 7, detect_divergence().
#     Was:  if (any(nds) && !all(nds)) ...
#     Now:  an explicit NA check that stops with the offending factor named.
#     Why:  belt and braces. Guard [1] should make this unreachable; if it ever
#           fires, the message says which factor and which stream.
#
# [5] Section 7, M3().
#     Was:  mv[... exposure_name == ex ...]
#     Now:  mv$exposure_name passed through canon() when mv is loaded.
#     Why:  the multivariable output may use either label, same problem as [1].
#
# [6] Sections 7 and 8, output labels.  NEW.
#     Now:  classification carries factor_label, mapping eGFR back to KD for the
#           manuscript tables.
#     Why:  internal keys and printed labels should not be the same object.
# -----------------------------------------------------------------------------
#
# Outputs
#   B5_out_scales.csv          table S-B5c panel a   exposure scales, both ancestries
#   B5_out_rank_freq.csv       table S-B5c panel b   population-stream ordering inputs
#   B5_out_xanc.csv            table S-B5c panel c   cross-ancestry contrasts, per SD
#   B5_out_xanc_natural.csv    appendix 3.4          the same on a natural-unit scale
#   B5_out_subtype.csv         table S-B5b, M2       within-ancestry subtype contrasts
#   B5_out_perSD.csv           Results, appendix 3.3 per-SD estimates with 95% CI
#   B5_out_orderings.csv       appendix 3.1          within-stream positions
#   B5_out_classification.csv  table S-B5b           the classification itself
#   B5_out_controls.txt        table S-B5b footnote  the checks that must fail
#
# Before sourcing, put paf_mat in the environment:
#   paf_mat <- build_matrix(dat, CAUSE, 2023, "paf", rf_order) * 100
# =============================================================================

library(data.table); library(dplyr)

ROOT     <- DATA_ROOT
EUR_DIR  <- file.path(ROOT, "rawdat/MR/EUR/res")
EAS_DIR  <- file.path(ROOT, "rawdat/MR/EAS/res")
MVMR_DIR <- file.path(ROOT, "rawdat/MR/EUR/MVMR_res")                          # [7]
OUT      <- file.path(ROOT, "B5_outputs"); dir.create(OUT, showWarnings = FALSE)

EXPOSURES <- c("SBP", "LDL", "BMI", "FG", "eGFR")   # canonical keys
SUBTYPES  <- c("IS", "ICH", "SAH")                  # total stroke excluded
PRIMARY   <- list(exposure = "LDL", outcome = "IS")

# --- [1] canonical factor keys -----------------------------------------------
CANON <- c(SBP = "SBP", LDL = "LDL", BMI = "BMI",
           FG = "FG", FPG = "FG",
           eGFR = "eGFR", KD = "eGFR")
DISPLAY <- c(SBP = "SBP", LDL = "LDL", BMI = "BMI", FG = "FG", eGFR = "KD")

canon <- function(x) {
  out <- unname(CANON[as.character(x)])
  if (any(is.na(out)))
    stop("unmapped factor name(s): ", paste(unique(x[is.na(out)]), collapse = ", "),
         "\nAdd them to CANON rather than letting the merge produce NA.")
  out
}

N_GWAS <- data.table(
  ancestry      = rep(c("EUR", "EAS"), each = 5),
  exposure_name = rep(EXPOSURES, 2),
  n_gwas = c(97656, 70814, 99998, 200622, 76511,
             196273, 165202, 191458, 111635, 154633)
)

# =============================================================================
# 0. EVIDENCE REGISTER
#    Inputs from analyses outside this script, each with its source. Everything
#    downstream is derived; nothing else is hard-coded.
# =============================================================================
REGISTER <- list(
  
  tmrel = data.table(
    exposure_name = c("SBP", "LDL", "BMI", "FG", "eGFR"),
    lo   = c(110, 0.7, 20, 4.90, 60),
    hi   = c(115, 1.3, 25, 5.30, Inf),
    unit = c("mmHg", "mmol/L", "kg/m2", "mmol/L", "mL/min/1.73m2"),
    source = "GBD 2019 risk factors methods; GBD 2021 CVD risk factors (BPRF update)"
  ),
  
  window_gradient = data.table(
    exposure_name = c("SBP", "LDL", "BMI", "FG", "eGFR"),
    covered    = c(TRUE, FALSE, TRUE, TRUE, FALSE),
    hr         = c(0.94, NA, 1.80, NA, NA),
    lci        = c(0.87, NA, 1.48, NA, NA),
    uci        = c(1.01, NA, 2.18, NA, NA),
    nonlin_fdr = c(NA, NA, 0.0093, NA, NA),
    note = c("110 vs referent 115",
             "cohort minimum 80 mg/dL, above the window",
             "20 vs 25",
             "window only partly covered; no within-window contrast estimated",
             "ACR not measured; window not testable"),
    source = "CHARLS restricted cubic splines, Model 2"
  ),
  
  cohort = data.table(
    exposure_name = c("SBP", "BMI", "FG", "LDL", "eGFR"),
    shr = c(1.37, 1.22, 1.10, 1.02, 1.02),
    lci = c(1.26, 1.14, 1.04, 0.94, 0.93),
    uci = c(1.48, 1.30, 1.17, 1.11, 1.12),
    source = "CHARLS Fine-Gray Model 2, per 1 SD"
  ),
  
  subtype_share = data.table(
    outcome_name = c("IS", "ICH", "SAH"),
    share_lo = c(0.39, NA, NA), share_hi = c(0.42, NA, NA),
    source = "GBD 2023, C1 member states, share of stroke DALYs"
  ),
  
  pop_sd = data.table(
    ancestry      = "EAS",
    exposure_name = EXPOSURES,
    sd_raw  = c(21.1872, 34.31764, 3.845321, 32.56558, 14.01594),
    sd_unit = c("mmHg", "mg/dL", "kg/m2", "mg/dL", "mL/min/1.73m2"),
    source  = "CHARLS baseline, n=8012"
  )
)
for (nm in c("tmrel", "window_gradient", "cohort", "pop_sd"))
  REGISTER[[nm]][, exposure_name := canon(exposure_name)]

CONVERT <- data.table(exposure_name = c("LDL", "FG"),
                      from = "mg/dL", to = "mmol/L", divisor = c(38.67, 18.016))

# =============================================================================
# 1. Exposure scales, both ancestries through one function
# =============================================================================
STANDARDISED_BAND <- c(0.85, 1.20)
NATURAL_RANGE <- list(SBP = c(13, 24), LDL = c(0.70, 1.15), BMI = c(3.2, 6.0),
                      FG = c(0.40, 1.10), eGFR = c(11, 24))
NATURAL_UNIT <- c(SBP = "mmHg", LDL = "mmol/L", BMI = "kg/m2",
                  FG = "mmol/L", eGFR = "mL/min/1.73m2")
STEP <- c(SBP = 10, LDL = 1, BMI = 1, FG = 1, eGFR = 10)

recover_sd <- function(snp_dir, exposure, n_gwas) {
  f <- list.files(snp_dir, sprintf("^%s_on_.*_snps\\.csv$", exposure), full.names = TRUE)
  if (!length(f)) return(data.table(trait_sd = NA_real_, sd_p25 = NA_real_,
                                    sd_p75 = NA_real_, af_source = "no file", n_snp = 0L))
  d <- fread(f[1]); if ("mr_keep" %in% names(d)) d <- d[mr_keep == TRUE]
  af <- if ("eaf.exposure" %in% names(d) && any(!is.na(d$eaf.exposure)))
    list(f = d$eaf.exposure, src = "eaf.exposure")
  else if ("eaf.outcome" %in% names(d) && any(!is.na(d$eaf.outcome)))
    list(f = d$eaf.outcome, src = "eaf.outcome")
  else list(f = rep(NA_real_, nrow(d)), src = "unavailable")
  ok <- !is.na(af$f) & af$f > 0.01 & af$f < 0.99 &
    !is.na(d$se.exposure) & !is.na(d$beta.exposure)
  if (!any(ok)) return(data.table(trait_sd = NA_real_, sd_p25 = NA_real_,
                                  sd_p75 = NA_real_, af_source = af$src, n_snp = 0L))
  v <- sqrt(2 * af$f[ok] * (1 - af$f[ok]) *
              (n_gwas * d$se.exposure[ok]^2 + d$beta.exposure[ok]^2))
  data.table(trait_sd = median(v, na.rm = TRUE),
             sd_p25 = unname(quantile(v, .25, na.rm = TRUE)),
             sd_p75 = unname(quantile(v, .75, na.rm = TRUE)),
             af_source = af$src, n_snp = sum(ok))
}

scales <- N_GWAS[, {
  cbind(.SD, recover_sd(if (ancestry == "EUR") file.path(EUR_DIR, "snp_tables")
                        else file.path(EAS_DIR, "snp_tables"), exposure_name, n_gwas))
}, by = .(ancestry, exposure_name)]

scales[, unit := mapply(function(ex, s) {
  if (is.na(s)) return("unresolved")
  if (s >= STANDARDISED_BAND[1] && s <= STANDARDISED_BAND[2]) return("SD")
  r <- NATURAL_RANGE[[ex]]
  if (s >= r[1] && s <= r[2]) return(unname(NATURAL_UNIT[ex]))
  "unresolved"
}, exposure_name, trait_sd)]
scales[, mult_to_sd := fifelse(unit == "SD", 1,
                               fifelse(unit == "unresolved", NA_real_, trait_sd))]
scales[, factor_label := unname(DISPLAY[exposure_name])]
fwrite(scales, file.path(OUT, "B5_out_scales.csv"))
cat("\n== Panel a: exposure scales ==\n"); print(as.data.frame(scales))

# =============================================================================
# 2. Estimates on a per-SD scale, with intervals
# =============================================================================
res <- rbind(fread(file.path(EUR_DIR, "01_primary_IVW_results_EUR.csv"))[, ancestry := "EUR"],
             fread(file.path(EAS_DIR, "01_primary_IVW_results_EAS.csv"))[, ancestry := "EAS"],
             fill = TRUE)
res[, exposure_name := canon(exposure_name)]                                   # [1]
res <- merge(res, scales[, .(ancestry, exposure_name, unit, trait_sd, mult_to_sd)],
             by = c("ancestry", "exposure_name"))

res[, `:=`(b_sd = b * mult_to_sd, se_sd = se * mult_to_sd)]
res[, `:=`(or_sd = exp(b_sd), lci_sd = exp(b_sd - 1.96 * se_sd),
           uci_sd = exp(b_sd + 1.96 * se_sd))]

perSD <- res[outcome_name %in% SUBTYPES,
             .(ancestry, factor_label = unname(DISPLAY[exposure_name]),
               exposure = exposure_name, outcome = outcome_name, nsnp,
               or = round(or_sd, 3), lci = round(lci_sd, 3), uci = round(uci_sd, 3),
               reported = sprintf("%.3f (%.3f-%.3f)", or_sd, lci_sd, uci_sd),
               p = signif(pval, 2), p_fdr = signif(p_fdr, 2))][order(ancestry, outcome, exposure)]
fwrite(perSD, file.path(OUT, "B5_out_perSD.csv"))
cat("\n== Per-SD estimates ==\n"); print(as.data.frame(perSD))

# =============================================================================
# 3. Cross-ancestry contrasts, per SD
# =============================================================================
xanc_test <- function(d, bcol, secol, label) {
  w <- dcast(d[outcome_name %in% SUBTYPES & !is.na(get(bcol))],
             exposure_name + outcome_name ~ ancestry, value.var = c(bcol, secol))
  need <- c(paste0(bcol, c("_EUR", "_EAS")), paste0(secol, c("_EUR", "_EAS")))
  if (!all(need %in% names(w))) return(NULL)
  setnames(w, need, c("bE", "bA", "sE", "sA"))
  w <- w[!is.na(bE) & !is.na(bA) & !is.na(sE) & !is.na(sA)]
  if (!nrow(w)) return(NULL)
  w[, `:=`(diff = bE - bA, se_diff = sqrt(sE^2 + sA^2))]
  w[, Q := (diff / se_diff)^2]
  w[, `:=`(
    scale = label,
    OR_EUR = round(exp(bE), 3), OR_EAS = round(exp(bA), 3),
    p_het = pchisq(Q, 1, lower.tail = FALSE),
    I2 = round(pmax(0, (Q - 1) / Q) * 100, 1),
    prec_ratio = round(sA / sE, 2),
    power = round(pnorm(abs(fifelse(sA > sE, bE, bA)) / pmax(sA, sE) - 1.96), 2),
    mde_or = round(exp(2.802 * pmax(sA, sE)), 2)
  )]
  prim <- w$exposure_name == PRIMARY$exposure & w$outcome_name == PRIMARY$outcome
  w[, role := fifelse(prim, "primary", "exploratory")]
  w[, p_bh := NA_real_]
  if (any(!prim)) w[!prim, p_bh := p.adjust(p_het, "BH")]
  w[, .(scale, factor_label = unname(DISPLAY[exposure_name]),
        exposure = exposure_name, outcome = outcome_name, role,
        OR_EUR, OR_EAS, Q = round(Q, 3), p_het = signif(p_het, 2),
        p_bh = signif(p_bh, 2), I2, prec_ratio, power, mde_or)]
}

xanc <- xanc_test(res, "b_sd", "se_sd", "per 1 SD")
setorder(xanc, -role, exposure, outcome)
fwrite(xanc, file.path(OUT, "B5_out_xanc.csv"))
cat("\n== Panel c: cross-ancestry contrasts, per SD ==\n"); print(as.data.frame(xanc))

# =============================================================================
# 4. Cross-ancestry contrasts on a natural-unit scale
#    Secondary and partial. A standardised GWAS needs an external population SD
#    to return to natural units; where that is unavailable on either side the
#    exposure is skipped rather than half-converted. Intervals here understate
#    uncertainty, the external SD carrying no error term.
# =============================================================================
popsd <- merge(REGISTER$pop_sd, CONVERT, by = "exposure_name", all.x = TRUE)
popsd[, pop_sd := fifelse(is.na(divisor), sd_raw, sd_raw / divisor)]
popsd <- popsd[, .(ancestry, exposure_name, pop_sd,
                   pop_sd_unit = fifelse(is.na(divisor), sd_unit, to))]

res2 <- merge(res, popsd, by = c("ancestry", "exposure_name"), all.x = TRUE)
res2[, step := unname(STEP[exposure_name])]
res2[, `:=`(b_nat = NA_real_, se_nat = NA_real_)]
res2[unit == "SD" & !is.na(pop_sd),
     `:=`(b_nat = b_sd / pop_sd * step, se_nat = se_sd / pop_sd * step)]
res2[!unit %in% c("SD", "unresolved"),
     `:=`(b_nat = b * step, se_nat = se * step)]

ok_nat <- res2[outcome_name %in% SUBTYPES & !is.na(b_nat),
               .N, by = .(exposure_name, outcome_name)][N == 2]
xanc_nat <- if (nrow(ok_nat))
  xanc_test(res2[paste(exposure_name, outcome_name) %in%
                   paste(ok_nat$exposure_name, ok_nat$outcome_name)],
            "b_nat", "se_nat", "per natural unit") else NULL
if (!is.null(xanc_nat)) {
  xanc_nat[, `:=`(step = unname(STEP[exposure]), unit = unname(NATURAL_UNIT[exposure]))]
  fwrite(xanc_nat, file.path(OUT, "B5_out_xanc_natural.csv"))
  cat("\n== Cross-ancestry contrasts, natural units ==\n")
  print(as.data.frame(xanc_nat))
}
cat("\nSkipped on the natural-unit scale (external SD unavailable on one side): ",
    paste(setdiff(EXPOSURES, unique(xanc_nat$exposure)), collapse = ", "), "\n")

# =============================================================================
# 5. Within-ancestry subtype contrasts - scale-free, needed by mechanism M2
# =============================================================================
subtype <- rbindlist(lapply(unique(res$ancestry), function(anc) {
  d <- res[ancestry == anc & outcome_name %in% SUBTYPES]
  cb <- t(combn(SUBTYPES, 2))
  rbindlist(lapply(EXPOSURES, function(ex) rbindlist(lapply(seq_len(nrow(cb)), function(i) {
    a <- d[exposure_name == ex & outcome_name == cb[i, 1]]
    z <- d[exposure_name == ex & outcome_name == cb[i, 2]]
    if (!nrow(a) || !nrow(z)) return(NULL)
    Q <- ((a$b - z$b) / sqrt(a$se^2 + z$se^2))^2       # raw scale, units cancel
    data.table(ancestry = anc, exposure = ex,
               factor_label = unname(DISPLAY[ex]),
               contrast = paste(cb[i, 1], "vs", cb[i, 2]),
               Q = round(Q, 3), p = pchisq(Q, 1, lower.tail = FALSE))
  }))))
}))
subtype[, p_bh := p.adjust(p, "BH"), by = ancestry]
subtype[, `:=`(p = signif(p, 3), p_bh = signif(p_bh, 3))]
fwrite(subtype, file.path(OUT, "B5_out_subtype.csv"))
cat("\n== Subtype contrasts ==\n"); print(as.data.frame(subtype))

# =============================================================================
# 6. Within-stream orderings
# =============================================================================
# --- genetic -----------------------------------------------------------------
gen <- res[ancestry %in% c("EUR", "EAS") & outcome_name %in% SUBTYPES,
           .(cells = sum(b > 0 & p_fdr < 0.05, na.rm = TRUE),
             min_fdr = suppressWarnings(min(p_fdr[b > 0 & p_fdr < 0.05], na.rm = TRUE))),
           by = exposure_name]
gen[is.infinite(min_fdr), min_fdr := NA_real_]

# --- [8] multivariable input -------------------------------------------------
# The dedicated MVMR output carries no `method` column: every row is a
# multivariable IVW estimate. Column names are normalised across plausible
# spellings, and pval and p_fdr are back-filled from b and se when absent.
mv_file <- file.path(MVMR_DIR, "01_MVMR_IVW_EUR.csv")
if (!file.exists(mv_file))
  stop("MVMR file not found at ", mv_file,
       "\nCriterion 2 of the genetic ordering and mechanism M3 both depend on it.")
mv <- fread(mv_file)
setnames(mv, tolower(names(mv)))
ALIAS <- c(exposure = "exposure_name", exposure_id = "exposure_name",
           outcome = "outcome_name", outcome_id = "outcome_name",
           beta = "b", b_mvmr = "b", se_mvmr = "se",
           p = "pval", pvalue = "pval", p_value = "pval",
           fdr = "p_fdr", pval_fdr = "p_fdr", p_adj = "p_fdr")
for (from in names(ALIAS))
  if (from %in% names(mv) && !(ALIAS[[from]] %in% names(mv)))
    setnames(mv, from, ALIAS[[from]])

need_mv <- c("exposure_name", "outcome_name", "b", "se")
if (!all(need_mv %in% names(mv)))
  stop("MVMR file lacks required columns. Needed: ", paste(need_mv, collapse = ", "),
       "\nPresent: ", paste(names(mv), collapse = ", "))

mv[, exposure_name := canon(exposure_name)]                                    # [5]
if (!"pval" %in% names(mv))  mv[, pval := 2 * pnorm(-abs(b / se))]
if (!"p_fdr" %in% names(mv)) mv[, p_fdr := p.adjust(pval, "BH"), by = outcome_name]
if ("method" %in% names(mv)) mv <- mv[grepl("IVW|inverse", method, ignore.case = TRUE)]

# --- [9] criterion 2: do the cells counted under criterion 1 persist? --------
sig_cells <- res[ancestry == "EUR" & outcome_name %in% SUBTYPES &
                   b > 0 & p_fdr < 0.05, .(exposure_name, outcome_name)]
if (!nrow(sig_cells))
  stop("No European univariable cell passed criterion 1; the ordering cannot be built.")

mvp <- merge(mv, sig_cells, by = c("exposure_name", "outcome_name"))
if (!nrow(mvp))
  stop("No multivariable row matched any cell counted under criterion 1.\n",
       "Exposure values in the MVMR file: ", paste(unique(mv$exposure_name), collapse = ", "),
       "\nOutcome values: ", paste(unique(mv$outcome_name), collapse = ", "),
       "\nCells sought: ", paste(sprintf("%s-%s", sig_cells$exposure_name,
                                         sig_cells$outcome_name), collapse = ", "))
mvp <- mvp[, .(mvmr_persists = any(b > 0 & p_fdr < 0.05, na.rm = TRUE)), by = exposure_name]

gen <- merge(gen, mvp, by = "exposure_name", all.x = TRUE)
gen[is.na(mvmr_persists), mvmr_persists := FALSE]
cat("\n== Criterion 2: persistence in multivariable analysis ==\n")
print(as.data.frame(gen[, .(exposure_name, cells, mvmr_persists, min_fdr)]))

gen[, nd := cells == 0]                                                        # [2]
n_gen_nonnull <- sum(!gen$nd)
setorder(gen, nd, -cells, -mvmr_persists, min_fdr)
gen[, position := seq_len(.N)]
gen[nd == TRUE, position := n_gen_nonnull + 1L]
gen[, stream := "genetic"]

# --- population --------------------------------------------------------------
if (!exists("paf_mat")) stop("paf_mat not in the environment. Run:\n",
                             "  paf_mat <- build_matrix(dat, CAUSE, 2023, \"paf\", rf_order) * 100\n",
                             "before sourcing this script.")
rank_mat <- apply(paf_mat, 2, function(x) rank(-x, ties.method = "min"))
rank_freq <- as.data.table(t(apply(rank_mat, 1, function(r)
  table(factor(r, levels = seq_len(nrow(paf_mat)))))))
setnames(rank_freq, paste0("rank", seq_len(nrow(paf_mat))))
rank_freq[, `:=`(factor_label = rownames(paf_mat),
                 median_paf = round(apply(paf_mat, 1, median), 1),
                 paf_lo = round(apply(paf_mat, 1, min), 1),
                 paf_hi = round(apply(paf_mat, 1, max), 1),
                 modal_rank = apply(rank_mat, 1, function(r)
                   as.integer(names(which.max(table(r))))))]
setorder(rank_freq, -median_paf)
rank_freq[, position := seq_len(.N)]
fwrite(rank_freq, file.path(OUT, "B5_out_rank_freq.csv"))

pairwise <- rbindlist(lapply(seq_len(nrow(rank_freq) - 1), function(i) {
  a <- rank_freq$factor_label[i]; b <- rank_freq$factor_label[i + 1]
  data.table(higher = a, lower = b,
             states_consistent = sum(rank_mat[a, ] < rank_mat[b, ]),
             states_total = ncol(rank_mat))
}))
cat("\n== Panel b: population ordering and its stability ==\n")
print(as.data.frame(rank_freq)); print(as.data.frame(pairwise))

pop <- rank_freq[, .(exposure_name = canon(factor_label),                      # [1]
                     position, nd = FALSE, stream = "population")]

# --- individual --------------------------------------------------------------
ind <- copy(REGISTER$cohort)
ind[, nd := lci <= 1 & uci >= 1]                                               # [3]
n_ind_nonnull <- sum(!ind$nd)
setorder(ind, nd, -shr)
ind[, position := seq_len(.N)]
ind[nd == TRUE, position := n_ind_nonnull + 1L]
ind[, stream := "individual"]

orderings <- rbind(gen[, .(stream, exposure_name, position, nd)],
                   pop[, .(stream, exposure_name, position, nd)],
                   ind[, .(stream, exposure_name, position, nd)])

# --- [1] guard: every stream must carry the same keys before any merge -------
pos <- dcast(orderings, exposure_name ~ stream, value.var = "position")
if (anyNA(pos)) {
  miss <- melt(pos, id.vars = "exposure_name")[is.na(value)]
  stop("factor missing from a stream after harmonisation:\n",
       paste(sprintf("  %s absent from the %s stream", miss$exposure_name, miss$variable),
             collapse = "\n"))
}
if (!setequal(pos$exposure_name, EXPOSURES))
  stop("factor set does not match EXPOSURES: ",
       paste(symdiff <- union(setdiff(pos$exposure_name, EXPOSURES),
                              setdiff(EXPOSURES, pos$exposure_name)), collapse = ", "))

ndf <- dcast(orderings, exposure_name ~ stream, value.var = "nd")
setnames(ndf, setdiff(names(ndf), "exposure_name"),
         paste0("nd_", setdiff(names(ndf), "exposure_name")))
ev <- merge(pos, ndf, by = "exposure_name")
stopifnot(!anyNA(ev))

fwrite(merge(ev, data.table(exposure_name = names(DISPLAY),
                            factor_label = unname(DISPLAY)), by = "exposure_name"),
       file.path(OUT, "B5_out_orderings.csv"))
cat("\n== Within-stream positions ==\n"); print(as.data.frame(ev))

# =============================================================================
# 7. Divergence detection and mechanism tests
#    Each mechanism returns (admitted, reason). A mechanism admits only if it
#    predicts both the presence and the direction of the divergence observed and
#    its precondition holds.
# =============================================================================
detect_divergence <- function(ex, r) {
  nds <- c(r$nd_genetic, r$nd_population, r$nd_individual)
  p   <- c(r$genetic, r$population, r$individual)
  if (anyNA(nds) || anyNA(p))                                                  # [4]
    stop("missing stream value for ", ex,
         " - check CANON and the guard in section 6")
  if (any(nds) && !all(nds)) return(list(TRUE, "not-detected criterion"))
  if (max(p) - min(p) >= 2)
    return(list(TRUE, sprintf("positional, %d places", max(p) - min(p))))
  list(FALSE, "none")
}

M1 <- function(ex, r) {
  t1 <- REGISTER$tmrel[exposure_name == ex]
  g  <- REGISTER$window_gradient[exposure_name == ex]
  if (!nrow(g) || !g$covered)
    return(list(FALSE, "precondition fails: TMREL window not covered by the cohort"))
  if (is.na(g$hr) || (g$lci <= 1 & g$uci >= 1))
    return(list(FALSE, sprintf("precondition fails: no gradient inside the window (%s)", g$note)))
  if (!(r$population > r$individual))
    return(list(FALSE, "direction fails: the population stream does not rank the factor lower than the individual stream"))
  list(TRUE, sprintf("TMREL %g-%g %s within cohort range; within-window HR %.2f (%.2f-%.2f)",
                     t1$lo, t1$hi, t1$unit, g$hr, g$lci, g$uci))
}

M2 <- function(ex, r) {
  gsub <- res[ancestry == "EUR" & exposure_name == ex & outcome_name %in% SUBTYPES]
  hits <- gsub[b > 0 & p_fdr < 0.05, outcome_name]
  if (length(hits) != 1)
    return(list(FALSE, "precondition fails: effect not confined to a single subtype"))
  sh <- REGISTER$subtype_share[outcome_name == hits]
  if (!nrow(sh) || is.na(sh$share_hi) || sh$share_hi >= 0.5)
    return(list(FALSE, "precondition fails: that subtype is not a minority of the aggregated outcome"))
  if (!isTRUE(r$nd_individual))
    return(list(FALSE, "direction fails: the aggregating stream did not flag the factor"))
  cont <- subtype[ancestry == "EUR" & exposure == ex]
  list(TRUE, sprintf("genetic effect confined to %s (%.0f-%.0f%% of stroke DALYs); %s",
                     hits, sh$share_lo * 100, sh$share_hi * 100,
                     paste(sprintf("%s Q=%.2f p=%.3g", cont$contrast, cont$Q, cont$p),
                           collapse = "; ")))
}

M3 <- function(ex, r) {
  if (is.null(mv)) return(list(FALSE, "precondition fails: multivariable results unavailable"))
  u <- res[ancestry == "EUR" & exposure_name == ex & outcome_name == "IS"]
  m <- mv[exposure_name == ex & outcome_name == "IS"]                          # [10]
  if (!nrow(u) || !nrow(m))
    return(list(FALSE, "precondition fails: both estimates not available"))
  if (!(u$p_fdr < 0.05 && m$p_fdr[1] >= 0.05))
    return(list(FALSE, "direction fails: no attenuation from univariable to multivariable"))
  list(TRUE, sprintf("univariable %.3f (FDR %.3g) attenuates to %.3f (FDR %.3g)",
                     exp(u$b), u$p_fdr, exp(m$b[1]), m$p_fdr[1]))
}

M4 <- function(ex, r)
  list(FALSE, "qualifies the genetic position only; never admits a between-stream divergence")

X_BLOCK <- function(ex, r) {
  obs_pos <- !isTRUE(r$nd_population)
  if (obs_pos && isTRUE(r$nd_individual) && isTRUE(r$nd_genetic))
    return("X2 blocks prevalence-driven attribution")
  if (obs_pos && !isTRUE(r$nd_individual) && isTRUE(r$nd_genetic))
    return("X1 blocks residual confounding and reverse causation")
  NA_character_
}

classify <- rbindlist(lapply(ev$exposure_name, function(ex) {
  r  <- as.list(ev[exposure_name == ex])
  dv <- detect_divergence(ex, r)
  if (!dv[[1]])
    return(data.table(exposure_name = ex, divergence = "none",
                      admitted = "", declined = "", blocked = NA_character_, class = "C"))
  tests <- list(M1 = M1(ex, r), M2 = M2(ex, r), M3 = M3(ex, r), M4 = M4(ex, r))
  adm <- names(tests)[vapply(tests, function(t) isTRUE(t[[1]]), logical(1))]
  data.table(
    exposure_name = ex, divergence = dv[[2]],
    admitted = if (length(adm))
      paste(sprintf("%s: %s", adm, vapply(tests[adm], `[[`, character(1), 2)),
            collapse = " | ") else "",
    declined = paste(sprintf("%s: %s", setdiff(names(tests), adm),
                             vapply(tests[setdiff(names(tests), adm)], `[[`, character(1), 2)),
                     collapse = " | "),
    blocked = X_BLOCK(ex, r),
    class = if (length(adm)) "SED" else "UD")
}))

classification <- merge(ev, classify, by = "exposure_name")
classification[, factor_label := unname(DISPLAY[exposure_name])]               # [6]
setorder(classification, population)
setcolorder(classification, c("factor_label", "exposure_name",
                              "genetic", "population", "individual"))
fwrite(classification, file.path(OUT, "B5_out_classification.csv"))
cat("\n== Table S-B5b: classification ==\n")
print(as.data.frame(classification[, .(factor_label, genetic, population, individual,
                                       divergence, class, blocked)]))

# =============================================================================
# 8. Control checks - the tests that must fail
# =============================================================================
ctl <- c(
  sprintf("M1 declined for SBP  : %s", M1("SBP",  as.list(ev[exposure_name == "SBP"]))[[2]]),
  sprintf("M1 declined for KD   : %s", M1("eGFR", as.list(ev[exposure_name == "eGFR"]))[[2]]),
  sprintf("UD is non-empty      : %d of %d factors unclassified",
          sum(classification$class == "UD"), nrow(classification)),
  sprintf("Cross-ancestry power : %d of %d contrasts below 0.80",
          sum(xanc$power < 0.80), nrow(xanc)),
  sprintf("Contrasts surviving correction: %d",
          sum(xanc$p_bh < 0.05, na.rm = TRUE)),
  sprintf("Genetic ordering, criterion 2 exercised: %s",
          if (length(unique(gen$mvmr_persists)) > 1)
            "yes, it separates at least two factors" else
              "no, all factors share the same value; the tie fell to criterion 3"),
  sprintf("Subtype contrasts surviving correction in European ancestry: %s",
          paste(subtype[ancestry == "EUR" & p_bh < 0.05,
                        sprintf("%s %s (p=%.3g)", factor_label, contrast, p_bh)],
                collapse = "; "))
)
writeLines(ctl, file.path(OUT, "B5_out_controls.txt"))
cat("\n== Control checks ==\n"); cat(paste(ctl, collapse = "\n"), "\n")

# =============================================================================
# 9. What is still not machine-checked
#    The evidence register in section 0 is the boundary. TMREL windows, spline
#    gradients, cohort hazard ratios and subtype shares are typed in from other
#    analyses: traceable, but not recomputed here. If any changes, the
#    classification may change and this script must be re-run. Everything
#    downstream of the register is derived.
# =============================================================================