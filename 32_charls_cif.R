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
# Figure S18 — Cumulative incidence functions of incident stroke by clinical
#              cut-off, with non-stroke death as the competing event (CHARLS)
#
# Panel order follows the effect-strength order used in figure 4 and tables
# S15-S18: SBP, BMI, FG, LDL, eGFR.
#
# Formatting conventions (manuscript-wide):
#   - decimal separator: middle dot U+00B7
#   - p values:  p < 0.0001            -> "<0\u00B70001"
#                0.0001 <= p < 0.01    -> four decimals
#                p >= 0.01             -> two decimals
#   - en dash U+2013 for ranges, minus U+2212 for negative numbers
#   - base font 8 pt, panel letters 11 pt bold, target width 180 mm
#
# Output: FigureS18.tif        (600 dpi LZW, rendered by ragg)
#         FigureS18.tif        (600 dpi, LZW)
#         FigureS18_300dpi.tif (300 dpi, LZW)
# =============================================================================

library(cmprsk)
library(ggplot2)
library(patchwork)
library(grid)
library(ragg)

# ---- 0. Global constants ----------------------------------------------------

DEC       <- "\u00B7"    # middle dot
EN_DASH   <- "\u2013"
MINUS     <- "\u2212"
BASE_SIZE <- 8           # pt
FONT      <- "Helvetica"   # metric-compatible Helvetica substitute

COL_LOW   <- "#2C5985"   # same blue as figure 5
COL_HIGH  <- "#A63603"

FIG_W_MM  <- 180         # hard upper limit, double column
FIG_H_MM  <- 118         # 2 rows; well inside the 240 mm limit

# ---- Output location: set this to wherever the figure files should go ----
OUT_DIR   <- file.path(RAWDAT, "CHARLS/CIF")         # e.g. file.path(DATA_ROOT, "figures")
out_path  <- function(f) {
  if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)
  file.path(OUT_DIR, f)
}

# ---- 1. Formatting helpers --------------------------------------------------

# Replace the decimal point with the middle dot, and the ASCII hyphen used by
# format() for negative numbers with a true minus sign.
dec <- function(x) {
  x <- gsub("\\.", DEC, x)
  gsub("(?<![0-9A-Za-z])-(?=[0-9])", MINUS, x, perl = TRUE)
}

# Manuscript p-value rule. Returns a display string, already dot-converted.
fmt_p <- function(p) {
  if (is.na(p)) return(NA_character_)
  if (p < 0.0001)  return(paste0("<0", DEC, "0001"))
  if (p < 0.01)    return(dec(sprintf("%.4f", p)))
  dec(sprintf("%.2f", p))
}

p_label <- function(p) paste0("Gray's p", if (p < 0.0001) "" else "=", fmt_p(p))

# Safety net: walk a finished gtable and convert any decimal point left in a
# text grob (e.g. a tick label produced by a scale we did not format by hand).
dotify <- function(g) {
  if (inherits(g, "gTree") || inherits(g, "gtable")) {
    g$children <- lapply(g$children, dotify)
  }
  if (inherits(g, "text") && !is.null(g$label)) {
    g$label <- dec(as.character(g$label))
  }
  g
}

# ---- 2. Exposure definitions (effect-strength order) ------------------------

expo_params <- list(
  SBP = list(
    var   = "sbp_high",
    title = "Systolic blood pressure (mmHg)",
    low   = "<140",   high = "\u2265140"
  ),
  BMI = list(
    var   = "bmi_high",
    title = paste0("Body-mass index (kg/m", "\u00B2", ")"),
    low   = "<28",    high = "\u226528"
  ),
  FG = list(
    var   = "fg_high",
    title = "Fasting plasma glucose (mg/dL)",
    low   = "<126",   high = "\u2265126"
  ),
  LDL = list(
    var   = "ldl_high",
    title = "LDL cholesterol (mg/dL)",
    low   = "<160",   high = "\u2265160"
  ),
  eGFR = list(
    var   = "egfr_low",
    title = paste0("eGFR (mL/min per 1", DEC, "73 m", "\u00B2", ")"),
    # for eGFR the higher-risk stratum is the one BELOW the cut-off
    low   = "\u226560", high = "<60"
  )
)

# ---- 3. Extraction ----------------------------------------------------------

extract_cif_stratified <- function(expo_name) {
  
  params <- expo_params[[expo_name]]
  
  # Build the grouping as an explicit factor so that stratum order is fixed and
  # never depends on alphabetical sorting inside cuminc().
  raw <- charls[[params$var]]
  grp <- factor(ifelse(as.numeric(raw) == 1, params$high, params$low),
                levels = c(params$low, params$high))
  
  keep <- !is.na(grp) & !is.na(charls$time_yr) & !is.na(charls$new_stroke)
  grp  <- grp[keep]
  
  cif <- cuminc(ftime   = charls$time_yr[keep],
                fstatus = charls$new_stroke[keep],
                group   = grp,
                cencode = 0)          # 0 censored, 1 stroke, 2 non-stroke death
  
  p_gray <- cif$Tests["1", "pv"]      # Gray's test for cause 1 (stroke)
  
  # Element names are "<group label> <cause>"; parse the label back out rather
  # than assuming the list order.
  nm     <- grep(paste0(" ", 1, "$"), names(cif), value = TRUE)
  labels <- sub(" 1$", "", nm)
  
  cif_df <- do.call(rbind, Map(function(el, lab) {
    data.frame(time  = cif[[el]]$time,
               est   = cif[[el]]$est,
               group = lab,
               stringsAsFactors = FALSE)
  }, nm, labels))
  
  n_by <- table(grp)
  keys <- setNames(sprintf("%s (n=%d)", levels(grp), as.integer(n_by)),
                   levels(grp))
  
  cif_df$group <- factor(keys[cif_df$group], levels = unname(keys))
  
  list(data    = cif_df,
       p_label = p_label(p_gray),
       p_raw   = p_gray,
       title   = params$title,
       keys    = unname(keys))
}

# ---- 4. Panel ---------------------------------------------------------------

plot_cif_strat <- function(expo_name) {
  
  res    <- extract_cif_stratified(expo_name)
  cif_df <- res$data
  
  # Extend every curve to the end of follow-up so the steps do not stop short.
  t_max <- max(charls$time_yr, na.rm = TRUE)
  tail_rows <- do.call(rbind, lapply(levels(cif_df$group), function(g) {
    sub <- cif_df[cif_df$group == g, ]
    data.frame(time = t_max, est = sub$est[which.max(sub$time)], group = g)
  }))
  tail_rows$group <- factor(tail_rows$group, levels = levels(cif_df$group))
  cif_df <- rbind(cif_df, tail_rows)
  
  ggplot(cif_df, aes(x = time, y = est, colour = group)) +
    geom_step(linewidth = 0.5) +
    scale_x_continuous(
      breaks = seq(0, 7, 1),
      limits = c(0, t_max),
      expand = expansion(mult = c(0.02, 0.02))
    ) +
    scale_y_continuous(
      breaks = seq(0, 0.15, 0.05),
      limits = c(0, 0.15),
      labels = function(x) dec(sprintf("%.2f", x)),
      expand = expansion(mult = c(0.01, 0.02))
    ) +
    scale_colour_manual(values = setNames(c(COL_LOW, COL_HIGH), res$keys)) +
    annotate("text", x = t_max, y = 0.15, label = res$p_label,
             hjust = 1, vjust = 1, size = BASE_SIZE / .pt,
             colour = "grey25", family = FONT) +
    labs(title = res$title, x = "Time (years)",
         y = "Cumulative incidence", colour = NULL) +
    theme_classic(base_size = BASE_SIZE, base_family = FONT) +
    theme(
      plot.title           = element_text(size = BASE_SIZE, face = "plain",
                                          hjust = 0, margin = margin(b = 3)),
      axis.title           = element_text(size = BASE_SIZE),
      axis.text            = element_text(size = BASE_SIZE, colour = "black"),
      axis.line            = element_line(linewidth = 0.3, colour = "black"),
      axis.ticks           = element_line(linewidth = 0.3, colour = "black"),
      # ggplot2 >= 3.5 deprecates a numeric legend.position; if you are on 3.5+
      # replace the next line with:
      #   legend.position = "inside", legend.position.inside = c(0.03, 0.86),
      legend.position      = c(0.03, 0.86),
      legend.justification = c(0, 1),
      legend.background    = element_blank(),
      legend.key           = element_blank(),
      legend.key.height    = unit(3.2, "mm"),
      legend.key.width     = unit(4.5, "mm"),
      legend.text          = element_text(size = BASE_SIZE - 1),
      legend.margin        = margin(0, 0, 0, 0),
      plot.margin          = margin(2, 4, 2, 2, "mm")
    )
}

# ---- 5. Assembly ------------------------------------------------------------

plots <- lapply(names(expo_params), plot_cif_strat)
names(plots) <- names(expo_params)

figS18 <- wrap_plots(plots, ncol = 3) +
  plot_spacer() +                      # keeps the 3x2 grid square
  plot_layout(ncol = 3) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(plot.margin = margin(2, 2, 2, 2, "mm"))
  ) &
  theme(plot.tag          = element_text(size = 11, face = "bold",
                                         family = FONT),
        plot.tag.position = c(0, 1))

# ---- 6. Output --------------------------------------------------------------

mm2in <- function(x) x / 25.4

draw <- function() {
  g <- patchwork::patchworkGrob(figS18)
  grid.newpage()
  grid.draw(dotify(g))               # final decimal-separator sweep
}

# 原先此处另出一份 cairo_pdf 矢量母版, 已移除: cairo 把每个字形的步进四舍五入
# 到整数 point (0.353 mm), 再用 TJ 微调补偿, 同一单词内字距不匀 (最坏 1 pt)。
# 投稿件本就是下面两份 TIF, ragg 在目标分辨率排字, 600 dpi 误差 <= 0.042 mm。

# 600 dpi submission master, LZW
agg_tiff(out_path("FigureS18.tif"),
         width = FIG_W_MM, height = FIG_H_MM, units = "mm",
         res = 600, compression = "lzw", background = "white")
draw()
dev.off()

# 300 dpi copy for the review PDF
agg_tiff(out_path("FigureS18_300dpi.tif"),
         width = FIG_W_MM, height = FIG_H_MM, units = "mm",
         res = 300, compression = "lzw", background = "white")
draw()
dev.off()

# ---- 7. Machine-readable record of the plotted p values ---------------------
# Write the exact Gray's p values so the figure, the legend, and the text can
# be checked against one source.

p_tab <- do.call(rbind, lapply(names(expo_params), function(e) {
  r <- extract_cif_stratified(e)
  data.frame(exposure = e, p_gray = r$p_raw, displayed = r$p_label)
}))
write.csv(p_tab, out_path("FigureS18_grays_p.csv"), row.names = FALSE)
print(p_tab)