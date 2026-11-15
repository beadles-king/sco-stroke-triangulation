# =============================================================================
# Reproducibility repository - Wang et al., transportability-aware triangulation
# of metabolic stroke evidence across the Shanghai Cooperation Organisation.
# Paths are resolved by config.R; see README.md.
# =============================================================================
if (!exists("RAWDAT")) {
  if (file.exists("config.R")) source("config.R") else
    stop("Run this script from the code/ directory (the folder containing config.R).")
}

# ===========================================================================
# Figure 1a: structural heterogeneity vs geographical coherence across
# thirteen intergovernmental blocs.
#
# Changes in this version:
#   - CIS and the Arab League added. Both were previously absent without an
#     inclusion rule; the CIS in particular is the closest comparator to the
#     SCO and its omission was the most conspicuous.
#   - ASEAN now has 11 members: Timor-Leste acceded on Oct 26, 2025, before
#     the Jan 1, 2026 census date. Largest component 8 -> 9, SDI range
#     0.370 -> 0.385, component count unchanged at 3.
#   - Point area no longer encodes membership size. Membership is not one of
#     the selection dimensions, and encoding it made the OIC (57) and the
#     Commonwealth (56) the visual focus of a figure about the SCO. Area now
#     encodes the number of ancestry groups, which is the second axis of
#     structural heterogeneity and was otherwise absent from the panel.
#
# Requires: ggplot2, ggrepel
# ===========================================================================
setwd(file.path(RAWDAT, "大改"))
library(ggplot2)
library(ggrepel)

# Timor-Leste is a boundary case: an Austronesian-speaking majority east of
# Wallace's Line with substantial Papuan admixture, matched by neither the
# 1000 Genomes east Asian panel nor an Oceanian one. Coded EAS on the rule
# that the predominant population is mapped to the nearest available panel.
# Set to "OCE" to reproduce the sensitivity analysis: ASEAN then has two
# ancestry groups, one without a matched panel. Neither coding alters the
# comparison with the SCO (five groups, four super-regions).
TIMOR_LESTE_GROUP <- "EAS"

# "ancestry" (default) or "none". Membership size is deliberately not offered.
SIZE_BY <- "ancestry"

blocs <- data.frame(
  bloc = c("Nordic", "GCC", "EAEU", "CIS", "EU", "ASEAN", "SAARC",
           "LAS", "AU", "SCO", "CW", "BRICS", "OIC"),
  n_members    = c( 5,  6,  5,  9, 27, 11,  8, 22, 54, 10, 56, 10, 57),
  n_superreg   = c( 1,  1,  1,  1,  2,  2,  2,  2,  3,  4,  5,  6,  7),
  n_ancestry   = c( 1,  1,  2,  3,  1,  1,  1,  2,  2,  5,  6,  6,  7),
  n_nopanel    = c( 0,  0,  1,  1,  0,  0,  0,  0,  0,  1,  1,  0,  1),
  n_components = c( 3,  2,  2,  2,  5,  3,  3,  4,  7,  1, 39,  8, 11),
  largest_n    = c( 3,  5,  4,  8, 22,  9,  6, 18, 48, 10, 12,  3, 43),
  sdi_range    = c(0.028, 0.096, 0.197, 0.271, 0.147, 0.385, 0.397,
                   0.648, 0.541, 0.310, 0.546, 0.468, 0.648),
  stringsAsFactors = FALSE
)

if (TIMOR_LESTE_GROUP == "OCE") {
  i <- blocs$bloc == "ASEAN"
  blocs$n_ancestry[i] <- 2
  blocs$n_nopanel[i]  <- 1
}

blocs$largest_pct <- round(100 * blocs$largest_n / blocs$n_members)
blocs$is_study    <- blocs$bloc == "SCO"

## ---- axis -------------------------------------------------------------
## y is -log10(components): one component maps to 0 at the top, more
## components map further down. Breaks are relabelled with the raw counts.

blocs$yv <- -log10(blocs$n_components)
brk      <- c(1, 2, 3, 5, 10, 20, 40)
y_breaks <- -log10(brk)
y_limits <- c(min(blocs$yv) - 0.22, 0.10)

## ---- dodge coincident points ------------------------------------------
## Three blocs share (1 super-region, 2 components) and two share (2, 3).
## Without a dodge only one of each group is visible.

key <- paste(blocs$n_superreg, blocs$n_components)
blocs$xv <- blocs$n_superreg
for (k in unique(key[duplicated(key)])) {
  i <- which(key == k)
  blocs$xv[i] <- blocs$n_superreg[i] + seq(-0.13, 0.13, length.out = length(i))
}

## ---- label offsets ----------------------------------------------------
## nudge_y is in log10 units; anything above ~0.15 runs a label off the panel.

blocs$nx <- 0; blocs$ny <- 0
off <- function(b, x, y) { i <- blocs$bloc == b; blocs$nx[i] <<- x; blocs$ny[i] <<- y }

off("SCO",             0.52,  0.00)
off("Commonwealth",   -0.90,  0.00)
off("BRICS",           0.34,  0.07)
off("OIC",            -0.58, -0.08)
off("African Union",  -0.85,  0.05)
off("EU-27",          -0.70,  0.06)
off("Arab League",     0.42, -0.05)
off("Nordic Council", -0.42, -0.08)
off("GCC",            -0.34,  0.08)
off("EAEU",            0.10,  0.10)
off("CIS",             0.44, -0.02)
off("ASEAN",          -0.50,  0.07)
off("SAARC",           0.44, -0.07)

theme_lancet <- function(base_size = 8, base_family = "") {
  theme_bw(base_size = base_size, base_family = base_family) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.25, colour = "grey90"),
      panel.border     = element_blank(),
      axis.line        = element_line(linewidth = 0.4, colour = "black"),
      axis.ticks       = element_line(linewidth = 0.4, colour = "black"),
      axis.text        = element_text(colour = "black", size = base_size),
      axis.title       = element_text(colour = "black", size = base_size),
      legend.position  = "none",
      plot.margin      = margin(4, 10, 4, 4)
    )
}

p <- ggplot(blocs, aes(x = xv, y = yv))

if (SIZE_BY == "ancestry") {
  p <- p +
    geom_point(aes(size = n_ancestry),
               shape = 21, fill = "white", colour = "grey35", stroke = 0.6) +
    geom_point(data = subset(blocs, is_study), aes(size = n_ancestry),
               shape = 23, fill = "#1F5C99", colour = "#1F5C99", stroke = 0.6) +
    scale_size_continuous(range = c(1.6, 4.8), breaks = c(1, 3, 5, 7))
} else {
  p <- p +
    geom_point(shape = 21, size = 2.6, fill = "white",
               colour = "grey35", stroke = 0.6) +
    geom_point(data = subset(blocs, is_study), shape = 23, size = 3.2,
               fill = "#1F5C99", colour = "#1F5C99", stroke = 0.6)
}

p <- p +
  geom_text_repel(
    aes(label = bloc, fontface = ifelse(is_study, "bold", "plain")),
    size = 2.6,
    nudge_x = blocs$nx, nudge_y = blocs$ny,
    box.padding = 0.38, point.padding = 0.22,
    min.segment.length = 0.30, segment.size = 0.25,
    segment.colour = "grey55",
    max.overlaps = Inf, force = 1.4, seed = 42,
    xlim = c(0.25, 7.95), ylim = y_limits) +
  scale_x_continuous(breaks = 1:7, limits = c(0.25, 7.95),
                     expand = expansion(mult = 0.02)) +
  scale_y_continuous(breaks = y_breaks, labels = brk, limits = y_limits,
                     expand = expansion(mult = 0.03)) +
  labs(x = "GBD super-regions spanned",
       y = "Connected components formed by shared land borders") +
  coord_cartesian(clip = "off") +
  theme_lancet()

## ---- guard: every bloc must be labelled -------------------------------
built <- ggplot_build(p)
lab_layer <- which(vapply(p$layers, function(l)
  inherits(l$geom, "GeomTextRepel"), logical(1)))
n_lab <- nrow(built$data[[lab_layer]])
if (n_lab != nrow(blocs)) {
  warning(sprintf("only %d of %d labels retained - adjust off() offsets",
                  n_lab, nrow(blocs)))
} else {
  message("all ", n_lab, " labels retained")
}

if (!exists("EXPORT_PANEL_A")) EXPORT_PANEL_A <- TRUE
if (EXPORT_PANEL_A) {
  ggsave("fig1a_bloc_selection.pdf", p, width = 175, height = 110, units = "mm")
  if (requireNamespace("ragg", quietly = TRUE)) {
    ggsave("fig1a_bloc_selection.tiff", p, device = ragg::agg_tiff,
           width = 175, height = 110, units = "mm", dpi = 600, background = "white")
  } else {
    ggsave("fig1a_bloc_selection.tiff", p, width = 175, height = 110, units = "mm",
           dpi = 600, compression = "lzw")
  }
}
print(p)

## Draft legend for panel a:
## Structural heterogeneity and geographical coherence of thirteen
## intergovernmental blocs. Point area is proportional to the number of
## continental ancestry groups represented. The vertical axis shows the number
## of connected components formed by treating member states as vertices and
## shared international land borders between member states as edges, on a
## reversed logarithmic scale, so that blocs higher in the panel are more
## geographically coherent. The Shanghai Cooperation Organisation (filled
## diamond) is the only bloc forming a single component. Points sharing
## identical coordinates are offset horizontally for legibility. Full values,
## inclusion criteria, and assignment rules are given in the appendix.