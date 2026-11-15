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
# Figure 1b - Three evidence streams, three orderings
#
# Companion panel to figure 1a. Both panels are now 175 mm wide, so figure 1
# occupies the full double-column width and stacks to 175 x 210 mm.
#
# -----------------------------------------------------------------------------
# CHANGES FROM v9
#
# [25] The class validator accepts a legacy "SD" label from an older pipeline
#      run, rewrites it to SED, and warns. This is a migration aid, not a
#      permanent tolerance: regenerate B5_out_classification.csv from the
#      pipeline and the branch stops firing.
#
# -----------------------------------------------------------------------------
# CHANGES FROM v8
#
# [24] The structurally explicable category is abbreviated SED, not SD.
#      Why:  SD is already defined in the manuscript's abbreviation list as
#            standard deviation, and both meanings appeared in this panel -
#            "SHR per 1 SD" in the individual column and SD in the
#            classification strip. SED derives properly from the full term and
#            keeps "structurally", which is the word that distinguishes a
#            prespecified structural mechanism from a post-hoc story.
#      X_MECH moved 9.14 to 9.25 for the extra character.
#
# -----------------------------------------------------------------------------
# CHANGES FROM v7
#
# [22] The multivariable-persistence marker changed from a dagger to an
#      asterisk, and a run-time check added for characters outside Latin-1.
#      Was:  "4/6 \u2020" and "\u2020 persists in multivariable analysis"
#      Now:  "4/6 *" and "* persists in multivariable analysis"
#      Why:  the dagger is outside Latin-1, which is the default encoding of
#            grDevices::pdf(), so it was silently dropped and rendered as a full
#            stop in the exported PDF. The middle dot survives because it is
#            inside Latin-1. cairo_pdf would handle the dagger but distorts type
#            on this system, so the character is changed rather than the device.
#            An asterisk is ASCII and is the first footnote marker in
#            Lancet-family style in any case.
#
# [23] Export: ragg TIFF is the deliverable; PNG added for quick inspection;
#      PDF is optional and emits a warning naming any glyph the default device
#      will drop.
#
# -----------------------------------------------------------------------------
# CHANGES FROM v6
#
# [18] Paths. ROOT is the output directory for this revision; the scripts live
#      on a separate branch, so CODE_DIR is declared independently and FIG_DIR
#      is ROOT itself rather than a new subdirectory.
#
# [19] White background written into the theme.
#      Was:  theme_void() alone, with bg = "white" passed to the device.
#      Now:  plot.background and panel.background set to white in the theme.
#      Why:  theme_void() leaves both elements blank, so the exported TIFF
#            carried an alpha channel with a fully transparent canvas. Viewers
#            that ignore TIFF alpha render that as black, which is the dark
#            background. A device argument cannot override a theme element, so
#            the fix has to be in the theme. The same theme is applied to the
#            assembled figure through plot_annotation.
#
# [20] Assembly conditions split, and the working directory switched while
#      sourcing panel a.
#      Was:  one combined if() whose message named both possible causes.
#      Now:  the missing file and the missing package report separately, and a
#            missing file lists any fig1a*.R found under CODE_DIR or ROOT.
#      Why:  the combined message could not say which of the two had failed,
#            and the path it suggested was in the wrong tree. Sourcing from
#            another directory also breaks any relative path inside panel a's
#            script, so the working directory is moved for the duration.
#
# [21] Individual subtitle shortened to "SHR per 1 SD" on one line; the
#      abbreviation is already defined in the manuscript's abbreviation list.
#
# -----------------------------------------------------------------------------
# CHANGES FROM v5
#
# [16] Classification strip reordered.
#      Was:  by population rank - SBP, LDL, KD, FG, BMI - so that a structurally
#            explicable factor sat below two unexplained ones.
#      Now:  by genetic position, ties broken by population rank - SBP, LDL,
#            BMI, KD, FG - which groups the classes as C, SD, SD, UD, UD.
#      The same ordering now drives the horizontal offset of tied not-detected
#      nodes, so one variable governs both.
#
# [17] Strip header and subtitle moved right by 0.25 units, and the individual
#      subtitle re-split.
#      Why:  "Subdistribution hazard ratio" centred on the individual column
#            crossed the divider rule, and "Pattern" began immediately to its
#            right, so the two ran together. The subtitle now breaks after
#            "hazard" and the strip header starts clear of the rule.
#
# -----------------------------------------------------------------------------
# CHANGES FROM v4
#
# [13] FG's population value returned to the left of its node, and the KD label
#      returned to dy 0.40 so that all five labels share one offset.
#      Was:  the value sat above the FG node, to avoid breaking the KD line.
#      Why:  above the node it falls between the KD and FG nodes and reads as
#            though it belonged to KD. Ambiguous attribution is a worse fault
#            than a short break in a line, and breaks already occur elsewhere
#            for the same reason - the KD value crosses the descending LDL
#            segment. Label padding cut from 0.25 mm to 0.15 mm so that the
#            breaks are as small as the type allows.
#
# [14] Figure 1 is assembled from two live ggplot objects rather than by reading
#      the exported panel a back in as a raster.
#      Was:  magick::image_read() on fig1a_bloc_selection.tiff, drawn with
#            annotation_raster() into a new bitmap.
#      Why:  that resamples an already-rasterised image onto a second pixel
#            grid, and the mismatch between the two grids drops pixels out of
#            thin strokes - the mottled letterforms in panel a's labels. Keeping
#            both panels as vector until a single final rasterisation removes
#            the problem at source.
#
# [15] Rasterisation moved to ragg::agg_tiff at 600 dpi, with PDF exported
#      alongside.
#      Why:  ragg renders type far more cleanly than the default TIFF device,
#            and this is line art, for which 300 dpi is the halftone floor
#            rather than an appropriate resolution. Submit the PDF if the
#            journal accepts vector artwork.
#
# -----------------------------------------------------------------------------
# CHANGES FROM v3
#
# [10] Column geometry compressed and shifted left.
#      Was:  X_GEN 1.30, X_POP 4.00, X_IND 6.70, X_RULE 7.80, X_MAX 10.5
#      Now:  X_GEN 1.05, X_POP 3.65, X_IND 6.25, X_RULE 7.30, X_MAX 10.05,
#            with the classification strip tightened to match.
#      Why:  the column subtitles were set narrow and tall. Widening them
#            requires the blocks to sit closer together, and closing the gap
#            would have pushed the individual subtitle into "Pattern" unless
#            the whole slopegraph moved left and the strip narrowed.
#
# [11] Node values: size 2.40 to 2.75, colour #4D4D4D to #333333.
#      Subtitles: size 2.55 to 2.70, population split rebalanced to two lines
#      of similar length.
#
# [12] FG label returned to below its population node, as for the other four.
#      Its value moved above the node instead.
#      Why:  three lines converge on the FG node, and both the left and right
#            sides of it are crossed - the left by the KD line at y 3.9-4.3,
#            the right by the BMI line at y 4.4 - so a value label on either
#            side breaks a line. Above the node is the only clear space. The KD
#            label was raised by 0.06 to open the gap between the two.
#
# -----------------------------------------------------------------------------
# CHANGES FROM v2
#
# [6] Column headers were clipped by the panel edge and overlapped their
#     subtitles.
#     Was:  annotate(..., y = HEAD_Y, vjust = 0) with HEAD_Y = -0.72 and
#           scale_y_reverse(limits = c(Y_ND + 1.20, -1.02))
#     Now:  vjust = 1 throughout the header block, HEAD_Y = -0.95,
#           SUB_Y = -0.60, upper limit -1.02, panel height 105 mm.
#     Why:  under scale_y_reverse, vjust = 0 anchors the bottom of the text and
#           the string grows upward on the device, so a header placed at the
#           top of the range ran off the panel and into its own subtitle.
#           vjust = 1 anchors the top and the string grows downward.
#
# [7] Type sizes raised throughout: subtitles 2.20 to 2.55, node values 2.10 to
#     2.40, position numbers and footnote to 2.50, factor labels 2.60 to 2.75.
#     Why:  the smallest type was at or below the 6 pt floor for print.
#
# [8] Rules darkened: grid #E5E5E5 to #D0D0D0 at 0.28 to 0.45 linewidth, the
#     strip divider to #BFBFBF at 0.65, the band fill to #EBEBEB, connecting
#     lines 0.62 to 0.80, nodes slightly larger, muted text #595959 to #4D4D4D.
#
# [9] Subtitles now start with a capital and are shortened so that adjacent
#     columns do not collide at the larger size. FG carries its population
#     value to the left of its node, where it does not break the two lines that
#     cross there.
#
# -----------------------------------------------------------------------------
# CHANGES FROM v1
#
# [1] Validation block.
#     Was:  read B5_out_classification.csv and infer the not-detected set as
#           ref[, nd_ref := pos_ref == max(pos_ref), by = stream]
#     Now:  read B5_out_orderings.csv, which carries nd_genetic, nd_population
#           and nd_individual explicitly.
#     Why:  the population stream has no not-detected factors and its positions
#           run 1 to 5, so inferring the flag from the lowest position marked
#           BMI as not detected and the run stopped with "BMI / population:
#           figure 5, output 5". Not-detected is a separate flag, not a function
#           of position, and in the population stream it can never be set - the
#           point the panel itself makes.
#
# [2] Width, 120 mm to 175 mm; height 100 mm.
#     Why:  at 120 mm the panel had room for the orderings and nothing else.
#           Three things were dropped that are worth showing: the stream-native
#           value at each node, the range of within-state ranks behind each
#           population position, and the reason LDL sits above BMI in the
#           genetic stream when both count one significant cell.
#
# [3] Node values.  NEW.
#     Each node carries the quantity its own stream was ordered by. The three
#     quantities are visibly different kinds of number, which is the intended
#     reading: they cannot be compared across columns.
#
# [4] Range of within-state ranks, shown as text beside the population value
#     rather than as a bar at the node.
#     Why:  four of the five ranges overlap in y, so bars drawn at a single x
#           pile on one another and resolve nothing.
#
# [5] Not-detected node values are right-aligned on the left-hand node of each
#     tied pair, so the two labels do not run together.
# -----------------------------------------------------------------------------
#
# Encoding
#   y            within-stream position; internal to each stream, with no
#                quantitative equivalence between columns
#   band below   factors flagged not detected in that stream
#   columns      ordered by decreasing independence from observational
#                confounding: genetic, population, individual
#   right strip  classification from the triangulation decision rules, in
#                population-stream order
#
# Positions and classes are those sealed in B5 and are validated against the
# pipeline output, so the figure cannot drift from the classification it depicts.
#
# NOTE for figure 1a: change its export to
#   ggsave("fig1a_bloc_selection.pdf",  p, width = 175, height = 110, units = "mm")
#   ggsave("fig1a_bloc_selection.tiff", p, width = 175, height = 110, units = "mm",
#          dpi = 300, compression = "lzw")
# Leave base_size at 8. The canvas grows, the type stays at its print size, and
# ggrepel re-lays the crowded GCC / EAEU / CIS labels with the extra room.
# =============================================================================

library(ggplot2)
library(data.table)

# [18] ROOT is where this revision's outputs live; the scripts sit elsewhere.
ROOT      <- file.path(RAWDAT, "大改")
# CODE_DIR is defined in config.R
ORDER_CSV <- file.path(ROOT, "B5_outputs/B5_out_orderings.csv")
CLASS_CSV <- file.path(ROOT, "B5_outputs/B5_out_classification.csv")
FIG_DIR   <- ROOT

# --- geometry, in plot units -------------------------------------------------
X_GEN <- 1.05; X_POP <- 3.65; X_IND <- 6.25
X_RULE <- 7.30; X_SWATCH <- 7.55; X_NAME <- 7.92; X_CLASS <- 8.72; X_MECH <- 9.25
X_MIN <- -0.35; X_MAX <- 10.05
Y_ND  <- 6.30       # the not-detected band
NUDGE <- 0.30       # horizontal offset for tied not-detected nodes
VDX   <- 0.26       # offset of a node's value label

# --- palette -----------------------------------------------------------------
# Avoids the orange and green used for clusters C1 and C2 in figure 2, so that a
# reader moving between figures cannot read these lines as cluster membership.
PAL <- c(SBP = "#1F5C99",   # matches the SCO highlight in panel a
         LDL = "#6A4C93",
         BMI = "#C1272D",
         KD  = "#00707C",
         FG  = "#8B6914")
INK <- "#000000"; MUTED <- "#4D4D4D"; GRID <- "#D0D0D0"; BAND <- "#EBEBEB"
RULE_COL <- "#BFBFBF"

# --- the sealed positions ----------------------------------------------------
nodes <- data.table(
  factor = rep(c("SBP", "LDL", "BMI", "FG", "KD"), each = 3),
  stream = rep(c("genetic", "population", "individual"), 5),
  pos = c(1, 1, 1,
          2, 2, NA,
          3, 5, 2,
          NA, 4, 3,
          NA, 3, NA)
)
nodes[, nd := is.na(pos)]

# Each node's value is the quantity its own stream was ordered by: counted cells
# out of six for the genetic stream, median attributable fraction with the range
# of within-state ranks for the population stream, and the mutually adjusted
# subdistribution hazard ratio per 1 SD for the individual stream.
vals <- data.table(
  factor     = c("SBP", "LDL", "BMI", "FG", "KD"),
  genetic    = c("4/6 *", "1/6 *", "1/6", "0/6", "0/6"),
  population = c("63\u00b71 (1)", "17\u00b70 (2\u20133)", "8\u00b77 (4\u20135)",
                 "9\u00b76 (3\u20135)", "11\u00b74 (2\u20134)"),
  individual = c("1\u00b737", "1\u00b702", "1\u00b722", "1\u00b710", "1\u00b702")
)
vals <- melt(vals, id.vars = "factor", variable.name = "stream", value.name = "val")

# [16] Strip order: genetic position first, ties broken by population rank.
# SBP, LDL and BMI hold genetic positions 1 to 3; KD and FG are both flagged not
# detected there and are separated by their population positions, 3 and 4. The
# resulting sequence groups the classes as C, SD, SD, UD, UD.
classes <- data.table(
  factor     = c("SBP", "LDL", "BMI", "KD", "FG"),
  strip_rank = 1:5,
  class      = c("C", "SED", "SED", "UD", "UD"),
  mech       = c("", "M2", "M1, M3", "", "")
)

# Label offsets from the population node. That column is the only one in which
# all five positions differ, so it is the only collision-free place to label the
# lines. FG is nudged left so its label does not break the BMI line where they
# cross.
LABOFF <- data.table(factor = c("SBP", "LDL", "KD", "FG", "BMI"),
                     dx = c(0, 0, 0, 0, 0),
                     dy = c(0.40, 0.40, 0.40, 0.40, 0.40))

# [13] Nodes whose value sits to the left of the node. The default right-hand
# side of the FG population node is crossed by the BMI line; the left side is
# crossed by KD, so a short break in one line or the other is unavoidable there
# and the left side is chosen because it keeps the value unambiguously on FG's
# row.
VAL_LEFT_EXTRA <- data.table(factor = "FG", stream = "population")

# --- [1] validation against the sealed output --------------------------------
if (file.exists(ORDER_CSV)) {
  ord <- fread(ORDER_CSV)
  ref <- melt(ord[, .(factor = factor_label, genetic, population, individual)],
              id.vars = "factor", variable.name = "stream", value.name = "pos_ref")
  ndr <- melt(ord[, .(factor = factor_label, genetic = nd_genetic,
                      population = nd_population, individual = nd_individual)],
              id.vars = "factor", variable.name = "stream", value.name = "nd_ref")
  ref <- merge(ref, ndr, by = c("factor", "stream"))
  cmp <- merge(nodes, ref, by = c("factor", "stream"))
  if (nrow(cmp) != nrow(nodes))
    stop("factor or stream labels do not match between the figure and ", basename(ORDER_CSV))
  bad <- cmp[(nd != nd_ref) | (!nd & pos != pos_ref)]
  if (nrow(bad))
    stop("figure positions disagree with the sealed classification:\n",
         paste(sprintf("  %s / %s: figure pos %s nd %s, output pos %s nd %s",
                       bad$factor, bad$stream, bad$pos, bad$nd,
                       bad$pos_ref, bad$nd_ref), collapse = "\n"))
  message("positions validated against ", basename(ORDER_CSV))
} else {
  warning("B5_out_orderings.csv not found; positions not validated")
}

if (file.exists(CLASS_CSV)) {
  cls <- fread(CLASS_CSV)[, .(factor = factor_label, class_ref = class)]
  # [25] the pipeline emitted SD for the structurally explicable category before
  # the abbreviation changed; accept the legacy label once, and say so
  if (any(cls$class_ref == "SD")) {
    warning("B5_out_classification.csv still labels the structurally explicable ",
            "category SD. Change class = \"SD\" to \"SED\" in the pipeline and ",
            "re-run, so that the sealed output carries the final label.")
    cls[class_ref == "SD", class_ref := "SED"]
  }
  badc <- merge(classes, cls, by = "factor")[class != class_ref]
  if (nrow(badc))
    stop("figure classes disagree with the sealed classification: ",
         paste(sprintf("%s figure %s, output %s", badc$factor, badc$class, badc$class_ref),
               collapse = "; "))
  message("classes validated against ", basename(CLASS_CSV))
}

# --- resolve coordinates -----------------------------------------------------
nodes[, x0 := fifelse(stream == "genetic", X_GEN,
                      fifelse(stream == "population", X_POP, X_IND))]
nodes[, y := fifelse(nd, Y_ND, as.numeric(pos))]
nodes <- merge(nodes, classes[, .(factor, strip_rank)], by = "factor")

# Tied not-detected nodes are separated horizontally in the same order as the
# classification strip, so that one variable governs both.
nodes[, x := x0]
nodes[nd == TRUE, x := x0 + fifelse(rank(strip_rank, ties.method = "first") == 1,
                                    -NUDGE, NUDGE), by = stream]
nodes <- merge(nodes, vals, by = c("factor", "stream"))

# [5] the left-hand node of a tied pair carries its value to the left, [13] as
# does any node listed in VAL_LEFT_EXTRA
nodes[, val_left := (nd & (x < x0)) |
        paste(factor, stream) %in% paste(VAL_LEFT_EXTRA$factor, VAL_LEFT_EXTRA$stream)]
nodes[, val_x := fifelse(val_left, x - VDX, x + VDX)]
nodes[, val_h := fifelse(val_left, 1, 0)]

segs <- nodes[order(factor, x)][
  , .(x = head(x, -1), y = head(y, -1), xend = tail(x, -1), yend = tail(y, -1)),
  by = factor]

lab   <- merge(nodes[stream == "population"], LABOFF, by = "factor")
strip <- copy(classes)[, y := strip_rank]

HEAD_Y <- -0.95; SUB_Y <- -0.60

# --- plot --------------------------------------------------------------------
p <- ggplot() +
  
  annotate("rect", xmin = X_MIN + 0.55, xmax = X_RULE - 0.25,
           ymin = Y_ND - 0.42, ymax = Y_ND + 0.42, fill = BAND, colour = NA) +
  annotate("segment", x = X_MIN + 0.55, xend = X_RULE - 0.25, y = 1:5, yend = 1:5,
           colour = GRID, linewidth = 0.45) +
  annotate("segment", x = X_RULE, xend = X_RULE, y = -0.30, yend = Y_ND + 0.55,
           colour = RULE_COL, linewidth = 0.65) +
  
  geom_segment(data = segs, aes(x = x, y = y, xend = xend, yend = yend, colour = factor),
               linewidth = 0.80, lineend = "round") +
  geom_point(data = nodes[nd == FALSE], aes(x = x, y = y, colour = factor),
             size = 2.25, shape = 16) +
  geom_point(data = nodes[nd == TRUE], aes(x = x, y = y, colour = factor),
             size = 2.15, shape = 21, fill = "white", stroke = 0.75) +
  
  # [3] stream-native value at each node
  geom_label(data = nodes, aes(x = val_x, y = y, label = val, hjust = val_h),
             fill = "white", label.size = 0, label.padding = unit(0.15, "mm"),
             size = 2.75, colour = "#333333") +
  
  geom_label(data = lab, aes(x = x + dx, y = y + dy, label = factor, colour = factor),
             fill = "white", label.size = 0, label.padding = unit(0.35, "mm"),
             size = 2.75, fontface = "bold", show.legend = FALSE) +
  
  annotate("text", x = c(X_GEN, X_POP, X_IND), y = HEAD_Y, hjust = 0.5, vjust = 1,
           label = c("Genetic", "Population", "Individual"),
           size = 2.95, fontface = "bold", colour = INK) +
  annotate("text", x = c(X_GEN, X_POP, X_IND), y = SUB_Y, hjust = 0.5, vjust = 1,
           label = c("Directionally consistent,\nsignificant cells (of six)",
                     "Median attributable fraction\n% and range of state ranks",
                     "SHR per 1 SD"),
           size = 2.70, colour = MUTED, lineheight = 1.10) +
  
  annotate("text", x = X_MIN + 0.42, y = 1:5, label = 1:5, hjust = 1,
           size = 2.55, colour = MUTED) +
  annotate("text", x = X_MIN + 0.42, y = Y_ND, label = "ND", hjust = 1,
           size = 2.55, colour = MUTED) +
  annotate("text", x = X_MIN - 0.05, y = 3.4, label = "Within-stream position",
           angle = 90, size = 2.65, colour = INK) +
  
  annotate("text", x = X_MIN + 0.55, y = Y_ND + 0.82, hjust = 0, vjust = 0, size = 2.50,
           colour = MUTED,
           label = paste("ND, not detected. The population stream cannot return this",
                         "category.   * persists in multivariable analysis.")) +
  
  annotate("text", x = X_SWATCH + 0.25, y = HEAD_Y, hjust = 0, vjust = 1,
           label = "Pattern", size = 2.95, fontface = "bold", colour = INK) +
  annotate("text", x = X_SWATCH + 0.25, y = SUB_Y, hjust = 0, vjust = 1,
           label = "Genetic order, ties\nby population rank",
           size = 2.70, colour = MUTED, lineheight = 1.10) +
  geom_point(data = strip, aes(x = X_SWATCH, y = y, colour = factor),
             size = 1.9, shape = 15) +
  geom_text(data = strip, aes(x = X_NAME, y = y, label = factor, colour = factor),
            hjust = 0, size = 2.75, fontface = "bold", show.legend = FALSE) +
  geom_text(data = strip, aes(x = X_CLASS, y = y, label = class),
            hjust = 0, size = 2.75, colour = INK) +
  geom_text(data = strip[mech != ""], aes(x = X_MECH, y = y, label = mech),
            hjust = 0, size = 2.40, colour = MUTED) +
  
  scale_colour_manual(values = PAL, guide = "none") +
  scale_y_reverse(limits = c(Y_ND + 1.20, HEAD_Y - 0.10), expand = c(0, 0)) +
  scale_x_continuous(limits = c(X_MIN - 0.35, X_MAX), expand = c(0, 0)) +
  theme_void(base_size = 8) +
  # [19] theme_void() leaves the canvas transparent; a transparent TIFF renders
  # black in viewers that ignore alpha, so white is set in the theme itself
  theme(plot.margin      = margin(2, 1, 2, 2, "mm"),
        plot.background  = element_rect(fill = "white", colour = NA),
        panel.background = element_rect(fill = "white", colour = NA))

# --- [22] guard: characters the default PDF device cannot encode -------------
# grDevices::pdf() defaults to Latin-1. Anything above U+00FF is dropped
# silently, which is how the dagger became a full stop in an earlier proof.
check_latin1 <- function(...) {
  txt <- unlist(list(...))
  bad <- unique(unlist(lapply(txt, function(z)
    grep("[^\\x01-\\xFF]", strsplit(z, "")[[1]], value = TRUE, useBytes = FALSE))))
  bad <- bad[!is.na(bad) & utf8ToInt(paste(bad, collapse = "")) > 255]
  if (length(bad))
    warning("characters outside Latin-1 in figure text: ",
            paste(sprintf("%s (U+%04X)", bad, vapply(bad, utf8ToInt, integer(1))),
                  collapse = ", "),
            "\nThese will be dropped by grDevices::pdf(). ragg renders them correctly.")
  invisible(bad)
}
check_latin1(vals$val, strip$factor, strip$class, strip$mech, levels(nodes$stream))

# --- [15][23] export ---------------------------------------------------------
# ragg renders type far more cleanly than the default TIFF device. PDF is the
# preferable submission format for line art; the TIFF is a fallback.
save_panel <- function(plot, stem, h_mm, pdf = FALSE) {
  if (!requireNamespace("ragg", quietly = TRUE))
    stop("ragg is required. install.packages(\"ragg\")")
  ggsave(file.path(FIG_DIR, paste0(stem, ".tiff")), plot, device = ragg::agg_tiff,
         width = 175, height = h_mm, units = "mm", dpi = 600, background = "white")
  ggsave(file.path(FIG_DIR, paste0(stem, ".png")), plot, device = ragg::agg_png,
         width = 175, height = h_mm, units = "mm", dpi = 300, background = "white")
  if (pdf) {
    warning("the default PDF device is Latin-1; check every non-ASCII glyph in ",
            "the output before using it")
    ggsave(file.path(FIG_DIR, paste0(stem, ".pdf")), plot,
           width = 175, height = h_mm, units = "mm", bg = "white")
  }
}

save_panel(p, "fig1b_three_orderings", 105)

# --- [14] assemble figure 1 --------------------------------------------------
# Panel a is rebuilt as a ggplot object rather than read back in as a raster.
# Reading the exported TIFF and drawing it with annotation_raster resamples an
# already-rasterised image onto a second pixel grid; the mismatch between grids
# is what breaks up the letterforms in panel a's labels.
#
# This requires 60_fig1a_bloc_selection.R to guard its own exports, so that
# sourcing it here does not re-write panel a:
#
#   if (!exists("EXPORT_PANEL_A")) EXPORT_PANEL_A <- TRUE
#   if (EXPORT_PANEL_A) { ggsave(...) ; ggsave(...) }
#
FIG1A_SRC <- file.path(CODE_DIR, "60_fig1a_bloc_selection.R")

if (!file.exists(FIG1A_SRC)) {
  cand <- unique(unlist(lapply(c(CODE_DIR, ROOT), function(d)
    list.files(d, pattern = "^fig1a.*\\.R$", recursive = TRUE, full.names = TRUE))))
  message("figure 1 not assembled: no file at\n  ", FIG1A_SRC,
          if (length(cand))
            paste0("\nSet FIG1A_SRC to one of:\n  ", paste(cand, collapse = "\n  "))
          else "\nNo fig1a*.R found under CODE_DIR or ROOT.")
  
} else if (!requireNamespace("patchwork", quietly = TRUE)) {
  message("figure 1 not assembled: run install.packages(\"patchwork\")")
  
} else {
  library(patchwork)
  e <- new.env()
  assign("EXPORT_PANEL_A", FALSE, envir = e)
  # [20] panel a's script may use relative paths, so source it from its own
  # directory and restore the working directory afterwards
  owd <- getwd(); on.exit(setwd(owd), add = TRUE)
  setwd(dirname(FIG1A_SRC))
  source(FIG1A_SRC, local = e, encoding = "UTF-8")
  setwd(owd)
  
  if (!exists("p", envir = e))
    stop("60_fig1a_bloc_selection.R did not leave its plot in an object named p")
  
  fig1 <- (get("p", envir = e) / p) + plot_layout(heights = c(110, 105)) +
    plot_annotation(
      tag_levels = "A",
      theme = theme(plot.background = element_rect(fill = "white", colour = NA))) &
    theme(plot.tag = element_text(face = "bold", size = 12))
  save_panel(fig1, "figure1_combined", 215)
}