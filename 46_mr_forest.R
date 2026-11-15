# =============================================================================
# Reproducibility repository - Wang et al., transportability-aware triangulation
# of metabolic stroke evidence across the Shanghai Cooperation Organisation.
# Paths are resolved by config.R; see README.md.
# =============================================================================
if (!exists("RAWDAT")) {
  if (file.exists("config.R")) source("config.R") else
    stop("Run this script from the code/ directory (the folder containing config.R).")
}

###############################################################################
# Figure 5: 三族裔 MR 森林图（IVW）— v3
# 风格与 Figure 4 一致：表格内嵌绘图列 + 三线表边框
#
# 设计要点：
#   - 三族裔暴露单位不同，OR 不可跨族裔比较 → 各自独立坐标轴、独立面板
#   - 坐标轴按点估计范围设定；超出量程的 95% CI 自动截断并加箭头
#   - 祖系名作为面板标题居中；坐标轴标签与方向箭头仅在最底部面板出现
#   - 列宽用不间断空格(U+00A0)补齐，保证三面板列位置一致（ASCII 空格会被裁掉）
###############################################################################

library(forestploter)
library(grid)

# ---- 出图: 一律走 ragg -------------------------------------------------------
# cairo_pdf / pdf() 会把每个字形的步进四舍五入到整数 point (0.353 mm) 或不嵌
# 字体, 前者导致同一单词内字距不匀 (最坏 1 pt), 后者栅格化时字面被替换。
# ragg 在目标分辨率的像素网格上排字, 600 dpi 下取整误差 <= 0.042 mm, 且全部
# 图形出自同一引擎、同一份 Helvetica, 与其余图观感一致。
if (!requireNamespace("ragg", quietly = TRUE)) install.packages("ragg")

save_panel <- function(p, stem, w_mm, h_mm) {
  for (r in c(600, 300)) {
    f <- if (r == 600) paste0(stem, ".tif") else paste0(stem, "_300dpi.tif")
    ragg::agg_tiff(f, width = w_mm, height = h_mm, units = "mm",
                   res = r, compression = "lzw", background = "white",
                   scaling = 1)
    print(p); dev.off()
    cat(sprintf("已保存: %s (%.1f x %.1f mm @ %d dpi)\n", f, w_mm, h_mm, r))
  }
}


# 小数点符号：最终统一为 Lancet 中点时改成 "\u00B7"
DEC <- "\u00B7"
dec <- function(x) if (DEC == ".") x else gsub(".", DEC, x, fixed = TRUE)

# 递归把 gtable 内所有文字标签的小数点替换为 DEC（用于坐标轴刻度）
dotify <- function(x) {
  if (inherits(x, "gtable")) {
    x$grobs <- lapply(x$grobs, dotify)
  } else if (inherits(x, "gTree") && length(x$children)) {
    for (nm in names(x$children)) x$children[[nm]] <- dotify(x$children[[nm]])
  } else if (inherits(x, "text") && !is.null(x$label)) {
    x$label <- gsub(".", DEC, as.character(x$label), fixed = TRUE)
  }
  x
}

# ============================================================================
# 绘图函数
# ============================================================================
make_mr_forest <- function(df,
                           xlim, ticks_at,
                           title        = NULL,
                           outcome_lab  = "Outcome",
                           term_lab     = "Exposure",
                           n_lab        = "nSNP",
                           est_lab      = "OR",
                           ci_lab       = "95% CI",
                           p_lab        = "p",
                           fdr_lab      = "FDR p",
                           xlab         = NULL,
                           arrow_lab    = NULL,
                           ci_colour    = "#2C5985",
                           base_size    = 8,
                           plot_width   = 60,
                           est_pad      = 2,
                           num_pad      = 2.5,   # OR/CI/p/FDR 之间的列间距（字符数）
                           group_width  = 26,   # 结局列字符宽度
                           p_width      = 7) {  # p / FDR p 列宽（"<0.0001" 为 7 字符）
  
  # 用不间断空格右补位：ASCII 空格在渲染时会被裁掉，导致各面板列宽不一致
  padr <- function(x, w) {
    n <- pmax(0, w - nchar(x))
    paste0(x, strrep("\u00A0", n))
  }
  pad <- strrep("\u00A0", est_pad)
  npad <- strrep("\u00A0", num_pad)
  
  est_txt <- ifelse(is.na(df$est), "", dec(sprintf("%.3f", df$est)))
  ci_txt  <- ifelse(is.na(df$est), "",
                    dec(sprintf("%.3f\u2013%.3f", df$low, df$hi)))
  
  tab <- data.frame(
    g     = padr(df$group, group_width),
    t     = df$term,
    n     = ifelse(is.na(df$n), "", as.character(df$n)),
    blank = strrep(" ", plot_width),
    e     = paste0(pad,  est_txt),
    c     = paste0(npad, ci_txt),
    p     = paste0(npad, padr(dec(ifelse(is.na(df$pval), "", df$pval)), p_width)),
    f     = paste0(npad, padr(dec(ifelse(is.na(df$fdr),  "", df$fdr)),  p_width)),
    stringsAsFactors = FALSE, check.names = FALSE
  )
  names(tab) <- c(padr(outcome_lab, group_width), term_lab, n_lab, " ",
                  paste0(pad,  est_lab),
                  paste0(npad, ci_lab),
                  paste0(npad, padr(p_lab,   p_width)),
                  paste0(npad, padr(fdr_lab, p_width)))
  
  tm <- forest_theme(
    base_size        = base_size,
    ci_pch           = 15,
    ci_col           = ci_colour,
    ci_lwd           = 1.1,
    ci_Theight       = 0.18,
    refline_gp       = gpar(lwd = 0.8, lty = "solid", col = "grey45"),
    xaxis_gp         = gpar(fontsize = base_size - 0.5, cex = 1),
    title_gp         = gpar(fontsize = base_size + 1, fontface = "bold"),
    footnote_gp      = gpar(fontsize = base_size - 1, col = "grey25"),
    arrow_type       = "closed",
    arrow_label_just = "start",
    arrow_gp         = gpar(lwd = 0.8, col = "grey25", fontsize = base_size - 1),
    core    = list(bg_params = list(fill = "white"),
                   padding    = unit(c(2.5, 1.8), "mm")),
    colhead = list(padding    = unit(c(2.5, 1.8), "mm")),
  )
  
  p <- forest(
    tab,
    est       = df$est,
    lower     = df$low,
    upper     = df$hi,
    sizes     = 0.4,
    ci_column = 4,
    ref_line  = 1,
    x_trans   = "log",
    xlim      = xlim,
    ticks_at  = ticks_at,
    title     = title,
    xlab      = xlab,
    arrow_lab = arrow_lab,
    theme     = tm
  )
  
  # 坐标轴刻度标签由 forestploter 内部生成，无法通过参数改小数点，
  # 这里递归遍历 gtable 把所有 text grob 里的 "." 换成 DEC
  if (DEC != ".") p <- dotify(p)
  
  p <- add_border(p, part = "header", where = "top",    gp = gpar(lwd = 1.2))
  p <- add_border(p, part = "header", where = "bottom", gp = gpar(lwd = 0.7))
  p <- add_border(p, row = nrow(tab), where = "bottom", gp = gpar(lwd = 1.2))
  p
}

# ============================================================================
# 数据
# ============================================================================
blank_row <- data.frame(group = "", term = "", n = NA, est = NA, low = NA,
                        hi = NA, pval = NA, fdr = NA, stringsAsFactors = FALSE)

mk <- function(outcome, term, n, est, low, hi, pval, fdr) {
  data.frame(group = c(outcome, rep("", length(term) - 1)),
             term = term, n = n, est = est, low = low, hi = hi,
             pval = pval, fdr = fdr, stringsAsFactors = FALSE)
}

exp5 <- c("SBP", "LDL", "BMI", "FG", "eGFR")

# ---- A. European ----
eur <- rbind(
  mk("Ischaemic stroke", exp5, c(20, 62, 37, 85, 41),
     c(1.047, 1.093, 1.030, 1.044, 1.001),
     c(1.033, 1.034, 1.008, 0.941, 0.995),
     c(1.061, 1.155, 1.054, 1.159, 1.006),
     c("<0.0001", "0.0017", "0.0082", "0.42", "0.86"),
     c("<0.0001", "0.0042", "0.01", "0.52", "0.86")),
  blank_row,
  mk("Intracerebral haemorrhage", exp5, c(20, 62, 37, 85, 41),
     c(1.034, 0.989, 1.022, 0.806, 0.998),
     c(0.998, 0.852, 0.958, 0.603, 0.983),
     c(1.070, 1.148, 1.091, 1.077, 1.014),
     c("0.06", "0.88", "0.51", "0.14", "0.81"),
     c("0.31", "0.88", "0.84", "0.36", "0.88")),
  blank_row,
  mk("Subarachnoid haemorrhage", exp5, c(20, 62, 37, 84, 41),
     c(1.045, 0.872, 1.022, 0.636, 1.015),
     c(1.004, 0.760, 0.958, 0.447, 0.996),
     c(1.089, 1.002, 1.090, 0.906, 1.034),
     c("0.03", "0.05", "0.51", "0.01", "0.11"),
     c("0.08", "0.09", "0.51", "0.06", "0.14"))
)

# ---- B. East Asian ----
eas <- rbind(
  mk("Ischaemic stroke", exp5, c(32, 61, 40, 31, 53),
     c(2.939, 1.106, 1.204, 0.974, 0.976),
     c(2.138, 0.987, 0.993, 0.831, 0.664),
     c(4.039, 1.240, 1.459, 1.142, 1.434),
     c("<0.0001", "0.08", "0.06", "0.75", "0.90"),
     c("<0.0001", "0.14", "0.14", "0.90", "0.90")),
  blank_row,
  mk("Intracerebral haemorrhage", exp5, c(32, 61, 40, 31, 54),
     c(4.767, 0.907, 1.378, 0.956, 0.971),
     c(2.986, 0.784, 1.045, 0.741, 0.537),
     c(7.608, 1.049, 1.817, 1.234, 1.757),
     c("<0.0001", "0.19", "0.02", "0.73", "0.92"),
     c("<0.0001", "0.32", "0.06", "0.91", "0.92")),
  blank_row,
  mk("Subarachnoid haemorrhage", exp5, c(32, 60, 40, 31, 54),
     c(4.811, 1.010, 2.100, 0.572, 0.811),
     c(1.960, 0.685, 0.942, 0.348, 0.251),
     c(11.812, 1.488, 4.679, 0.941, 2.618),
     c("0.0006", "0.96", "0.07", "0.03", "0.73"),
     c("0.0030", "0.96", "0.12", "0.07", "0.91"))
)

# ---- C. South Asian ----
sas <- mk("Ischaemic stroke", exp5, c(10, 5, 12, 13, 6),
          c(1.362, 0.842, 1.390, 1.098, 1.087),
          c(0.920, 0.536, 0.960, 0.805, 0.626),
          c(2.016, 1.322, 2.012, 1.497, 1.888),
          c("0.12", "0.45", "0.08", "0.56", "0.77"),
          c("0.31", "0.70", "0.31", "0.70", "0.77"))

# ============================================================================
# 绘制
#   量程按点估计范围设定；超出量程的 CI 端点被截断并画箭头
#   European   : 点估计 0.636–1.093（FG→SAH 下限 0.447 截断）
#   East Asian : 点估计 0.572–4.811（SBP→ICH/SAH 上限、FG/eGFR→SAH 下限截断）
#   South Asian: 点估计 0.842–1.390（无截断）
# ============================================================================
p_eur <- make_mr_forest(eur,
                        xlim     = c(0.60, 1.30),
                        ticks_at = c(0.6, 0.8, 1.0, 1.2),
                        title    = "European")

p_eas <- make_mr_forest(eas,
                        xlim     = c(0.50, 6.0),
                        ticks_at = c(0.5, 1, 2, 4),
                        title    = "East Asian")

p_sas <- make_mr_forest(sas,
                        xlim      = c(0.50, 2.1),
                        ticks_at  = c(0.5, 0.75, 1, 1.5, 2),
                        title     = "South Asian",
                        xlab      = "OR (log scale)",
                        arrow_lab = c("Lower risk", "Higher risk"))

# 强制三面板列宽一致（以 European 面板为准）
for (nm in c("eas", "sas")) {
  p <- get(paste0("p_", nm))
  if (length(p$widths) == length(p_eur$widths)) {
    p$widths <- p_eur$widths
    assign(paste0("p_", nm), p)
  } else {
    warning(nm, " 的列数与 eur 不一致，未强制对齐")
  }
}

# ============================================================================
# 保存（三面板分别出图，之后拼接为 A/B/C）
# ============================================================================
setwd(file.path(RAWDAT, "MR"))
for (nm in c("eur", "eas", "sas")) {
  p  <- get(paste0("p_", nm))
  wh <- get_wh(p, unit = "in")
  cat(sprintf("%s: %.1f x %.1f mm\n", nm, wh[1] * 25.4, wh[2] * 25.4))
  for (r in c(600, 300)) {
    f <- sprintf(if (r == 600) "Figure5_%s.tif" else "Figure5_%s_300dpi.tif", nm)
    ragg::agg_tiff(f, width = wh[1], height = wh[2], units = "in",
                   res = r, compression = "lzw", background = "white",
                   scaling = 1)
    plot(p); dev.off()
    cat("已保存:", f, "\n")
  }
}
cat("已保存: Figure5_eur / Figure5_eas / Figure5_sas\n")