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
# Figure 4: CHARLS 五项代谢暴露与卒中发病的 Fine-Gray SHR
# 森林图通用函数 v4
#
# 相对 v3 的唯一结构性改动：
#   数值不再写死在脚本里，改为从 24_charls_fine_gray.R 输出的 FG_models1to3.csv 读取。
#   v3 把 SHR/CI/p 硬编码，主分析重跑后图不会跟着更新（哑变量化后
#   Model 2 的 FG CI 下界 1·03→1·04、BMI CI 下界 1·15→1·14，
#   Model 3 的 SBP 1·26→1·27、BMI 1·16→1·15 等都不会反映到图上）。
#
#   p 值标签也改为按全文规则由原始 p 计算，不再手工键入：
#     p < 0·0001            -> "<0·0001"
#     0·0001 <= p < 0·01    -> 四位小数
#     p >= 0·01             -> 两位小数
#
# 依赖：forestploter, grid, ragg
# 前置：先跑 24_charls_fine_gray.R（生成 FG_models1to3.csv）
###############################################################################

library(forestploter)
library(grid)

RES_DIR  <- file.path(RAWDAT, "CHARLS")
RES_FILE <- file.path(RES_DIR, "FG_models1to3.csv")

DEC <- "\u00B7"

# ============================================================================
# 格式化工具（与 25_charls_rcs.R / 21_charls_model_utils.R 的规则一致）
# ============================================================================
dec <- function(x) gsub("\\.", DEC, as.character(x))

fmt_p <- function(p) {
  vapply(p, function(pi) {
    if (is.na(pi))   return(NA_character_)
    if (pi < 0.0001) return(paste0("<0", DEC, "0001"))
    if (pi < 0.01)   return(dec(sprintf("%.4f", pi)))
    dec(sprintf("%.2f", pi))
  }, character(1))
}

# 递归把 gtable 内所有 text grob 的小数点换成中点（坐标轴刻度的唯一途径）
dotify <- function(g) {
  if (!is.null(g$children)) {
    ch <- lapply(g$children, dotify)
    class(ch) <- class(g$children)          # 保住 gList 类
    g$children <- ch
  }
  if (!is.null(g$grobs)) g$grobs <- lapply(g$grobs, dotify)
  if (inherits(g, "text") && !is.null(g$label)) {
    g$label <- gsub("(?<=\\d)\\.(?=\\d)", DEC, as.character(g$label), perl = TRUE)
  }
  g
}

# ============================================================================
# 通用绘图函数（与 v3 完全一致，未改）
# ============================================================================
make_forest <- function(df,
                        group_lab   = "Model",
                        term_lab    = "Exposure",
                        est_lab     = "SHR",
                        ci_lab      = "95% CI",
                        p_lab       = "p value",
                        ratio       = TRUE,
                        xlim        = NULL,
                        ticks_at    = NULL,
                        xlab        = NULL,
                        arrow_lab   = NULL,
                        footnote    = NULL,
                        ci_colour   = "#2C5985",
                        base_size   = 8,
                        plot_width  = 42,
                        ticks_digits = 1,
                        est_pad     = 2,
                        three_line  = TRUE) {
  
  fmt <- function(x, digits = 2) ifelse(is.na(x), "",
                                        gsub("\\.", DEC,
                                             sprintf(paste0("%.", digits, "f"), x)))
  pad <- strrep(" ", est_pad)
  
  tab <- data.frame(
    g     = df$group,
    t     = df$term,
    blank = strrep(" ", plot_width),
    e     = paste0(pad, fmt(df$est)),
    c     = ifelse(is.na(df$low), "",
                   paste0(fmt(df$low), "\u2013", fmt(df$hi))),
    p     = ifelse(is.na(df$pval), "", df$pval),
    stringsAsFactors = FALSE, check.names = FALSE
  )
  names(tab) <- c(group_lab, term_lab, " ",
                  paste0(pad, est_lab), ci_lab, p_lab)
  
  ref <- if (ratio) 1 else 0
  if (is.null(xlab)) xlab <- if (ratio) paste0(est_lab, " per 1 SD") else "Effect estimate"
  
  tm <- forest_theme(
    base_size         = base_size,
    ci_pch            = 15,
    ci_col            = ci_colour,
    ci_lwd            = 1.1,
    ci_Theight        = 0.18,
    refline_gp        = gpar(lwd = 0.8, lty = "solid", col = "grey45"),
    xaxis_gp          = gpar(fontsize = base_size - 0.5, cex = 1),
    footnote_gp       = gpar(fontsize = base_size - 1, col = "grey25"),
    arrow_type        = "closed",
    arrow_label_just  = "start",
    arrow_gp          = gpar(lwd = 0.8, col = "grey25", fontsize = base_size - 1),
    core              = list(bg_params = list(fill = "white")),
    padding           = unit(c(2.5, 3), "mm")
  )
  
  p <- forest(
    tab,
    est        = df$est,
    lower      = df$low,
    upper      = df$hi,
    sizes      = 0.4,
    ci_column  = 3,
    ref_line   = ref,
    x_trans    = if (ratio) "log" else "none",
    xlim       = xlim,
    ticks_at   = ticks_at,
    ticks_digits = ticks_digits,
    xlab       = xlab,
    arrow_lab  = arrow_lab,
    footnote   = footnote,
    theme      = tm
  )
  
  if (three_line) {
    p <- add_border(p, part = "header", where = "top",    gp = gpar(lwd = 1.2))
    p <- add_border(p, part = "header", where = "bottom", gp = gpar(lwd = 0.7))
    p <- add_border(p, row = nrow(tab), where = "bottom", gp = gpar(lwd = 1.2))
  }
  dotify(p)
}

# ============================================================================
# Figure 4 数据：从主分析输出读取
# ============================================================================
stopifnot(file.exists(RES_FILE))
res <- read.csv(RES_FILE, stringsAsFactors = FALSE)

# 展示顺序：按效应强度 SBP、BMI、FG、LDL、eGFR
# （与 Table S15–S18、Figure S19 一致）
plot_order <- c("SBP_sd", "BMI_sd", "FG_sd", "LDL_sd", "eGFR_sd")
term_label <- c(SBP_sd = "SBP", BMI_sd = "BMI", FG_sd = "FG",
                LDL_sd = "LDL", eGFR_sd = "eGFR")

blank_row <- data.frame(group = "", term = "", est = NA, low = NA, hi = NA,
                        pval = NA, stringsAsFactors = FALSE)

build_block <- function(model_name) {
  d <- res[res$model == model_name, ]
  d <- d[match(plot_order, d$variable), ]
  stopifnot(!any(is.na(d$SHR)))
  data.frame(
    group = c(model_name, rep("", nrow(d) - 1)),
    term  = unname(term_label[d$variable]),
    est   = d$SHR,
    low   = d$CI_low,
    hi    = d$CI_high,
    pval  = fmt_p(d$P),
    stringsAsFactors = FALSE
  )
}

models <- c("Model 1", "Model 2", "Model 3")
fig4 <- build_block(models[1])
for (m in models[-1]) fig4 <- rbind(fig4, blank_row, build_block(m))

# ---- 核对：把即将画上图的数值打印出来，便于与表格逐行比对 ----
cat("\n===== Figure 4 绘图数值（应与 Table S15 完全一致）=====\n")
chk <- fig4[fig4$term != "", ]
chk$display <- sprintf("%s (%s\u2013%s)",
                       gsub("\\.", DEC, sprintf("%.2f", chk$est)),
                       gsub("\\.", DEC, sprintf("%.2f", chk$low)),
                       gsub("\\.", DEC, sprintf("%.2f", chk$hi)))
print(chk[, c("group", "term", "display", "pval")], row.names = FALSE)

# 各模型分析人数（写进图注）
cat("\n===== 各模型分析人数（Figure 4 图注需引用）=====\n")
print(unique(res[, c("model", "n", "n_stroke")]), row.names = FALSE)

# 轴范围自适应核查：确保没有 CI 被裁切
XLIM <- c(0.85, 1.60)
if (min(chk$low) < XLIM[1] || max(chk$hi) > XLIM[2]) {
  warning(sprintf("有 CI 超出 xlim：实际范围 %.3f-%.3f，当前 xlim %.2f-%.2f",
                  min(chk$low), max(chk$hi), XLIM[1], XLIM[2]), call. = FALSE)
}

# ============================================================================
# 绘制并保存
# ============================================================================
p4 <- make_forest(
  fig4,
  est_lab    = "SHR",
  xlim       = XLIM,
  ticks_at   = seq(0.9, 1.5, by = 0.1),
  xlab       = "SHR per 1 SD (log scale)",
  arrow_lab  = c("Lower risk", "Higher risk"),
  plot_width = 110,
  est_pad    = 2
)

wh <- get_wh(p4, unit = "in")
cat(sprintf("\n图面尺寸: %.1f x %.1f mm\n", wh[1] * 25.4, wh[2] * 25.4))
if (wh[1] * 25.4 > 180) warning("宽度超过 180 mm 上限，请调小 plot_width。", call. = FALSE)
if (wh[2] * 25.4 > 240) warning("高度超过 240 mm 上限。", call. = FALSE)

setwd(RES_DIR)

# 出图统一走 ragg：cairo_pdf 会把字形步进四舍五入到整数 point 再用 TJ 补偿，
# 同一单词内字距不匀；png() 在 macOS 上走 quartz，是第三种排版引擎。
# ragg 在目标分辨率排字，600 dpi 下取整误差 <= 0.042 mm。
if (!requireNamespace("ragg", quietly = TRUE)) install.packages("ragg")

for (r in c(600, 300)) {
  f <- if (r == 600) "Figure4.tif" else "Figure4_300dpi.tif"
  ragg::agg_tiff(f, width = wh[1], height = wh[2], units = "in",
                 res = r, compression = "lzw", background = "white",
                 scaling = 1)
  plot(p4)
  dev.off()
  cat(sprintf("已保存: %s (%.1f x %.1f mm @ %d dpi)\n",
              f, wh[1] * 25.4, wh[2] * 25.4, r))
}

# 供表格与正文取用的显示字符串
write.csv(chk[, c("group", "term", "est", "low", "hi", "display", "pval")],
          file.path(RES_DIR, "Figure4_display_values.csv"), row.names = FALSE)
cat("\n已写出 Figure4_display_values.csv（Table S15 与正文可直接取用）\n")