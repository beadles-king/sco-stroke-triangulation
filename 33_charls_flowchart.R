# =============================================================================
# Reproducibility repository - Wang et al., transportability-aware triangulation
# of metabolic stroke evidence across the Shanghai Cooperation Organisation.
# Paths are resolved by config.R; see README.md.
# =============================================================================
if (!exists("RAWDAT")) {
  if (file.exists("config.R")) source("config.R") else
    stop("Run this script from the code/ directory (the folder containing config.R).")
}

# ============================================================================
# Figure S17：CHARLS 纳入/排除流程图（STROBE 风格）
#
# 设计说明
#   - 全图由下面的 FLOW 定义驱动：每个框就是一串"行"，每行 = 左侧文字 +
#     可选的右对齐数字。改文字、改数字、增删行都只动这一段，版面（框高、
#     箭头、总画布高度）由脚本自动重算，不需要手工调坐标。
#   - 单位一律 mm。画布宽固定 180 mm（供应图上限），高度由内容算出后
#     断言 <= 240 mm。
#   - 出图走 ragg::agg_tiff，600 dpi 与 300 dpi 各渲染一次（禁用
#     cairo_pdf / pdf() / ggsave / png()，见交接文档 v3 陷阱 1）。
#   - 字体 Helvetica，正文 7 pt，框标题 7 pt 粗体。
#
# 数值来源
#   主链           20_charls_prepare_data.R 的逐步计数
#   主分析分支     ledger（Models 1-2: 7896/565；Model 3: 7710/556）
#   两波均值分支   sens2_flow.csv + sens2_n_ledger.csv
#   胱抑素 C 分支  Sens_eGFR_n_ledger  ← 441 / 435 两个事件数尚待核验，
#                  见交接文档 v6 §7 第 2 条，核对后若有变动只改 FLOW 内的字符串
# ============================================================================

library(ggplot2)
library(ragg)

OUT_STEM <- "FigS17_flowchart"
OUT_DIR  <- file.path(RAWDAT, "CHARLS")

# ---- 版面常数（mm / pt）----------------------------------------------------
W_TOTAL  <- 180      # 画布宽，供应图上限
H_MAX    <- 240      # 画布高上限
FS       <- 7        # 正文字号 pt
FS_HEAD  <- 7        # 框标题字号 pt
LINE_H   <- 3.4      # 行高 mm
PAD_X    <- 2.4      # 框内左右留白 mm
PAD_Y    <- 2.0      # 框内上下留白 mm
GAP_V    <- 3.5      # 纵向框间距 mm
GAP_BR   <- 7.0      # 合格队列框与分支之间的间距 mm
MARGIN_T <- 3.0      # 画布上边距 mm
MARGIN_B <- 3.0      # 画布下边距 mm

MAIN_X   <- 2;   MAIN_W <- 92     # 主链框
EXCL_X   <- 102; EXCL_W <- 76     # 排除框
COL_N    <- 3
COL_GAP  <- 7
COL_X0   <- 2
COL_W    <- (W_TOTAL - 2 * COL_X0 - (COL_N - 1) * COL_GAP) / COL_N
COL_X    <- COL_X0 + (seq_len(COL_N) - 1) * (COL_W + COL_GAP)

# ---- 配色（克制，与全文聚类稳定性系列同源）--------------------------------
COL_MAIN   <- "#2C5985"
COL_LINE   <- "#4D4D4D"
FILL_PLAIN <- "#FFFFFF"
FILL_EXCL  <- "#FAFAFA"
FILL_KEY   <- "#EAF0F6"
FILL_HEAD  <- "#F2F2F2"

pt2mm <- function(pt) pt / .pt

# ============================================================================
# 内容定义 —— 需要改字改数，只动这一段
# ============================================================================
ln <- function(text, n = NA_character_, face = "plain")
  list(text = text, n = n, face = face)

MAIN <- list(
  A = list(ln("Harmonised CHARLS, wave 1 (2011)", "n=25 586", "bold")),
  B = list(ln("Alive and stroke-free with follow-up", "n=16 378", "bold")),
  C = list(ln("Complete anthropometric and age data", "n=12 326", "bold")),
  D = list(ln("Valid blood-based measures", "n=8061", "bold")),
  E = list(ln("Eligible cohort", "n=8012", "bold"),
           ln("573 incident strokes; 682 competing deaths; 6757 event-free"))
)

EXCL <- list(
  a = list(ln("Excluded", "n=9208", "bold"),
           ln("Not alive and responding at wave 1, prevalent"),
           ln("stroke, or lost to follow-up after wave 2 and"),
           ln("never re-contacted")),
  b = list(ln("Excluded", "n=4052", "bold"),
           ln("Missing blood pressure",  "3370"),
           ln("Missing body-mass index",  "201"),
           ln("Missing age",               "59"),
           ln("Aged <45 years",           "422")),
  c = list(ln("Excluded", "n=4265", "bold"),
           ln("No wave-1 venous blood subsample",         "3209"),
           ln("Non-fasting blood sample",                  "872"),
           ln("Missing fasting glucose, LDL, or creatinine","184")),
  d = list(ln("Excluded", "n=49", "bold"),
           ln("Physiologically implausible value for any of"),
           ln("the five exposures"))
)

# 分支顺序：主分析（左）→ 两波均值（中）→ 胱抑素 C（右）
BRANCH <- list(
  list(head   = list(ln("PRIMARY ANALYSIS &", face = "bold"),
                     ln("SENSITIVITY Cause-specific Cox", face = "bold")),
       excl   = list(ln("No further restriction", face = "bold"),
                     ln("Cause-specific Cox models use"),
                     ln("the same samples")),
       excl_fill = FILL_PLAIN,
       sample = list(ln("All eligible participants", "n=8012", "bold")),
       m12    = list(ln("Models 1\u20132", "n=7896; 565 events", "bold")),
       m3     = list(ln("Model 3", "n=7710; 556 events", "bold"))),
  
  list(head   = list(ln("SENSITIVITY", face = "bold"),
                     ln("Two-wave mean exposures", face = "bold")),
       excl   = list(ln("Excluded", "n=3728", "bold"),
                     ln("No stroke-free wave-3 interview", "939"),
                     ln("No wave-4 follow-up",             "368"),
                     ln("Missing wave-3 SBP",              "849"),
                     ln("Missing wave-3 BMI",               "47"),
                     ln("No fasting wave-3 blood data",   "1500"),
                     ln("Implausible wave-3 value",         "25")),
       sample = list(ln("Two-wave mean sample", "n=4284", "bold"),
                     ln("265 incident strokes; 176 competing deaths")),
       m12    = list(ln("Models 1\u20132", "n=4233; 263 events", "bold")),
       m3     = list(ln("Model 3", "n=4134; 257 events", "bold"))),
  
  list(head   = list(ln("SENSITIVITY", face = "bold"),
                     ln("Cystatin C", face = "bold")),
       excl   = list(ln("Excluded", "n=1903", "bold"),
                     ln("No cystatin C measurement", "1902"),
                     ln("Implausible cystatin C eGFR",  "1")),
       sample = list(ln("Cystatin C subsample", "n=6109", "bold")),
       m12    = list(ln("Models 1\u20132", "n=6014; 441 events", "bold")),
       m3     = list(ln("Model 3", "n=5876; 435 events", "bold")))
)

# ---- 算术自检：任何一步对不上就停，不出图 ----------------------------------
chk <- function(prev, excluded, nxt, what) {
  if (prev - excluded != nxt)
    stop(sprintf("%s 对不上：%d - %d = %d，但下一步写的是 %d",
                 what, prev, excluded, prev - excluded, nxt))
  invisible(TRUE)
}
chk(25586, 9208, 16378, "主链第 1 步")
chk(16378, 4052, 12326, "主链第 2 步")
chk(12326, 4265,  8061, "主链第 3 步")
chk( 8061,   49,  8012, "主链第 4 步")
stopifnot(3370 + 201 + 59 + 422 == 4052)
stopifnot(3209 + 872 + 184 == 4265)
stopifnot(573 + 682 + 6757 == 8012)   # Eligible cohort 框内那句话的三个数
chk(8012, 3728, 4284, "两波均值分支")
stopifnot(939 + 368 + 849 + 47 + 1500 + 25 == 3728)
chk(8012, 1903, 6109, "胱抑素 C 分支")
stopifnot(1902 + 1 == 1903)
cat("算术自检通过\n")

# ============================================================================
# 版面引擎
# ============================================================================
rects <- list(); texts <- list(); segs <- list()

box_h <- function(lines) length(lines) * LINE_H + 2 * PAD_Y

# top 为距画布顶端的距离（mm），内部一律用"向下为正"，绘图时再翻转
add_box <- function(lines, x, w, top, fill = FILL_PLAIN,
                    border = COL_LINE, lwd = 0.4) {
  h <- box_h(lines)
  rects[[length(rects) + 1]] <<- data.frame(
    xmin = x, xmax = x + w, ytop = top, ybot = top + h,
    fill = fill, border = border, lwd = lwd, stringsAsFactors = FALSE)
  for (i in seq_along(lines)) {
    L  <- lines[[i]]
    yc <- top + PAD_Y + (i - 0.5) * LINE_H
    texts[[length(texts) + 1]] <<- data.frame(
      x = x + PAD_X, y = yc, label = L$text, hjust = 0,
      face = L$face, stringsAsFactors = FALSE)
    if (!is.na(L$n))
      texts[[length(texts) + 1]] <<- data.frame(
        x = x + w - PAD_X, y = yc, label = L$n, hjust = 1,
        face = L$face, stringsAsFactors = FALSE)
  }
  invisible(top + h)
}

add_seg <- function(x1, y1, x2, y2, arrow = TRUE) {
  segs[[length(segs) + 1]] <<- data.frame(
    x = x1, y = y1, xend = x2, yend = y2, arrow = arrow)
  invisible(NULL)
}

# ---- 主链 ------------------------------------------------------------------
mid_x <- MAIN_X + MAIN_W / 2
cur   <- MARGIN_T
main_ids <- c("A", "B", "C", "D")
excl_ids <- c("a", "b", "c", "d")
tops <- numeric(0)

for (k in seq_along(main_ids)) {
  bx  <- MAIN[[main_ids[k]]]
  bot <- add_box(bx, MAIN_X, MAIN_W, cur, fill = FILL_PLAIN)
  tops <- c(tops, cur)
  
  ex   <- EXCL[[excl_ids[k]]]
  eh   <- box_h(ex)
  etop <- bot + GAP_V
  add_box(ex, EXCL_X, EXCL_W, etop, fill = FILL_EXCL, lwd = 0.3)
  
  nxt <- max(etop + eh + GAP_V, bot + 14)
  add_seg(mid_x, bot, mid_x, nxt)                       # 竖直主干
  add_seg(mid_x, etop + eh / 2, EXCL_X, etop + eh / 2)  # 横向指向排除框
  cur <- nxt
}

e_top <- cur
e_bot <- add_box(MAIN$E, MAIN_X, MAIN_W, e_top, fill = FILL_KEY,
                 border = COL_MAIN, lwd = 0.8)

# ---- 分支 ------------------------------------------------------------------
head_top <- e_bot + GAP_BR + 4          # 4 mm 留给分叉横线
head_h   <- max(vapply(BRANCH, function(b) box_h(b$head), 0))
excl_h   <- vapply(BRANCH, function(b) if (is.null(b$excl)) 0 else box_h(b$excl), 0)
samp_h   <- vapply(BRANCH, function(b) box_h(b$sample), 0)

# 样本框底边对齐：中间区域高度取三列所需的最大值
mid_need <- max(excl_h + ifelse(excl_h > 0, GAP_V, 0) + samp_h)
samp_bot <- head_top + head_h + GAP_V + mid_need
m12_top  <- samp_bot + GAP_V
m12_h    <- max(vapply(BRANCH, function(b) box_h(b$m12), 0))
m3_top   <- m12_top + m12_h + GAP_V
m3_h     <- max(vapply(BRANCH, function(b) box_h(b$m3), 0))
H_TOTAL  <- m3_top + m3_h + MARGIN_B

if (H_TOTAL > H_MAX)
  stop(sprintf("画布高 %.1f mm 超过 %d mm 上限，请精简文字或缩小 LINE_H", H_TOTAL, H_MAX))

fork_y <- e_bot + GAP_BR / 2
add_seg(mid_x, e_bot, mid_x, fork_y, arrow = FALSE)

for (j in seq_along(BRANCH)) {
  b  <- BRANCH[[j]]
  cx <- COL_X[j] + COL_W / 2
  
  add_seg(mid_x, fork_y, cx, fork_y, arrow = FALSE)
  add_seg(cx, fork_y, cx, head_top)
  
  add_box(b$head, COL_X[j], COL_W, head_top, fill = FILL_HEAD, lwd = 0.4)
  
  s_top <- samp_bot - box_h(b$sample)
  if (!is.null(b$excl)) {
    etop <- head_top + head_h + GAP_V
    add_box(b$excl, COL_X[j], COL_W, etop,
            fill = if (is.null(b$excl_fill)) FILL_EXCL else b$excl_fill, lwd = 0.3)
    add_seg(cx, head_top + head_h, cx, etop, arrow = FALSE)
    add_seg(cx, etop + box_h(b$excl), cx, s_top)
  } else {
    add_seg(cx, head_top + head_h, cx, s_top)
  }
  add_box(b$sample, COL_X[j], COL_W, s_top, fill = FILL_KEY,
          border = COL_MAIN, lwd = 0.6)
  
  add_seg(cx, samp_bot, cx, m12_top)
  add_box(b$m12, COL_X[j], COL_W, m12_top, fill = FILL_PLAIN)
  add_seg(cx, m12_top + box_h(b$m12), cx, m3_top)
  add_box(b$m3,  COL_X[j], COL_W, m3_top,  fill = FILL_PLAIN)
}

# ============================================================================
# 绘制（内部 y 向下为正，此处翻转成绘图坐标）
# ============================================================================
R <- do.call(rbind, rects); T <- do.call(rbind, texts); S <- do.call(rbind, segs)
R$ymin <- H_TOTAL - R$ybot; R$ymax <- H_TOTAL - R$ytop
T$yy   <- H_TOTAL - T$y
S$y1   <- H_TOTAL - S$y;    S$y2   <- H_TOTAL - S$yend

p <- ggplot() +
  geom_segment(data = subset(S, !arrow),
               aes(x = x, y = y1, xend = xend, yend = y2),
               colour = COL_LINE, linewidth = 0.25) +
  geom_segment(data = subset(S, arrow),
               aes(x = x, y = y1, xend = xend, yend = y2),
               colour = COL_LINE, linewidth = 0.25,
               arrow = arrow(length = unit(1.4, "mm"), type = "closed")) +
  geom_rect(data = R,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = R$fill, colour = R$border, linewidth = R$lwd) +
  geom_text(data = T,
            aes(x = x, y = yy, label = label, hjust = hjust),
            vjust = 0.5, size = pt2mm(FS), family = "Helvetica",
            fontface = T$face, colour = "black") +
  coord_cartesian(xlim = c(0, W_TOTAL), ylim = c(0, H_TOTAL), expand = FALSE) +
  theme_void()

# ---- 出图：600 dpi 与 300 dpi 各渲染一次 -----------------------------------
save_panel <- function(plot, stem, w_mm, h_mm) {
  for (dpi in c(600, 300)) {
    f <- file.path(OUT_DIR,
                   if (dpi == 600) paste0(stem, ".tif") else paste0(stem, "_300dpi.tif"))
    agg_tiff(f, width = w_mm, height = h_mm, units = "mm",
             res = dpi, background = "white", compression = "lzw")
    print(plot); dev.off()
    cat(sprintf("已写出 %s（%.1f x %.1f mm, %d dpi）\n", f, w_mm, h_mm, dpi))
  }
}

save_panel(p, OUT_STEM, W_TOTAL, H_TOTAL)
cat(sprintf("\n最终画布：%.1f x %.1f mm（上限 180 x 240）\n", W_TOTAL, H_TOTAL))