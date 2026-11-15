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
# GBD SCO SDI-ASDR 前沿分析 —— 发表版
#
# 相对上一版的改动只有两类，建模与效率计算逻辑一字未动：
#   (1) p 值精度：原 fmt_p 在 p<0.001 处封顶，现按全文规则输出
#       p<0.0001 -> "<0·0001"；0.0001<=p<0.01 -> 四位小数；p>=0.01 -> 两位小数
#   (2) 出图规格：8 pt 字号、180 mm 幅面、中点小数、en dash、真减号、
#       面板字母 A/B 11 pt 粗体、cairo PDF + 600/300 dpi LZW TIFF
#
# 另将 cor_label_xy 参数真正用起来（原版计算了 GAM 文本却从未 annotate），
# 并把 R 与 p 的注释按规则格式化。
###############################################################################

setwd(RAWDAT)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(cowplot)
library(mgcv)
library(grid)
library(ragg)

# ---- 输出位置：改这一行即可 ----
OUT_DIR  <- file.path(RAWDAT, "前沿分析_GAM")
out_path <- function(f) {
  if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)
  file.path(OUT_DIR, f)
}

# ============================================================================
# 全文格式约定
# ============================================================================
DEC        <- "\u00B7"   # 中点
EN_DASH    <- "\u2013"   # –
MINUS      <- "\u2212"   # −
BASE_SIZE  <- 8
FONT       <- "Helvetica"        # 与 Figure S18/S19 一致；Lancet 接受 Helvetica/Arial
LABEL_SIZE <- 8 / .pt            # 国名标签 7 pt，等于下限
FIG_W_MM   <- 180                # 双栏硬上限
FIG_H_MM   <- 220                # A、B 上下排列；侧排见 PANEL_NCOL
PANEL_NCOL <- 1                  # 1 = A 上 B 下（推荐）；2 = 左右并排
COR_METHOD <- "pearson"          # 图内 R 的算法，按需改 "spearman"
ANN_POS    <- c(0.015, 0.03)     # 注释位置，面板相对坐标 (x, y)，左下角为 (0,0)

mm2in <- function(x) x / 25.4

# 小数点 -> 中点；ASCII 连字符 -> 真减号
dec <- function(x) {
  x <- gsub("\\.", DEC, as.character(x))
  gsub("(?<![0-9A-Za-z])-(?=[0-9])", MINUS, x, perl = TRUE)
}

# 全文 p 规则（本次的核心修改）
fmt_p <- function(p) {
  vapply(p, function(pi) {
    if (is.na(pi))   return(NA_character_)
    if (pi < 0.0001) return(paste0("<0", DEC, "0001"))
    if (pi < 0.01)   return(dec(sprintf("%.4f", pi)))
    dec(sprintf("%.2f", pi))
  }, character(1))
}

# "p<0·0001" 不带等号，其余带等号
p_inline <- function(prefix, p) {
  paste0(prefix, if (!is.na(p) && p < 0.0001) "" else "=", fmt_p(p))
}

# 坐标轴数字：五位及以上加细空格千分位，四位以内不加（与正文 8012 一致）
fmt_num <- function(v) {
  out <- ifelse(abs(v) >= 10000,
                formatC(v, format = "d", big.mark = "\u2009"),
                formatC(v, format = "fg", digits = 7))
  dec(trimws(out))
}

# 兜底：把成图 gtable 内所有 text grob 的小数号与负号统一
dotify <- function(g) {
  if (inherits(g, "gTree") || inherits(g, "gtable")) {
    g$children <- lapply(g$children, dotify)
  }
  if (inherits(g, "text") && !is.null(g$label)) g$label <- dec(g$label)
  g
}

theme_frontier <- function() {
  theme_bw(base_size = BASE_SIZE, base_family = FONT) +
    theme(
      plot.title        = element_text(size = BASE_SIZE, face = "plain",
                                       hjust = 0, margin = margin(b = 3)),
      axis.title        = element_text(size = BASE_SIZE),
      axis.text         = element_text(size = BASE_SIZE, colour = "black"),
      panel.border      = element_rect(linewidth = 0.3, colour = "black"),
      panel.grid.major  = element_line(linewidth = 0.2, colour = "grey92"),
      panel.grid.minor  = element_blank(),
      axis.ticks        = element_line(linewidth = 0.3, colour = "black"),
      legend.title      = element_text(size = BASE_SIZE),
      legend.text       = element_text(size = BASE_SIZE - 1),
      legend.key        = element_blank(),
      legend.key.height = unit(3.2, "mm"),
      legend.key.width  = unit(5, "mm"),
      legend.position   = "bottom",
      legend.box        = "horizontal",
      legend.spacing.x   = unit(24, "mm"),
      legend.box.margin  = margin(1, 0, 1, 0, "mm"),
      legend.margin     = margin(0, 0, 0, 0),
      legend.box.spacing = unit(2, "mm"),
      plot.margin       = margin(2, 3, 2, 2, "mm")
    )
}

save_figure <- function(gg, stem, w_mm = FIG_W_MM, h_mm = FIG_H_MM) {
  draw <- function() {
    grid.newpage()
    grid.draw(dotify(ggplotGrob(gg)))
  }
  # 原先此处另出一份 cairo_pdf 矢量母版, 已移除: cairo 把每个字形的步进四舍五入
  # 到整数 point (0.353 mm), 再用 TJ 微调补偿, 同一单词内字距不匀 (最坏 1 pt)。
  # 投稿件本就是下面两份 TIF, ragg 在目标分辨率排字, 600 dpi 误差 <= 0.042 mm。
  agg_tiff(out_path(paste0(stem, ".tif")), width = w_mm, height = h_mm,
           units = "mm", res = 600, compression = "lzw", background = "white")
  draw(); dev.off()
  agg_tiff(out_path(paste0(stem, "_300dpi.tif")), width = w_mm, height = h_mm,
           units = "mm", res = 300, compression = "lzw", background = "white")
  draw(); dev.off()
  invisible(NULL)
}

SCO <- c("Belarus", "China", "India", "Iran (Islamic Republic of)",
         "Kazakhstan", "Kyrgyzstan", "Pakistan", "Russian Federation",
         "Tajikistan", "Uzbekistan")

# 简化SCO国家名用于标签
sco_short <- c("Belarus" = "BLR", "China" = "CHN", "India" = "IND",
               "Iran (Islamic Republic of)" = "IRN",
               "Kazakhstan" = "KAZ", "Kyrgyzstan" = "KGZ",
               "Pakistan" = "PAK", "Russian Federation" = "RUS",
               "Tajikistan" = "TJK", "Uzbekistan" = "UZB")

# ============================================================================
# 通用分析与绘图函数
# ============================================================================

###GAM建模
# ============================================================================
# 辅助函数: 拟合GAM并提取结果   —— 建模部分未改动
# ============================================================================
fit_gam_and_summarize <- function(data) {
  # data: 需包含 SDI 和 ASDR 两列, 且已去除NA
  data <- data[!is.na(data$SDI) & !is.na(data$ASDR), ]
  
  # 1. 拟合GAM(带默认thin-plate样条)
  gam_fit <- gam(ASDR ~ s(SDI), data = data, method = "REML")
  
  # 2. 拟合线性对照模型(用于似然比检验)
  lm_fit <- gam(ASDR ~ SDI, data = data, method = "REML")
  
  # 3. 提取GAM summary
  gam_sum <- summary(gam_fit)
  edf <- gam_sum$s.table[1, "edf"]         # 有效自由度
  F_stat <- gam_sum$s.table[1, "F"]        # F统计量
  p_smooth <- gam_sum$s.table[1, "p-value"] # 平滑项显著性P值
  r2_adj <- gam_sum$r.sq                   # 调整R²
  dev_exp <- gam_sum$dev.expl              # 偏差解释比例
  
  # 4. 非线性检验(线性 vs GAM的似然比检验)
  gam_ml <- gam(ASDR ~ s(SDI), data = data, method = "ML")
  lm_ml  <- gam(ASDR ~ SDI,    data = data, method = "ML")
  lrt <- anova(lm_ml, gam_ml, test = "Chisq")
  p_nonlinear <- lrt$`Pr(>Chi)`[2]
  
  list(
    gam_fit = gam_fit,
    edf = edf,
    F_stat = F_stat,
    p_smooth = p_smooth,
    r2_adj = r2_adj,
    dev_exp = dev_exp,
    p_nonlinear = p_nonlinear,
    n = nrow(data)
  )
}

# ============================================================================
# 辅助函数: 格式化GAM结果为图上注释文本
#   原 fmt_p 在 p<0.001 处封顶，是本次要修的地方；现调用全局 fmt_p
# ============================================================================
format_gam_text <- function(gam_res, compact = TRUE) {
  if (compact) {
    # 紧凑版: 适合图内annotate使用
    sprintf("GAM: edf=%s, %s\nAdj. R\u00b2=%s",
            dec(sprintf("%.2f", gam_res$edf)),
            p_inline("p", gam_res$p_smooth),
            dec(sprintf("%.3f", gam_res$r2_adj)))
  } else {
    # 完整版: 适合论文表格/结果文字
    sprintf(paste0("GAM smooth: edf=%s, F=%s, %s | Adj. R\u00b2=%s | ",
                   "Deviance explained=%s%% | Nonlinearity test %s | n=%d"),
            dec(sprintf("%.2f", gam_res$edf)),
            dec(sprintf("%.2f", gam_res$F_stat)),
            p_inline("p", gam_res$p_smooth),
            dec(sprintf("%.3f", gam_res$r2_adj)),
            dec(sprintf("%.1f", gam_res$dev_exp * 100)),
            p_inline("p", gam_res$p_nonlinear),
            gam_res$n)
  }
}

# 相关系数注释：R 两位小数（与正文 r=0·68、R=−0·13 一致）
cor_annotation <- function(d) {
  d <- d[!is.na(d$SDI) & !is.na(d$ASDR), ]
  ct <- suppressWarnings(cor.test(d$SDI, d$ASDR, method = COR_METHOD))
  list(R = unname(ct$estimate), p = ct$p.value,
       label = paste0("R=", dec(sprintf("%.2f", unname(ct$estimate))), "\n",
                      p_inline("p", ct$p.value)))
}
# 两行注释：R（线性成分）+ 非线性 p（平滑项）
ann_label <- function(cor_res, gam_res) {
  paste0("R=", dec(sprintf("%.2f", cor_res$R)), "\n",
         p_inline("Nonlinearity p", gam_res$p_nonlinear))
}

###分析绘图函数
plot_frontier_analysis <- function(boostrap_DEA, frontier_SDI,
                                   cause_short, cor_label_xy) {
  # cor_label_xy: list(x=0.1, y_left=5000, y_right=7500) 注释位置（y 为 ASDR 值）
  
  # --- 1. 计算效率差(相对) ---
  boostrap_DEA <- boostrap_DEA %>%
    mutate(eff_diff = ASDR - frontier,
           eff_ratio = (ASDR - frontier) / frontier)
  
  dat_2001 <- boostrap_DEA %>% filter(year == 2001) %>%
    select(location, eff_ratio_2001 = eff_ratio)
  dat_2023 <- boostrap_DEA %>% filter(year == 2023)
  
  dat_2023 <- dat_2023 %>%
    left_join(dat_2001, by = "location") %>%
    mutate(trend = ifelse(eff_ratio < eff_ratio_2001,
                          "Closer to frontier", "Farther from frontier"))
  
  # --- 2. SCO 2023数据+排序 ---
  SCO_2023 <- dat_2023 %>% filter(location %in% SCO) %>%
    mutate(location_label = sco_short[location]) %>%
    arrange(eff_ratio) %>%
    mutate(emph = ifelse(row_number() <= 2 |
                           row_number() > n() - 2, "bold", "plain"))
  
  best   <- SCO_2023[1:2, ]
  worst  <- SCO_2023[(nrow(SCO_2023)-1):nrow(SCO_2023), ]
  middle <- SCO_2023[3:(nrow(SCO_2023)-2), ]
  
  # --- 3. GAM 与相关系数 ---
  gam_all  <- fit_gam_and_summarize(frontier_SDI)
  gam_2023 <- fit_gam_and_summarize(frontier_SDI[frontier_SDI$year == 2023, ])
  
  cor_all  <- cor_annotation(frontier_SDI)
  cor_2023 <- cor_annotation(frontier_SDI[frontier_SDI$year == 2023, ])
  
  # 控制台打印完整版(用于 Table S8 与正文引用)
  cat(sprintf("\n[%s - All years] %s\n", cause_short,
              format_gam_text(gam_all, compact = FALSE)))
  cat(sprintf("[%s - All years] R=%.4f, p=%.3e\n",
              cause_short, cor_all$R, cor_all$p))
  cat(sprintf("[%s - 2023]      %s\n", cause_short,
              format_gam_text(gam_2023, compact = FALSE)))
  cat(sprintf("[%s - 2023]      R=%.4f, p=%.3e\n",
              cause_short, cor_2023$R, cor_2023$p))
  
  y_lab <- paste0("ASDR of ", cause_short, " per 100 000")
  # 面板顶端（y 轴反转，ASDR 最小值即最上方），用于放 R/p 注释
  y_top_A <- min(c(boostrap_DEA$ASDR, boostrap_DEA$frontier), na.rm = TRUE)
  y_top_B <- min(c(dat_2023$ASDR,     boostrap_DEA$frontier), na.rm = TRUE)
  
  # ==========================================================================
  # 左图(A): SCO轨迹图 (2001–2023)
  # ==========================================================================
  SCO_traj <- boostrap_DEA %>%
    filter(location %in% SCO) %>%
    mutate(location_label = sco_short[location]) %>%
    arrange(location, year)
  
  SCO_start <- SCO_traj %>% filter(year == 2001) %>%
    left_join(SCO_2023[, c("location", "emph")], by = "location")
  SCO_end   <- SCO_traj %>% filter(year == 2023)
  
  p_left <- ggplot() +
    geom_point(data = boostrap_DEA %>% filter(!(location %in% SCO)),
               aes(SDI, ASDR),
               color = "grey85", size = 0.3, alpha = 0.5) +
    stat_smooth(data = boostrap_DEA,
                aes(SDI, frontier, linetype = "Frontier"),
                color = "black", method = "loess", span = 0.3,
                se = FALSE, fullrange = TRUE, linewidth = 0.5) +
    stat_smooth(data = boostrap_DEA,
                aes(SDI, ASDR, linetype = "Regression"),
                color = "black", method = "gam", formula = y ~ s(x),
                se = FALSE, fullrange = TRUE, linewidth = 0.5) +
    geom_path(data = SCO_traj,
              aes(SDI, ASDR, group = location, color = year),
              linewidth = 0.45, alpha = 0.85) +
    geom_point(data = SCO_start, aes(SDI, ASDR),
               shape = 21, fill = "white", color = "#1A5B1A",
               size = 1.3, stroke = 0.5) +
    geom_point(data = SCO_end, aes(SDI, ASDR),
               shape = 16, color = "#D2691E", size = 1.5) +
    geom_label_repel(data = SCO_start,
                     aes(SDI, ASDR, label = location_label, fontface = emph),
                     colour = "black",
                     size = LABEL_SIZE, family = FONT,
                     fill = scales::alpha("white", 0.85),
                     label.size = 0,
                     label.padding = unit(0.5, "mm"),
                     label.r = unit(0.3, "mm"),
                     max.overlaps = Inf, min.segment.length = 0,
                     box.padding = 0.6, point.padding = 0.25,
                     force = 6, force_pull = 0.4,
                     max.iter = 20000, max.time = 2, seed = 42,
                     segment.colour = "grey35", segment.size = 0.25) +
    scale_discrete_identity(aesthetics = "fontface") +
    annotation_custom(grid::textGrob(
      ann_label(cor_all, gam_all),
      x = unit(ANN_POS[1], "npc"), y = unit(ANN_POS[2], "npc"),
      hjust = 0, vjust = 0,
      gp = grid::gpar(fontsize = BASE_SIZE, fontfamily = FONT,
                      col = "grey25", lineheight = 1.05))) +
    scale_color_gradient(name = "Year",
                         low = "#1A5B1A", high = "#D2691E",
                         breaks = c(2001, 2010, 2020),
                         labels = function(v) dec(format(v)),
                         guide = guide_colourbar(title.position = "top",
                                                 barheight = unit(2.5, "mm"),
                                                 barwidth = unit(22, "mm"),
                                                 order = 1)) +
    scale_x_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1),
                       labels = function(v) dec(sprintf("%.1f", v))) +
    scale_y_reverse(labels = fmt_num) +
    scale_linetype_manual(name = "Line type",
                          values = c("Frontier" = "solid",
                                     "Regression" = "dashed"),
                          guide = guide_legend(title.position = "top",
                                               order = 2)) +
    labs(x = "SDI", y = y_lab,
         title = paste0("SCO trajectories, 2001", EN_DASH, "2023")) +
    theme_frontier()
  
  # ==========================================================================
  # 右图(B): 2023截面图
  # ==========================================================================
  trend_colors <- c("Closer to frontier" = "#2C5985",
                    "Farther from frontier" = "#A63603")
  
  p_right <- ggplot() +
    geom_point(data = dat_2023,
               aes(SDI, ASDR, color = trend),
               size = 1.1, alpha = 0.75) +
    stat_smooth(data = boostrap_DEA,
                aes(SDI, frontier, linetype = "Frontier"),
                color = "black", method = "loess", span = 0.3,
                se = FALSE, fullrange = TRUE, linewidth = 0.5) +
    stat_smooth(data = dat_2023,
                aes(SDI, ASDR, linetype = "Regression (2023)"),
                color = "black", method = "gam", formula = y ~ s(x),
                se = FALSE, fullrange = TRUE, linewidth = 0.5) +
    # 所有成员国统一黑圈；最接近/最远离前沿的两国用更粗的圈强调
    geom_point(data = SCO_2023, aes(SDI, ASDR),
               shape = 21, fill = NA, colour = "black",
               size = 2.0, stroke = 0.4) +
    geom_point(data = rbind(best, worst), aes(SDI, ASDR),
               shape = 21, fill = NA, colour = "black",
               size = 2.8, stroke = 0.9) +
    # 十个标签一次性 repel，彼此可见；白色描边保证压在散点上仍可读
    geom_label_repel(data = SCO_2023,
                     aes(SDI, ASDR, label = location_label, fontface = emph),
                     colour = "black",
                     size = LABEL_SIZE, family = FONT,
                     fill = scales::alpha("white", 0.85),
                     label.size = 0,                      # 无边框
                     label.padding = unit(0.5, "mm"),
                     label.r = unit(0.3, "mm"),
                     max.overlaps = Inf, min.segment.length = 0,
                     box.padding = 0.6, point.padding = 0.25,
                     force = 8, force_pull = 0.4,
                     max.iter = 20000, max.time = 2, seed = 42,
                     segment.colour = "grey35", segment.size = 0.25) +
    scale_discrete_identity(aesthetics = "fontface") +
    annotation_custom(grid::textGrob(
      ann_label(cor_2023, gam_2023),
      x = unit(ANN_POS[1], "npc"), y = unit(ANN_POS[2], "npc"),
      hjust = 0, vjust = 0,
      gp = grid::gpar(fontsize = BASE_SIZE, fontfamily = FONT,
                      col = "grey25", lineheight = 1.05))) +
    scale_color_manual(name = paste0("Trend (2001", EN_DASH, "2023)"),
                       values = trend_colors,
                       guide = guide_legend(title.position = "top", order = 1)) +
    scale_linetype_manual(name = "Line type",
                          values = c("Frontier" = "solid",
                                     "Regression (2023)" = "dashed"),
                          guide = guide_legend(title.position = "top", order = 2)) +
    scale_x_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1),
                       labels = function(v) dec(sprintf("%.1f", v))) +
    scale_y_reverse(labels = fmt_num) +
    labs(x = "SDI", y = y_lab, title = "2023 cross-section") +
    theme_frontier()
  
  # ==========================================================================
  # 拼图：面板字母由 plot_grid 统一加，标题里不再写 "A." / "B."
  # ==========================================================================
  combined <- plot_grid(p_left, p_right, ncol = PANEL_NCOL,
                        align = if (PANEL_NCOL == 1) "v" else "h",
                        axis = if (PANEL_NCOL == 1) "lr" else "tb",
                        labels = c("A", "B"),
                        label_size = 11, label_fontface = "bold",
                        label_fontfamily = FONT,
                        hjust = 0, vjust = 1)
  
  list(plot = combined, left = p_left, right = p_right,
       SCO_2023 = SCO_2023, best = best, worst = worst,
       gam_all = gam_all, gam_2023 = gam_2023,
       cor_all = cor_all, cor_2023 = cor_2023)
}

# 汇总一行 Table S8 记录
s8_row <- function(res, cause_label) {
  data.frame(
    cause              = cause_label,
    edf_all            = res$gam_all$edf,
    dev_expl_all       = res$gam_all$dev_exp * 100,
    R_all              = res$cor_all$R,
    p_R_all            = res$cor_all$p,
    p_R_all_d          = fmt_p(res$cor_all$p),
    p_nonlinear_all    = res$gam_all$p_nonlinear,
    p_nonlinear_all_d  = fmt_p(res$gam_all$p_nonlinear),
    edf_2023           = res$gam_2023$edf,
    dev_expl_2023      = res$gam_2023$dev_exp * 100,
    p_nonlinear_2023   = res$gam_2023$p_nonlinear,
    p_nonlinear_2023_d = fmt_p(res$gam_2023$p_nonlinear),
    R_2023             = res$cor_2023$R,
    p_R_2023           = res$cor_2023$p,
    p_R_2023_d         = fmt_p(res$cor_2023$p),
    stringsAsFactors   = FALSE
  )
}

# ============================================================================
# 运行: Stroke (主图)
# ============================================================================
ASDR <- read.csv('0123-allASDR.csv')
SDI_data <- read.csv('SDI_2023.csv') %>%
  rename(location = location_name, year = year_id, SDI = mean_value) %>%
  select(location, year, SDI)

# Stroke
ASDR_str <- ASDR %>% filter(cause_name == 'Stroke') %>%
  rename(location = location_name, ASDR = val) %>%
  select(location, year, ASDR)
frontier_SDI_str <- left_join(ASDR_str, SDI_data, by = c('location', 'year'))
boostrap_DEA_str <- read.csv('boostrap_DEA_str.csv')

res_str <- plot_frontier_analysis(
  boostrap_DEA_str, frontier_SDI_str,
  cause_short = "stroke",
  cor_label_xy = list(x = 0.01, y_left = 500, y_right = 400)
)
save_figure(res_str$plot, "FigureS7")

cat("Stroke 2023排序(从最接近到最远离前沿):\n")
print(res_str$SCO_2023[, c("location_label", "eff_ratio", "trend")])

s8 <- s8_row(res_str, "Total stroke")

###############################################################################
# Figure S8 拼图块 —— 三个亚型合并为一张，A=IS、B=ICH、C=SAH 左右排布
#
# 替换范围：原脚本 "运行: 3个亚型 (补充图)" 整段（subtypes 定义 + for 循环）
#           及其后的 write.csv(s8, ...)
# 建模与 plot_frontier_analysis() 本体不动。
#
# 版面：横向 250 × 170 mm（附录横排页），3 列 × 2 行
#       列 = 亚型（A/B/C），行 = 轨迹图 / 2023 截面
#       图例每行共享一套，不再重复六次
###############################################################################

# ---- 亚型定义：cause_short 改用英式拼写，与正文/图注一致 --------------------
subtypes <- list(
  list(name = "Ischemic stroke",            short = "IS",
       label = "ischaemic stroke",          header = "Ischaemic stroke",
       csv = "boostrap_DEA_isc.csv",
       cor_xy = list(x = 0.01, y_left = 0, y_right = 0)),
  list(name = "Intracerebral hemorrhage",   short = "ICH",
       label = "intracerebral haemorrhage", header = "Intracerebral haemorrhage",
       csv = "boostrap_DEA_ich.csv",
       cor_xy = list(x = 0.01, y_left = 0, y_right = 0)),
  list(name = "Subarachnoid hemorrhage",    short = "SAH",
       label = "subarachnoid haemorrhage",  header = "Subarachnoid haemorrhage",
       csv = "boostrap_DEA_sah.csv",
       cor_xy = list(x = 0.01, y_left = 0, y_right = 0))
)

# 拼图后单格只有约 78 mm 宽，国名标签回到 7 pt（下限）
LABEL_SIZE <- 7 / .pt

res_list <- list()

for (st in subtypes) {
  ASDR_sub <- ASDR %>% filter(cause_name == st$name) %>%
    rename(location = location_name, ASDR = val) %>%
    select(location, year, ASDR)
  frontier_SDI_sub <- left_join(ASDR_sub, SDI_data, by = c('location', 'year'))
  boostrap_DEA_sub <- read.csv(st$csv)
  
  res <- plot_frontier_analysis(
    boostrap_DEA_sub, frontier_SDI_sub,
    cause_short  = st$label,          # 英式拼写
    cor_label_xy = st$cor_xy
  )
  
  cat("\n", st$name, "2023排序:\n")
  print(res$SCO_2023[, c("location_label", "eff_ratio", "trend")])
  
  res_list[[st$short]] <- res
  s8 <- rbind(s8, s8_row(res, st$name))
}

write.csv(s8, out_path("TableS8_GAM_frontier.csv"), row.names = FALSE)
cat("\n===== Table S8 =====\n"); print(s8)


# ===========================================================================
# 拼接
# ===========================================================================

# 单格：去掉图例、标题、x 轴标题；最左列保留 y 轴标题并改为通用措辞
strip_panel <- function(p, show_y = FALSE, show_x = FALSE) {
  p <- p + theme(
    legend.position = "none",
    plot.title      = element_blank(),
    plot.margin     = margin(1.5, 2, 1.5, 1.5, "mm")
  )
  p <- p + labs(y = if (show_y) "ASDR per 100 000" else NULL,
                x = if (show_x) "SDI" else NULL)
  if (!show_y) p <- p + theme(axis.title.y = element_blank())
  if (!show_x) p <- p + theme(axis.title.x = element_blank())
  p
}

# 提取底部图例（兼容 ggplot2 3.5 前后的 cowplot 接口）
get_leg <- function(p) {
  g <- tryCatch(cowplot::get_plot_component(p, "guide-box-bottom"),
                error = function(e) cowplot::get_legend(p))
  if (inherits(g, "list")) g <- g[[1]]
  g
}

# 列表头：面板字母 + 亚型名
col_header <- function(letter, txt) {
  ggdraw() +
    draw_label(letter, x = 0.015, hjust = 0, fontface = "bold",
               fontfamily = FONT, size = 11) +
    draw_label(txt, x = 0.5, hjust = 0.5,
               fontfamily = FONT, size = BASE_SIZE)
}

# 行标签（竖排）
row_label <- function(txt) {
  ggdraw() + draw_label(txt, angle = 90, fontfamily = FONT, size = BASE_SIZE)
}

ord <- c("IS", "ICH", "SAH")

row1_panels <- lapply(seq_along(ord), function(i)
  strip_panel(res_list[[ord[i]]]$left,  show_y = (i == 1), show_x = FALSE))
row2_panels <- lapply(seq_along(ord), function(i)
  strip_panel(res_list[[ord[i]]]$right, show_y = (i == 1), show_x = TRUE))

leg_row1 <- get_leg(res_list[["IS"]]$left)    # Year + Line type
leg_row2 <- get_leg(res_list[["IS"]]$right)   # Trend + Line type

# 最左列因带 y 轴标题略宽
W <- c(1.10, 1, 1)

hdr <- plot_grid(plotlist = lapply(seq_along(ord), function(i)
  col_header(LETTERS[i], subtypes[[i]]$header)),
  nrow = 1, rel_widths = W)

row1 <- plot_grid(plotlist = row1_panels, nrow = 1,
                  align = "hv", axis = "tblr", rel_widths = W)
row2 <- plot_grid(plotlist = row2_panels, nrow = 1,
                  align = "hv", axis = "tblr", rel_widths = W)

# 左侧竖排行标签 + 面板
band1 <- plot_grid(row_label(paste0("SCO trajectories, 2001", EN_DASH, "2023")),
                   row1, nrow = 1, rel_widths = c(0.035, 1))
band2 <- plot_grid(row_label("2023 cross-section"),
                   row2, nrow = 1, rel_widths = c(0.035, 1))
hdr2  <- plot_grid(NULL, hdr, nrow = 1, rel_widths = c(0.035, 1))

figS8 <- plot_grid(
  hdr2, band1, leg_row1, band2, leg_row2,
  ncol = 1,
  rel_heights = c(0.06, 1, 0.20, 1, 0.16)
)

print(figS8)
save_figure(figS8, "FigureS8",
            w_mm = 250, h_mm = 170)      # 附录横排页

cat("\n========== 全部完成 ==========\n")