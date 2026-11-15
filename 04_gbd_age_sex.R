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
# GBD SCO 10国 中风 年龄-性别负担热图
#   Figure S5  总卒中     : 3 行(Incidence / Mortality / DALYs) × 2 列(Male / Female)
#   Figure S6  三个亚型   : 3 行(measure) × 6 列(subtype × sex，嵌套分面)
#
# 本版相对旧版的改动：
#   [幅面] 旧版 14 × 9 in = 355 × 229 mm，远超 180 mm；期刊等比缩到 180 mm 会把
#          8 pt 压到 4 pt。现改为按毫米出图，宽度 180 mm 硬上限。
#   [字号] base_size 8，图内最小 7 pt，缩放后仍达标。
#   [字体] 改 cairo_pdf 嵌入字体。旧版 pdf() 用 8 位编码，已把 30 个年龄组
#          连字符静默换成减号 U+2212（"20−24"）；现统一为 en dash "20–24"。
#   [数字] 图例 "1,000" 改 "1000"，10 000 用空格分位，与全文一致。
#   [配色] 顺序色标锚定到全文主色 #A63603。
#   [版面] 弃用 cowplot 手工拼装，改单个 ggplot + 分面：面板自动对齐，
#          行/列标题由 strip 承担，也绕开了 ggplot2 3.5+ 的
#          cowplot::get_legend 返回空图例的问题。
#          S6 的两层列表头用 ggh4x::facet_nested，亚型名跨两列合并居中；
#          原先用 facet_grid(measure ~ subtype + sex) 会逐列重复亚型名并
#          因列宽不足而截断（"acerebral haemorrh"）。
#          S6 每列仅约 22 mm，年龄轴改在色块分割线上标整岁年龄
#          （20 / 45 / 70 / 95），间距完全均匀、左右对称，标签可水平排版。
#   [删除] 旧版第 4 节与第 5 节的折线图（SuppFig_age_sex_facet_line /
#          SuppFig_*_age_sex_line）没有任何图注引用，且同样不合规，已移除。
#
# 注意: 旧版把总卒中热图存成 "Figure2_age_sex_heatmap"，但正文 Figure 2 是
#       代谢归因热图。此图是 Figure S5，文件名已更正。
###############################################################################

library(tidyverse)
library(scales)
library(grid)

# S6 的两层列表头需要跨列合并，ggplot2 原生 facet_grid 做不到（会把亚型名
# 逐列重复并截断），改用 ggh4x::facet_nested
if (!requireNamespace("ggh4x", quietly = TRUE)) install.packages("ggh4x")
library(ggh4x)

# ============================================================================
# 全文格式约定（同交接文档第 2 节）
# ============================================================================
BASE_SIZE <- 8
FONT      <- "Helvetica"
EN        <- "\u2013"                 # en dash
PAL       <- c("#FFF7BC", "#FEE391", "#FEC44F", "#FE9929",
               "#EC7014", "#CC4C02", "#A63603")   # 暖色顺序标，锚定 #A63603

en <- function(x) gsub("-", EN, x, fixed = TRUE)

# 出图: 直接以目标分辨率栅格化 (ragg), 不经 cairo。
# 说明: cairo_pdf 会把每个字形的步进四舍五入到整数 point (0.353 mm), 再用
#       TJ 微调补偿, 导致同一单词内字距不匀 (最坏 1 pt)。ragg 在目标分辨率的
#       像素网格上排字, 600 dpi 下取整误差 <= 0.042 mm, 肉眼不可见。
#       file 传入的扩展名会被忽略, 实际输出 <stem>.tif (600 dpi) 与
#       <stem>_300dpi.tif (300 dpi)。
if (!requireNamespace("ragg", quietly = TRUE)) install.packages("ragg")

save_panel <- function(p, file, w_mm, h_mm) {
  stem <- sub("\\.(pdf|tif|tiff|png)$", "", file)
  for (r in c(600, 300)) {
    f <- if (r == 600) paste0(stem, ".tif") else paste0(stem, "_300dpi.tif")
    ragg::agg_tiff(f, width = w_mm, height = h_mm, units = "mm",
                   res = r, compression = "lzw", background = "white",
                   scaling = 1)
    print(p); dev.off()
    cat(sprintf("  saved: %s (%.1f x %.1f mm @ %d dpi)\n", f, w_mm, h_mm, r))
  }
}

# ============================================================================
# 0. 数据加载 [逻辑未改]
# ============================================================================
setwd(RAWDAT)
agesex        <- read.csv('age&sex.csv', header = TRUE)
agesex_global <- read.csv('age&sex_global.csv', header = TRUE)
agesex        <- rbind(agesex, agesex_global)

age_levels_raw <- c("<20", "20-24", "25-29", "30-34", "35-39", "40-44",
                    "45-49", "50-54", "55-59", "60-64", "65-69", "70-74",
                    "75-79", "80-84", "85-89", "90-94", "95+")
age_levels <- en(age_levels_raw)

measure_map <- c("Incidence"                             = "Incidence",
                 "Deaths"                                = "Mortality",
                 "DALYs (Disability-Adjusted Life Years)" = "DALYs")

# GBD 用美式拼写，展示时改用与正文一致的英式拼写
subtype_map <- c("Ischemic stroke"          = "Ischaemic stroke",
                 "Intracerebral hemorrhage" = "Intracerebral haemorrhage",
                 "Subarachnoid hemorrhage"  = "Subarachnoid haemorrhage")

agesex <- agesex %>%
  mutate(
    location_name = case_when(
      location_name == "Russian Federation"        ~ "Russia",
      location_name == "Iran (Islamic Republic of)" ~ "Iran",
      TRUE ~ location_name
    ),
    age_label = factor(en(gsub(" years", "", age_name)), levels = age_levels),
    measure   = factor(unname(measure_map[measure_name]),
                       levels = c("Incidence", "Mortality", "DALYs")),
    sex       = factor(sex_name, levels = c("Male", "Female")),
    val_log   = log10(val + 1)                 # +1 避免 log(0)，与旧版一致
  )

# 热图 y 轴自下而上绘制，反转字母序使 A 在顶部，Global 置顶
country_alpha   <- sort(unique(agesex$location_name[agesex$location_name != "Global"]))
country_levels  <- c(rev(country_alpha), "Global")
agesex$location_name <- factor(agesex$location_name, levels = country_levels)

# Global 行加粗；向量化 face 在 ggplot2 3.5+ 会有一条 "not officially supported"
# 的提示，可忽略，输出正确
country_face <- ifelse(country_levels == "Global", "bold", "plain")

# ============================================================================
# 1. 共用主题与色标
# ============================================================================
theme_heat <- function(base = BASE_SIZE) {
  theme_minimal(base_size = base, base_family = FONT) +
    theme(
      text              = element_text(size = base, family = FONT),
      axis.title.x      = element_text(size = base, margin = margin(t = 1, unit = "mm")),
      axis.title.y      = element_blank(),
      axis.text.y       = element_text(size = base - 1, colour = "black",
                                       face = country_face),
      axis.ticks        = element_blank(),
      panel.grid        = element_blank(),
      panel.spacing.x   = unit(1, "mm"),
      panel.spacing.y   = unit(1.5, "mm"),
      strip.text        = element_text(size = base, face = "bold"),
      strip.placement   = "outside",
      legend.position   = "bottom",
      legend.title      = element_text(size = base - 1, vjust = 1),
      legend.text       = element_text(size = base - 1),
      legend.key.height = unit(2.5, "mm"),
      legend.key.width  = unit(14, "mm"),
      legend.margin     = margin(t = 1, unit = "mm"),
      plot.margin       = margin(2, 2, 2, 2, "mm")
    )
}

fill_scale <- scale_fill_gradientn(
  colours = PAL,
  name    = "Rate per 100 000 (log10 scale)",
  breaks  = log10(c(10, 100, 1000, 10000) + 1),
  labels  = c("10", "100", "1000", "10 000"),
  guide   = guide_colourbar(title.position = "top", title.hjust = 0.5)
)

# 年龄组是 5 岁一档，相邻色块的分割线正好落在整岁边界上：<20 与 20–24 之间
# 是 20 岁，35–39 与 40–44 之间是 40 岁，依此类推。x_mode = "boundary" 时把
# x 轴改成连续坐标，刻度画在分割线上、标注该处的年龄，间距因此可以做到完全均匀。
# 边界年龄 a 对应的位置 = (a - 20)/5 + 1.5
age_pos <- function(a) (a - 20) / 5 + 1.5

# x_breaks: 离散模式下要标出的年龄组；NULL 表示全标
# x_ages  : 边界模式下要标出的整岁年龄
build_heat <- function(dat, nested, x_angle = 45,
                       x_mode = c("discrete", "boundary"),
                       x_breaks = NULL, x_ages = NULL) {
  x_mode <- match.arg(x_mode)
  
  if (x_mode == "boundary") {
    dat$age_idx <- as.integer(dat$age_label)
    p <- ggplot(dat, aes(x = age_idx, y = location_name, fill = val_log)) +
      geom_tile(width = 1, colour = "white", linewidth = 0.15) +
      scale_x_continuous(breaks = age_pos(x_ages),
                         labels = as.character(x_ages),
                         limits = c(0.5, length(age_levels) + 0.5),
                         expand = expansion(0, 0)) +
      labs(x = "Age, years")
    x_theme <- theme(
      axis.text.x         = element_text(size = BASE_SIZE - 1, colour = "black"),
      axis.ticks.x        = element_line(colour = "grey40", linewidth = 0.3),
      axis.ticks.length.x = unit(1, "mm"))
  } else {
    brk <- if (is.null(x_breaks)) age_levels else x_breaks
    stopifnot(all(brk %in% age_levels))
    p <- ggplot(dat, aes(x = age_label, y = location_name, fill = val_log)) +
      geom_tile(colour = "white", linewidth = 0.15) +
      scale_x_discrete(breaks = brk) +
      labs(x = "Age group, years")
    x_theme <- theme(axis.text.x = if (x_angle == 90) {
      element_text(size = BASE_SIZE - 1, colour = "black",
                   angle = 90, hjust = 1, vjust = 0.5)
    } else {
      element_text(size = BASE_SIZE - 1, colour = "black",
                   angle = 45, hjust = 1)
    })
  }
  
  # theme_heat 里 axis.ticks = element_blank()，模式专属主题必须加在它之后
  p <- p + fill_scale +
    scale_y_discrete(expand = expansion(0, 0)) +
    theme_heat() + x_theme
  
  if (nested) {
    # facet_nested 会把相邻且相同的外层 strip 合并成一条跨列居中的标题，
    # 亚型名因此只出现一次、且有两列的宽度可用（不再被截断）
    p + ggh4x::facet_nested(
      measure ~ subtype + sex, switch = "y",
      nest_line = element_line(colour = "grey70", linewidth = 0.3),
      strip     = ggh4x::strip_nested(bleed = FALSE))
  } else {
    p + facet_grid(measure ~ sex, switch = "y")
  }
}

# ============================================================================
# 2. Figure S5: 总卒中
# ============================================================================
setwd(file.path(RAWDAT, "总体、年龄、性别负担/年龄、性别"))

cat("\n== Figure S5: total stroke ==\n")
dat_s5 <- agesex %>% filter(cause_name == "Stroke", !is.na(measure))

save_panel(build_heat(dat_s5, nested = FALSE, x_angle = 45,
                      x_mode = "discrete", x_breaks = NULL),
           "FigureS5", w_mm = 180, h_mm = 135)   # 单面板, 直接就是投稿版

# ============================================================================
# 3. Figure S6: 三个亚型 (3 measure × 6 列 = 3 subtype × 2 sex)
#    嵌套 strip 由 facet_grid(measure ~ subtype + sex) 自动生成:
#    上层跨列写亚型名，下层写性别
# ============================================================================
cat("\n== Figure S6: three subtypes ==\n")
dat_s6 <- agesex %>%
  filter(cause_name %in% names(subtype_map), !is.na(measure)) %>%
  mutate(subtype = factor(unname(subtype_map[cause_name]),
                          levels = unname(subtype_map)))

# 每列仅约 22 mm。刻度画在色块分割线上、标整岁年龄。
# 选点约束: 相邻两列之间的刻度间距（末刻度→列右缘 + 列间距 + 下一列左缘→首刻度）
# 必须 ≥ 列内相邻刻度间距，否则跨列的两个数字会挤在一起。
# 取 30 / 50 / 70 / 90: 列内间距 4 格，跨列间距 2 + 3 = 5 格（再加 1 mm 列间距），
# 且与左右边缘各留 3 格和 2 格，接近对称。
S6_AGE_TICKS <- c(30, 50, 70, 90)

# 校验（单位: 格）
{
  pp    <- age_pos(S6_AGE_TICKS)
  within <- min(diff(pp))
  cross  <- (length(age_levels) + 0.5 - max(pp)) + (min(pp) - 0.5)
  cat(sprintf("  年龄刻度校验: 列内 %.1f 格, 跨列 %.1f 格 -> %s\n",
              within, cross, if (cross >= within) "通过" else "不通过"))
  stopifnot(cross >= within,
            min(pp) >= 0.5, max(pp) <= length(age_levels) + 0.5)
}

save_panel(build_heat(dat_s6, nested = TRUE, x_mode = "boundary",
                      x_ages = S6_AGE_TICKS),
           "FigureS6", w_mm = 180, h_mm = 120)   # 单面板, 直接就是投稿版

cat("\n出片检查: pdffonts 应显示 emb=yes；pdftotext 后应无 U+2212，只剩 en dash\n")