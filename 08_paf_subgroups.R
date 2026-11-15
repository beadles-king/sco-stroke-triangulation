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
# Figures S12 / S13：代谢归因 PAF 的年龄与性别亚组偏移（ΔPAF）
#
# 相对 PAF_subgroup_age_sex_.R 的重构：
#
#  【D1】四个结局合成一张图，不再每个结局各出一张。
#        原做法每个结局一张 PDF，S12/S13 各需 4 张拼接，必然超出 180×240 mm，
#        只能拆图改编号。改用 facet_grid(outcome ~ subgroup, scales = "free_y",
#        space = "free_y")：ICH 的 4 行、SAH 的 2 行按行数分配高度，不被拉伸。
#
#  【D2】格内三行删到一行。原来每格是
#            +0.037 / (0.571) / Higher
#        第三行与 ΔPAF 的正负号、与填色完全等价，同一件事说了三遍；
#        第二行的总体 PAF 在 Figure 2A（总卒中）与 Figure S9（亚型）里已逐格
#        印过。删掉两行后格高降到三分之一，这才是原版超限的真正来源。
#
#  【D3】单位改百分点、一位小数。全文 PAF 一律百分数一位小数（v3 §2.3），
#        原版用 0–1 比例三位小数，读者要在 "PAF 53·5%" 与 "ΔPAF +0·037"
#        之间换算。附录第 217 段已同步改为百分点。
#
#  【D4】聚类成员改从 Table S11 取，并加断言。原脚本 cause_config 的成员是
#        照热图截图手抄的（其注释亦如此说明），脑出血多算了 Tajikistan、
#        蛛网膜下腔多算了 Iran，与 Table S11 和 Figure S10 冲突。
#
#  【D5】列顺序统一为总卒中的聚类顺序，与 Figure 2A/2B 一致。四个结局共用
#        一条 x 轴是合成单图的前提；各结局自身的聚类归属由每块顶部的注释条
#        表示，不再靠列序表达。
#
#  【D6】出图改 ragg 双分辨率 TIF，弃用 pdf()（见 v3 陷阱 1）。
#
# 输出：FigureS12_delta_paf_age / FigureS13_delta_paf_sex，各 600 与 300 dpi
###############################################################################

setwd(RAWDAT)

library(tidyverse)
library(ragg)
if (!requireNamespace("ggnewscale", quietly = TRUE)) install.packages("ggnewscale")
library(ggnewscale)

# ---- 全文格式常数 ----------------------------------------------------------
BASE_SIZE <- 8
FONT      <- "Helvetica"
DEC       <- "\u00b7"
COL_C1    <- "#E6550D"   # 与 Figure 2A/2B、Figure S10 同色
COL_C2    <- "#7FBFA0"
COL_NEG   <- "#2C5985"   # 全文主蓝
COL_POS   <- "#A63603"   # 全文橙红

dec <- function(x) gsub("\\.", DEC, as.character(x))

OUT_DIR <- file.path(RAWDAT, "PAF")

save_panel <- function(p, stem, w_mm, h_mm) {
  dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
  for (r in c(600, 300)) {
    f <- file.path(OUT_DIR,
                   if (r == 600) paste0(stem, ".tif") else paste0(stem, "_300dpi.tif"))
    agg_tiff(f, width = w_mm, height = h_mm, units = "mm",
             res = r, compression = "lzw", background = "white", scaling = 1)
    print(p); dev.off()
    cat(sprintf("  已保存 %s (%.1f x %.1f mm @ %d dpi)\n", f, w_mm, h_mm, r))
  }
}

# ============================================================================
# 0. 数据
# ============================================================================
unify_country <- function(df) df %>% mutate(location_name = case_when(
  location_name == "Russian Federation"         ~ "Russia",
  location_name == "Iran (Islamic Republic of)" ~ "Iran",
  TRUE ~ location_name))

rf_short <- c("High systolic blood pressure" = "SBP",
              "High LDL cholesterol"         = "LDL",
              "Kidney dysfunction"           = "KD",
              "High fasting plasma glucose"  = "FG",
              "High body-mass index"         = "BMI")

# GBD 原始导出共 20 列，PAF 在 val 列。三份文件必须同 measure、同 metric、
# 同年份，否则相减没有意义 —— 下面逐项断言。
MEASURE <- "DALYs (Disability-Adjusted Life Years)"
YEAR    <- 2023

# PAF_v2.csv 同时含 Percent（PAF）与 Rate（归因 ASDR，即 Figure 2A 下半块），
# 故此处是"先筛后验"，不能断言全部为 Percent。
read_gbd <- function(file, value_name, keys) {
  d <- read_csv(file, show_col_types = FALSE)
  need <- c("location_name", "sex_name", "age_name", "cause_name",
            "rei_name", "measure_name", "metric_name", "year", "val")
  miss <- setdiff(need, names(d))
  if (length(miss))
    stop(sprintf("%s 缺少列：%s", file, paste(miss, collapse = ", ")))
  
  d <- d %>% filter(metric_name == "Percent", year == YEAR)
  if (!nrow(d))
    stop(sprintf("%s 里没有 %d 年的 Percent 行。", file, YEAR))
  if (!all(d$measure_name == MEASURE))
    stop(sprintf("%s 的 measure 不是 %s，与主分析口径不符。", file, MEASURE))
  
  d <- d %>% unify_country() %>% mutate(rf = rf_short[rei_name])
  if (any(is.na(d$rf)))
    stop(sprintf("%s 里有无法映射的 rei_name：%s", file,
                 paste(unique(d$rei_name[is.na(d$rf)]), collapse = ", ")))
  
  # 主键重复会让后面的 inner_join 静默做多对多展开，必须先拦下
  k <- d %>% select(all_of(keys)) %>% duplicated()
  if (any(k))
    stop(sprintf("%s 的主键 (%s) 有 %d 处重复。",
                 file, paste(keys, collapse = " + "), sum(k)))
  
  cat(sprintf("  %-14s Percent 行 %3d，国家 %d，cause×rei 组合 %d\n",
              file, nrow(d), n_distinct(d$location_name),
              n_distinct(paste(d$cause_name, d$rf))))
  d %>% rename(!!value_name := val)
}

cat("读取 GBD 导出：\n")
dat_main <- read_gbd('PAF_v2.csv', "paf_overall",
                     c("location_name", "cause_name", "rf")) %>%
  select(location_name, cause_name, rf, paf_overall)
stopifnot(all(dat_main$paf_overall >= 0), all(dat_main$paf_overall <= 1))

# 性别亚组：age_name 应为 Age-standardized，与主分析同口径
dat_sex <- read_gbd('PAF_sex.csv', "paf_sub",
                    c("location_name", "sex_name", "cause_name", "rf"))
stopifnot(all(dat_sex$age_name == "Age-standardized"),
          setequal(dat_sex$sex_name, c("Male", "Female")))
dat_sex <- dat_sex %>% mutate(subgroup = sex_name)

# 年龄亚组：sex_name 应为 Both；age_name 是具体年龄段，按定义不作年龄标化
dat_age <- read_gbd('PAF_age.csv', "paf_sub",
                    c("location_name", "age_name", "cause_name", "rf"))
stopifnot(all(dat_age$sex_name == "Both"),
          setequal(dat_age$age_name, c("50-69 years", "70+ years")))
dat_age <- dat_age %>% mutate(subgroup = age_name)

cat(sprintf("口径：%s，%d 年，Percent\n", MEASURE, YEAR))

# ============================================================================
# 1. 结局配置：列顺序统一，聚类成员取自 Table S11【D4】【D5】
# ============================================================================
COUNTRY_ORDER <- c("India", "China", "Pakistan",              # C1（总卒中）
                   "Kazakhstan", "Kyrgyzstan", "Belarus", "Russia",
                   "Iran", "Tajikistan", "Uzbekistan")

# Table S11 的 C1 成员。含 China 的簇一律为 C1（命名惯例见 v5 §5）
C1_MEMBERS <- list(
  "Stroke"                   = c("China", "India", "Pakistan"),
  "Ischemic stroke"          = c("China", "India", "Pakistan", "Tajikistan", "Uzbekistan"),
  "Intracerebral hemorrhage" = c("China", "India", "Pakistan"),
  "Subarachnoid hemorrhage"  = c("China", "India", "Pakistan")
)

OUTCOME_LABEL <- c("Stroke"                   = "Total stroke",
                   "Ischemic stroke"          = "Ischaemic stroke",
                   "Intracerebral hemorrhage" = "Intracerebral haemorrhage",
                   "Subarachnoid hemorrhage"  = "Subarachnoid haemorrhage")

# 分面左侧的显示名。首版单行旋转 90°，25 个字符沿分块高度排不下，被截成
# "tracerebral haemorrhag"；改水平后不再截断，但占去约 28 mm 宽度，20 列每格
# 只剩 6·5 mm，"+14·2" 几乎顶到格边。现改为竖排 + 两行折行：最长的一行
# "Intracerebral"（13 字符，约 18 mm）沿分块高度排布，脑出血块 5 行约 33 mm、
# 蛛网膜下腔块 3 行约 20 mm，均放得下；宽度则降到两行行高约 7 mm。
STRIP_LABEL <- c("Total stroke"              = "Total\nstroke",
                 "Ischaemic stroke"          = "Ischaemic\nstroke",
                 "Intracerebral haemorrhage" = "Intracerebral\nhaemorrhage",
                 "Subarachnoid haemorrhage"  = "Subarachnoid\nhaemorrhage")

# 各结局纳入的风险因子（GBD 只对部分亚型估计部分因子）
RF_BY_OUTCOME <- list(
  "Stroke"                   = c("SBP", "LDL", "KD", "FG", "BMI"),
  "Ischemic stroke"          = c("SBP", "LDL", "KD", "FG", "BMI"),
  "Intracerebral hemorrhage" = c("SBP", "KD", "FG", "BMI"),
  "Subarachnoid hemorrhage"  = c("SBP", "BMI")
)

stopifnot(setequal(names(C1_MEMBERS), names(RF_BY_OUTCOME)),
          all(unlist(C1_MEMBERS) %in% COUNTRY_ORDER))
cat("聚类成员（对照 Table S11）：\n")
for (o in names(C1_MEMBERS))
  cat(sprintf("  %-26s C1 = %s\n", OUTCOME_LABEL[[o]],
              paste(sort(C1_MEMBERS[[o]]), collapse = ", ")))

# ============================================================================
# 2. 组装长表
# ============================================================================
build_delta <- function(dat_sub, sub_levels, sub_labels) {
  out <- lapply(names(RF_BY_OUTCOME), function(o) {
    rfs <- RF_BY_OUTCOME[[o]]
    ref <- dat_main %>% filter(cause_name == o, rf %in% rfs) %>%
      select(location_name, rf, paf_overall)
    dat_sub %>%
      filter(cause_name == o, rf %in% rfs, subgroup %in% sub_levels) %>%
      select(location_name, subgroup, rf, paf_sub) %>%
      inner_join(ref, by = c("location_name", "rf")) %>%
      mutate(outcome = o,
             # 【D3】百分点。GBD 的 PAF 是 0–1 比例，总体值同尺度，故 ×100
             delta_pp = (paf_sub - paf_overall) * 100)
  }) %>% bind_rows()
  
  miss <- out %>% filter(is.na(delta_pp))
  if (nrow(miss)) stop(sprintf("有 %d 格缺失，先查数据再出图。", nrow(miss)))
  
  # 每个结局 × 每个亚组 × 10 国 × 该结局的因子数，一格都不能少
  expect <- sum(vapply(RF_BY_OUTCOME, length, 1L)) * length(sub_levels) * 10
  if (nrow(out) != expect)
    stop(sprintf("组装出 %d 行，应为 %d 行，请查数据完整性。", nrow(out), expect))
  
  out %>% mutate(
    country  = factor(location_name, levels = COUNTRY_ORDER),
    subgroup = factor(unname(sub_labels[subgroup]), levels = unname(sub_labels)),
    outcome  = factor(unname(OUTCOME_LABEL[outcome]), levels = unname(OUTCOME_LABEL)),
    # y 轴在每块内独立，故用 "结局:因子" 作水平，显示时只留因子名
    row_key  = paste(outcome, rf, sep = ":")
  )
}

# 顶部聚类注释条：每个结局一行，随该结局自身的成员
build_cluster_strip <- function(delta_df) {
  lapply(names(C1_MEMBERS), function(o) {
    lab <- unname(OUTCOME_LABEL[o])
    expand_grid(country = factor(COUNTRY_ORDER, levels = COUNTRY_ORDER),
                subgroup = levels(delta_df$subgroup)) %>%
      mutate(outcome = factor(lab, levels = unname(OUTCOME_LABEL)),
             cluster = ifelse(as.character(country) %in% C1_MEMBERS[[o]], "C1", "C2"),
             row_key = paste(lab, "Cluster", sep = ":"))
  }) %>% bind_rows() %>%
    mutate(subgroup = factor(subgroup, levels = levels(delta_df$subgroup)))
}

# ============================================================================
# 3. 绘图
# ============================================================================
# 四舍五入后为零的格不带符号，避免出现 "−0·0"
fmt_dpaf <- function(v) {
  r <- round(v, 1)
  dec(ifelse(r == 0, sprintf("%.1f", 0), sprintf("%+.1f", r)))
}

draw_delta_fig <- function(delta_df, strip_df) {
  # 色标上限取 |ΔPAF| 的 95 百分位（向上取整到 0·5），超出者以极端色显示。
  # 若按最大值定标，少数离群格（如蛛网膜下腔的 +14·2）会把总卒中那一块
  # （±4 之内，正是附录引用的部分）压成近乎全白。每格都印着数字，信息不丢。
  lim  <- ceiling(quantile(abs(delta_df$delta_pp), 0.95) * 2) / 2
  n_sq <- sum(abs(delta_df$delta_pp) > lim)
  cat(sprintf("  色标 ±%.1f 百分点（95 百分位）；%d / %d 格超出范围，以极端色显示\n",
              lim, n_sq, nrow(delta_df)))
  
  # y 轴自上而下：Cluster 条在每块顶端，其后按 RF_BY_OUTCOME 的顺序。
  # 首版靠 factor levels 排序，但 scales = "free_y" 下每个分面的次序实际取决于
  # 图层出现顺序（聚类层先画，Cluster 因此落到了底部）。改为给 scale_y_discrete
  # 传 limits 函数：它拿到本分面实际出现的值，返回想要的次序，与图层顺序无关。
  lev <- unlist(lapply(names(RF_BY_OUTCOME), function(o) {
    lab <- unname(OUTCOME_LABEL[o])
    paste(lab, c("Cluster", RF_BY_OUTCOME[[o]]), sep = ":")
  }))
  y_order <- function(present) rev(lev[lev %in% present])
  
  # 格内文字按底色亮度切黑白（v5 §5：不要用固定阈值手切）
  txt_col <- ifelse(abs(delta_df$delta_pp) > lim * 0.62, "white", "grey10")
  
  ggplot() +
    # 聚类注释条
    geom_tile(data = strip_df,
              aes(country, row_key, fill = cluster),
              colour = "white", linewidth = 0.2) +
    scale_fill_manual(values = c(C1 = COL_C1, C2 = COL_C2), name = "Cluster",
                      guide = guide_legend(order = 1, keyheight = unit(3, "mm"),
                                           keywidth = unit(3, "mm"))) +
    ggnewscale::new_scale_fill() +
    # ΔPAF 热图
    geom_tile(data = delta_df,
              aes(country, row_key, fill = delta_pp),
              colour = "white", linewidth = 0.2) +
    geom_text(data = delta_df,
              aes(country, row_key,
                  label = fmt_dpaf(delta_pp)),
              # 7 pt 是全文图内字号下限；格宽约 7·5 mm，"+14·2" 约 5·7 mm，留有余量
              size = (BASE_SIZE - 1) / .pt, family = FONT, colour = txt_col) +
    scale_fill_gradient2(low = COL_NEG, mid = "white", high = COL_POS,
                         midpoint = 0, limits = c(-lim, lim),
                         oob = scales::squish,
                         name = "\u0394PAF, percentage points",
                         labels = function(x) dec(sprintf("%g", x)),
                         guide = guide_colourbar(order = 2, barheight = unit(3, "mm"),
                                                 barwidth = unit(34, "mm"),
                                                 title.position = "top")) +
    scale_y_discrete(limits = y_order,
                     labels = function(x) sub("^.*:", "", x), expand = c(0, 0)) +
    scale_x_discrete(position = "bottom", expand = c(0, 0)) +
    facet_grid(outcome ~ subgroup, scales = "free_y", space = "free_y", switch = "y",
               labeller = labeller(outcome = STRIP_LABEL)) +
    labs(x = NULL, y = NULL) +
    theme_minimal(base_size = BASE_SIZE, base_family = FONT) +
    theme(
      panel.grid        = element_blank(),
      panel.spacing.x   = unit(2, "mm"),
      panel.spacing.y   = unit(1.5, "mm"),
      axis.text.x       = element_text(size = BASE_SIZE - 1, angle = 45,
                                       hjust = 1, colour = "black"),
      axis.text.y       = element_text(size = BASE_SIZE - 1, colour = "black"),
      axis.ticks        = element_blank(),
      strip.text.x      = element_text(size = BASE_SIZE, face = "bold"),
      strip.text.y.left = element_text(size = BASE_SIZE, face = "bold", angle = 90,
                                       hjust = 0.5, lineheight = 0.9),
      strip.placement   = "outside",
      legend.position   = "bottom",
      legend.box        = "horizontal",
      legend.title      = element_text(size = BASE_SIZE - 1),
      legend.text       = element_text(size = BASE_SIZE - 1),
      plot.margin       = margin(2, 3, 2, 2, "mm")
    )
}


# ============================================================================
# 4. 出图 + 打印附录所需的均值
# ============================================================================
report_means <- function(delta_df, tag) {
  cat(sprintf("\n===== %s：总卒中各因子的平均 ΔPAF（百分点，附录第 217 段用）=====\n", tag))
  delta_df %>%
    filter(outcome == "Total stroke") %>%
    group_by(subgroup, rf) %>%
    summarise(mean_dpaf = mean(delta_pp), .groups = "drop") %>%
    mutate(mean_dpaf = sprintf("%+.1f", mean_dpaf)) %>%
    pivot_wider(names_from = rf, values_from = mean_dpaf) %>%
    as.data.frame() %>% print(row.names = FALSE)
  cat(sprintf("全体 mean |ΔPAF|（总卒中）：%s\n",
              paste(sprintf("%s = %.1f", levels(delta_df$subgroup),
                            tapply(abs(delta_df$delta_pp[delta_df$outcome == "Total stroke"]),
                                   delta_df$subgroup[delta_df$outcome == "Total stroke"], mean)),
                    collapse = " | ")))
}

# ---- Figure S12：年龄 ----
age_lab <- c("50-69 years" = "Age 50\u201369 years", "70+ years" = "Age 70 years or older")
d_age <- build_delta(dat_age, names(age_lab), age_lab)
report_means(d_age, "Figure S12（年龄）")
save_panel(draw_delta_fig(d_age, build_cluster_strip(d_age)),
           "FigureS12_delta_paf_age", 174, 168)

# ---- Figure S13：性别 ----
sex_lab <- c("Male" = "Men", "Female" = "Women")
d_sex <- build_delta(dat_sex, names(sex_lab), sex_lab)
report_means(d_sex, "Figure S13（性别）")
save_panel(draw_delta_fig(d_sex, build_cluster_strip(d_sex)),
           "FigureS13_delta_paf_sex", 174, 168)

cat("\n两张图各自单张出片，编号无需改动。\n")