# =============================================================================
# Reproducibility repository - Wang et al., transportability-aware triangulation
# of metabolic stroke evidence across the Shanghai Cooperation Organisation.
# Paths are resolved by config.R; see README.md.
# =============================================================================
if (!exists("RAWDAT")) {
  if (file.exists("config.R")) source("config.R") else
    stop("Run this script from the code/ directory (the folder containing config.R).")
}

setwd(RAWDAT)

###############################################################################
# C2  GBD 空间平滑敏感性分析  (PAF 风险结构聚类层, 仅针对 Mantel)   v4
#
# 针对的批评：Mantel r=0·50 所反映的地理结构，可能由 GBD 在数据稀疏国家借用
# 邻近地区信息 (ST-GPR、CODEm 的时空平滑) 所制造。
#
# 四层分析（按证据强度 = 报告顺序）
#   S0  全集基线                必须复现 r=0·50
#   S1  层级检验 + 偏 Mantel + 嵌套结构   ★主分析，直接检验平滑假说
#   S2  单一 super-region 内部的地理梯度   平滑对其成员一视同仁，故此梯度
#                                          不可能由 super-region 归属产生
#   S3  留一法刀切              无主观判据
#   S4  低质量死因数据国家排除  官方星级 <=2：仅印度、巴基斯坦
#
# 为何不做基于不确定区间的敏感性分析（附录须载明，勿仅默默省略）
#   PAF 的不确定性主要来自各国共用的成分（全球合并 RR、TMREL 分布），而本文
#   的聚类做的是按危险因素标准化后的国家间相对剖面，共用成分在该对比中被消
#   去；剩余的国家特异成分无法从已发表的汇总结果中分离。此外，UI 宽度与官方
#   死因数据星级不相关（见文末诊断），因此也不能用作数据质量的代理。
#   本脚本仅保留该相关系数的计算，不做任何 UI 扰动分析。
#
# 关键背景：GBD 的借力沿 region / super-region 层级进行，不沿陆地边界。
#   中国 East Asia；印巴 South Asia；白俄俄罗 Eastern Europe；
#   哈吉塔乌 Central Asia；伊朗 North Africa & Middle East。
#   印巴的估计借自孟加拉、尼泊尔、阿富汗等南亚国家，无一进入本分析。
#
# 解释上的告诫：GBD 的 super-region 本身即按流行病学特征划分，而这些特征与
# 祖源、环境等地理相关因素纠缠，因此"层级解释力强于距离"同时相容于伪影读法
# 与构念效度读法。S1 的嵌套结构与 S2 分别从两侧收窄这个歧义，但无法消除它；
# 附录必须如实说明。
#
# 依赖：PAF_v2.csv
# 输出：C2_output/
###############################################################################

library(vegan)
library(tidyverse)
library(geosphere)
library(cluster)
library(mclust)

SEED_C2    <- 20260823   # 独立于主脚本 SEED=2023，不干扰其随机数流
NPERM_BASE <- 99999      # S0，与主脚本一致
NPERM_SENS <- 99999       # S1/S2/S3/S4
K_RANGE    <- 2:5

OUT_DIR <- file.path(RAWDAT, "大改/C2_output"); dir.create(OUT_DIR, showWarnings = FALSE)
OUT_DIR_final <- file.path(RAWDAT, "大改")
DEC <- "\u00B7"; dec <- function(x) gsub("\\.", DEC, as.character(x))

# ============================================================================
# 1. 首都距离与 GBD 层级
# ============================================================================
sco_countries <- c("Belarus", "China", "India", "Iran", "Kazakhstan",
                   "Kyrgyzstan", "Pakistan", "Russia", "Tajikistan", "Uzbekistan")

# 坐标原样照抄主脚本 07_paf_spatial_tests.R，勿改
capitals <- data.frame(
  country = sco_countries,
  lon = c(27.5667, 116.4074, 77.2090, 51.3890, 71.4491,
          74.5698, 73.0479, 37.6173, 68.7870, 69.2401),
  lat = c(53.9000, 39.9042, 28.6139, 35.6892, 51.1605,
          42.8746, 33.6844, 55.7558, 38.5598, 41.2995)
)
rownames(capitals) <- capitals$country
geo_dist_of <- function(cty)
  as.dist(distm(capitals[cty, c("lon", "lat")], fun = distHaversine) / 1000)

# 星级出处：GBD 2016 Causes of Death Collaborators. Lancet 2017; 390: 1151-210.
gbd_meta <- tribble(
  ~country,      ~region,                      ~superregion,                    ~stars,
  "Belarus",     "Eastern Europe",             "CEE, EE, and Central Asia",     4L,
  "China",       "East Asia",                  "SEA, East Asia, and Oceania",   3L,
  "India",       "South Asia",                 "South Asia",                    2L,
  "Iran",        "North Africa & Middle East", "North Africa and Middle East",  3L,
  "Kazakhstan",  "Central Asia",               "CEE, EE, and Central Asia",     4L,
  "Kyrgyzstan",  "Central Asia",               "CEE, EE, and Central Asia",     4L,
  "Pakistan",    "South Asia",                 "South Asia",                    1L,
  "Russia",      "Eastern Europe",             "CEE, EE, and Central Asia",     5L,
  "Tajikistan",  "Central Asia",               "CEE, EE, and Central Asia",     3L,
  "Uzbekistan",  "Central Asia",               "CEE, EE, and Central Asia",     4L
)

# 层级距离：同 region = 0；同 super-region 不同 region = 1；不同 super-region = 2
hier_dist_of <- function(cty) {
  g <- gbd_meta[match(cty, gbd_meta$country), ]
  n <- length(cty); M <- matrix(0, n, n, dimnames = list(cty, cty))
  for (i in seq_len(n)) for (j in seq_len(n))
    M[i, j] <- if (g$region[i] == g$region[j]) 0L
  else if (g$superregion[i] == g$superregion[j]) 1L else 2L
  as.dist(M)
}

# region 距离（用于 S2 内部的偏 Mantel）：同 region = 0，否则 1
region_dist_of <- function(cty) {
  g <- gbd_meta[match(cty, gbd_meta$country), ]
  n <- length(cty); M <- matrix(0, n, n, dimnames = list(cty, cty))
  for (i in seq_len(n)) for (j in seq_len(n))
    M[i, j] <- as.integer(g$region[i] != g$region[j])
  as.dist(M)
}

# ============================================================================
# 2. 单次分析
# ============================================================================
mantel_suite <- function(feat_matrix, cty, nperm = NPERM_SENS, k_range = K_RANGE) {
  fm <- feat_matrix[, cty, drop = FALSE]
  feat_z <- scale(t(fm))                 # 按危险因素标准化，国家在行（同主脚本）
  feat_z[is.na(feat_z)] <- 0
  feat_dist <- dist(feat_z, method = "euclidean")
  mt <- vegan::mantel(geo_dist_of(cty), feat_dist,
                      method = "pearson", permutations = nperm)
  hc  <- hclust(feat_dist, method = "ward.D2")
  kr  <- k_range[k_range < length(cty)]
  sil <- sapply(kr, function(ki) mean(silhouette(cutree(hc, ki), feat_dist)[, "sil_width"]))
  best_k <- kr[which.max(sil)]
  list(n = length(cty), feat_dist = feat_dist,
       mantel_r = mt$statistic, mantel_p = mt$signif,
       k = best_k, sil = max(sil), clusters = cutree(hc, best_k))
}

as_row <- function(res, label, dropped = "\u2014", ari = NA_real_)
  tibble(analysis = label, dropped = dropped, n = res$n,
         mantel_r = res$mantel_r, mantel_p = res$mantel_p,
         k = res$k, sil = res$sil, ari_vs_base = ari)

# ============================================================================
# 3. 数据
# ============================================================================
raw <- read_csv("PAF_v2.csv", show_col_types = FALSE)

rf_short <- c("High systolic blood pressure" = "SBP",
              "High LDL cholesterol"         = "LDL",
              "Kidney dysfunction"           = "KD",
              "High fasting plasma glucose"  = "FG",
              "High body-mass index"         = "BMI")

pafui <- raw %>%
  filter(metric_name == "Percent") %>%
  mutate(location_name = case_when(
    location_name == "Russian Federation"         ~ "Russia",
    location_name == "Iran (Islamic Republic of)" ~ "Iran",
    TRUE ~ location_name),
    rei_short = rf_short[rei_name]) %>%
  filter(!is.na(rei_short)) %>%
  select(location_name, cause_name, rei_short, year,
         any_of(c("val", "upper", "lower")))

mat_of <- function(cause, col, rfs, yr = 2023) {
  m <- pafui %>%
    filter(cause_name == cause, year == yr, rei_short %in% rfs) %>%
    select(location_name, rei_short, all_of(col)) %>%
    pivot_wider(names_from = location_name, values_from = all_of(col)) %>%
    column_to_rownames("rei_short")
  as.matrix(m[rfs, sco_countries, drop = FALSE])
}

# ============================================================================
# 4. 主流程
# ============================================================================
run_C2 <- function(cause, rfs, tag) {
  
  cat("\n\n############ C2:", tag, "############\n")
  P <- mat_of(cause, "val", rfs)
  
  # ---------------- S0 基线 ----------------
  set.seed(SEED_C2)
  b <- mantel_suite(P, sco_countries, nperm = NPERM_BASE)
  cat(sprintf("[S0] n=10  Mantel r=%.3f  p=%.5f  |  k=%d  sil=%.3f\n",
              b$mantel_r, b$mantel_p, b$k, b$sil))
  # Mantel 的 r 是确定性统计量，只有 p 来自置换；r 必须精确复现。
  if (tag == "Stroke (Total)" && abs(b$mantel_r - 0.50) > 0.015)
    warning("S0 未复现已发表的 Mantel r=0\u00B750，请先核对输入矩阵。")
  base_cl <- b$clusters
  rows <- as_row(b, paste0("S0 full set [", tag, "]"), "\u2014", 1)
  
  # ---------------- S1 层级检验 + 偏 Mantel + 嵌套结构（主分析） ----------------
  set.seed(SEED_C2 + 4)
  dg <- geo_dist_of(sco_countries); dh <- hier_dist_of(sco_countries); df <- b$feat_dist
  
  m_geo   <- vegan::mantel(dg, df, method = "pearson", permutations = NPERM_SENS)
  m_hier  <- vegan::mantel(dh, df, method = "pearson", permutations = NPERM_SENS)
  m_gh    <- vegan::mantel(dg, dh, method = "pearson", permutations = NPERM_SENS)
  m_part  <- vegan::mantel.partial(dg, df, dh, method = "pearson", permutations = NPERM_SENS)
  m_parth <- vegan::mantel.partial(dh, df, dg, method = "pearson", permutations = NPERM_SENS)
  
  cat("\n[S1] 平滑假说的直接检验  \u2605 主分析\n")
  cat(sprintf("  地理 ~ PAF                     r=%.3f  p=%.5f\n", m_geo$statistic,  m_geo$signif))
  cat(sprintf("  GBD 层级 ~ PAF                 r=%.3f  p=%.5f\n", m_hier$statistic, m_hier$signif))
  cat(sprintf("  地理 ~ GBD 层级（共线性）      r=%.3f  p=%.5f\n", m_gh$statistic,   m_gh$signif))
  cat(sprintf("  地理 ~ PAF | 控制 GBD 层级     r=%.3f  p=%.5f\n", m_part$statistic,  m_part$signif))
  cat(sprintf("  GBD 层级 ~ PAF | 控制地理      r=%.3f  p=%.5f\n", m_parth$statistic, m_parth$signif))
  if (m_gh$statistic > 0.7)
    cat("  [注] 地理与层级高度共线，偏 Mantel 检验力受限。\n")
  
  # --- 嵌套结构：簇是否为若干完整 super-region 的并集？---
  # 平滑只能在 super-region 内借力，因此它能解释 super-region 内部的同质，
  # 但无法解释把两个 super-region 并入同一簇。若无任何 super-region 被拆开、
  # 且每个簇都跨越多个 super-region，则该划分不可能由借力单独产生。
  g <- gbd_meta[match(sco_countries, gbd_meta$country), ]
  nest <- tibble(country = sco_countries, cluster = base_cl[sco_countries],
                 region = g$region, superregion = g$superregion)
  tab_sr <- table(nest$superregion, nest$cluster)
  split_sr <- rownames(tab_sr)[rowSums(tab_sr > 0) > 1]
  cat("\n  簇 x super-region 交叉表：\n"); print(tab_sr)
  cat(sprintf("  被拆开的 super-region：%s\n",
              if (length(split_sr)) paste(split_sr, collapse = ", ") else "无"))
  for (ci in sort(unique(nest$cluster)))
    cat(sprintf("  C%d 跨越 %d 个 super-region：%s\n", ci,
                n_distinct(nest$superregion[nest$cluster == ci]),
                paste(sort(unique(nest$superregion[nest$cluster == ci])), collapse = " | ")))
  
  perm_ari <- function(lab, ref, nperm = NPERM_SENS) {
    obs <- mclust::adjustedRandIndex(lab, ref)
    nul <- replicate(nperm, mclust::adjustedRandIndex(sample(lab), ref))
    c(ari = obs, p = (sum(nul >= obs) + 1) / (nperm + 1))
  }
  a_reg <- perm_ari(base_cl[sco_countries], g$region)
  a_sup <- perm_ari(base_cl[sco_countries], g$superregion)
  cat(sprintf("\n  簇 vs GBD region       ARI=%.3f  p=%.5f\n", a_reg["ari"], a_reg["p"]))
  cat(sprintf("  簇 vs GBD super-region ARI=%.3f  p=%.5f\n", a_sup["ari"], a_sup["p"]))
  
  s1 <- tibble(outcome = tag,
               test = c("geo~PAF", "hier~PAF", "geo~hier", "geo~PAF | hier",
                        "hier~PAF | geo", "cluster vs region (ARI)",
                        "cluster vs superregion (ARI)"),
               statistic = c(m_geo$statistic, m_hier$statistic, m_gh$statistic,
                             m_part$statistic, m_parth$statistic,
                             a_reg["ari"], a_sup["ari"]),
               p = c(m_geo$signif, m_hier$signif, m_gh$signif,
                     m_part$signif, m_parth$signif, a_reg["p"], a_sup["p"]))
  write_csv(s1, file.path(OUT_DIR, paste0("S1_hierarchy_", make.names(tag), ".csv")))
  write_csv(nest, file.path(OUT_DIR, paste0("S1_nesting_", make.names(tag), ".csv")))
  
  # ---------------- S2 单一 super-region 内部的地理梯度 ----------------
  # 同一 super-region 内的成员在借力上被一视同仁，因此其内部若仍存在
  # "距离越近剖面越像"，该梯度不可能由 super-region 归属产生。
  # 注意：n 很小，检验力有限；报告方向与量级，不以显著与否下结论。
  sr_tab <- sort(table(g$superregion), decreasing = TRUE)
  sr_big <- names(sr_tab)[1]
  cty_in <- sco_countries[g$superregion == sr_big]
  cat(sprintf("\n[S2] 最大 super-region 内部的地理梯度：%s（n=%d，国家对 %d）\n",
              sr_big, length(cty_in), choose(length(cty_in), 2)))
  s2 <- NULL
  if (length(cty_in) >= 5) {
    set.seed(SEED_C2 + 5)
    fz <- scale(t(P[, cty_in, drop = FALSE])); fz[is.na(fz)] <- 0
    dfin <- dist(fz, method = "euclidean")
    mw  <- vegan::mantel(geo_dist_of(cty_in), dfin,
                         method = "pearson", permutations = NPERM_SENS)
    mwp <- vegan::mantel.partial(geo_dist_of(cty_in), dfin, region_dist_of(cty_in),
                                 method = "pearson", permutations = NPERM_SENS)
    cat(sprintf("  地理 ~ PAF（区内）             r=%.3f  p=%.5f\n", mw$statistic,  mw$signif))
    cat(sprintf("  地理 ~ PAF | 控制 region（区内）r=%.3f  p=%.5f\n", mwp$statistic, mwp$signif))
    cat("  成员：", paste(cty_in, collapse = ", "), "\n")
    cat("  [注] 检验力受 n 限制，此处报告方向与量级，不据显著性下结论。\n")
    s2 <- tibble(outcome = tag, superregion = sr_big, n = length(cty_in),
                 test = c("geo~PAF within", "geo~PAF | region within"),
                 statistic = c(mw$statistic, mwp$statistic),
                 p = c(mw$signif, mwp$signif))
    write_csv(s2, file.path(OUT_DIR, paste0("S2_within_sr_", make.names(tag), ".csv")))
  } else cat("  成员不足 5 国，跳过。\n")
  
  # ---------------- S3 留一法 ----------------
  set.seed(SEED_C2 + 2)
  loo <- map_dfr(sco_countries, function(x) {
    keep <- setdiff(sco_countries, x)
    s <- mantel_suite(P, keep)
    as_row(s, paste0("S3 drop ", x), x,
           mclust::adjustedRandIndex(s$clusters[keep], base_cl[keep]))
  })
  cat("\n[S3] 留一法刀切（按 Mantel r 升序）\n")
  print(loo %>% select(dropped, n, mantel_r, mantel_p, k, ari_vs_base) %>%
          mutate(across(where(is.numeric), ~round(.x, 4))) %>%
          arrange(mantel_r) %>% as.data.frame(), row.names = FALSE)
  
  # ---------------- S4 低质量死因数据国家排除 ----------------
  set.seed(SEED_C2 + 3)
  keep <- setdiff(sco_countries, c("India", "Pakistan"))
  s4 <- mantel_suite(P, keep)
  exr <- as_row(s4, "S4 CoD quality <=2 stars", "India, Pakistan",
                mclust::adjustedRandIndex(s4$clusters[keep], base_cl[keep]))
  cat("\n[S4] 排除官方星级 <=2 的国家（印度 2\u2605，巴基斯坦 1\u2605）\n")
  print(exr %>% select(analysis, n, dropped, mantel_r, mantel_p, k, ari_vs_base) %>%
          mutate(across(where(is.numeric), ~round(.x, 4))) %>% as.data.frame(),
        row.names = FALSE)
  cat("  注：移走的是 C1 三国中的两国，剩余划分必然重排，ARI 列不可解读；\n")
  cat("      其价值只在 Mantel 的连续梯度。n=8 时国家对数 28。\n")
  
  all_rows <- bind_rows(rows, loo, exr)
  write_csv(all_rows, file.path(OUT_DIR, paste0("C2_", make.names(tag), ".csv")))
  invisible(list(rows = all_rows, s1 = s1, s2 = s2, nest = nest, base = b))
}

# ============================================================================
# 5. 运行
# ============================================================================
rf_order <- c("SBP", "LDL", "KD", "FG", "BMI")

res_stroke <- run_C2("Stroke", rf_order, "Stroke (Total)")
res_is  <- run_C2("Ischemic stroke",          rf_order,                 "Ischemic stroke")
res_ich <- run_C2("Intracerebral hemorrhage", c("SBP","KD","FG","BMI"), "ICH")
res_sah <- run_C2("Subarachnoid hemorrhage",  c("SBP","BMI"),           "SAH")

# ============================================================================
# 6. 诊断：UI 宽度与死因数据星级的关系
#    唯一保留的 UI 相关计算。不做任何 UI 扰动分析，理由见文件头。
#    用途：为附录中"UI 宽度不能代理数据质量"一句提供数字。
# ============================================================================
if (all(c("upper", "lower") %in% names(pafui))) {
  eps <- 1e-6
  lgv <- function(x) { y <- pmin(pmax(x, eps), 1 - eps); log(y / (1 - y)) }
  d_ui <- pafui %>%
    filter(cause_name == "Stroke", year == 2023, rei_short %in% rf_order) %>%
    mutate(ui_sd = (lgv(upper) - lgv(lower)) / (2 * qnorm(0.975))) %>%
    group_by(location_name) %>%
    summarise(ui_sd = median(ui_sd), .groups = "drop") %>%
    left_join(gbd_meta %>% select(country, stars), by = c("location_name" = "country")) %>%
    arrange(desc(ui_sd))
  rho_sp <- suppressWarnings(cor(d_ui$ui_sd, d_ui$stars, method = "spearman"))
  cat(sprintf("\n[诊断] UI 宽度 vs 官方死因数据星级：Spearman rho = %.3f（n=%d）\n",
              rho_sp, nrow(d_ui)))
  cat("        两者不相关，故 UI 宽度不能用作数据质量的代理。\n")
  print(as.data.frame(d_ui %>% mutate(ui_sd = round(ui_sd, 3))), row.names = FALSE)
  write_csv(d_ui, file.path(OUT_DIR, "diag_ui_vs_stars.csv"))
} else {
  cat("\n[诊断] 输入缺少 upper/lower，跳过 UI-星级相关性计算。\n")
}

# ============================================================================
# 7. 附录汇总
# ============================================================================
bind_rows(res_stroke$rows, res_is$rows, res_ich$rows, res_sah$rows) %>%
  transmute(Analysis = analysis, `Countries excluded` = dropped, n = n,
            `Mantel r` = dec(sprintf("%.3f", mantel_r)),
            `Mantel p` = dec(sprintf("%.5f", mantel_p)),
            `Optimal k` = k,
            `ARI vs primary` = dec(sprintf("%.3f", ari_vs_base))) %>%
  write_csv(file.path(OUT_DIR_final, "C2_summary_table.csv"))

bind_rows(res_stroke$s1, res_is$s1, res_ich$s1, res_sah$s1) %>%
  write_csv(file.path(OUT_DIR_final, "C2_hierarchy_table.csv"))

bind_rows(res_stroke$s2, res_is$s2, res_ich$s2, res_sah$s2) %>%
  write_csv(file.path(OUT_DIR_final, "C2_within_superregion_table.csv"))

cat("\n输出：C2_summary_table.csv / C2_hierarchy_table.csv /",
    "C2_within_superregion_table.csv / S1_nesting_*.csv / diag_ui_vs_stars.csv\n")
cat(sprintf("种子 C2=%d；S0 置换 %d 次，S1-S4 置换 %d 次\n",
            SEED_C2, NPERM_BASE, NPERM_SENS))
cat("\n========== C2 完成 ==========\n")