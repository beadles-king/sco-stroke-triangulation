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
# SCO 10国中风亚型构成聚类分析（基于2023年数据）
# 验证亚型构成（IS/ICH/SAH）是否为代谢风险PAF异质性的机制基础
#
# 流程：
#   1. ILR 变换 + 层次聚类（Euclidean + Ward.D2，与PAF一致）
#   2. 与 PAF 聚类一致性检验：ARI + Mantel + Cophenetic + Baker's gamma
#   3. 堆叠柱状图（DALY + Incidence 并排）
#   4. 双树对比（Tanglegram）
#   5. 树状图 + 柱状图一体化展示
#
# 前置要求：
#   - 已运行 05_paf_attribution.R，环境中存在 main_res 对象
#   - main_res 包含 res_paf$hc / res_paf$clusters / res_paf$k / res_paf$dist
# ============================================================================

setwd(file.path(RAWDAT, "PAF/亚型占比"))

library(tidyverse)
library(compositions)   # ILR 变换
library(dendextend)     # 树状图比较
library(mclust)         # adjustedRandIndex
library(vegan)          # Mantel test
library(patchwork)      # 拼图
library(ggdendro)       # ggplot 版树状图
library(scales)
library(cowplot)        # plot_grid (备用)

# ---- 出图: 一律走 ragg -------------------------------------------------------
# pdf() 不嵌字体, 栅格化时字面会被替换成 Nimbus Sans; cairo_pdf 则把字形步进
# 四舍五入到整数 point, 同一单词内字距不匀。ragg 在目标分辨率排字, 600 dpi 下
# 取整误差 <= 0.042 mm, 且与其余图出自同一引擎、同一份 Helvetica。
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


# ============================================================================
# 0. 全局配置
# ============================================================================
sco_countries <- c("China", "India", "Pakistan", "Russia", "Belarus",
                   "Kazakhstan", "Kyrgyzstan", "Tajikistan", "Uzbekistan", "Iran")

subtype_levels <- c("Ischemic stroke",
                    "Intracerebral hemorrhage",
                    "Subarachnoid hemorrhage")

target_year <- 2023

country_order <- c("India", "China", "Pakistan",
                   "Kazakhstan", "Kyrgyzstan", "Belarus", "Russia",
                   "Iran", "Tajikistan", "Uzbekistan")

# 柱状图填充: IS=深青绿, ICH=橙红, SAH=浅橙
col_subtype <- c("IS" = "#1F8A70", "ICH" = "#E07560", "SAH" = "#FDD0A2")

# 两层配色: PAF=中等饱和度, 亚型聚类=深色(同色系)
# C1/SC1: 红橙系  -- 多ICH群(出血→暖色)
# C2/SC2: 青绿系  -- 多IS群 (缺血→冷色,避开海洋蓝)
# 色盲友好: 通过色相+亮度双通道编码
col_pafclu     <- c("C1" = "#E6550D", "C2" = "#7FBFA0", "C3" = "#9970AB")  # 中: PAF 注释带 + tanglegram C 树
col_subtypeclu <- c("SC1"= "#A81E13", "SC2"= "#1F4E9C", "SC3"= "#4A148C")


# ============================================================================
# 1. 数据准备：提取 2023 年 DALY rate 和 Incidence rate
# ============================================================================
gbd_raw <- read_csv(file.path(RAWDAT, "general.csv"))

# 诊断：先看哪些国家在原始数据中真的有亚型数据 ----------------------------
cat("\n=== 数据诊断：检查 location_name 实际拼写 ===\n")
loc_actual <- unique(gbd_raw$location_name)
loc_match  <- loc_actual[grepl(paste(c("Russia", "Iran", "China", "India",
                                       "Pakistan", "Belarus", "Kazakhstan",
                                       "Kyrgyzstan", "Tajikistan", "Uzbekistan"),
                                     collapse = "|"),
                               loc_actual, ignore.case = TRUE)]
print(loc_match)

# 自动适配 location_name 的别名（修复 Iran/Russia 缺失问题）
location_alias <- c(
  "China"      = "China",
  "India"      = "India",
  "Pakistan"   = "Pakistan",
  "Russia"     = grep("^Russia", loc_actual, value = TRUE)[1],
  "Belarus"    = "Belarus",
  "Kazakhstan" = "Kazakhstan",
  "Kyrgyzstan" = "Kyrgyzstan",
  "Tajikistan" = "Tajikistan",
  "Uzbekistan" = "Uzbekistan",
  "Iran"       = grep("^Iran", loc_actual, value = TRUE)[1]
)
cat("\n=== 自动匹配的国家名 ===\n")
print(location_alias)

# 检查是否有 NA（说明该国家在数据里完全找不到）
if (any(is.na(location_alias))) {
  warning("部分国家未匹配到：", paste(names(location_alias)[is.na(location_alias)],
                             collapse = ", "))
}

# 数据准备函数
prep_subtype_data <- function(gbd_raw, measure_target, year_target = 2023) {
  
  # 用别名表中的实际数据库名做筛选
  target_locations <- as.character(location_alias)
  target_locations <- target_locations[!is.na(target_locations)]
  
  result <- gbd_raw %>%
    filter(measure_name == measure_target,
           metric_name == "Rate",
           age_name == "Age-standardized",
           sex_name == "Both",
           location_name %in% target_locations,
           cause_name %in% subtype_levels,
           year == year_target) %>%
    # 把数据库的全名映射回标准短名
    mutate(country = names(location_alias)[match(location_name, location_alias)]) %>%
    select(country,
           subtype = cause_name,
           rate = val,
           rate_lower = lower,
           rate_upper = upper) %>%
    group_by(country) %>%
    mutate(prop = rate / sum(rate),
           prop_lower = rate_lower / sum(rate),
           prop_upper = rate_upper / sum(rate)) %>%
    ungroup()
  
  # 报告：哪些国家成功提取
  cat(sprintf("\n[%s] 提取国家数: %d / %d\n",
              measure_target,
              n_distinct(result$country),
              length(sco_countries)))
  missing <- setdiff(sco_countries, unique(result$country))
  if (length(missing) > 0) {
    cat("[缺失]", paste(missing, collapse = ", "), "\n")
  }
  result
}

subtype_daly <- prep_subtype_data(gbd_raw,
                                  "DALYs (Disability-Adjusted Life Years)",
                                  target_year)
subtype_inc  <- prep_subtype_data(gbd_raw, "Incidence", target_year)

# ============================================================================
# 2. 转宽表 + ILR 变换 + 层次聚类
# ============================================================================
to_wide <- function(df) {
  df %>%
    select(country, subtype, prop) %>%
    pivot_wider(names_from = subtype, values_from = prop) %>%
    rename(IS  = `Ischemic stroke`,
           ICH = `Intracerebral hemorrhage`,
           SAH = `Subarachnoid hemorrhage`) %>%
    column_to_rownames("country")
}

prop_daly_wide <- to_wide(subtype_daly)
prop_inc_wide  <- to_wide(subtype_inc)

stopifnot(all(abs(rowSums(prop_daly_wide) - 1) < 1e-6))
stopifnot(all(abs(rowSums(prop_inc_wide)  - 1) < 1e-6))

cat("\n=== DALY 占比矩阵 ===\n"); print(round(prop_daly_wide, 3))
cat("\n=== Incidence 占比矩阵 ===\n"); print(round(prop_inc_wide, 3))

# ============================================================================
# 3. 提取 PAF 聚类结果（来自 05_paf_attribution.R 的 main_res）
# ============================================================================
hc_paf           <- main_res$res_paf$hc
paf_clusters_num <- main_res$res_paf$clusters
k_paf            <- main_res$res_paf$k
d_paf            <- main_res$res_paf$dist   # PAF 的距离矩阵（用于 Mantel）

cat("\n=== PAF 聚类信息 ===\n")
cat("k_paf =", k_paf, "\n")
cat("PAF clusters:\n"); print(paf_clusters_num)

paf_cluster_assign <- paste0("C", paf_clusters_num)
names(paf_cluster_assign) <- names(paf_clusters_num)

# ============================================================================
# 4. ILR + 层次聚类（与PAF同方法：Euclidean + Ward.D2）
# ============================================================================
do_clustering <- function(prop_wide, k_force = k_paf, k_range = 2:5) {
  ilr_coords <- compositions::ilr(as.matrix(prop_wide))
  rownames(ilr_coords) <- rownames(prop_wide)
  
  d  <- dist(ilr_coords, method = "euclidean")
  hc <- hclust(d, method = "ward.D2")
  
  # silhouette 自动选最佳 k（与PAF代码一致的策略）
  sil_widths <- sapply(k_range, function(k) {
    cl <- cutree(hc, k = k)
    if (length(unique(cl)) < 2) return(NA_real_)
    mean(cluster::silhouette(cl, d)[, "sil_width"])
  })
  names(sil_widths) <- k_range
  k_best <- as.integer(names(sil_widths)[which.max(sil_widths)])
  
  cat(sprintf("\n[Silhouette] k_range = %s\n", paste(k_range, collapse = ",")))
  cat("  Silhouette widths:\n"); print(round(sil_widths, 3))
  cat(sprintf("  → k_best = %d (max sil width = %.3f)\n",
              k_best, max(sil_widths, na.rm = TRUE)))
  
  list(ilr      = ilr_coords,
       dist     = d,
       hc       = hc,
       sil_widths       = sil_widths,
       k_best           = k_best,
       cluster_k_paf    = cutree(hc, k = k_force),  # 与PAF同k(用于一致性检验)
       cluster_k_best   = cutree(hc, k = k_best))   # 数据驱动的最佳k(用于可视化)
}

clust_daly <- do_clustering(prop_daly_wide)
clust_inc  <- do_clustering(prop_inc_wide)

cat("\n=== 亚型聚类结果（k=k_paf=2，与PAF同k）===\n")
cat("DALY:\n");      print(clust_daly$cluster_k_paf)
cat("Incidence:\n"); print(clust_inc$cluster_k_paf)

cat("\n=== 亚型聚类结果（k=3，数据自然结构）===\n")
cat("DALY:\n");      print(clust_daly$cluster_k_best)
cat("Incidence:\n"); print(clust_inc$cluster_k_best)

# ============================================================================
# 5. 与 PAF 聚类的一致性检验（多指标）
# ============================================================================
test_consistency <- function(clust_obj, paf_clusters_num, d_paf,
                             label = "DALY", n_perm = 99999) {
  
  common <- intersect(names(clust_obj$cluster_k_paf), names(paf_clusters_num))
  cat(sprintf("\n========== %s vs PAF ==========\n", label))
  cat(sprintf("Common countries: %d  (%s)\n",
              length(common), paste(common, collapse = ", ")))
  
  sc <- clust_obj$cluster_k_paf[common]
  pc <- paf_clusters_num[common]
  
  # ---- 主指标: Mantel test (国家间距离矩阵相关性,最稳健) ----
  d_sub_m <- as.matrix(clust_obj$dist)[common, common]
  d_paf_m <- as.matrix(d_paf)[common, common]
  mt <- vegan::mantel(as.dist(d_sub_m), as.dist(d_paf_m),
                      method = "spearman", permutations = n_perm)
  cat(sprintf("[Mantel r, Spearman]   %.3f  (p = %.6f, %d perm)  ← primary\n",
              mt$statistic, mt$signif, n_perm))
  
  # ---- 辅助指标: ARI (传统硬分组指标) ----
  ari_obs <- mclust::adjustedRandIndex(sc, pc)
  ari_perm <- replicate(n_perm,
                        mclust::adjustedRandIndex(sc, sample(pc)))
  p_ari <- (sum(ari_perm >= ari_obs) + 1) / (n_perm + 1)
  cat(sprintf("[ARI, k=%d]             %.3f  (p = %.6f, %d perm)  ← supplementary\n",
              k_paf, ari_obs, p_ari, n_perm))
  
  # 交叉表(辅助理解)
  cat("Cross-tabulation:\n")
  print(table(PAF = pc, Subtype = sc))
  
  invisible(tibble(
    Comparison  = label,
    N_countries = length(common),
    Mantel_r    = mt$statistic, Mantel_p = mt$signif,
    ARI         = ari_obs,      ARI_p    = p_ari,
    N_perm      = n_perm
  ))
}

# 随机种子：只在此设一次，放在两次调用之外。
# 原来在函数内部 set.seed(seed)，两次调用会共用同一段随机数，DALY 与 ASIR
# 的置换分布完全相同、不再独立（v7 §4.3 同一问题）。
set.seed(42)
res_daly <- test_consistency(clust_daly, paf_clusters_num, d_paf, "DALY")
res_inc  <- test_consistency(clust_inc,  paf_clusters_num, d_paf, "Incidence")

# 不再 round：99999 次置换下 p 最小可达 1e-5，round(.,4) 会把它压成 0
summary_tbl <- bind_rows(res_daly, res_inc)
print(summary_tbl)
write_csv(summary_tbl, "consistency_test_results.csv")

# ============================================================================
# 6. 堆叠柱状图（DALY + Incidence 上下排列）
# ============================================================================
# plot_stacked <- function(df, paf_assign, country_order, title_text) {
#   pdat <- df %>%
#     mutate(country = factor(country, levels = country_order),
#            subtype = factor(subtype,
#                             levels = c("Subarachnoid hemorrhage",
#                                        "Intracerebral hemorrhage",
#                                        "Ischemic stroke"),
#                             labels = c("SAH", "ICH", "IS")),
#            paf_cluster = paf_assign[as.character(country)])
#   anno <- pdat %>% distinct(country, paf_cluster)
#   
#   p_anno <- ggplot(anno, aes(x = country, y = 1, fill = paf_cluster)) +
#     geom_tile(color = "white", linewidth = 0.5) +
#     geom_text(aes(label = paf_cluster), size = 3, fontface = "bold",
#               color = "white") +
#     scale_fill_manual(values = col_pafclu, guide = "none") +
#     scale_x_discrete(expand = c(0, 0)) +
#     scale_y_continuous(expand = c(0, 0)) +
#     theme_void() +
#     theme(plot.margin = margin(0, 0, 0, 0))
#   
#   p_main <- ggplot(pdat, aes(x = country, y = prop, fill = subtype)) +
#     geom_col(width = 0.78, color = "white", linewidth = 0.4) +
#     geom_text(aes(label = ifelse(prop >= 0.03,
#                                  percent(prop, accuracy = 1), "")),
#               position = position_stack(vjust = 0.5),
#               size = 3, color = "white", fontface = "bold") +
#     scale_y_continuous(labels = percent_format(),
#                        expand = expansion(mult = c(0, 0.02))) +
#     scale_fill_manual(values = col_subtype,
#                       breaks = c("IS", "ICH", "SAH"), name = "Subtype") +
#     labs(x = NULL, y = "Proportion", title = title_text) +
#     theme_minimal(base_size = 11) +
#     theme(axis.text.x = element_text(angle = 35, hjust = 1),
#           panel.grid.major.x = element_blank(),
#           plot.title = element_text(face = "bold", size = 12),
#           legend.position = "right")
#   
#   combined <- (p_anno / p_main) + plot_layout(heights = c(0.06, 1))
#   wrap_elements(full = combined)
# }
# 
# p_daly <- plot_stacked(subtype_daly, paf_cluster_assign, country_order,
#                        "Stroke subtype composition by DALY rate (2023)")
# p_inc  <- plot_stacked(subtype_inc, paf_cluster_assign, country_order,
#                        "Stroke subtype composition by Incidence rate (2023)")

# p_full <- p_daly / p_inc
# ggsave("subtype_stacked_DALY_Incidence_2023.pdf", p_full,
#        width = 11, height = 9)
# ggsave("subtype_stacked_DALY_Incidence_2023.png", p_full,
#        width = 11, height = 9, dpi = 300)

# ============================================================================
# 7. 双树对比 (Tanglegram)
# ============================================================================
plot_tanglegram <- function(hc_subtype, hc_paf_obj, label_subtype,
                            fname_base, k = k_paf,
                            scale = 1, width = 14 * scale, height = 7 * scale,
                            show_legend = TRUE) {
  
  # 原画布 14x7 in（裁剪后内容仍有 340x159 mm）从未对齐过"两面板拼入
  # 180x240 mm 供应图"的目标，字号/线宽是按更大的展示尺寸调的。scale 把
  # 画布物理尺寸与下面所有硬编码的 cex/lwd 同比例缩小，做原生重绘而非
  # 事后缩放位图，相对比例不变，只是整体印刷尺寸变小。
  
  # 小数号统一为中点 U+00B7，与 plot_dendro_with_bars() 的 dec_lab() 保持一致。
  # tanglegram() 内部走 base graphics 的默认 axis()，不经过 ggplot 的 labels=
  # 参数，因此不能像 dec_lab() 那样在 scale 层面替换；改用 options(OutDec=...)
  # 在数值转字符时统一替换小数点，覆盖 axis() 标签的格式化路径。仅在本函数
  # 作用域内生效，退出时用 on.exit 还原，不影响脚本其余部分。
  old_out_dec <- getOption("OutDec")
  options(OutDec = "\u00B7")
  on.exit(options(OutDec = old_out_dec), add = TRUE)
  
  common <- intersect(hc_subtype$labels, hc_paf_obj$labels)
  
  d1 <- as.dendrogram(hc_subtype) %>% prune(setdiff(labels(.), common))
  d2 <- as.dendrogram(hc_paf_obj) %>% prune(setdiff(labels(.), common))
  
  # 颜色定义(局部变量,避免作用域问题)
  sc_colors_local <- col_subtypeclu[1:k]
  c_colors_local  <- col_pafclu[1:k]
  
  # 左侧(亚型聚类): SC 配色(深)
  d1 <- d1 %>%
    set("branches_k_color", value = sc_colors_local, k = k) %>%
    set("labels_col",       value = sc_colors_local, k = k) %>%
    set("branches_lwd", 2.6 * scale) %>%
    set("labels_cex", 1.75 * scale) # 1.25 → 1.5
  
  # 右侧(PAF聚类): C 配色(中)
  d2 <- d2 %>%
    set("branches_k_color", value = c_colors_local, k = k) %>%
    set("labels_col",       value = c_colors_local, k = k) %>%
    set("branches_lwd", 2.6 * scale) %>%
    set("labels_cex", 1.75 * scale) # 1.25 → 1.5
  
  dl <- dendlist(d1, d2) %>% untangle(method = "step2side")
  
  # 内部绘图函数: tanglegram + (可选)在外 margin 添加图例
  draw_one <- function() {
    # oma: 上下各留 2 行字符空间(拉开树状图与图例距离)
    #      右侧留 12 行字符空间放图例；不画图例时右边距缩到 2 行即可
    par(oma = c(2, 0, 2, if (show_legend) 12 else 2) * scale,
        mar = c(4, 2, 4, 2) * scale,
        cex.main = 1.6 * scale,
        cex.axis = 1.3 * scale,
        cex.lab  = 1.3 * scale)
    
    tanglegram(dl,
               main_left  = label_subtype,
               main_right = "Metabolic PAF",
               common_subtrees_color_branches = FALSE,
               highlight_distinct_edges = FALSE,
               color_lines = "grey60",
               margin_inner = 9 * scale,
               lwd = 2.6 * scale,
               columns_width = c(5, 2, 5),
               cex_main = 1.6 * scale,
               cex_main_left = 1.6 * scale,
               cex_main_right = 1.6 * scale)
    
    if (!show_legend) return(invisible())
    
    # 在画布外的右侧 margin 区域绘制图例
    par(xpd = NA)
    
    legend_labels <- c(paste0("C", 1:k), paste0("SC", 1:k))
    legend_colors <- c(c_colors_local, sc_colors_local)
    
    legend(x = grconvertX(0.88, from = "ndc", to = "user"),
           y = grconvertY(0.5,  from = "ndc", to = "user"),
           legend = legend_labels,
           col = legend_colors,
           text.col = legend_colors,
           lwd = 3 * scale, lty = 1,
           seg.len = 2.0 * scale,
           bty = "n",
           cex = 1.25 * scale,
           y.intersp = 1.7,
           text.font = 2,
           xjust = 0, yjust = 0.5)
  }
  
  # ---- 600 / 300 dpi TIFF 输出 (Figure S15 的两个面板, 拼图阶段合成) ----
  for (r in c(600, 300)) {
    f <- if (r == 600) paste0(fname_base, ".tif")
    else          paste0(fname_base, "_300dpi.tif")
    ragg::agg_tiff(f, width = width, height = height, units = "in",
                   res = r, compression = "lzw", background = "white",
                   scaling = 1)
    draw_one()
    dev.off()
    cat("Saved:", f, "\n")
  }
}

plot_tanglegram(clust_daly$hc, hc_paf,
                "Subtype (ASDR %, ILR)", "tanglegram_ASDR", scale = 0.5,
                show_legend = FALSE)
plot_tanglegram(clust_inc$hc, hc_paf,
                "Subtype (ASIR %, ILR)", "tanglegram_ASIR", scale = 0.5,
                show_legend = FALSE)

# ---- 独立图例(Figure S15 两面板共用一份, 拼图阶段居中放置) ----
# C1/C2/SC1/SC2 两面板配色相同, 不必各画一份; 单独出图省得拼图时要去
# 裁"图例和树状图之间的空白"——右侧根分支横线和图例色块在同一列范围内
# 交叠(只是行不同), 裁图没法保证不切到线或字, 干脆各自不画, 单出一张。
save_tanglegram_legend <- function(k = k_paf, scale = 0.5,
                                   fname_base = "tanglegram_legend",
                                   width = 2.2, height = 2.2) {
  legend_labels <- c(paste0("C", 1:k), paste0("SC", 1:k))
  legend_colors <- c(col_pafclu[1:k], col_subtypeclu[1:k])
  
  for (r in c(600, 300)) {
    f <- if (r == 600) paste0(fname_base, ".tif")
    else          paste0(fname_base, "_300dpi.tif")
    ragg::agg_tiff(f, width = width, height = height, units = "in",
                   res = r, compression = "lzw", background = "white",
                   scaling = 1)
    par(mar = c(0, 0, 0, 0))
    plot.new()
    legend(x = "center",
           legend = legend_labels,
           col = legend_colors,
           text.col = legend_colors,
           lwd = 3 * scale, lty = 1,
           seg.len = 2.0 * scale,
           bty = "n",
           cex = 1.25 * scale,
           y.intersp = 1.7,
           text.font = 2)
    dev.off()
    cat("Saved:", f, "\n")
  }
}
save_tanglegram_legend()

# ============================================================================
# 8. 树状图 + 柱状图一体化展示（推荐主图）
#    左: 横向树状图   右: 横向堆叠柱状图   共享y轴(国家)
# ============================================================================
plot_dendro_with_bars <- function(hc, prop_df, paf_assign,
                                  k_cut = k_paf,
                                  title_text = "",
                                  bar_colors = col_subtype) {
  
  # ---- 1. dendrogram 数据 ----
  dend_data <- ggdendro::dendro_data(hc, type = "rectangle")
  cluster_assign <- cutree(hc, k = k_cut)
  ordered_countries <- hc$labels[hc$order]
  
  cluster_along_order <- cluster_assign[ordered_countries]
  cluster_relabel <- as.integer(factor(cluster_along_order,
                                       levels = unique(cluster_along_order)))
  names(cluster_relabel) <- ordered_countries
  
  country_to_sc <- paste0("SC", cluster_relabel[ordered_countries])
  names(country_to_sc) <- ordered_countries
  
  subtype_clu_colors <- col_subtypeclu[1:k_cut]
  
  # ---- 2. 树状图(无文字标签) ----
  # 小数号统一为中点 U+00B7
  dec_lab <- function(v) gsub("\\.", "\u00B7",
                              scales::label_number(accuracy = 0.1)(v))
  p_dend <- ggplot() +
    geom_segment(data = dend_data$segments,
                 aes(x = x, y = y, xend = xend, yend = yend),
                 linewidth = 0.6, color = "grey30") +
    coord_flip() +
    scale_x_reverse(expand = expansion(add = c(0.6, 0.6))) +
    scale_y_reverse(expand = expansion(mult = c(0.05, 0.02)),
                    labels = dec_lab) +
    labs(x = NULL, y = "Distance") +
    theme_minimal(base_size = 8) +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor = element_blank(),
          plot.margin = margin(3, 0, 3, 3))
  
  # ---- 3. 柱状图(国家名作为 y 轴标签) ----
  pdat <- prop_df %>%
    mutate(country = factor(country, levels = rev(ordered_countries)),
           subtype = factor(subtype,
                            levels = c("Subarachnoid hemorrhage",
                                       "Intracerebral hemorrhage",
                                       "Ischemic stroke"),
                            labels = c("SAH", "ICH", "IS")),
           paf_cluster = paf_assign[as.character(country)])
  
  axis_text_colors <- subtype_clu_colors[country_to_sc[rev(ordered_countries)]]
  
  p_bar <- ggplot(pdat, aes(x = country, y = prop, fill = subtype)) +
    geom_col(width = 0.78, color = "white", linewidth = 0.3) +
    geom_text(aes(label = ifelse(prop >= 0.03,
                                 percent(prop, accuracy = 1), ""),
                  colour = subtype),
              position = position_stack(vjust = 0.5),
              size = 2.1, fontface = "bold") +
    scale_colour_manual(values = c("IS" = "white", "ICH" = "white",
                                   "SAH" = "grey20"), guide = "none") +
    scale_y_continuous(labels = percent_format(),
                       expand = expansion(mult = c(0, 0.02))) +
    scale_fill_manual(values = bar_colors,
                      breaks = c("IS", "ICH", "SAH"), name = "Subtype") +
    coord_flip() +
    labs(x = NULL, y = "Proportion") +
    theme_minimal(base_size = 8) +
    theme(axis.text.y  = element_text(face = "bold", size = 7.5,
                                      color = axis_text_colors, hjust = 0.5),
          axis.ticks.y = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor = element_blank(),
          plot.margin = margin(3, 3, 3, 3),
          legend.position = "none")
  
  # ---- 4. PAF 注释带 ----
  paf_anno_df <- tibble(country = factor(rev(ordered_countries),
                                         levels = rev(ordered_countries)),
                        paf = paf_assign[rev(ordered_countries)])
  
  p_paf_anno <- ggplot(paf_anno_df, aes(x = 1, y = country, fill = paf)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = paf), size = 2.1, fontface = "bold", color = "white") +
    scale_fill_manual(values = col_pafclu, name = "PAF cluster") +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    labs(x = NULL, y = NULL) +
    theme_void() +
    theme(legend.position = "none",
          axis.text = element_blank())
  
  # ---- 5. 三个图例：统一主题，统一用填充色块 ----
  leg_theme <- theme_minimal(base_size = 8) +
    theme(legend.position   = "bottom",
          legend.direction  = "horizontal",
          legend.title      = element_text(size = 7.5, face = "bold"),
          legend.text       = element_text(size = 7),
          legend.key.size   = unit(0.4, "cm"),
          legend.spacing.x  = unit(0.15, "cm"),
          legend.margin     = margin(0, 4, 0, 4))
  
  # 图例 1: Subtype
  p_leg1 <- ggplot(pdat, aes(x = country, y = prop, fill = subtype)) +
    geom_col() +
    scale_fill_manual(values = bar_colors,
                      breaks = c("IS", "ICH", "SAH"), name = "Subtype") +
    leg_theme
  
  # 图例 2: PAF cluster
  p_leg2 <- ggplot(paf_anno_df, aes(x = 1, y = country, fill = paf)) +
    geom_tile() +
    scale_fill_manual(values = col_pafclu, name = "PAF cluster") +
    leg_theme
  
  # 图例 3: Subtype cluster（改为填充色块，与前两个等大）
  sc_legend_df <- tibble(
    cluster = factor(paste0("SC", 1:k_cut),
                     levels = paste0("SC", 1:k_cut))
  )
  p_leg3 <- ggplot(sc_legend_df, aes(x = 1, y = 1, fill = cluster)) +
    geom_tile() +
    scale_fill_manual(values = subtype_clu_colors, name = "Subtype cluster") +
    leg_theme
  
  legend1 <- cowplot::get_legend(p_leg1)
  legend2 <- cowplot::get_legend(p_leg2)
  legend3 <- cowplot::get_legend(p_leg3)
  
  # ---- 6. 三等分排列，各自居中；不再手调 rel_widths 以免截断 ----
  # 顺序与主图左中右对应：国名颜色(SC) → 柱状(Subtype) → 右侧色带(PAF)
  legend_row <- cowplot::plot_grid(legend3, legend1, legend2, nrow = 1)
  
  
  # ---- 7. 主图三栏拼接(无图例) ----
  main_row <- p_dend + p_bar + p_paf_anno +
    plot_layout(widths = c(0.9, 1.9, 0.13))
  
  # ---- 8. 主图 + 图例条 上下拼接 ----
  cowplot::plot_grid(
    main_row, legend_row,
    ncol = 1,
    rel_heights = c(1, 0.12)
  )
}


# DALY 一体化版
# 都强制 k=2，与 PAF 聚类对齐（且 Incidence k=2/k=3 silhouette 差距<0.01）
p_combo_daly <- plot_dendro_with_bars(clust_daly$hc, subtype_daly,
                                      paf_cluster_assign,
                                      k_cut = k_paf,
                                      title_text = "ASDR 2023")
# 正文 Figure 3: 单幅整合图 (左树状图 / 中堆积条 / 右色条), 无面板字母,
# R 输出即投稿版
save_panel(p_combo_daly, "Figure3", w_mm = 180, h_mm = 110)


# Incidence 一体化版
p_combo_inc <- plot_dendro_with_bars(clust_inc$hc, subtype_inc,
                                     paf_cluster_assign,
                                     k_cut = k_paf,
                                     title_text = "ASIR 2023")
# Figure S14: 同结构的 ASIR 版本, 同样是单幅, R 输出即投稿版
save_panel(p_combo_inc, "FigureS14", w_mm = 180, h_mm = 110)

# # ============================================================================
# # 9. ILR 空间二维散点图（可选辅助）
# # ============================================================================
# make_2d_scatter <- function(clust_obj, paf_assign, title_text) {
#   df <- as.data.frame(clust_obj$ilr) %>%
#     rownames_to_column("country") %>%
#     mutate(paf      = paf_assign[country],
#            subtype_cl = factor(clust_obj$cluster_k_paf[country]))
#   
#   ggplot(df, aes(x = V1, y = V2)) +
#     geom_point(aes(color = paf, shape = subtype_cl), size = 5, alpha = 0.85) +
#     geom_text(aes(label = country), vjust = -1.3, size = 3.3) +
#     scale_color_manual(values = col_pafclu, name = "PAF cluster") +
#     scale_shape_manual(values = c("1" = 16, "2" = 17, "3" = 15),
#                        name = "Subtype cluster (k=k_paf)") +
#     labs(title = title_text,
#          x = "ILR axis 1 (IS vs ICH/SAH balance)",
#          y = "ILR axis 2 (ICH vs SAH balance)") +
#     theme_minimal(base_size = 11) +
#     theme(plot.title = element_text(face = "bold"),
#           legend.position = "bottom")
# }
# 
# p_scatter_daly <- make_2d_scatter(clust_daly, paf_cluster_assign,
#                                   "ILR space: DALY-based")
# p_scatter_inc  <- make_2d_scatter(clust_inc, paf_cluster_assign,
#                                   "ILR space: Incidence-based")
# 
# p_scatter <- (p_scatter_daly + p_scatter_inc) +
#   plot_layout(guides = "collect") & theme(legend.position = "bottom")
# ggsave("subtype_ILR_scatter.pdf", p_scatter, width = 12, height = 6)
# 
# # ============================================================================
# # 10. 最终结果输出
# # ============================================================================
# cat("\n\n=========== Final summary ===========\n")
# print(summary_tbl)
# cat("\nFiles saved:\n",
#     "  - subtype_stacked_DALY_Incidence_2023.[pdf|png]\n",
#     "  - tanglegram_DALY.pdf, tanglegram_Incidence.pdf\n",
#     "  - dendro_with_bars_DALY.[pdf|png]\n",
#     "  - dendro_with_bars_Incidence.[pdf|png]\n",
#     "  - dendro_with_bars_BOTH.[pdf|png]   ← 推荐主图\n",
#     "  - subtype_ILR_scatter.pdf\n",
#     "  - consistency_test_results.csv\n")