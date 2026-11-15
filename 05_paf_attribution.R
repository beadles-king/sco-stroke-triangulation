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
# 改进版：2023年代谢归因 热图 + 层次聚类 (Ward.D2)
# 主图：Stroke (全中风) — 双面板：PAF% + 归因ASDR，共享聚类结果
# 附图：IS / ICH / SAH 亚型
#
# 聚类策略：
#   - 基于归因占比(PAF%) 做层次聚类 (主分析, 识别"风险结构模式")
#   - 归因ASDR 热图共享同样的列顺序, 便于对照
#   - 敏感性分析: 基于 PAF × log(ASDR) 复合指标聚类
#
# 数据：PAF.csv (已有占比), 归因ASDR.csv (待补, 同样格式)
###############################################################################
setwd(RAWDAT)

library(tidyverse)
library(ComplexHeatmap)   # 行/列dendrogram + annotation
library(circlize)         # colorRamp2 连续色带
library(cluster)          # silhouette
library(factoextra)       # fviz_nbclust
library(grid)
library(dendextend)
if (!requireNamespace("ragg", quietly = TRUE)) install.packages("ragg")

# ---- 全文格式约定 ----
DEC <- "\u00B7"
dec <- function(x) gsub("\\.", DEC, as.character(x))
# 半值向上舍入。先乘到目标位数，加 0·5 再取整；1e-9 用来抵消
# 二进制表示对半值的微小下偏（如 0·00045 实际存为 0·00044999999999999999）。
# p 落在 1/100000 网格上，最近的非半值也在 0·1 之外，这个容差是安全的。
round_half_up <- function(x, digits) {
  z <- x * 10^digits
  floor(z + 0.5 + 1e-9) / 10^digits
}

fmt_p <- function(p) {
  if (is.na(p)) return(NA_character_)
  if (p < 0.0001) return(paste0("<0", DEC, "0001"))
  if (p < 0.01)   return(dec(sprintf("%.4f", round_half_up(p, 4))))
  dec(sprintf("%.2f", round_half_up(p, 2)))
}

# 底色偏暗时用白字，保证格内数字在深蓝/深红上仍可读
txt_col <- function(fill) {
  rgb <- col2rgb(fill)[, 1] / 255
  lum <- 0.299 * rgb[1] + 0.587 * rgb[2] + 0.114 * rgb[3]
  if (lum < 0.55) "white" else "black"
}

# ============================================================================
# 0. 数据加载
# ============================================================================
PAF  <- read_csv('PAF_v2.csv')  %>%
  filter(metric_name == 'Percent') %>%
  select(location_name, cause_name, rei_name, year, val) %>%
  rename(paf = val)

# 预期归因 ASDR 的文件（请确认文件名与列名；如不同请修改）
ASDR <- read_csv('PAF_v2.csv') %>%
  filter(metric_name == 'Rate') %>%
  select(location_name, cause_name, rei_name, year, val) %>%
  rename(asdr = val)

# 合并
dat <- PAF %>%
  inner_join(ASDR, by = c("location_name", "cause_name", "rei_name", "year"))

# 国家名统一（与01_gbd_burden_trends.R一致）
dat <- dat %>%
  mutate(location_name = case_when(
    location_name == "Russian Federation" ~ "Russia",
    location_name == "Iran (Islamic Republic of)" ~ "Iran",
    TRUE ~ location_name
  ))

# RF简写（热图更紧凑）
rf_short <- c(
  "High systolic blood pressure" = "SBP",
  "High LDL cholesterol"         = "LDL",
  "Kidney dysfunction"           = "KD",
  "High fasting plasma glucose"  = "FG",
  "High body-mass index"         = "BMI"
)
dat$rei_short <- rf_short[dat$rei_name]

# RF 固定顺序（按生物学意义, 血管-代谢-肾-糖-体质）
rf_order <- c("SBP", "LDL", "KD", "FG", "BMI")

# ============================================================================
# 1. 通用函数：准备宽矩阵
# ============================================================================
build_matrix <- function(data, cause, yr, value_col, rfs) {
  m <- data %>%
    filter(cause_name == cause, year == yr, rei_short %in% rfs) %>%
    select(location_name, rei_short, all_of(value_col)) %>%
    pivot_wider(names_from = location_name, values_from = all_of(value_col)) %>%
    column_to_rownames("rei_short")
  m <- m[rfs, ]  # 固定行序
  as.matrix(m)
}

# ============================================================================
# 2. 基于 PAF 的层次聚类 (列聚类：10国)
#    对每个 RF 做 z-score，让每个 RF 平权
# ============================================================================
cluster_countries_by_paf <- function(paf_mat, k_range = 2:5) {
  # 按行 (RF) 标准化，国家作为样本
  paf_z <- t(scale(t(paf_mat)))   # 行标准化
  d <- dist(t(paf_z), method = "euclidean")   # 国家间距离
  hc <- hclust(d, method = "ward.D2")
  
  # 轮廓系数选最优 k
  sil <- sapply(k_range, function(ki) {
    cl <- cutree(hc, k = ki)
    mean(silhouette(cl, d)[, "sil_width"])
  })
  names(sil) <- paste0("k=", k_range)
  cat("轮廓系数:\n"); print(round(sil, 3))
  best_k <- k_range[which.max(sil)]
  cat("最优 k =", best_k, "\n")
  
  clusters <- cutree(hc, k = best_k)
  list(hc = hc, dist = d, k = best_k, clusters = clusters,
       sil_scores = sil, mat_z = paf_z)
}

# ============================================================================
# 3. 主图绘制函数 (ComplexHeatmap)
# ============================================================================
draw_dual_heatmap <- function(paf_mat, asdr_mat, hc_col, clusters, k,
                              title_prefix = "Stroke", save_path = NULL,
                              # ---- 版面参数; 默认值 = 旧行为, 保证 Figure 2A 不变 ----
                              body_w_cm = 11,      # 热图体宽度 (cm)
                              cell_mm   = NULL,    # 每行格高 (mm); NULL = 旧的定高 4.5 cm
                              dend_mm   = 15,      # 列树状图高度 (mm)
                              # 多面板拼图时 Cluster 与 PAF 图例各张完全相同,
                              # 只需在第一张保留; ASDR 色标各自归一, 每张都要留
                              show_cluster_legend = TRUE,
                              show_paf_legend     = TRUE) {
  paf_pct <- 100 * paf_mat   # 百分数尺度，仅用于绘图；聚类仍用原比例矩阵
  
  # === 标签统一 (中国所在cluster = C1) ===
  if (k == 2 && "China" %in% names(clusters)) {
    if (clusters["China"] != 1) {
      clusters <- ifelse(clusters == 1, 2, 1)
      names(clusters) <- names(clusters)
    }
  }
  
  # === 新增: dendrogram 左右翻转, 让 C1 国家在左侧 ===
  if (k == 2) {
    c1_countries <- names(clusters)[clusters == 1]
    current_order <- hc_col$labels[hc_col$order]
    # 检查 C1 国家在当前 dendrogram 顺序中是不是在左半边
    c1_positions <- which(current_order %in% c1_countries)
    # 如果 C1 国家的平均位置 > n/2 (即在右半边), 翻转
    if (mean(c1_positions) > length(current_order) / 2) {
      desired_order <- c(c1_countries, setdiff(current_order, c1_countries))
      hc_col <- as.hclust(rotate(as.dendrogram(hc_col), order = desired_order))
    }
  } 
  col_order <- hc_col$labels[hc_col$order]
  
  # ---- 热图体尺寸: 定高改为随危险因素个数按格高伸缩 ----
  n_rf   <- nrow(paf_mat)
  body_w <- unit(body_w_cm, "cm")
  body_h <- if (is.null(cell_mm)) unit(4.5, "cm") else unit(n_rf * cell_mm, "mm")
  
  # 固定开销 (15.5 + dend_mm) mm 为实测值: 聚类条 + 两热图间隙 + 45 度国名 +
  # 内边距 + 列树状图。用于预估裁切后的内容高度, 便于核对 240 mm 预算
  if (!is.null(cell_mm)) {
    cat(sprintf("  [版面] RF = %d, 格高 %.1f mm -> 内容高约 %.1f mm\n",
                n_rf, cell_mm, 2 * n_rf * cell_mm + 15.5 + dend_mm))
  }
  
  # 聚类分组注释
  cluster_colors <- setNames(
    c("#E6550D", "#7FBFA0", "#9970AB", "#F4A261", "#8AB17D")[1:k],
    paste0("C", 1:k)
  )
  col_anno <- HeatmapAnnotation(
    Cluster = paste0("C", clusters[colnames(paf_mat)]),  # 关键：用 colnames 而非 col_order
    col = list(Cluster = cluster_colors),
    show_annotation_name = FALSE,
    show_legend = show_cluster_legend,                 # 可关
    annotation_legend_param = list(title = "Cluster")
  )
  
  # === Panel A: PAF ===
  # rank 注释（每国内RF排序）
  rank_mat_paf <- apply(paf_mat, 2, function(x) rank(-x, ties.method = "min"))
  
  # h1 (PAF面板)
  h1 <- Heatmap(
    paf_pct,
    name = "PAF",
    col = colorRamp2(c(0, 35, 70), c("white", "#FDBF6F", "#B2182B")),
    cluster_columns = hc_col,
    cluster_rows = FALSE,
    column_dend_height = unit(dend_mm, "mm"),          # 参数化
    top_annotation = col_anno,
    # column_title 删除
    row_names_side = "left",
    column_names_rot = 45,
    cell_fun = function(j, i, x, y, w, h, fill) {
      grid.text(sprintf("%d\n%s", rank_mat_paf[i, j],
                        dec(sprintf("%.1f", paf_pct[i, j]))),
                x, y, gp = gpar(fontsize = 7, col = txt_col(fill)))
    },
    row_names_gp = gpar(fontsize = 8),                   # 新增
    column_names_gp = gpar(fontsize = 8),                # 新增
    width = body_w, height = body_h,                      # 参数化
    show_heatmap_legend = show_paf_legend,                # 可关
    heatmap_legend_param = list(
      title = "PAF, %",
      title_gp = gpar(fontsize = 8, fontface = "bold"),
      labels_gp = gpar(fontsize = 7),
      at = c(0, 35, 70),
      labels = c("0", "35", "70"),
      direction = "vertical"
    )
  )
  
  # === Panel B: Attributable ASDR ===
  rank_mat_asdr <- apply(asdr_mat, 2, function(x) rank(-x, ties.method = "min"))
  asdr_max <- max(asdr_mat, na.rm = TRUE)
  
  # h2 (ASDR面板)  
  h2 <- Heatmap(
    asdr_mat,
    name = "ASDR",
    col = colorRamp2(c(0, asdr_max/2, asdr_max),
                     c("white", "#A6BDDB", "#08306B")),
    cluster_columns = hc_col,
    cluster_rows = FALSE,
    show_column_dend = FALSE,
    # column_title 删除
    row_names_side = "left",
    column_names_rot = 45,
    cell_fun = function(j, i, x, y, w, h, fill) {
      grid.text(sprintf("%d\n%s", rank_mat_asdr[i, j],
                        dec(sprintf("%.0f", asdr_mat[i, j]))),
                x, y, gp = gpar(fontsize = 7, col = txt_col(fill)))
    },
    row_names_gp = gpar(fontsize = 8),                   # 新增
    column_names_gp = gpar(fontsize = 8),                # 新增
    width = body_w, height = body_h,                      # 参数化
    heatmap_legend_param = list(
      title = "ASDR\nper 100 000",             # 改这里
      title_gp = gpar(fontsize = 8, fontface = "bold"),
      labels_gp = gpar(fontsize = 7),
      at = pretty(c(0, asdr_max), 4),
      labels = dec(format(pretty(c(0, asdr_max), 4), trim = TRUE)),
      direction = "vertical"
    )
  )
  
  # 上下拼接
  ht_list <- h1 %v% h2
  
  if (!is.null(save_path)) {
    # 出图改走 ragg: cairo_pdf 会把每个字形步进四舍五入到整数 point (0.353 mm),
    # 再用 TJ 微调补偿, 同一单词内字距不匀 (最坏 1 pt)。ragg 在目标分辨率的像素
    # 网格上排字, 600 dpi 下取整误差 <= 0.042 mm。画布 180 x 190.5 mm 足以容纳
    # 各亚型的内容 (最大 122 x 91 mm), 多余留白由拼图阶段的包围盒裁切吃掉。
    stem <- sub("\\.(pdf|tif|tiff)$", "", save_path)
    for (r in c(600, 300)) {
      f <- if (r == 600) paste0(stem, ".tif") else paste0(stem, "_300dpi.tif")
      ragg::agg_tiff(f, width = 180, height = 190.5, units = "mm",
                     res = r, compression = "lzw", background = "white",
                     scaling = 1)
      draw(ht_list, merge_legend = FALSE,
           heatmap_legend_side = "right", annotation_legend_side = "right")
      dev.off()
      cat("已保存:", f, sprintf("(%d dpi)\n", r))
    }
  }
  ht_list
}

# ============================================================================
# 4. 敏感性分析：基于 PAF × log(ASDR+1) 的复合指标聚类
# ============================================================================
cluster_countries_composite <- function(paf_mat, asdr_mat, k_range = 2:5) {
  # 复合指标：归因占比 × log(归因ASDR)，反映结构×规模
  composite <- paf_mat * log1p(asdr_mat)
  comp_z <- t(scale(t(composite)))
  d <- dist(t(comp_z), method = "euclidean")
  hc <- hclust(d, method = "ward.D2")
  
  sil <- sapply(k_range, function(ki) {
    cl <- cutree(hc, k = ki)
    mean(silhouette(cl, d)[, "sil_width"])
  })
  best_k <- k_range[which.max(sil)]
  cat("\n[敏感性分析·复合指标聚类] 最优 k =", best_k, "\n")
  list(hc = hc, k = best_k, clusters = cutree(hc, k = best_k))
}

# ============================================================================
# 5. 一体化分析函数
# ============================================================================
run_cause <- function(data, cause_name_in, rfs, save_prefix, ...) {
  cat("\n\n==========", cause_name_in, "==========\n")
  
  # 宽矩阵
  paf_m  <- build_matrix(data, cause_name_in, 2023, "paf",  rfs)
  asdr_m <- build_matrix(data, cause_name_in, 2023, "asdr", rfs)
  
  cat("矩阵维度: RF =", nrow(paf_m), "| 国家 =", ncol(paf_m), "\n")
  
  # 主聚类 (基于PAF)
  res_paf <- cluster_countries_by_paf(paf_m)
  
  # 敏感性聚类 (复合指标)
  res_comp <- cluster_countries_composite(paf_m, asdr_m)
  
  # 一致性: 两种聚类的 Rand index
  if (length(res_paf$clusters) == length(res_comp$clusters)) {
    # 简化实现: 调整兰德指数
    ari <- if (requireNamespace("mclust", quietly = TRUE)) {
      mclust::adjustedRandIndex(res_paf$clusters, res_comp$clusters)
    } else {
      NA
    }
    cat("PAF聚类 vs 复合指标聚类 ARI =", round(ari, 3), "\n")
  }
  
  # 打印分组
  cat("\n[主聚类·基于PAF] k =", res_paf$k, "\n")
  for (ci in sort(unique(res_paf$clusters))) {
    cat("  C", ci, ":", paste(names(res_paf$clusters)[res_paf$clusters == ci],
                              collapse = ", "), "\n")
  }
  
  # 画图
  ht <- draw_dual_heatmap(
    paf_m, asdr_m, res_paf$hc, res_paf$clusters, res_paf$k,
    title_prefix = cause_name_in,
    save_path = paste0(save_prefix, "_heatmap_dual.pdf"),
    ...                                   # 透传版面参数
  )
  
  # 导出聚类表
  cl_df <- data.frame(
    country = names(res_paf$clusters),
    cluster_paf = res_paf$clusters,
    cluster_composite = res_comp$clusters[names(res_paf$clusters)]
  )
  write.csv(cl_df, paste0(save_prefix, "_clusters.csv"), row.names = FALSE)
  
  invisible(list(paf_mat = paf_m, asdr_mat = asdr_m,
                 res_paf = res_paf, res_comp = res_comp))
}

# ============================================================================
# 6. 运行：主分析 + 三亚型
# ============================================================================
setwd(file.path(RAWDAT, "PAF"))

# Figure S9 的三张亚型面板要竖排进 180 x 240 mm。
#   宽度: 热图体 10 cm -> 每列 10 mm, 色块 10 x 6 mm, 格内左右各留约 2.2 mm,
#         与 Figure 2A 的疏朗度接近; 面板内容宽约 122 mm (介于单双栏之间)
#   高度: 格高 6 mm (两行 7 pt 文字约 5.9 mm)。但 ICH 与 SAH 的热图体较矮,
#         其内容高度被右侧图例栈 (Cluster + PAF + ASDR, 约 78 mm) 撑住,
#         三张合计 247 mm 会超出 240 mm 预算。
#         Cluster 色与 PAF 色标 (0/35/70) 三张完全相同, 只在 (A) 保留一次;
#         ASDR 色标各自归一到本亚型最大值, 三张都必须保留。
#         关掉重复图例后: IS 90.5 / ICH 73.5 / SAH 49.5 mm, 合计 213.5 mm
S9_LAYOUT <- list(body_w_cm = 10, cell_mm = 6, dend_mm = 10)

# 主图: 全中风 (5个RF) —— 正文 Figure 2A, 不传版面参数, 保持原尺寸
main_res <- run_cause(dat, "Stroke",
                      rfs = rf_order,
                      save_prefix = "Fig3_Stroke")

# 附图: 三亚型 —— 正文 Figure S9
# 注: 原脚本在本节和末尾"空间检验准备"处各调用了一次 run_cause, 后一次会用
#     同名 save_prefix 覆盖前一次的面板文件。此处合并为一次调用并保留返回值,
#     供后续空间检验使用, 避免面板被写回旧尺寸。
is_res  <- do.call(run_cause, c(list(dat, "Ischemic stroke",
                                     rfs = rf_order, save_prefix = "SuppFig_IS"), S9_LAYOUT))

# (B) 与 (C) 关掉与 (A) 重复的 Cluster / PAF 图例
S9_NOLEG <- c(S9_LAYOUT, list(show_cluster_legend = FALSE,
                              show_paf_legend     = FALSE))

ich_res <- do.call(run_cause, c(list(dat, "Intracerebral hemorrhage",
                                     rfs = c("SBP", "KD", "FG", "BMI"),
                                     save_prefix = "SuppFig_ICH"), S9_NOLEG))

sah_res <- do.call(run_cause, c(list(dat, "Subarachnoid hemorrhage",
                                     rfs = c("SBP", "BMI"),
                                     save_prefix = "SuppFig_SAH"), S9_NOLEG))

cat("\n========== 全部完成 ==========\n")


#========== PAF排名稳健性检验 ==========
CAUSE   <- "Stroke"          # 与第 306 行调用处传入的字符串保持一致
paf_mat <- build_matrix(dat, CAUSE, 2023, "paf", rf_order) * 100
rank_mat_paf <- apply(paf_mat, 2, function(x) rank(-x, ties.method = "min"))
rank_freq <- as.data.frame(t(apply(rank_mat_paf, 1, function(r)
  table(factor(r, levels = 1:5)))))
names(rank_freq) <- paste0("rank", 1:5)
rank_freq$median_paf <- round(apply(paf_mat, 1, median), 1)
rank_freq$modal_rank <- apply(rank_mat_paf, 1, function(r)
  as.integer(names(which.max(table(r)))))
rank_freq$paf_range  <- sprintf("%.1f-%.1f", apply(paf_mat,1,min), apply(paf_mat,1,max))
rank_freq <- rank_freq[rf_order, ]
write.csv(rank_freq, file.path(RAWDAT, "PAF/B5_population_rank_frequency.csv"))


#========== k与GBD超区域交互 ==========
# 需跑完PAF_v2_空间检验_C2敏感性.R
hc  <- hclust(res_stroke$base$feat_dist, method = "ward.D2")
for (ki in 2:5) print(split(sco_countries, cutree(hc, ki)[sco_countries]))

