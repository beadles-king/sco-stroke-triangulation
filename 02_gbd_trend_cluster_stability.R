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
# Bootstrap 聚类稳定性补充模块 (基于 01_gbd_burden_trends.R 的输出对象)
#
# 用法: 在 01_gbd_burden_trends.R 跑完后, 直接 source 本文件:
#       source('02_gbd_trend_cluster_stability.R')
#
# 依赖对象: res_stroke (主分析), res_list (3 亚型: IS / ICH / SAH)
#
# 输出 (写到当前工作目录):
#   stability_jaccard_<name>.csv / .pdf      # cluster 级 Jaccard
#   stability_coassign_<name>.csv / .pdf     # 国家级 co-assignment
#   其中 <name> ∈ {stroke, IS, ICH, SAH}
#
# 本版改动仅限绘图与出图部分（分析逻辑、随机种子、B 次数一律未动）：
#   - 小数号一律中点 ·（格内数字、y 轴刻度标签）
#   - 字号：基础 8 pt，图内最小 7 pt；幅面 ≤180 mm，无需期刊缩放
#   - 字体 Helvetica，改用 cairo_pdf 嵌字体
#     （原 pdf() 不嵌字体且用 8 位编码，已把 "Co-cluster"/"co-assignment"
#       里的连字符静默换成了减号 U+2212）
#   - 去掉图内标题（面板字母与图注在拼图阶段添加）
#   - 配色改为全文主色 #2C5985 / #A63603 一族；三条阈值线改统一灰色虚线，
#     不再用红/橙/绿三色——阈值已由 y 轴刻度标出，三色会与方法配色抢编码
#   - 热图格内数字按底色亮度自动切黑白字
#   - 方法顺序改为 Bootstrap → Jitter → Subset，与附录行文一致
###############################################################################

# ============================================================================
# 0. 依赖
# ============================================================================
if (!requireNamespace("fpc", quietly = TRUE))      install.packages("fpc")
if (!requireNamespace("reshape2", quietly = TRUE)) install.packages("reshape2")

library(fpc)
library(reshape2)
library(scales)
library(grid)
# tidyverse / ggplot2 / cluster 已由原脚本载入

# ============================================================================
# 配置
# ============================================================================
B_BOOT  <- 1000          # bootstrap 重复次数
SEED    <- 20260509      # 可复现
K_RANGE <- 2:4           # 同主脚本

# ---- 全文格式约定（同交接文档第 2 节）----
DEC       <- "\u00B7"    # 中点
BASE_SIZE <- 8
FONT      <- "Helvetica"
COL_MAIN  <- "#2C5985"   # 主蓝
COL_HIGH  <- "#A63603"   # 橙红
COL_NEUT  <- "#7F7F7F"   # 中性灰

# 小数点 -> 中点
dec <- function(x) gsub("\\.", DEC, as.character(x))

# 两位小数 + 中点，供坐标轴与格内数字统一调用
fmt2 <- function(x) dec(sprintf("%.2f", x))

# 深色底自动转白字
txt_col <- function(fill) {
  rgb <- col2rgb(fill)[, 1] / 255
  lum <- 0.299 * rgb[1] + 0.587 * rgb[2] + 0.114 * rgb[3]
  if (lum < 0.55) "white" else "black"
}

# 统一主题
theme_pub <- function(base = BASE_SIZE) {
  theme_minimal(base_size = base, base_family = FONT) +
    theme(
      text             = element_text(size = base, family = FONT),
      axis.title       = element_text(size = base),
      axis.text        = element_text(size = base - 1, colour = "black"),
      legend.text      = element_text(size = base - 1),
      legend.title     = element_text(size = base - 1),
      strip.text       = element_text(size = base, face = "bold"),
      panel.grid.minor = element_blank(),
      plot.title       = element_blank(),      # 面板标题在拼图阶段处理
      plot.margin      = margin(2, 2, 2, 2, "mm")
    )
}

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
# 1. 函数: 用 fpc::clusterboot 运行 cluster 级 Jaccard 稳定性  [未改动]
# ============================================================================
run_stability <- function(feat_scaled, k_range = K_RANGE,
                          B = B_BOOT, seed = SEED) {
  out <- list()
  
  for (ki in k_range) {
    cat("    k =", ki, "...\n")
    
    # 用三种 bootstrap 方法对同一聚类设置评估
    cb <- clusterboot(
      feat_scaled,
      B             = B,
      bootmethod    = c("boot", "jitter", "subset"),
      clustermethod = hclustCBI,
      method        = "ward.D2",          # 与主脚本一致
      k             = ki,
      seed          = seed,
      jittertuning  = 0.05,                # SD 的 5% 噪声
      subtuning     = ceiling(0.8 * nrow(feat_scaled)),  # 抽 80% 国家
      count         = FALSE,
      showplots     = FALSE
    )
    
    df <- data.frame(
      k                = ki,
      cluster          = seq_len(ki),
      n_members        = as.integer(table(cb$partition)),
      members          = sapply(seq_len(ki), function(c) {
        paste(rownames(feat_scaled)[cb$partition == c], collapse = ", ")
      }),
      jaccard_boot     = round(cb$bootmean,    3),
      jaccard_jitter   = round(cb$jittermean,  3),
      jaccard_subset   = round(cb$subsetmean,  3),
      stringsAsFactors = FALSE
    )
    out[[paste0("k", ki)]] <- df
  }
  
  do.call(rbind, out)
}

# ============================================================================
# 2. 函数: 国家级 co-assignment matrix (用 jitter 方法)  [未改动]
# ============================================================================
run_coassignment <- function(feat_scaled, k, B = B_BOOT,
                             jitter_sd = 0.05, seed = SEED) {
  set.seed(seed)
  n         <- nrow(feat_scaled)
  countries <- rownames(feat_scaled)
  M <- matrix(0, n, n, dimnames = list(countries, countries))
  
  for (b in seq_len(B)) {
    noise <- matrix(rnorm(length(feat_scaled), 0, jitter_sd), nrow = n)
    X     <- feat_scaled + noise
    
    cl <- cutree(hclust(dist(X), method = "ward.D2"), k = k)
    
    # 累计同簇出现次数 (向量化)
    same <- outer(cl, cl, "==")
    M    <- M + same
  }
  
  M / B
}

# ============================================================================
# 3. 函数: Jaccard 条形图  → 面板 B
# ============================================================================
plot_jaccard <- function(stab_df) {
  long <- stab_df %>%
    pivot_longer(
      cols      = c("jaccard_boot", "jaccard_jitter", "jaccard_subset"),
      names_to  = "method",
      values_to = "jaccard"
    ) %>%
    mutate(method = recode(method,
                           jaccard_boot   = "Bootstrap",
                           jaccard_jitter = "Jitter",
                           jaccard_subset = "Subset"),
           method = factor(method, levels = c("Bootstrap", "Jitter", "Subset")))
  
  ggplot(long, aes(x = factor(cluster), y = jaccard, fill = method)) +
    geom_col(position = position_dodge(0.8), width = 0.72) +
    # Hennig (2007) 阈值: 0·60 / 0·75 / 0·85，统一灰虚线，数值由 y 轴刻度给出
    geom_hline(yintercept = c(0.60, 0.75, 0.85), linetype = "dashed",
               colour = COL_NEUT, linewidth = 0.3) +
    facet_wrap(~ paste0("k = ", k), scales = "free_x", nrow = 1) +
    scale_fill_manual(values = c("Bootstrap" = COL_MAIN,
                                 "Jitter"    = COL_HIGH,
                                 "Subset"    = COL_NEUT)) +
    scale_y_continuous(limits = c(0, 1.02),
                       breaks = c(0, 0.50, 0.60, 0.75, 0.85, 1.00),
                       labels = fmt2,
                       expand = expansion(mult = c(0, 0))) +
    labs(x = "Cluster", y = "Jaccard similarity", fill = NULL) +
    theme_pub() +
    theme(legend.position      = "top",
          legend.justification = "left",
          legend.key.size      = unit(3, "mm"),
          legend.margin        = margin(0, 0, 1, 0, "mm"),
          panel.grid.major.x   = element_blank())
}

# ============================================================================
# 4. 函数: Co-assignment heatmap (按 cluster 排序)  → 面板 A
# ============================================================================
plot_coassign <- function(M, clusters) {
  ord    <- order(clusters, names(clusters))
  M      <- M[ord, ord]
  cl_ord <- clusters[ord]
  
  long <- melt(M)
  colnames(long) <- c("c1", "c2", "freq")
  long$c1 <- factor(long$c1, levels = rownames(M))
  long$c2 <- factor(long$c2, levels = rownames(M))
  
  # 格内数字按底色亮度切黑白，阈值与全文一致 (0·55)
  ramp <- colorRampPalette(c("#FFFFFF", COL_MAIN))
  long$lab_col <- vapply(long$freq,
                         function(f) txt_col(ramp(101)[round(f * 100) + 1]),
                         character(1))
  
  # cluster 边界标记
  bdy <- cumsum(table(cl_ord))
  bdy <- bdy[-length(bdy)] + 0.5
  
  ggplot(long, aes(c1, c2, fill = freq)) +
    geom_tile(colour = "white", linewidth = 0.3) +
    geom_text(aes(label = fmt2(freq)), colour = long$lab_col,
              size = 7 / .pt) +                       # 7 pt
    geom_hline(yintercept = bdy, colour = "black", linewidth = 0.4) +
    geom_vline(xintercept = bdy, colour = "black", linewidth = 0.4) +
    scale_fill_gradient(low = "#FFFFFF", high = COL_MAIN, limits = c(0, 1),
                        name   = "Co-cluster\nfrequency",
                        breaks = c(0, 0.25, 0.50, 0.75, 1.00),
                        labels = fmt2) +
    coord_equal() +
    labs(x = NULL, y = NULL) +
    theme_pub() +
    theme(axis.text.x     = element_text(angle = 45, hjust = 1),
          panel.grid      = element_blank(),
          legend.position = "right",
          legend.key.width  = unit(3, "mm"),
          legend.key.height = unit(6, "mm"))
}

# ============================================================================
# 5. 主流程: 对一个 res 对象跑完全套
# ============================================================================
setwd(file.path(RAWDAT, "总体、年龄、性别负担/总体/总体聚类稳定性"))

run_full_stability <- function(res, name, label,
                               k_use = NULL,        # 指定要用的 k (覆盖 res$k)
                               B = B_BOOT) {
  # 如未指定,沿用 res$k (主脚本自动选的)
  if (is.null(k_use)) k_use <- res$k
  clusters_use <- cutree(res$hc, k = k_use)
  
  cat("\n========== Bootstrap stability:", label,
      "(reporting k =", k_use, ") ==========\n")
  
  # ----- (a) cluster 级 Jaccard (k = 2:4 全部跑,作为 k 选择稳健性证据) -----
  cat("  Step 1/2: cluster-level Jaccard (k =",
      paste(K_RANGE, collapse = ", "), ", B =", B, ")...\n")
  stab <- run_stability(res$feat_scaled, K_RANGE, B, SEED)
  
  csv_jac <- paste0("stability_jaccard_", name, ".csv")
  write.csv(stab, csv_jac, row.names = FALSE)
  cat("  saved:", csv_jac, "\n")
  print(stab[, c("k","cluster","n_members",
                 "jaccard_jitter","jaccard_subset","jaccard_boot")])
  
  # 面板 B：宽 180 mm（双栏满幅），高 62 mm
  save_panel(plot_jaccard(stab),
             paste0("stability_jaccard_", name),
             w_mm = 176, h_mm = 62)
  
  # ----- (b) 国家级 co-assignment (用 k_use,即论文报告的 k) -----
  cat("  Step 2/2: country-level co-assignment matrix (k =", k_use, ")...\n")
  M <- run_coassignment(res$feat_scaled, k = k_use, B = B, seed = SEED)
  
  csv_coa <- paste0("stability_coassign_", name, ".csv")
  write.csv(round(M, 3), csv_coa)
  
  # 面板 A：coord_equal 会自留边，多余留白由拼图阶段的包围盒裁切吃掉
  save_panel(plot_coassign(M, clusters_use),
             paste0("stability_coassign_", name),
             w_mm = 152, h_mm = 126)
  
  invisible(list(jaccard = stab, coassign = M, k_used = k_use))
}

# ============================================================================
# 6. 执行: 主分析 + 3 亚型 (使用论文报告的 k: 3, 4, 3, 4)
# ============================================================================
cat("\n##############################################################\n")
cat("# Bootstrap 聚类稳定性 (B =", B_BOOT, ")\n")
cat("# 总体负担轨迹聚类 (Stroke + 3 subtypes)\n")
cat("##############################################################\n")

# 论文报告的 k 值: Stroke=3, IS=4, ICH=3, SAH=4
K_REPORT <- list(stroke = 3, IS = 4, ICH = 3, SAH = 4)

stab_stroke <- run_full_stability(res_stroke, "stroke", "Stroke (overall)",
                                  k_use = K_REPORT$stroke)

stab_subtypes <- list()
subtype_label <- list(IS  = "Ischaemic stroke",
                      ICH = "Intracerebral haemorrhage",
                      SAH = "Subarachnoid haemorrhage")
for (sshort in names(res_list)) {
  stab_subtypes[[sshort]] <- run_full_stability(
    res_list[[sshort]], sshort, subtype_label[[sshort]],
    k_use = K_REPORT[[sshort]])
}

cat("\n========== Bootstrap stability 全部完成 ==========\n")
cat("生成文件 (在当前工作目录):\n")
cat("  - stability_jaccard_<name>.{csv,pdf}    # cluster 级 Jaccard 表+面板 B\n")
cat("  - stability_coassign_<name>.{csv,pdf}   # co-assignment 矩阵+面板 A\n")
cat("  其中 <name> ∈ {stroke, IS, ICH, SAH}\n")
cat("\n出片检查: pdffonts 应显示 emb=yes；pdftotext 后应无 U+2212\n")