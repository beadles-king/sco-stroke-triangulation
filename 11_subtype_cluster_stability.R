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
# 11_subtype_cluster_stability.R  v3 (绘图对齐 02_gbd_trend_cluster_stability.R / Figure S4)
#
# 本版改动仅限绘图与出图部分（run_stability / run_coassignment 的分析逻辑、
# 随机种子、B 次数、jitter 比例一律未动，仍是 v2 BUG-FIXED 的口径）：
#   - 小数号一律中点 ·（格内数字、y 轴刻度标签），与 02_gbd_trend_cluster_stability.R 的
#     dec()/fmt2() 一致
#   - 字号：基础 8 pt，图内最小 7 pt；幅面 ≤180 mm，无需期刊缩放
#   - 出图改用 ragg::agg_tiff 直接栅格化到目标 dpi（600/300 两档），不再用
#     ggsave 的 png/pdf —— 与全文其余图同一引擎、同一份 Helvetica
#   - 去掉图内标题（面板字母与图注在拼图阶段添加）
#   - 配色改为全文主色 COL_MAIN 单色序列（原 viridis 多色 magma 与全文配色
#     体系不统一）；三条阈值线改统一灰色虚线，不再用红/橙/绿三色
#   - 热图格内数字按底色亮度自动切黑白字
#   - 方法顺序改为 Bootstrap → Jitter → Subset，与附录行文一致
#   - Figure S16 幅面：仅出 ASDR 两面板（A 热图 + B 条形图），与 Figure S4
#     的处理方式一致（S4 也只展示总卒中、不展示三个亚型）——四面板会超出
#     180×240 mm 供应图上限。ASIR 的稳定性数据已在 Table S13 展示，不重复
#     出图；本脚本仍计算 ASIR 的 CSV（供 Table S13 使用）但不再为它出图。
#
# 用法 (在 10_subtype_composition.R 跑完后):
#   source('11_subtype_cluster_stability.R')
#
# 输出 (写到当前工作目录):
#   stability_subtype_jaccard_<name>.csv          # ASDR 和 ASIR 都出
#   stability_subtype_coassign_<name>.csv         # ASDR 和 ASIR 都出
#   stability_subtype_coassign_ASDR(.tif/_300dpi.tif)   # Figure S16 面板 A
#   stability_subtype_jaccard_ASDR(.tif/_300dpi.tif)    # Figure S16 面板 B
###############################################################################
setwd(file.path(RAWDAT, "PAF/亚型占比聚类稳定性"))
# ============================================================================
# 0. 依赖
# ============================================================================
if (!requireNamespace("fpc",      quietly = TRUE)) install.packages("fpc")
if (!requireNamespace("reshape2", quietly = TRUE)) install.packages("reshape2")
if (!requireNamespace("ragg",     quietly = TRUE)) install.packages("ragg")
if (!requireNamespace("scales",   quietly = TRUE)) install.packages("scales")

library(fpc)
library(reshape2)
library(scales)
library(ggplot2)
library(dplyr)
library(tidyr)
library(grid)

# ============================================================================
# 配置
# ============================================================================
B_BOOT       <- 1000
SEED         <- 20260509
K_RANGE      <- 2:4
JITTER_RATIO <- 0.05   # 噪声 SD = JITTER_RATIO × 每列 SD

K_REPORT_SUB <- list(ASDR = NULL, ASIR = NULL)   # NULL = 用 k_paf

# ---- 全文格式约定（同 02_gbd_trend_cluster_stability.R / 交接文档第 2 节）----
DEC       <- "\u00B7"
BASE_SIZE <- 8
FONT      <- "Helvetica"
COL_MAIN  <- "#2C5985"
COL_HIGH  <- "#A63603"
COL_NEUT  <- "#7F7F7F"

dec  <- function(x) gsub("\\.", DEC, as.character(x))
fmt2 <- function(x) dec(sprintf("%.2f", x))

txt_col <- function(fill) {
  rgb <- col2rgb(fill)[, 1] / 255
  lum <- 0.299 * rgb[1] + 0.587 * rgb[2] + 0.114 * rgb[3]
  if (lum < 0.55) "white" else "black"
}

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
      plot.title       = element_blank(),
      plot.margin      = margin(2, 2, 2, 2, "mm")
    )
}

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
# 1. 准备特征矩阵  [未改动]
# ============================================================================
prepare_feat_subtype <- function(clust_obj) {
  feat <- as.matrix(clust_obj$ilr)
  colnames(feat) <- paste0("ILR", seq_len(ncol(feat)))
  feat
}

infer_k_used <- function(clust_obj) length(unique(clust_obj$cluster_k_paf))

# ============================================================================
# 2. cluster 级 Jaccard 稳定性  [未改动]
# ============================================================================
run_stability <- function(feat, k_range = K_RANGE,
                          B = B_BOOT, seed = SEED) {
  out <- list()
  for (ki in k_range) {
    cat("    k =", ki, "...\n")
    cb <- clusterboot(
      feat,
      B             = B,
      bootmethod    = c("boot", "jitter", "subset"),
      clustermethod = hclustCBI,
      scaling       = FALSE,
      method        = "ward.D2",
      k             = ki,
      seed          = seed,
      jittertuning  = JITTER_RATIO,
      subtuning     = ceiling(0.8 * nrow(feat)),
      count         = FALSE,
      showplots     = FALSE
    )
    df <- data.frame(
      k                = ki,
      cluster          = seq_len(ki),
      n_members        = as.integer(table(cb$partition)),
      members          = sapply(seq_len(ki), function(c) {
        paste(rownames(feat)[cb$partition == c], collapse = ", ")
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
# 3. co-assignment matrix  [未改动]
# ============================================================================
run_coassignment <- function(feat, k, B = B_BOOT,
                             jitter_ratio = JITTER_RATIO, seed = SEED) {
  set.seed(seed)
  n         <- nrow(feat)
  countries <- rownames(feat)
  
  col_sds   <- apply(feat, 2, sd)
  noise_sds <- jitter_ratio * col_sds
  
  M <- matrix(0, n, n, dimnames = list(countries, countries))
  for (b in seq_len(B)) {
    noise <- sapply(noise_sds, function(s) rnorm(n, 0, s))
    X     <- feat + noise
    cl    <- cutree(hclust(dist(X), method = "ward.D2"), k = k)
    M     <- M + outer(cl, cl, "==")
  }
  M / B
}

# ============================================================================
# 4. Jaccard 条形图  → Figure S16 面板 B
# ============================================================================
plot_jaccard <- function(stab_df) {
  long <- stab_df %>%
    pivot_longer(cols = c("jaccard_boot", "jaccard_jitter", "jaccard_subset"),
                 names_to = "method", values_to = "jaccard") %>%
    mutate(method = recode(method,
                           jaccard_boot   = "Bootstrap",
                           jaccard_jitter = "Jitter",
                           jaccard_subset = "Subset"),
           method = factor(method, levels = c("Bootstrap", "Jitter", "Subset")))
  
  ggplot(long, aes(x = factor(cluster), y = jaccard, fill = method)) +
    geom_col(position = position_dodge(0.8), width = 0.72) +
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
# 5. Co-assignment heatmap  → Figure S16 面板 A
# ============================================================================
plot_coassign <- function(M, clusters) {
  ord    <- order(clusters, names(clusters))
  M      <- M[ord, ord]
  cl_ord <- clusters[ord]
  
  long <- melt(M)
  colnames(long) <- c("c1", "c2", "freq")
  long$c1 <- factor(long$c1, levels = rownames(M))
  long$c2 <- factor(long$c2, levels = rownames(M))
  
  ramp <- colorRampPalette(c("#FFFFFF", COL_MAIN))
  long$lab_col <- vapply(long$freq,
                         function(f) txt_col(ramp(101)[round(f * 100) + 1]),
                         character(1))
  
  bdy <- cumsum(table(cl_ord))
  bdy <- bdy[-length(bdy)] + 0.5
  
  ggplot(long, aes(c1, c2, fill = freq)) +
    geom_tile(colour = "white", linewidth = 0.3) +
    geom_text(aes(label = fmt2(freq)), colour = long$lab_col,
              size = 7 / .pt) +
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
# 6. 主流程 + 一致性检查
#    ASDR: 计算 + 出图 (Figure S16 A/B)
#    ASIR: 只计算 CSV (供 Table S13), 不出图 —— 与 Figure S4 只展示总卒中、
#          不展示三个亚型是同一处理方式, 四面板会超出 180x240 mm 供应图上限
# ============================================================================
run_full_stability_subtype <- function(clust_obj, name, label,
                                       k_use = NULL, B = B_BOOT,
                                       make_figure = TRUE) {
  feat <- prepare_feat_subtype(clust_obj)
  if (is.null(k_use)) k_use <- infer_k_used(clust_obj)
  
  clusters_true <- cutree(clust_obj$hc, k = k_use)
  
  cat("\n========== Subtype stability:", label,
      "(n =", nrow(feat), ", p =", ncol(feat),
      ", reporting k =", k_use, ") ==========\n")
  
  cat("  真实 cutree(hc, k =", k_use, ") 分组:\n")
  for (ci in sort(unique(clusters_true))) {
    cat("    C", ci, ":",
        paste(names(clusters_true)[clusters_true == ci], collapse = ", "), "\n")
  }
  
  cat("  Step 1/2: cluster-level Jaccard ...\n")
  stab <- run_stability(feat, K_RANGE, B, SEED)
  
  stab_kuse <- stab[stab$k == k_use, ]
  cat("  稳定性脚本 original partition (k =", k_use, "):\n")
  for (i in seq_len(nrow(stab_kuse))) {
    cat("    C", stab_kuse$cluster[i], ":", stab_kuse$members[i], "\n")
  }
  if (!identical(sort(stab_kuse$n_members), sort(as.integer(table(clusters_true))))) {
    cat("  ⚠ 警告: 真实 cutree 与稳定性 partition 大小不一致, 请检查!\n")
  } else {
    cat("  ✓ 真实 cutree 与稳定性 partition 一致\n")
  }
  
  csv_jac <- paste0("stability_subtype_jaccard_", name, ".csv")
  write.csv(stab, csv_jac, row.names = FALSE)
  cat("  saved:", csv_jac, "\n")
  
  cat("  Step 2/2: country-level co-assignment matrix (k =", k_use, ")...\n")
  M <- run_coassignment(feat, k = k_use, B = B, seed = SEED)
  
  csv_coa <- paste0("stability_subtype_coassign_", name, ".csv")
  write.csv(round(M, 3), csv_coa)
  cat("  saved:", csv_coa, "\n")
  
  if (make_figure) {
    # 面板尺寸与 Figure S4 (02_gbd_trend_cluster_stability.R) 完全一致: 同为 10 国、
    # 同为 k=2:4 的 Jaccard 条形图, 直接复用同一套幅面, 无需重新试参数。
    save_panel(plot_coassign(M, clusters_true),
               paste0("stability_subtype_coassign_", name),
               w_mm = 152, h_mm = 126)
    save_panel(plot_jaccard(stab),
               paste0("stability_subtype_jaccard_", name),
               w_mm = 176, h_mm = 62)
  }
  
  invisible(list(jaccard = stab, coassign = M, k_used = k_use))
}

cat("\n##############################################################\n")
cat("# 亚型占比模式聚类稳定性 v3 (B =", B_BOOT, ")\n")
cat("##############################################################\n")

stab_subtype <- list()
stab_subtype$ASDR <- run_full_stability_subtype(
  clust_daly, "ASDR", "ASDR-based subtype composition",
  k_use = K_REPORT_SUB$ASDR, make_figure = TRUE)
stab_subtype$ASIR <- run_full_stability_subtype(
  clust_inc,  "ASIR", "ASIR-based subtype composition",
  k_use = K_REPORT_SUB$ASIR, make_figure = FALSE)

cat("\n========== 亚型占比稳定性分析完成 ==========\n")
cat("Figure S16 用 stability_subtype_coassign_ASDR / stability_subtype_jaccard_ASDR\n")
cat("ASIR 的 CSV 已生成供 Table S13 核对, 不出图\n")