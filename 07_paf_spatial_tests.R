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
# 空间聚集性检验：检验风险结构聚类(基于2023年PAF)是否具有地理空间聚集性
# 沿用 Figure 1 趋势聚类同样的方法学框架
#
# 方法：
#   1. Join Count Statistic (主方法, Monte Carlo 99999次置换)
#   2. Mantel Test (地理距离 vs PAF特征距离, 99999次置换)
#   3. Moran's I on PC1/PC2 (敏感性, Monte Carlo 99999次)
#
# 空间权重：SCO 10国 Queen邻接 (与 spatial_cluster_test.R 完全一致)
#
# 运行前提：05_paf_attribution.R 已运行, main_res 等聚类结果对象已生成
#
# ---------------------------------------------------------------------------
# 2026-08 修订
#   【P1】加 set.seed(SEED)。此前脚本没有设种子，三个检验的 p 与联结数期望
#         每跑一次都不同，Table S12 是一次无法复现的运行快照（表里总卒中
#         Mantel p=0·0058、联结数期望 8·02，重跑得 0·0079 / 7·93）。
#         种子放在四次调用之外、只设一次：四个结局在同一随机数流中依次
#         运行，整体可复现，且彼此不共用同一段随机数。
#   【P2】置换次数 999 / 9999 → 99999。9999 次下 p̂ 在 0·006 附近的标准误
#         约 0·0008，而表里按四位小数报告，位数超出了蒙特卡洛精度。
#         联结数原本只有 999 次（函数默认参数），一并提高。
#   【P3】联结数的置换计数向量化。原写法在 99999 次下要跑约 450 万次
#         R 层循环；改为预先取邻接索引再向量比较，结果完全等价。
#
# 本脚本写出的 spatial_PAF_cluster_test_results.csv 是 Table S12、
# Figure 2B 与 Figure S10 三处标注的唯一数据源。重跑本脚本后，这三处
# 必须同步更新，不要让图和表各算各的。
###############################################################################

library(spdep)
library(vegan)
library(tidyverse)
library(geosphere)

SEED  <- 2023      # 见【P1】
NPERM <- 99999     # 三个检验统一的置换次数，见【P2】

# ============================================================================
# 1. SCO 10国 Queen 邻接矩阵 (与趋势聚类检验完全一致)
# ============================================================================
sco_countries <- c("Belarus", "China", "India", "Iran", "Kazakhstan",
                   "Kyrgyzstan", "Pakistan", "Russia", "Tajikistan", "Uzbekistan")

adj_pairs <- list(
  c("China", "Russia"), c("China", "Kazakhstan"), c("China", "Kyrgyzstan"),
  c("China", "Tajikistan"), c("China", "Pakistan"), c("China", "India"),
  c("Russia", "Kazakhstan"), c("Russia", "Belarus"),
  c("Kazakhstan", "Kyrgyzstan"), c("Kazakhstan", "Uzbekistan"),
  c("Kyrgyzstan", "Uzbekistan"), c("Kyrgyzstan", "Tajikistan"),
  c("Tajikistan", "Uzbekistan"),
  c("Iran", "Pakistan"),
  c("Pakistan", "India")
)

n <- length(sco_countries)
W <- matrix(0, n, n, dimnames = list(sco_countries, sco_countries))
for (p in adj_pairs) {
  W[p[1], p[2]] <- 1
  W[p[2], p[1]] <- 1
}

nb   <- mat2listw(W, style = "W", zero.policy = FALSE)
nb_B <- mat2listw(W, style = "B", zero.policy = FALSE)

# 邻接国家对的索引 (i<j)，供联结数置换向量化使用，见【P3】
adj_idx <- which(upper.tri(W) & W == 1, arr.ind = TRUE)
stopifnot(nrow(adj_idx) == length(adj_pairs))

# ============================================================================
# 2. 首都坐标 + 地理距离矩阵
# ============================================================================
capitals <- data.frame(
  country = sco_countries,
  lon = c(27.5667, 116.4074, 77.2090, 51.3890, 71.4491,
          74.5698, 73.0479, 37.6173, 68.7870, 69.2401),
  lat = c(53.9000, 39.9042, 28.6139, 35.6892, 51.1605,
          42.8746, 33.6844, 55.7558, 38.5598, 41.2995)
)
rownames(capitals) <- capitals$country
geo_dist <- distm(capitals[sco_countries, c("lon", "lat")],
                  fun = distHaversine) / 1000   # km

# ============================================================================
# 3. 核心检验函数 (与趋势聚类版本统一)
# ============================================================================
spatial_test_cluster <- function(cluster_vec, feat_matrix, outcome_name,
                                 nperm = NPERM) {
  # 保证与W顺序一致
  cluster_vec <- cluster_vec[sco_countries]
  cluster_fac <- factor(cluster_vec)
  k <- length(levels(cluster_fac))
  
  cat("\n========================================\n")
  cat("Outcome:", outcome_name, "| k =", k, "\n")
  cat("========================================\n")
  print(table(cluster_fac))
  
  # ---- 3.1 Join Count (Monte Carlo) ----
  # 向量化写法，与原双重循环等价：只数 i<j 且相邻、且同簇的国家对
  count_same_joins <- function(cvec) {
    sum(cvec[adj_idx[, 1]] == cvec[adj_idx[, 2]])
  }
  cvec_int  <- as.integer(cluster_fac)
  obs_same  <- count_same_joins(cvec_int)
  perm_same <- replicate(nperm, count_same_joins(sample(cvec_int)))
  p_mc <- (sum(perm_same >= obs_same) + 1) / (nperm + 1)
  
  cat(sprintf("\n[1] Monte Carlo Join Count:\n"))
  cat(sprintf("  观察值: %d  |  置换均值: %.2f  |  p (单侧, 聚集): %.4f\n",
              obs_same, mean(perm_same), p_mc))
  
  # ---- 3.2 Mantel ----
  # 风险结构聚类: PAF 矩阵是 RF×国家, 需要转置成 国家×RF, 然后按列(RF)标准化
  feat_for_dist <- scale(t(feat_matrix[, sco_countries]))   # 10×5, 国家在行
  feat_dist <- as.matrix(dist(feat_for_dist, method = "euclidean"))   # 10×10
  mantel_res <- mantel(as.dist(geo_dist),
                       as.dist(feat_dist),
                       method = "pearson",
                       permutations = nperm)
  cat("\n[2] Mantel Test (地理距离 vs PAF特征距离):\n")
  cat(sprintf("  Mantel r = %.3f  |  p = %.4f\n",
              mantel_res$statistic, mantel_res$signif))
  
  # ---- 3.3 Moran's I on PCs ----
  # PCA 也需要 国家×RF 的矩阵 (国家作为样本)
  pca_input <- scale(t(feat_matrix[, sco_countries]))   # 10×5
  pca <- prcomp(pca_input, scale. = FALSE)
  var_exp <- summary(pca)$importance[2, 1:min(2, ncol(pca$x))] * 100
  
  moran_pc1 <- moran.mc(pca$x[, 1], listw = nb, nsim = nperm, zero.policy = FALSE)
  moran_pc2 <- if (ncol(pca$x) >= 2) {
    moran.mc(pca$x[, 2], listw = nb, nsim = nperm, zero.policy = FALSE)
  } else NULL
  
  cat("\n[3] Moran's I on Principal Components:\n")
  cat(sprintf("  PC1 (%.1f%% var): I = %.3f, p = %.4f\n",
              var_exp[1], moran_pc1$statistic, moran_pc1$p.value))
  if (!is.null(moran_pc2)) {
    cat(sprintf("  PC2 (%.1f%% var): I = %.3f, p = %.4f\n",
                var_exp[2], moran_pc2$statistic, moran_pc2$p.value))
  }
  
  list(
    outcome = outcome_name, k = k, cluster = cluster_vec,
    join_count_obs = obs_same,
    join_count_perm_mean = mean(perm_same),
    join_count_p = p_mc,
    mantel_r = mantel_res$statistic,
    mantel_p = mantel_res$signif,
    moran_pc1 = moran_pc1$statistic,
    moran_pc1_p = moran_pc1$p.value,
    moran_pc2 = if (!is.null(moran_pc2)) moran_pc2$statistic else NA,
    moran_pc2_p = if (!is.null(moran_pc2)) moran_pc2$p.value else NA,
    var_exp_pc1 = var_exp[1],
    var_exp_pc2 = if (length(var_exp) >= 2) var_exp[2] else NA
  )
}

# ============================================================================
# 4. 提取聚类结果 + 特征矩阵
# 注：05_paf_attribution.R 中 run_cause 返回的是 invisible 对象，需要显式接收
# 如果之前没有保存，请重新运行：
#   main_res <- run_cause(dat, "Stroke", rfs = rf_order, save_prefix = "Fig3_Stroke")
#   is_res   <- run_cause(dat, "Ischemic stroke", rfs = rf_order, save_prefix = "SuppFig_IS")
#   ich_res  <- run_cause(dat, "Intracerebral hemorrhage", rfs = c("SBP","KD","FG","BMI"), save_prefix = "SuppFig_ICH")
#   sah_res  <- run_cause(dat, "Subarachnoid hemorrhage", rfs = c("SBP","BMI"), save_prefix = "SuppFig_SAH")
# ============================================================================

results_paf <- list()

# 【P1】整块运行只设一次种子。四个结局在同一随机数流中依次运行：
# 结果可复现，且四个检验不共用同一段随机数。
# 若只想重跑其中某一个结局，必须从这里整块重跑，否则随机数流对不上。
set.seed(SEED)

results_paf[["Stroke"]] <- spatial_test_cluster(
  cluster_vec  = main_res$res_paf$clusters,
  feat_matrix  = main_res$paf_mat,
  outcome_name = "Stroke (Total)"
)

results_paf[["IS"]] <- spatial_test_cluster(
  cluster_vec  = is_res$res_paf$clusters,
  feat_matrix  = is_res$paf_mat,
  outcome_name = "Ischemic stroke"
)

results_paf[["ICH"]] <- spatial_test_cluster(
  cluster_vec  = ich_res$res_paf$clusters,
  feat_matrix  = ich_res$paf_mat,
  outcome_name = "Intracerebral hemorrhage"
)

# SAH 只有2个RF，PCA 仅 PC1 有意义
results_paf[["SAH"]] <- spatial_test_cluster(
  cluster_vec  = sah_res$res_paf$clusters,
  feat_matrix  = sah_res$paf_mat,
  outcome_name = "Subarachnoid hemorrhage"
)

# ============================================================================
# 5. 汇总表
# ============================================================================
summary_paf <- do.call(rbind, lapply(results_paf, function(r) {
  data.frame(
    Outcome           = r$outcome,
    k                 = r$k,
    JoinCount_Obs     = r$join_count_obs,
    JoinCount_ExpMean = round(r$join_count_perm_mean, 2),
    JoinCount_P       = round(r$join_count_p, 5),
    Mantel_r          = round(r$mantel_r, 3),
    Mantel_P          = round(r$mantel_p, 5),
    Moran_PC1         = round(r$moran_pc1, 3),
    Moran_PC1_P       = round(r$moran_pc1_p, 5),
    Moran_PC2         = round(r$moran_pc2, 3),
    Moran_PC2_P       = round(r$moran_pc2_p, 5)
  )
}))

cat("\n\n========== 风险结构聚类·空间检验汇总 ==========\n")
print(summary_paf, row.names = FALSE)

write.csv(summary_paf,
          file.path(RAWDAT, "PAF/spatial_PAF_cluster_test_results.csv"),
          row.names = FALSE)
cat("\n结果已保存: spatial_PAF_cluster_test_results.csv\n")
cat(sprintf("种子 = %d，置换次数 = %d，p 的最小可能值 = %.5f\n",
            SEED, NPERM, 1 / (NPERM + 1)))
cat("提醒：Table S12、Figure 2B 与 Figure S10 的标注都源自这份 CSV，\n")
cat("      重跑后请三处同步更新。\n")