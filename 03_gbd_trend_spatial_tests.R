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
# 空间聚集性检验：检验 feature-based clustering 是否具有地理空间聚集性
# 方法：
#   1. Join Count Statistic (主方法, 针对类别型cluster)
#   2. Mantel Test (地理距离 vs 特征距离, 连续变量交叉验证)
#   3. Moran's I on PC1/PC2 (敏感性分析)
#
# 空间权重：SCO 10国 Queen邻接（共享陆地边界）
# 运行前提：res_stroke, res (for subtypes) 对象已在 01_gbd_burden_trends.R 中生成
###############################################################################

library(spdep)      # Join Count, Moran's I
library(vegan)      # Mantel
library(tidyverse)
library(geosphere)  # 首都大圆距离

# ============================================================================
# 1. 定义 SCO 10国 Queen 邻接矩阵（基于真实陆地边界）
# ============================================================================
sco_countries <- c("Belarus", "China", "India", "Iran", "Kazakhstan",
                   "Kyrgyzstan", "Pakistan", "Russia", "Tajikistan", "Uzbekistan")

# Queen邻接对 (无向)
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

# 构建 10x10 邻接矩阵
n <- length(sco_countries)
W <- matrix(0, n, n, dimnames = list(sco_countries, sco_countries))
for (p in adj_pairs) {
  W[p[1], p[2]] <- 1
  W[p[2], p[1]] <- 1
}

cat("邻接矩阵度数分布:\n")
print(rowSums(W))
cat("总连接数（单向）:", sum(W)/2, "\n\n")

# 转为 spdep 需要的 listw 对象 (行标准化)
nb <- mat2listw(W, style = "W", zero.policy = FALSE)  # 都连通，不必 zero.policy
nb_B <- mat2listw(W, style = "B", zero.policy = FALSE) # Join Count 用二元权重

# ============================================================================
# 2. 首都经纬度 & 地理距离矩阵（Mantel 用）
# ============================================================================
capitals <- data.frame(
  country = sco_countries,
  lon = c(27.5667, 116.4074, 77.2090, 51.3890, 71.4491,
          74.5698, 73.0479, 37.6173, 68.7870, 69.2401),
  lat = c(53.9000, 39.9042, 28.6139, 35.6892, 51.1605,
          42.8746, 33.6844, 55.7558, 38.5598, 41.2995)
)
rownames(capitals) <- capitals$country
geo_dist <- distm(capitals[sco_countries, c("lon", "lat")], fun = distHaversine) / 1000  # km
rownames(geo_dist) <- colnames(geo_dist) <- sco_countries

# ============================================================================
# 3. 核心检验函数
# ============================================================================
spatial_test_cluster <- function(cluster_vec, feat_matrix, outcome_name,
                                 nperm = 99999) {
  # 保证cluster顺序与W一致
  cluster_vec <- cluster_vec[sco_countries]
  cluster_fac <- factor(cluster_vec)
  k <- length(levels(cluster_fac))
  
  cat("\n========================================\n")
  cat("Outcome:", outcome_name, "| k =", k, "\n")
  cat("========================================\n")
  print(table(cluster_fac))
  
  # ---- 3.1 Join Count (主方法) ----
  # joincount.multi 给出每个类内 (BB) / 类间 (BW) 的计数
  jc <- joincount.multi(cluster_fac, listw = nb_B, zero.policy = FALSE)
  cat("\n[1] Join Count Statistic (渐近检验):\n")
  print(round(jc, 3))
  
  # Monte Carlo 置换检验 (对 Jtot = 所有同类连接数)
  # 注：以下 observed_same 从未被使用，且其索引写法依赖 nb_B 的内部结构，
  #     容易在 spdep 版本变动后静默算错。obs_same 才是实际使用的观察值。
  # observed_same <- sum(cluster_fac[unlist(lapply(nb_B$neighbours,
  #                                                function(x) x[x > seq_along(nb_B$neighbours)[1]]))] ==
  #                        cluster_fac[rep(seq_along(nb_B$neighbours),
  #                                        sapply(nb_B$neighbours, function(x) sum(x > seq_along(nb_B$neighbours)[1])))])
  
  # 更稳健的方式：手工置换（向量化：预先取出上三角的邻接对，避免 99999 次双重循环）
  edges <- which(upper.tri(W) & W == 1, arr.ind = TRUE)
  stopifnot(nrow(edges) == sum(W) / 2)      # 15 条无向边
  ei <- edges[, 1]; ej <- edges[, 2]
  count_same_joins <- function(cvec) sum(cvec[ei] == cvec[ej])
  
  cint      <- as.integer(cluster_fac)
  obs_same  <- count_same_joins(cint)
  perm_same <- replicate(nperm, count_same_joins(sample(cint)))
  p_mc      <- (sum(perm_same >= obs_same) + 1) / (nperm + 1)
  
  cat(sprintf("\n[1b] Monte Carlo Join Count (同类连接数, %d 次置换):\n", nperm))
  cat(sprintf("  观察值: %d  |  置换均值: %.2f  |  p (单侧, 聚集): %.6f\n",
              obs_same, mean(perm_same), p_mc))
  
  # ---- 3.2 Mantel 检验 ----
  # 特征距离矩阵 (标准化后欧氏)
  feat_dist <- as.matrix(dist(scale(feat_matrix[sco_countries, ]), method = "euclidean"))
  mantel_res <- mantel(as.dist(geo_dist[sco_countries, sco_countries]),
                       as.dist(feat_dist), method = "pearson",
                       permutations = nperm)
  cat("\n[2] Mantel Test (地理距离 vs 特征距离):\n")
  cat(sprintf("  Mantel r = %.3f  |  p = %.6f  |  permutations = %d\n",
              mantel_res$statistic, mantel_res$signif, nperm))
  
  # ---- 3.3 Moran's I on PC1/PC2 ----
  pca <- prcomp(scale(feat_matrix[sco_countries, ]), scale. = FALSE)
  var_exp <- summary(pca)$importance[2, 1:2] * 100
  
  moran_pc1 <- moran.mc(pca$x[, 1], listw = nb, nsim = nperm, zero.policy = FALSE)
  moran_pc2 <- moran.mc(pca$x[, 2], listw = nb, nsim = nperm, zero.policy = FALSE)
  
  cat(sprintf("\n[3] Moran's I on Principal Components (Monte Carlo, %d sim):\n", nperm))
  cat(sprintf("  PC1 (%.1f%% var): I = %.3f, p = %.4f\n",
              var_exp[1], moran_pc1$statistic, moran_pc1$p.value))
  cat(sprintf("  PC2 (%.1f%% var): I = %.3f, p = %.4f\n",
              var_exp[2], moran_pc2$statistic, moran_pc2$p.value))
  
  # 结果汇总
  list(
    outcome = outcome_name,
    k = k,
    nperm = nperm,
    cluster = cluster_vec,
    join_count_obs = obs_same,
    join_count_perm_mean = mean(perm_same),
    join_count_p = p_mc,
    mantel_r = mantel_res$statistic,
    mantel_p = mantel_res$signif,
    moran_pc1 = moran_pc1$statistic,
    moran_pc1_p = moran_pc1$p.value,
    moran_pc2 = moran_pc2$statistic,
    moran_pc2_p = moran_pc2$p.value,
    var_exp_pc1 = var_exp[1],
    var_exp_pc2 = var_exp[2]
  )
}

# ============================================================================
# 4. 对四个outcome分别运行 (需要 res_stroke, 以及循环中的 res 对象)
# ============================================================================
# --- 依赖的全局对象断言（v6 §4.1：凡依赖全局对象的，开头必须断言）---
stopifnot(exists("res_stroke"), exists("res_list"))
stopifnot(setequal(names(res_stroke$clusters), sco_countries))
stopifnot(all(c("IS", "ICH", "SAH") %in% names(res_list)))

# --- 随机种子：只在此设一次，放在四次调用之外 ---
# 放进函数里每次重置也能复现，但那样四个检验会共用同一段随机数、不再独立。
# 代价：四个 outcome 必须整块重跑，不能只重跑其中一个。
set.seed(2023)

# 先跑主分析 (Stroke)
results_all <- list()

results_all[["Stroke"]] <- spatial_test_cluster(
  cluster_vec  = res_stroke$clusters,
  feat_matrix  = res_stroke$feat_matrix,
  outcome_name = "Stroke (Total)"
)

# 三个亚型 —— 需要改造原循环以保存 res 对象
# 建议修改原 01_gbd_burden_trends.R 的亚型循环部分，新增保存：
#   res_list <- list()
#   for (si in seq_along(subtypes)) {
#     ...
#     res <- feature_based_cluster(general_year, stype)
#     res_list[[subtype_short[si]]] <- res
#     ...
#   }
# 然后在此脚本中：
subtype_short <- c("IS", "ICH", "SAH")
subtype_full  <- c("Ischemic stroke", "Intracerebral hemorrhage", "Subarachnoid hemorrhage")

for (si in seq_along(subtype_short)) {
  ss <- subtype_short[si]
  if (!is.null(res_list[[ss]])) {
    results_all[[ss]] <- spatial_test_cluster(
      cluster_vec  = res_list[[ss]]$clusters,
      feat_matrix  = res_list[[ss]]$feat_matrix,
      outcome_name = subtype_full[si]
    )
  }
}

# ============================================================================
# 5. 汇总表
# ============================================================================
summary_tbl <- do.call(rbind, lapply(results_all, function(r) {
  data.frame(
    Outcome = r$outcome,
    k = r$k,
    # p 值一律全精度写盘：99999 次置换下最小可达 1e-5，round(.,4) 会把它压成 0
    JoinCount_Obs = r$join_count_obs,
    JoinCount_ExpMean = r$join_count_perm_mean,
    JoinCount_P = r$join_count_p,
    Mantel_r = r$mantel_r,
    Mantel_P = r$mantel_p,
    Moran_PC1 = r$moran_pc1,
    Moran_PC1_P = r$moran_pc1_p,
    Moran_PC2 = r$moran_pc2,
    Moran_PC2_P = r$moran_pc2_p,
    N_perm = r$nperm
  )
}))
cat("\n\n========== 最终汇总 ==========\n")
print(summary_tbl, row.names = FALSE)

write.csv(summary_tbl, file.path(RAWDAT, "总体、年龄、性别负担/总体/spatial_clustering_test_results.csv"), row.names = FALSE)
cat("\n结果已保存: spatial_clustering_test_results.csv\n")