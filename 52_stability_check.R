# =============================================================================
# Reproducibility repository - Wang et al., transportability-aware triangulation
# of metabolic stroke evidence across the Shanghai Cooperation Organisation.
# Paths are resolved by config.R; see README.md.
# =============================================================================
if (!exists("RAWDAT")) {
  if (file.exists("config.R")) source("config.R") else
    stop("Run this script from the code/ directory (the folder containing config.R).")
}

# ############################################################################
# 检查：PAF 稳定性脚本的加噪写法能否与亚型脚本统一
#
# 背景：
#   11_subtype_cluster_stability.R  noise_sds <- 0.05 * apply(feat, 2, sd)   # 每列 SD 的 5%
#   06_paf_cluster_stability.R    rnorm(..., sd = 0.05)                    # 绝对值 0.05
# 二者只有在 feat 的每一列 SD 恰为 1 时才等价。PAF 的输入是 t(mat_z)，
# 若 mat_z 由 scale() 或 (x-mean)/sd 生成，则列 SD 应为 1，两种写法数值等价，
# 改代码不会改结果、图不必重出；否则必须保留代码、改附录与图注表述。
#
# 运行时机：source 过 05_paf_attribution.R 与 06_paf_cluster_stability.R 之后
# ############################################################################

check_paf_noise_equivalence <- function(res_obj = NULL) {

  # ---- 自动定位 res 对象（也可直接把对象传进来）----
  if (is.null(res_obj)) {
    cand <- Filter(function(nm) {
      o <- get(nm, envir = .GlobalEnv)
      is.list(o) && !is.null(o$res_paf) && !is.null(o$res_paf$mat_z)
    }, ls(envir = .GlobalEnv))
    if (length(cand) == 0)
      stop("没找到含 $res_paf$mat_z 的对象，请手工传入：",
           "check_paf_noise_equivalence(你的对象名)")
    cat("使用对象：", cand[1], "\n")
    res_obj <- get(cand[1], envir = .GlobalEnv)
  }

  feat <- prepare_feat_paf(res_obj)      # Country x RF

  cat("\n===== D5 检查 =====\n")
  cat("特征矩阵维度（国家 x 风险因子）：", paste(dim(feat), collapse = " x "), "\n")
  cat("列名：", paste(colnames(feat), collapse = ", "), "\n\n")

  col_sd   <- apply(feat, 2, sd)
  col_mean <- colMeans(feat)
  print(round(rbind(SD = col_sd, Mean = col_mean), 6))

  all_one <- all(abs(col_sd - 1) < 1e-9)
  max_dev <- max(abs(col_sd - 1))

  cat(sprintf("\n各列 SD 与 1 的最大偏离：%.3e\n", max_dev))
  if (all_one) {
    cat("→ 每列 SD 恰为 1。绝对 0·05 与 0·05 x SD 完全等价：\n",
        "   可以把 06_paf_cluster_stability.R 改成 noise_sds <- 0.05 * apply(feat, 2, sd)，\n",
        "   结果一字不变，Figure S11 不用重出，图注与附录也不必改数值。\n", sep = "")
  } else if (max_dev < 1e-6) {
    cat("→ 每列 SD 在数值误差内等于 1（可能用了 n 分母而非 n-1）。\n",
        "   改写法后 co-assignment 频率理论上可能在小数第三位有别，\n",
        "   建议改完重跑一次，用下面的对照确认 Figure S11 面板 A 仍全为 0/1。\n", sep = "")
  } else {
    cat("→ 每列 SD 不为 1。两种写法不等价：\n",
        "   不要改代码；改附录方法段与 Figure S11 图注，把 PAF 那一层写成\n",
        "   \"gaussian noise with an absolute standard deviation of 0·05 on the\n",
        "    unit-scaled risk-factor matrix\"，与亚型层的\"5% of each feature's SD\"分开表述。\n",
        sep = "")
  }

  invisible(data.frame(variable = names(col_sd), sd = as.numeric(col_sd),
                       mean = as.numeric(col_mean), row.names = NULL))
}

# ---- 运行 ----
# 若知道对象名（例如 res_paf_stroke），直接传入更稳妥：
d5 <- check_paf_noise_equivalence(main_res)
# d5 <- check_paf_noise_equivalence()
write.csv(d5, file.path(RAWDAT, "CHARLS/D5_paf_feature_sd.csv"),
          row.names = FALSE)


# ---- 附：若上面判定为"不等价"，用这段做一次实证对照 --------------------
# 两种写法各跑一次 co-assignment，直接比较矩阵，比理论推断更可靠。
compare_noise_schemes <- function(feat, k = 2, B = 1000, seed = 20260509) {
  run <- function(sd_vec) {
    set.seed(seed); n <- nrow(feat)
    co <- matrix(0, n, n, dimnames = list(rownames(feat), rownames(feat)))
    for (b in seq_len(B)) {
      X  <- feat + matrix(rnorm(length(feat), 0, 1), nrow = n) *
              matrix(sd_vec, nrow = n, ncol = ncol(feat), byrow = TRUE)
      cl <- cutree(hclust(dist(X), method = "ward.D2"), k = k)
      co <- co + outer(cl, cl, "==")
    }
    co / B
  }
  a <- run(rep(0.05, ncol(feat)))              # 现行：绝对 0·05
  b <- run(0.05 * apply(feat, 2, sd))          # 拟改：每列 SD 的 5%
  cat(sprintf("两种写法 co-assignment 的最大绝对差：%.4f\n", max(abs(a - b))))
  list(absolute = a, proportional = b, max_diff = max(abs(a - b)))
}
# cmp <- compare_noise_schemes(prepare_feat_paf(res_paf_stroke))
