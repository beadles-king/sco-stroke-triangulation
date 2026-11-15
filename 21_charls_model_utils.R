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
# CHARLS 竞争风险模型 —— 共用设计矩阵工具  (v3)
#
# File: 21_charls_model_utils.R (in this repository)
# 由 24_charls_fine_gray.R / 26_charls_subgroups.R / 30_charls_cystatin_c.R 三个脚本 source。
#
# 解决两个贯穿全部 Fine-Gray 分析的问题：
#
#   F1  marriage（1=已婚/同居, 2=分居/离异, 3=丧偶, 4=未婚）是无序名义变量，
#       education（1/2/3）是有序三分类，但此前都以单一数值列进入
#       crr(cov1=)，被拟合成一个线性系数。本工具统一转哑变量
#       （参考水平 = 编码最小值），与 Table S4 的分类呈现一致。
#
#   F2  Model 3 追加的三个药物变量（bpmed 40 / diamed 76 / lipidmed 165 例
#       缺失）使 Model 3 的分析人群小于 Models 1–2，此前未报告。本工具显式
#       做成列删除，并返回删除前/后的人数与事件数。
#
# ---------------------------------------------------------------------------
# v3 的唯一改动【p 值精度】
#   cmprsk::summary.crr 对 p 值那一列做了 signif(, digits)，digits 默认
#   max(options()$digits - 5, 2) = 2，即所有 p 只保留两位有效数字：
#       coef.table[, 5] <- signif(2 * (1 - pnorm(abs(z))), digits)
#   全文 p 值规则要求 0·0001 <= p < 0·01 保留四位小数，两位有效数字不够用
#   （Model 3 的 FG 只能得到 0·045，无法判定应显示 0·04 还是 0·05）。
#   另外 1 - pnorm(|z|) 在 |z| > 8·3 时下溢为 0，Model 1 的 SBP/BMI 因此
#   显示成 0·0e+00。
#
#   现改为直接由全精度的 coef 与 se 计算：P = 2 * pnorm(-|coef/se|)。
#   下侧尾概率不下溢，可算到 1e-300 量级，且不经过任何 signif()。
#   coef 与 se(coef) 两列本就是全精度，故 SHR 与 CI 完全不变；
#   本次改动只让 p 值变精确，不改动任何点估计或区间。
# ---------------------------------------------------------------------------
#
# v2 相对 v1 的两处加固（哑变量化后才会出现的问题）：
#   - X 建好后剔除零方差列（层内某个婚姻/教育水平为空时会产生全 0 列，
#     否则 crr 设计矩阵奇异）
#   - 对观测数少于 MIN_LEVEL_N 的分类水平发警告（层内极稀疏水平会导致
#     系数发散、CI 爆炸，需人工决定是否合并水平）
#
# 用法：
#   source("21_charls_model_utils.R")
#   des <- build_design(charls, expo = EXPO_SD, covs = COV_M2, label = "M2")
#   fit <- fit_crr(des)
#   extract_fg(fit, "SBP_sd"); des$n; des$n_stroke
# ============================================================================

# ---- 变量集定义（与 Table S4 一致）----------------------------------------
EXPO_RAW <- c("SBP", "FG", "LDL", "eGFR", "BMI")
EXPO_SD  <- c("SBP_sd", "FG_sd", "LDL_sd", "eGFR_sd", "BMI_sd")

COV_M2 <- c("age", "female", "marriage", "education",
            "rural", "smoke", "drink", "heartd", "depression")
COV_M3 <- c(COV_M2, "bpmed", "diamed", "lipidmed")

# 需要哑变量化的分类协变量（其余均为连续或 0/1 二值）
FACTOR_COVS  <- c("marriage", "education")
MIN_LEVEL_N  <- 10    # 分类水平观测数低于此值时发警告

DEC <- "\u00B7"

# ---- p 值：全精度计算与全文规则格式化 --------------------------------------
# 绕开 summary.crr 的 signif() 截断与 1-pnorm 的下溢
z_to_p <- function(coef, se) 2 * stats::pnorm(-abs(coef / se))

# 全文规则：p<0·0001 -> "<0·0001"；0·0001<=p<0·01 -> 四位小数；p>=0·01 -> 两位小数
fmt_p <- function(p) {
  vapply(p, function(pi) {
    if (is.na(pi))   return(NA_character_)
    if (pi < 0.0001) return(paste0("<0", DEC, "0001"))
    if (pi < 0.01)   return(gsub("\\.", DEC, sprintf("%.4f", pi)))
    gsub("\\.", DEC, sprintf("%.2f", pi))
  }, character(1))
}

# 效应量显示串：CHARLS 部分两位小数（MR 部分三位，另行处理）
fmt_est <- function(est, low, high, digits = 2) {
  f <- function(x) gsub("\\.", DEC, sprintf(paste0("%.", digits, "f"), x))
  sprintf("%s (%s\u2013%s)", f(est), f(low), f(high))
}

# ---- 设计矩阵构造 ----------------------------------------------------------
# 返回 list：
#   X            数值设计矩阵（无截距），可直接传给 crr(cov1=)
#   d            成列删除后的数据框（含 time / status / 原始列）
#   n_before / n            删除前 / 后行数（后者即"进模型人数"）
#   ev_before / n_stroke    删除前 / 后事件数
#   n_death      删除后竞争死亡数
#   dropped_cov  因零方差或与结局完全分离而剔除的协变量
#   dropped_col  因零方差而剔除的设计矩阵列（通常是空的哑变量水平）
#   dummy_map    哑变量列名 → 原变量与参考水平
build_design <- function(data,
                         expo,
                         covs,
                         factor_covs = FACTOR_COVS,
                         time_var    = "time_yr",
                         status_var  = "new_stroke",
                         extra_cols  = NULL,
                         check_separation = TRUE,
                         verbose     = TRUE,
                         label       = "") {
  
  need <- c(time_var, status_var, expo, covs, extra_cols)
  miss <- setdiff(need, names(data))
  if (length(miss) > 0) {
    stop("数据中缺少以下列：", paste(miss, collapse = ", "))
  }
  
  n_before  <- nrow(data)
  ev_before <- sum(data[[status_var]] == 1, na.rm = TRUE)
  
  # --- 成列删除：时间、结局、全部暴露、全部协变量 ---
  d <- data[stats::complete.cases(data[, need, drop = FALSE]), , drop = FALSE]
  if (nrow(d) == 0) stop("[", label, "] 成列删除后没有剩余观测")
  
  # --- 剔除该分析集内取值单一的协变量（如男性层内的 female）---
  keep <- covs[vapply(covs, function(v) length(unique(d[[v]])) > 1, logical(1))]
  dropped <- setdiff(covs, keep)
  
  # --- 剔除与结局完全分离的二分类协变量（否则 crr 不收敛）---
  if (check_separation && length(keep) > 0) {
    binary_like <- keep[vapply(keep, function(v) {
      !(v %in% factor_covs) && length(unique(d[[v]])) == 2
    }, logical(1))]
    sep <- binary_like[vapply(binary_like, function(v) {
      tab <- table(d[[v]], d[[status_var]] == 1)
      any(dim(tab) < 2) || min(tab) == 0
    }, logical(1))]
    if (length(sep) > 0) {
      keep    <- setdiff(keep, sep)
      dropped <- c(dropped, sep)
    }
  }
  
  # --- 分类协变量转因子（参考水平 = 编码最小值）---
  fc <- intersect(factor_covs, keep)
  for (v in fc) {
    d[[v]] <- factor(d[[v]], levels = sort(unique(d[[v]])))
    tb   <- table(d[[v]])
    thin <- names(tb)[tb < MIN_LEVEL_N]
    if (length(thin) > 0) {
      warning(sprintf(
        "[%s] %s 的水平 %s 观测数不足 %d（实际 %s），估计可能不稳定，考虑合并水平。",
        label, v, paste(thin, collapse = "/"), MIN_LEVEL_N,
        paste(tb[thin], collapse = "/")), call. = FALSE)
    }
  }
  
  # --- model.matrix 生成哑变量，去掉截距列 ---
  rhs <- paste(c(expo, keep), collapse = " + ")
  X   <- stats::model.matrix(stats::as.formula(paste("~", rhs)), data = d)
  X   <- X[, setdiff(colnames(X), "(Intercept)"), drop = FALSE]
  
  # --- 剔除零方差列（空的哑变量水平会产生全 0 列，使设计矩阵奇异）---
  nzv         <- apply(X, 2, function(col) length(unique(col)) > 1)
  dropped_col <- colnames(X)[!nzv]
  X <- X[, nzv, drop = FALSE]
  
  dummy_map <- NULL
  if (length(fc) > 0) {
    dummy_map <- do.call(rbind, lapply(fc, function(v) {
      cols <- grep(paste0("^", v), colnames(X), value = TRUE)
      if (length(cols) == 0) return(NULL)
      data.frame(column = cols, source = v,
                 reference = levels(d[[v]])[1], row.names = NULL)
    }))
  }
  
  if (verbose) {
    cat(sprintf("  [%s] n: %d -> %d (剔除 %d)；stroke: %d -> %d；设计矩阵 %d 列\n",
                label, n_before, nrow(d), n_before - nrow(d),
                ev_before, sum(d[[status_var]] == 1), ncol(X)))
    if (length(dropped) > 0)
      cat(sprintf("      剔除协变量（零方差/完全分离）：%s\n",
                  paste(dropped, collapse = ", ")))
    if (length(dropped_col) > 0)
      cat(sprintf("      剔除零方差设计列：%s\n",
                  paste(dropped_col, collapse = ", ")))
  }
  
  list(X = X, d = d,
       n_before = n_before, n = nrow(d),
       ev_before = ev_before,
       n_stroke = sum(d[[status_var]] == 1),
       n_death  = sum(d[[status_var]] == 2),
       dropped_cov = dropped, dropped_col = dropped_col,
       dummy_map   = dummy_map,
       time_var = time_var, status_var = status_var)
}

# ---- 便捷拟合（extra_X 用于追加交互项等额外列）----------------------------
fit_crr <- function(des, extra_X = NULL) {
  X <- if (is.null(extra_X)) des$X else cbind(des$X, extra_X)
  cmprsk::crr(ftime    = des$d[[des$time_var]],
              fstatus  = des$d[[des$status_var]],
              cov1     = X,
              failcode = 1, cencode = 0)
}

# ---- 结果提取 --------------------------------------------------------------
# 【v3】P 由 coef/se 直接算，不再读 summary 表里被 signif(,2) 截断的那一列
extract_fg <- function(fg_fit, var_name) {
  s   <- summary(fg_fit)$coef
  idx <- which(rownames(s) == var_name)
  if (length(idx) == 0) stop("模型中没有变量：", var_name)
  b  <- s[idx, "coef"]
  se <- s[idx, "se(coef)"]
  p  <- z_to_p(b, se)
  data.frame(
    variable = var_name,
    SHR      = exp(b),
    CI_low   = exp(b - 1.96 * se),
    CI_high  = exp(b + 1.96 * se),
    coef     = b,
    se       = se,
    P        = p,
    P_exact  = sprintf("%.3e", p),   # 供人工核对
    P_disp   = fmt_p(p),             # 按全文规则的显示串
    row.names = NULL
  )
}

# 一次提取多个变量
extract_fg_many <- function(fg_fit, var_names) {
  do.call(rbind, lapply(var_names, function(v) extract_fg(fg_fit, v)))
}

# ---- 分析人数记账表（供 Table S15/S19/S20/S21/S22 表注引用）---------------
n_ledger <- function(...) {
  rows <- list(...)
  do.call(rbind, lapply(names(rows), function(nm) {
    des <- rows[[nm]]
    data.frame(analysis = nm,
               n_before = des$n_before, n_analysis = des$n,
               events_before = des$ev_before, events_analysis = des$n_stroke,
               n_dropped = des$n_before - des$n,
               row.names = NULL)
  }))
}