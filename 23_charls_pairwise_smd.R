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
# Supplementary Table S1 (updated): Pairwise comparisons with SMD
# - P-values (Bonferroni-corrected) + Standardized Mean Differences
# ============================================================================

library(tidyverse)

# ============================================================================
# SMD 计算函数
# 连续变量: SMD = |mean_1 - mean_2| / sqrt((SD_1^2 + SD_2^2)/2)
# 分类变量: 使用 Yang & Dalton (2012) 多类别扩展公式
# ============================================================================

# ---- 连续变量 SMD ----
smd_cont <- function(x1, x2) {
  x1 <- x1[!is.na(x1)]; x2 <- x2[!is.na(x2)]
  m1 <- mean(x1); m2 <- mean(x2)
  s1 <- sd(x1);   s2 <- sd(x2)
  # Pooled SD: 两组方差的平均的开方 (Cohen's d for unequal n, 流行病学常用)
  pooled_sd <- sqrt((s1^2 + s2^2) / 2)
  if (pooled_sd == 0 || is.na(pooled_sd)) return(NA_real_)
  abs(m1 - m2) / pooled_sd
}

# ---- 分类变量 SMD (Yang & Dalton 2012 多类别公式) ----
# 适用于 2类别 (binary) 及 多类别 (polytomous) 变量
smd_cat <- function(x1, x2) {
  x1 <- x1[!is.na(x1)]; x2 <- x2[!is.na(x2)]
  if (length(x1) == 0 || length(x2) == 0) return(NA_real_)
  # 获取所有水平的比例
  all_levels <- union(levels(factor(x1)), levels(factor(x2)))
  p1 <- prop.table(table(factor(x1, levels = all_levels)))
  p2 <- prop.table(table(factor(x2, levels = all_levels)))
  # 对于 K 个水平, 取前 K-1 个 (最后一个冗余)
  K <- length(all_levels)
  if (K < 2) return(NA_real_)
  p1 <- as.numeric(p1)[1:(K-1)]
  p2 <- as.numeric(p2)[1:(K-1)]
  T_diff <- matrix(p1 - p2, ncol = 1)
  # 协方差矩阵 S = (S1 + S2) / 2
  build_S <- function(p) {
    K1 <- length(p)
    S <- matrix(0, K1, K1)
    for (i in 1:K1) for (j in 1:K1) {
      S[i, j] <- if (i == j) p[i] * (1 - p[i]) else -p[i] * p[j]
    }
    S
  }
  S <- (build_S(p1) + build_S(p2)) / 2
  # 若矩阵奇异 (极端不平衡)，用广义逆
  S_inv <- tryCatch(solve(S), error = function(e) MASS::ginv(S))
  smd_sq <- t(T_diff) %*% S_inv %*% T_diff
  sqrt(as.numeric(smd_sq))
}

# ---- 格式化 ----
fmt_p <- function(p) {
  if (is.na(p))  return("—")
  if (p < 0.0001) return("<0.0001")
  sprintf("%.4f", p)
}
fmt_smd <- function(s) {
  if (is.na(s)) return("—")
  # Cohen (1988): small = 0.20, medium = 0.50, large = 0.80
  stars <- if (s >= 0.80)      "***"
  else if (s >= 0.50) "**"
  else if (s >= 0.20) "*"
  else                 ""
  sprintf("%.3f%s", s, stars)
}

# ============================================================================
# Pairwise 函数 (同时返回 P 和 SMD)
# ============================================================================

pairwise_cont <- function(df, var) {
  pairs <- list(
    c("No event", "Stroke"),
    c("No event", "Death"),
    c("Stroke",   "Death")
  )
  # P 值 (Welch's t-test + Bonferroni)
  pvals <- sapply(pairs, function(pr) {
    x <- df[[var]][df$group == pr[1]]
    y <- df[[var]][df$group == pr[2]]
    t.test(x, y, var.equal = FALSE)$p.value
  })
  pvals_adj <- pmin(pvals * 3, 1)
  # SMD
  smds <- sapply(pairs, function(pr) {
    x <- df[[var]][df$group == pr[1]]
    y <- df[[var]][df$group == pr[2]]
    smd_cont(x, y)
  })
  tibble(
    Variable          = var,
    `No event vs Stroke P`      = fmt_p(pvals_adj[1]),
    `No event vs Stroke SMD`    = fmt_smd(smds[1]),
    `No event vs Death P`      = fmt_p(pvals_adj[2]),
    `No event vs Death SMD`    = fmt_smd(smds[2]),
    `Stroke vs Death P`      = fmt_p(pvals_adj[3]),
    `Stroke vs Death SMD`    = fmt_smd(smds[3])
  )
}

pairwise_cat <- function(df, var) {
  pairs <- list(
    c("No event", "Stroke"),
    c("No event", "Death"),
    c("Stroke",   "Death")
  )
  pvals <- sapply(pairs, function(pr) {
    sub <- df %>% filter(group %in% pr) %>% mutate(group = droplevels(group))
    tab <- table(sub[[var]], sub$group)
    if (any(dim(tab) < 2)) return(NA_real_)
    chi <- suppressWarnings(chisq.test(tab))
    if (any(chi$expected < 5)) {
      if (nrow(tab) > 2 || ncol(tab) > 2) {
        fisher.test(tab, simulate.p.value = TRUE, B = 1e5)$p.value
      } else {
        fisher.test(tab)$p.value
      }
    } else chi$p.value
  })
  pvals_adj <- pmin(pvals * 3, 1)
  smds <- sapply(pairs, function(pr) {
    x <- df[[var]][df$group == pr[1]]
    y <- df[[var]][df$group == pr[2]]
    smd_cat(x, y)
  })
  tibble(
    Variable          = var,
    `No event vs Stroke P`      = fmt_p(pvals_adj[1]),
    `No event vs Stroke SMD`    = fmt_smd(smds[1]),
    `No event vs Death P`      = fmt_p(pvals_adj[2]),
    `No event vs Death SMD`    = fmt_smd(smds[2]),
    `Stroke vs Death P`      = fmt_p(pvals_adj[3]),
    `Stroke vs Death SMD`    = fmt_smd(smds[3])
  )
}

# ============================================================================
# 组装表格
# ============================================================================
add_panel <- function(title) tibble(
  Variable = title,
  `No event vs Stroke P` = "", `No event vs Stroke SMD` = "",
  `No event vs Death P` = "", `No event vs Death SMD` = "",
  `Stroke vs Death P` = "", `Stroke vs Death SMD` = ""
)

tableS1 <- bind_rows(
  add_panel("Metabolic exposures"),
  map_dfr(c("SBP", "FG", "LDL", "eGFR", "BMI"), ~pairwise_cont(dat, .x)),
  add_panel("Demographics"),
  pairwise_cont(dat, "age"),
  map_dfr(c("gender", "marriage", "education", "rural"), ~pairwise_cat(dat, .x)),
  add_panel("Lifestyle and comorbidities"),
  map_dfr(c("smoke", "drink", "heart disease", "depression"), ~pairwise_cat(dat, .x)),
  add_panel("Follow-up time"),
  pairwise_cont(dat, "time_yr")
)

# 列名加上 group sizes
n_ne <- sum(dat$group == "No event")
n_st <- sum(dat$group == "Stroke")
n_dt <- sum(dat$group == "Death")

names(tableS1) <- c(
  "Characteristic",
  "No event vs Stroke, P", "No event vs Stroke, SMD",
  "No event vs Death, P", "No event vs Death, SMD",
  "Stroke vs Death, P", "Stroke vs Death, SMD"
)

# ---- 输出 ----
print(tableS1, n = Inf)

write.csv(tableS1, file.path(RAWDAT, "CHARLS/TableS1_pairwise_with_SMD.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

