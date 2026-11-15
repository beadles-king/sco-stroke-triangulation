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

# ============================================================================
# CHARLS Baseline Table 1: Metabolic Risk Factors & Incident Stroke
# Three groups: No event / Stroke / Death (non-stroke)
# Design: Mean ± SD + ANOVA for continuous; n (%) + Chi-square for categorical
# P-value strategy: 3-group overall test (ANOVA / Chi-square / Fisher's exact)
# ============================================================================

# ---- 0. Packages ----
#install.packages('tableone')
library(tidyverse)
library(rstatix)
library(tableone)


# ---- SMD helper functions (for Table 1 overall SMD) ----
smd_cont_pair <- function(x1, x2) {
  x1 <- x1[!is.na(x1)]; x2 <- x2[!is.na(x2)]
  if (length(x1) < 2 || length(x2) < 2) return(NA_real_)
  pooled <- sqrt((sd(x1)^2 + sd(x2)^2) / 2)
  if (pooled == 0 || is.na(pooled)) return(NA_real_)
  abs(mean(x1) - mean(x2)) / pooled
}

smd_cat_pair <- function(x1, x2) {
  x1 <- x1[!is.na(x1)]; x2 <- x2[!is.na(x2)]
  if (length(x1) == 0 || length(x2) == 0) return(NA_real_)
  all_levels <- union(levels(factor(x1)), levels(factor(x2)))
  K <- length(all_levels)
  if (K < 2) return(NA_real_)
  p1 <- as.numeric(prop.table(table(factor(x1, levels = all_levels))))[1:(K-1)]
  p2 <- as.numeric(prop.table(table(factor(x2, levels = all_levels))))[1:(K-1)]
  T_diff <- matrix(p1 - p2, ncol = 1)
  build_S <- function(p) {
    K1 <- length(p); S <- matrix(0, K1, K1)
    for (i in 1:K1) for (j in 1:K1)
      S[i, j] <- if (i == j) p[i] * (1 - p[i]) else -p[i] * p[j]
    S
  }
  S <- (build_S(p1) + build_S(p2)) / 2
  S_inv <- tryCatch(solve(S), error = function(e) MASS::ginv(S))
  sqrt(as.numeric(t(T_diff) %*% S_inv %*% T_diff))
}

# 多组整体 SMD = max pairwise SMD (tableone 默认做法)
max_smd <- function(df, var, type = c("cont", "cat")) {
  type <- match.arg(type)
  pairs <- list(
    c("No event", "Stroke"),
    c("No event", "Death"),
    c("Stroke",   "Death")
  )
  smds <- sapply(pairs, function(pr) {
    x1 <- df[[var]][df$group == pr[1]]
    x2 <- df[[var]][df$group == pr[2]]
    if (type == "cont") smd_cont_pair(x1, x2) else smd_cat_pair(x1, x2)
  })
  max(smds, na.rm = TRUE)
}

fmt_smd_star <- function(s) {
  if (is.na(s) || is.infinite(s)) return("—")
  stars <- if (s >= 0.80)      "***"
  else if (s >= 0.50) "**"
  else if (s >= 0.20) "*"
  else                 ""
  sprintf("%.3f%s", s, stars)
}


# ---- 1. Load data ----
# 请修改为您的实际路径
dat <- read.csv("charls_dat.csv", stringsAsFactors = FALSE)
# 如果您的数据是RData: load("your_charls_data.RData")

# ---- 2. Build 3-group outcome variable ----
# No event: new_stroke == 0 且未发生竞争死亡 (在您的数据中应有death变量；
#   若您的 new_stroke 已编码为 0=无事件, 1=中风, 则需要从原始数据中合并death信息)
# 若 new_stroke 中已区分: 请告知以便调整
# 这里假设您有一个 death 列（1=非中风死亡, 0=存活）; 若无, 请根据实际字段调整

# 方案A: 若您的数据中已有 death 列
# dat <- dat %>%
#   mutate(group = case_when(
#     new_stroke == 1 ~ "Stroke",
#     new_stroke == 0 & death == 1 ~ "Death",
#     TRUE ~ "No event"
#   ))

# 方案B: 若 new_stroke 已编码为三分类 (0=No event, 1=Stroke, 2=Death)
# 请根据您的Fine-Gray建模时使用的event编码确认
dat <- dat %>%
  mutate(group = factor(
    case_when(
      new_stroke == 1 ~ "Stroke",
      new_stroke == 2 ~ "Death",
      TRUE ~ "No event"
    ),
    levels = c("No event", "Stroke", "Death")
  ))

# 验证样本量 (应为 n=8079, n1=577, n2=696)
cat("Group sample sizes:\n")
print(table(dat$group, useNA = "ifany"))

# ---- 3. Variable setup ----
# Label categorical variable levels (must match your original coding)
dat <- dat %>%
  mutate(
    gender     = factor(gender,    levels = c(1, 2), labels = c("Male", "Female")),
    # marriage: 1=Married/partnered, 2=Separated/divorced, 3=Widowed, 4=Never married
    marriage   = factor(marriage,  levels = c(1, 2, 3, 4),
                        labels = c("Married/partnered", "Separated/divorced",
                                   "Widowed", "Never married")),
    # education: 1=Less than lower secondary, 2=Upper secondary & vocational, 3=Tertiary
    education  = factor(education, levels = c(1, 2, 3),
                        labels = c("Less than lower secondary",
                                   "Upper secondary/vocational", "Tertiary")),
    rural      = factor(rural,     levels = c(0, 1), labels = c("Urban", "Rural")),
    smoke      = factor(smoke,     levels = c(0, 1), labels = c("No", "Yes")),
    drink      = factor(drink,     levels = c(0, 1), labels = c("No", "Yes")),
    `heart disease`     = factor(heartd,    levels = c(0, 1), labels = c("No", "Yes")),
    depression = factor(depression,levels = c(0, 1), labels = c("No", "Yes")),
    bpmed      = factor(bpmed,     levels = c(0, 1), labels = c("No", "Yes")),
    diamed     = factor(diamed,    levels = c(0, 1), labels = c("No", "Yes")),
    lipidmed   = factor(lipidmed,  levels = c(0, 1), labels = c("No", "Yes"))
  )

# NOTE: 请根据您原始数据中 rural 的编码 (0/1 还是 1/2) 调整
# 您给的示例数据中 rural=1 多见，请确认 0=Urban, 1=Rural 的编码方向

# Variable groups
cont_vars <- c("SBP", "FG", "LDL", "eGFR", "BMI",        # Panel A: Exposures
               "age", "time_yr")                          # 其他连续
cat_vars  <- c("gender", "marriage", "education", "rural",# Panel B: Demographics
               "smoke", "drink", "heart disease", "depression")  # Panel C: Lifestyle & comorbidities

# ---- 4. Helper functions ----

# Continuous: Mean ± SD by group + ANOVA p
summ_cont <- function(df, var) {
  overall <- df %>%
    summarise(mean = mean(.data[[var]], na.rm = TRUE),
              sd   = sd(.data[[var]], na.rm = TRUE),
              n_missing = sum(is.na(.data[[var]])))
  by_grp <- df %>%
    group_by(group) %>%
    summarise(mean = mean(.data[[var]], na.rm = TRUE),
              sd   = sd(.data[[var]], na.rm = TRUE),
              .groups = "drop")
  # ANOVA
  fml <- as.formula(paste(var, "~ group"))
  aov_res <- aov(fml, data = df)
  p_val   <- summary(aov_res)[[1]][["Pr(>F)"]][1]
  
  # Format output
  fmt_ms <- function(m, s) sprintf("%.2f ± %.2f", m, s)
  tibble(
    Variable = var,
    Overall         = fmt_ms(overall$mean, overall$sd),
    `No event`  = fmt_ms(by_grp$mean[by_grp$group == "No event"],
                         by_grp$sd[by_grp$group == "No event"]),
    `Stroke`    = fmt_ms(by_grp$mean[by_grp$group == "Stroke"],
                         by_grp$sd[by_grp$group == "Stroke"]),
    `Death`     = fmt_ms(by_grp$mean[by_grp$group == "Death"],
                         by_grp$sd[by_grp$group == "Death"]),
    SMD = fmt_smd_star(max_smd(df, var, "cont")),
    P = ifelse(p_val < 0.001, "<0.001", sprintf("%.3f", p_val)),
    Missing_n = overall$n_missing
  )
}

# Categorical: n (%) by group + Chi-square (Fisher if expected < 5)
summ_cat <- function(df, var) {
  # Overall header row with variable name & p-value
  tab <- table(df[[var]], df$group)
  tab_overall <- table(df[[var]])
  # Test choice: Fisher if any expected cell < 5
  chi <- suppressWarnings(chisq.test(tab))
  use_fisher <- any(chi$expected < 5)
  if (use_fisher) {
    p_val <- fisher.test(tab, simulate.p.value = TRUE, B = 1e5)$p.value
    test_note <- "Fisher's exact"
  } else {
    p_val <- chi$p.value
    test_note <- "Chi-square"
  }
  p_fmt <- ifelse(p_val < 0.001, "<0.001", sprintf("%.3f", p_val))
  n_missing <- sum(is.na(df[[var]]))
  
  # Header row
  header <- tibble(
    Variable = var,
    Overall = "", `No event` = "", `Stroke` = "", `Death` = "",
    P = p_fmt, Missing_n = n_missing,
    SMD = fmt_smd_star(max_smd(df, var, "cat"))
  )
  # Level rows: n (%)
  fmt_np <- function(n, total) sprintf("%d (%.1f%%)", n, 100 * n / total)
  levs <- levels(df[[var]])
  level_rows <- map_dfr(levs, function(lev) {
    tibble(
      Variable = paste0("  ", lev),
      Overall     = fmt_np(sum(df[[var]] == lev, na.rm = TRUE), sum(!is.na(df[[var]]))),
      `No event`  = fmt_np(tab[lev, "No event"],  sum(tab[, "No event"])),
      `Stroke`    = fmt_np(tab[lev, "Stroke"],    sum(tab[, "Stroke"])),
      `Death`     = fmt_np(tab[lev, "Death"],     sum(tab[, "Death"])),
      P = "",       
      SMD = "",
      Missing_n = NA_integer_
    )
  })
  bind_rows(header, level_rows)
}

# ---- 5. Build Table 1 ----
# Panel A: Metabolic exposures
panelA <- map_dfr(c("SBP", "FG", "LDL", "eGFR", "BMI"), ~summ_cont(dat, .x))
# Panel B: Demographics
panelB_cont <- summ_cont(dat, "age")
panelB_cat  <- map_dfr(c("gender", "marriage", "education", "rural"), ~summ_cat(dat, .x))
panelC <- map_dfr(c("smoke", "drink", "heart disease", "depression"), ~summ_cat(dat, .x))
panelD_med <- map_dfr(c("bpmed", "diamed", "lipidmed"), ~summ_cat(dat, .x))
# Panel D: Follow-up (supplementary)
panelD <- summ_cont(dat, "time_yr")

# Combine with panel headers
add_header <- function(title) tibble(
  Variable = title, Overall = "", `No event` = "",
  `Stroke` = "", `Death` = "", P = "", SMD = "", Missing_n = NA_integer_
)

table1 <- bind_rows(
  add_header("Metabolic exposures"),
  panelA,
  add_header("Demographics"),
  panelB_cont,
  panelB_cat,
  add_header("Lifestyle and comorbidities"),
  panelC,
  add_header("Medication use"),    # ← 新增 panel
  panelD_med,                      # ← 新增
  add_header("Follow-up time"),
  panelD
)

# Rename header with n for each group
n_total <- nrow(dat)
n_ne <- sum(dat$group == "No event")
n_st <- sum(dat$group == "Stroke")
n_dt <- sum(dat$group == "Death")

names(table1) <- c(
  "Characteristic",
  sprintf("Overall (N=%d)", n_total),
  sprintf("No event (n=%d)", n_ne),
  sprintf("Stroke (n=%d)", n_st),
  sprintf("Death (n=%d)", n_dt),
  "P-value",
  "SMD",
  "Missing, n"
)

# ---- 6. View & export ----
print(table1, n = Inf, width = Inf)

# Export to CSV
write.csv(table1, file.path(RAWDAT, "CHARLS/Table1_CHARLS_baseline.csv"), row.names = FALSE, fileEncoding = "UTF-8")

# ============================================================================
# OPTIONAL: Sensitivity check using gtsummary (另一种官方推荐的做法, 可对照验证)
# ============================================================================
# library(gtsummary)
# tbl <- dat %>%
#   select(group, SBP, FG, LDL, eGFR, BMI,
#          age, gender, marriage, education, rural,
#          smoke, drink, heartd, depression, time_yr) %>%
#   tbl_summary(
#     by = group,
#     statistic = list(all_continuous() ~ "{mean} ± {sd}",
#                      all_categorical() ~ "{n} ({p}%)"),
#     digits = all_continuous() ~ 2,
#     missing = "no"
#   ) %>%
#   add_overall() %>%
#   add_p(test = list(all_continuous() ~ "aov",
#                     all_categorical() ~ "chisq.test"),
#         pvalue_fun = ~style_pvalue(.x, digits = 3))
# tbl