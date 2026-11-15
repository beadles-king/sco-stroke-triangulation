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
# 敏感性分析 2（回归稀释）：Fine-Gray 竞争风险模型 —— 修订版（Table S21）
# 暴露 = W1–W3 均值；基线 = W3；结局 = W3→W4
#
# 三处修改：
#  【F1】marriage / education 经 build_design() 转哑变量，不再作为线性项。
#  【F2】每个模型报告自身的进模型人数；旧版 cat(nrow(charls)) 报的是
#        协变量成列删除【之前】的行数，与实际拟合人群不符。
#  【尺度】显式报告本子样本的暴露 SD 与主队列 SD 的差异。本分析的暴露是
#        两波均值，方差本就小于单次测量，"每 1 SD" 与 Table S15 不同尺度，
#        Table S21 表注必须写明，否则读者会直接比较 SHR 大小。
#
# 另：比例亚分布风险检验改用 Model 2 协变量集，与主分析 24_charls_fine_gray.R 一致
#     （旧版用的是 Model 3 的 cov_mat_m3）。可用 PH_COVSET 切回。
#
# 依赖：21_charls_model_utils.R；输入 charls_dat_敏感.csv（由 28_charls_dilution_prepare.R 生成）
# ============================================================================

setwd(RAWDAT)
library(tidyverse)
library(cmprsk)

source(file.path(CODE_DIR, "21_21_charls_model_utils.R"))

OUT_DIR    <- file.path(RAWDAT, "CHARLS")
PH_COVSET  <- "M2"      # "M2" 与主分析一致；"M3" 复现旧版
MAIN_DATA  <- 'charls_dat.csv'   # 仅用于打印主队列 SD 作对照

# ============================================================================
# 1. 数据准备
# ============================================================================
charls_dil <- read.csv('charls_dat_敏感.csv')
charls_dil$female <- as.integer(charls_dil$gender == 2)
for (i in seq_along(EXPO_RAW))
  charls_dil[[EXPO_SD[i]]] <- as.numeric(scale(charls_dil[[EXPO_RAW[i]]]))

cat(sprintf("\n敏感性分析 2 数据：N = %d | 卒中 %d | 竞争死亡 %d | 无事件 %d\n",
            nrow(charls_dil), sum(charls_dil$new_stroke == 1),
            sum(charls_dil$new_stroke == 2), sum(charls_dil$new_stroke == 0)))
cat("（以上为暴露完整者；进模型人数见下方各模型的记账）\n")

# ---- 标准化尺度对照：本子样本 vs 主队列 ----
sd_sens <- vapply(EXPO_RAW, function(v) sd(charls_dil[[v]], na.rm = TRUE), numeric(1))
sd_main <- tryCatch({
  m <- read.csv(MAIN_DATA)
  vapply(EXPO_RAW, function(v) sd(m[[v]], na.rm = TRUE), numeric(1))
}, error = function(e) rep(NA_real_, length(EXPO_RAW)))

sd_cmp <- rbind(`两波均值(本分析)` = sd_sens,
                `主队列 W1 单次`   = sd_main,
                `比值`             = sd_sens / sd_main)
cat("\n===== 标准化基准 SD 对照 =====\n")
print(round(sd_cmp, 3))
cat("→ 比值明显 < 1 的暴露，其\"每 1 SD\" SHR 不能与 Table S15 直接比较大小。\n")

# ============================================================================
# 2. Models 1–3
# ============================================================================
cat("\n===== 各模型分析人群 =====\n")

m1 <- lapply(EXPO_SD, function(e) {
  des <- build_design(charls_dil, expo = e, covs = COV_M2, label = paste("M1", e))
  list(des = des,
       res = extract_fg(fit_crr(des), e) %>%
         mutate(model = "Model 1", n = des$n, n_stroke = des$n_stroke))
})
names(m1) <- EXPO_SD

des_m2 <- build_design(charls_dil, expo = EXPO_SD, covs = COV_M2, label = "Model 2")
des_m3 <- build_design(charls_dil, expo = EXPO_SD, covs = COV_M3, label = "Model 3")
fg_m2  <- fit_crr(des_m2)
fg_m3  <- fit_crr(des_m3)

res_m2 <- bind_rows(lapply(EXPO_SD, function(e) extract_fg(fg_m2, e))) %>%
  mutate(model = "Model 2", n = des_m2$n, n_stroke = des_m2$n_stroke)
res_m3 <- bind_rows(lapply(EXPO_SD, function(e) extract_fg(fg_m3, e))) %>%
  mutate(model = "Model 3", n = des_m3$n, n_stroke = des_m3$n_stroke)

sens2_table <- bind_rows(bind_rows(lapply(m1, `[[`, "res")), res_m2, res_m3) %>%
  select(model, variable, n, n_stroke, SHR, CI_low, CI_high, P)

cat("\n===== Table S21：两波平均 Fine-Gray（每 1 SD）=====\n")
print(sens2_table, digits = 3, row.names = FALSE)

ledger <- n_ledger(`Model 1 (SBP)` = m1$SBP_sd$des,
                   `Model 2`       = des_m2,
                   `Model 3`       = des_m3)
cat("\n===== 分析人数记账 =====\n"); print(ledger, row.names = FALSE)

# ============================================================================
# 3. Model 2 vs Model 3
# ============================================================================
compare_m2m3 <- res_m2 %>%
  select(variable, SHR_m2 = SHR, CI_low_m2 = CI_low,
         CI_high_m2 = CI_high, P_m2 = P) %>%
  left_join(res_m3 %>% select(variable, SHR_m3 = SHR, CI_low_m3 = CI_low,
                              CI_high_m3 = CI_high, P_m3 = P),
            by = "variable") %>%
  mutate(SHR_change_pct = (SHR_m3 - SHR_m2) / SHR_m2 * 100,
         n_m2 = des_m2$n, n_m3 = des_m3$n)
cat("\n===== Model 2 vs Model 3 =====\n")
print(compare_m2m3, digits = 3, row.names = FALSE)

# ============================================================================
# 4. 比例亚分布风险假设检验（暴露 × log 时间）
# ============================================================================
des_ph   <- if (PH_COVSET == "M3") des_m3 else des_m2
d_ph     <- des_ph$d
log_time <- log(d_ph$time_yr)
stopifnot(all(is.finite(log_time)))

t_cols <- sapply(EXPO_SD, function(e) d_ph[[e]] * log_time)
colnames(t_cols) <- sub("_sd$", "_t", EXPO_SD)

fg_ph <- fit_crr(des_ph, extra_X = t_cols)
s_ph  <- summary(fg_ph)$coef

ph_test <- bind_rows(lapply(EXPO_SD, function(e) {
  t_var <- sub("_sd$", "_t", e)
  b_m <- s_ph[e,     "coef"]
  b_t <- s_ph[t_var, "coef"];  se_t <- s_ph[t_var, "se(coef)"]
  p_t <- z_to_p(b_t, se_t)
  data.frame(variable   = e,
             SHR_at_1yr = exp(b_m),
             coef_logT  = b_t,
             se_logT    = se_t,
             P_interact = p_t,
             P_exact    = sprintf("%.3e", p_t),
             P_disp     = fmt_p(p_t),
             row.names  = NULL)
}))
cat(sprintf("\n===== 比例亚分布风险假设检验（%s 协变量集，n = %d）=====\n",
            PH_COVSET, des_ph$n))
print(ph_test, digits = 4, row.names = FALSE)

# ============================================================================
# 5. 输出
# ============================================================================
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
write.csv(sens2_table,  paste0(OUT_DIR, "TableS21_sens2_models.csv"),  row.names = FALSE)
write.csv(compare_m2m3, paste0(OUT_DIR, "sens2_m2_vs_m3.csv"),         row.names = FALSE)
write.csv(ph_test,      paste0(OUT_DIR, "sens2_ph_test.csv"),          row.names = FALSE)
write.csv(ledger,       paste0(OUT_DIR, "sens2_n_ledger.csv"),         row.names = FALSE)
write.csv(as.data.frame(sd_cmp), paste0(OUT_DIR, "sens2_sd_compare.csv"))
if (!is.null(des_m2$dummy_map))
  write.csv(des_m2$dummy_map, paste0(OUT_DIR, "sens2_dummy_map.csv"), row.names = FALSE)

cat("\n结果已保存至：", OUT_DIR, "\n")