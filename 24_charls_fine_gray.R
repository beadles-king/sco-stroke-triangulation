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
# CHARLS 主分析：Fine-Gray 竞争风险模型（Models 1–3）+ 比例亚分布风险检验
# —— 修订版（对应 Table S15、Figure 4、Table 3）
#
# 相对旧版的两处实质修改：
#
#   【F1】marriage / education 改为哑变量。旧版把 marriage（1=已婚/同居,
#         2=分居/离异, 3=丧偶, 4=未婚，纯名义编码）和 education（1/2/3）
#         各以一个数值列进入 crr(cov1=)，等于假设相邻编码间的对数亚分布
#         风险增量相等——对婚姻状况没有依据。现经 build_design() 以
#         model.matrix 展开，参考水平为编码最小值（marriage=1 已婚/同居，
#         education=1 初中以下），与 Table S4 的分类呈现一致。
#
#   【F2】Models 1–2 与 Model 3 的分析人群不同。Model 3 追加 bpmed /
#         diamed / lipidmed 三个协变量（Table 2 中分别缺失 40 / 76 / 165
#         例），成列删除后人数少于 Models 1–2。旧版三个模型并排呈现却不
#         报告任何 n，正文把 M2→M3 的衰减完全归因于药物调整，其中有一部
#         分实为换了人群。现每个模型各自报告 n 与事件数，并输出记账表。
#
# 注：*_sd 仍在【全队列 8012 人】上 scale()，与旧版一致，因此"每 1 SD"
#     的绝对尺度不变；成列删除只影响进模型的人，不影响标准化基准。
#
# 依赖：21_charls_model_utils.R
# ============================================================================

setwd(RAWDAT)

library(readr)
library(tidyverse)
library(survival)
library(cmprsk)
library(rms)
library(ggplot2)

source(file.path(CODE_DIR, "21_21_charls_model_utils.R"))

OUT_DIR <- file.path(RAWDAT, "CHARLS")

# ============================================================================
# 1. 数据准备
# ============================================================================
charls <- read.csv('charls_dat.csv')

charls <- charls %>%
  mutate(
    # 标准化 5 个连续暴露（全队列基准，Fine-Gray 用）
    SBP_sd  = scale(SBP)[, 1],
    BMI_sd  = scale(BMI)[, 1],
    FG_sd   = scale(FG)[, 1],
    eGFR_sd = scale(eGFR)[, 1],
    LDL_sd  = scale(LDL)[, 1],
    # 二分类暴露（CIF 分层用，基于临床 cutoff）
    sbp_high = as.integer(SBP  >= 140),
    fg_high  = as.integer(FG   >= 126),
    ldl_high = as.integer(LDL  >= 160),
    bmi_high = as.integer(BMI  >= 28),      # 中国肥胖标准
    egfr_low = as.integer(eGFR <  60),
    # 协变量整理
    female   = as.integer(gender == 2),
    status_f = factor(new_stroke, levels = c(0, 1, 2))
  )

cat(sprintf("\n原始队列：n = %d，stroke = %d，competing death = %d\n",
            nrow(charls), sum(charls$new_stroke == 1),
            sum(charls$new_stroke == 2)))
cat("各协变量缺失数：\n")
print(sapply(charls[, COV_M3], function(x) sum(is.na(x))))

# ============================================================================
# 2. Model 1：各暴露逐个纳入（+ 人口学协变量）
# ============================================================================
cat("\n===== Model 1（逐个暴露）=====\n")
m1_runs <- lapply(EXPO_SD, function(e) {
  des <- build_design(charls, expo = e, covs = COV_M2, label = paste("M1", e))
  list(des = des,
       res = extract_fg(fit_crr(des), e) %>%
         mutate(model = "Model 1", n = des$n, n_stroke = des$n_stroke))
})
names(m1_runs) <- EXPO_SD
table2 <- bind_rows(lapply(m1_runs, `[[`, "res"))

# ============================================================================
# 3. Model 2：5 暴露同时纳入（主模型）
#    Model 3：5 暴露 + 人口学 + 三种药物
# ============================================================================
cat("\n===== Model 2 / Model 3（5 暴露互相调整）=====\n")
des_m2 <- build_design(charls, expo = EXPO_SD, covs = COV_M2, label = "Model 2")
des_m3 <- build_design(charls, expo = EXPO_SD, covs = COV_M3, label = "Model 3")

fg_m2 <- fit_crr(des_m2)
fg_m3 <- fit_crr(des_m3)

table3 <- bind_rows(lapply(EXPO_SD, function(e) extract_fg(fg_m2, e))) %>%
  mutate(model = "Model 2", n = des_m2$n, n_stroke = des_m2$n_stroke)
table4 <- bind_rows(lapply(EXPO_SD, function(e) extract_fg(fg_m3, e))) %>%
  mutate(model = "Model 3", n = des_m3$n, n_stroke = des_m3$n_stroke)

main_table <- bind_rows(table2, table3, table4) %>%
  select(model, variable, n, n_stroke, SHR, CI_low, CI_high, P)

cat("\n===== Models 1–3 汇总（Table S15 / Figure 4）=====\n")
print(main_table, digits = 3, row.names = FALSE)

# 哑变量参考水平说明（写进附录方法段）
cat("\n===== 分类协变量哑变量映射 =====\n")
print(des_m2$dummy_map, row.names = FALSE)

# ============================================================================
# 4. 分析人数记账（F2 的直接产出，供 Table S15 表注引用）
# ============================================================================
ledger <- n_ledger(`Model 1 (SBP)` = m1_runs$SBP_sd$des,
                   `Model 2`       = des_m2,
                   `Model 3`       = des_m3)
cat("\n===== 分析人数记账 =====\n")
print(ledger, row.names = FALSE)
cat(sprintf("\nModel 2 与 Model 3 的人数差：%d 人 / %d 事件\n",
            des_m2$n - des_m3$n, des_m2$n_stroke - des_m3$n_stroke))
cat("→ 正文若称 Model 3 的衰减源于药物调整，需说明其中有多少来自人群变化。\n")

# ============================================================================
# 5. Model 2 vs Model 3 的 SHR 对比
#    注意：两者人群不同，此对比是描述性的，不能解读为纯粹的药物调整效应
# ============================================================================
compare_m2m3 <- table3 %>%
  select(variable, SHR_m2 = SHR, CI_low_m2 = CI_low,
         CI_high_m2 = CI_high, P_m2 = P) %>%
  left_join(table4 %>% select(variable, SHR_m3 = SHR, CI_low_m3 = CI_low,
                              CI_high_m3 = CI_high, P_m3 = P),
            by = "variable") %>%
  mutate(SHR_change_pct = (SHR_m3 - SHR_m2) / SHR_m2 * 100,
         n_m2 = des_m2$n, n_m3 = des_m3$n)
cat("\n===== Model 2 vs Model 3 =====\n")
print(compare_m2m3, digits = 3, row.names = FALSE)

# --- 敏感性：把 Model 2 限制在 Model 3 的人群上重跑，分离"换人群"与"调药物" ---
charls_m3set <- des_m3$d
des_m2_on_m3 <- build_design(charls_m3set, expo = EXPO_SD, covs = COV_M2,
                             label = "Model 2 (restricted to M3 sample)")
fg_m2_on_m3  <- fit_crr(des_m2_on_m3)
table3r <- bind_rows(lapply(EXPO_SD, function(e) extract_fg(fg_m2_on_m3, e))) %>%
  mutate(model = "Model 2 (M3 sample)",
         n = des_m2_on_m3$n, n_stroke = des_m2_on_m3$n_stroke)
cat("\n===== Model 2 限制在 Model 3 人群（分离人群变化与药物调整）=====\n")
print(table3r %>% select(variable, n, n_stroke, SHR, CI_low, CI_high, P),
      digits = 3, row.names = FALSE)
cat("→ 与上面的 Model 2 相比，差异即为\"换人群\"部分；\n")
cat("  与 Model 3 相比，差异才是真正的\"药物调整\"部分。\n")

# ============================================================================
# 6. 比例亚分布风险假设检验（暴露 × log 时间）
#    设计矩阵沿用 Model 2（含哑变量），再追加 5 个 *_t 交互列
# ============================================================================
d_ph     <- des_m2$d
log_time <- log(d_ph$time_yr)
stopifnot(all(is.finite(log_time)))

t_cols <- sapply(EXPO_SD, function(e) d_ph[[e]] * log_time)
colnames(t_cols) <- sub("_sd$", "_t", EXPO_SD)

fg_ph <- fit_crr(des_m2, extra_X = t_cols)
s_ph  <- summary(fg_ph)$coef

ph_test <- bind_rows(lapply(EXPO_SD, function(e) {
  t_var <- sub("_sd$", "_t", e)
  b_m <- s_ph[e,     "coef"];  se_m <- s_ph[e,     "se(coef)"]
  b_t <- s_ph[t_var, "coef"];  se_t <- s_ph[t_var, "se(coef)"]
  p_t <- z_to_p(b_t, se_t)                    # ← 全精度，取代 s_ph[t_var, "p-value"]
  data.frame(variable   = e,
             SHR_at_1yr = exp(b_m),
             coef_logT  = b_t,
             se_logT    = se_t,
             P_interact = p_t,
             P_exact    = sprintf("%.3e", p_t),
             P_disp     = fmt_p(p_t),
             row.names  = NULL)
}))
cat("\n===== 比例亚分布风险假设检验 =====\n")
cat(sprintf("（n = %d，事件 %d）\n", des_m2$n, des_m2$n_stroke))
print(ph_test, digits = 4, row.names = FALSE)

# 时间特异性 SHR（Table S17）：系数本就是全精度，此段无需改动
timepoints <- 1:7
time_shr <- bind_rows(lapply(EXPO_SD, function(e) {
  t_var <- sub("_sd$", "_t", e)
  b1 <- s_ph[e, "coef"]; b2 <- s_ph[t_var, "coef"]
  data.frame(variable = e, time_yr = timepoints,
             SHR = exp(b1 + b2 * log(timepoints)), row.names = NULL)
}))
cat("\n===== 时间特异性 SHR（Table S17）=====\n")
print(time_shr %>%
        pivot_wider(names_from = time_yr, values_from = SHR,
                    names_prefix = "yr"),
      digits = 3, row.names = FALSE)

# ============================================================================
# 7. 输出
# ============================================================================
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
write.csv(main_table,   paste0(OUT_DIR, "FG_models1to3.csv"),      row.names = FALSE)
write.csv(compare_m2m3, paste0(OUT_DIR, "FG_m2_vs_m3.csv"),        row.names = FALSE)
write.csv(table3r,      paste0(OUT_DIR, "FG_m2_on_m3sample.csv"),  row.names = FALSE)
write.csv(ph_test,      paste0(OUT_DIR, "FG_ph_test.csv"),         row.names = FALSE)
write.csv(time_shr,     paste0(OUT_DIR, "FG_time_specific_SHR.csv"), row.names = FALSE)
write.csv(ledger,       paste0(OUT_DIR, "FG_n_ledger.csv"),        row.names = FALSE)
if (!is.null(des_m2$dummy_map))
  write.csv(des_m2$dummy_map, paste0(OUT_DIR, "FG_dummy_map.csv"), row.names = FALSE)

cat("\n结果已保存至：", OUT_DIR, "\n")