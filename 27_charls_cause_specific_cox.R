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
# 敏感性分析 1：Cause-specific Cox 回归（Table S20）
# 竞争死亡 (new_stroke == 2) 重编码为删失
#
# 两处修改：
#  【F1】marriage / education 改为 factor()，与主分析的哑变量设定一致。
#        旧版在公式里直接写 marriage + education，被当成线性项拟合。
#  【F2】每个模型报告自身的分析人数与事件数（coxph 的 n / nevent），
#        Model 3 因追加三个药物协变量而人群更小。
#
# 注：本脚本用公式接口，不需要 build_design()，但 factor() 的参考水平
#     与 21_charls_model_utils.R 一致（编码最小值：marriage=1 已婚/同居，
#     education=1 初中以下），保证与 Fine-Gray 主分析可对照。
# ============================================================================

setwd(RAWDAT)
library(tidyverse)
library(survival)

OUT_DIR <- file.path(RAWDAT, "CHARLS")

charls <- read.csv('charls_dat.csv')

charls <- charls %>%
  mutate(
    female   = as.integer(gender == 2),
    SBP_sd   = as.numeric(scale(SBP)),
    BMI_sd   = as.numeric(scale(BMI)),
    FG_sd    = as.numeric(scale(FG)),
    eGFR_sd  = as.numeric(scale(eGFR)),
    LDL_sd   = as.numeric(scale(LDL)),
    # 竞争死亡 → 删失
    event_cs = ifelse(new_stroke == 1, 1, 0),
    # 【F1】分类协变量转因子，参考水平 = 编码最小值
    marriage  = factor(marriage,  levels = sort(unique(marriage))),
    education = factor(education, levels = sort(unique(education)))
  )

cat("\n===== marriage / education 水平分布 =====\n")
print(table(charls$marriage,  useNA = "ifany"))
print(table(charls$education, useNA = "ifany"))

expo_sd     <- c("SBP_sd", "FG_sd", "LDL_sd", "eGFR_sd", "BMI_sd")
cov_form    <- "age + female + marriage + education + rural + smoke + drink + heartd + depression"
cov_form_m3 <- paste(cov_form, "+ bpmed + diamed + lipidmed")

# ---- 提取函数：同时带回该模型的 n 与事件数 --------------------------------
extract_cs <- function(fit, vars, model_label) {
  s <- summary(fit)
  bind_rows(lapply(vars, function(v) {
    idx <- which(rownames(s$coefficients) == v)
    data.frame(model    = model_label,
               variable = v,
               n        = fit$n,
               n_event  = fit$nevent,
               HR       = s$conf.int[idx, 1],
               CI_low   = s$conf.int[idx, 3],
               CI_high  = s$conf.int[idx, 4],
               P        = s$coefficients[idx, 5],
               row.names = NULL)
  }))
}

# ============================================================================
# Model 1：各暴露逐个纳入
# ============================================================================
cs_model1 <- bind_rows(lapply(expo_sd, function(e) {
  f   <- as.formula(paste0("Surv(time_yr, event_cs) ~ ", e, " + ", cov_form))
  fit <- coxph(f, data = charls)
  extract_cs(fit, e, "Model 1")
}))
cat("\n===== Cause-specific Cox: Model 1 =====\n")
print(cs_model1, digits = 3, row.names = FALSE)

# ============================================================================
# Model 2：5 暴露同时纳入
# ============================================================================
f_m2   <- as.formula(paste0("Surv(time_yr, event_cs) ~ ",
                            paste(expo_sd, collapse = " + "), " + ", cov_form))
fit_m2 <- coxph(f_m2, data = charls)
cs_model2 <- extract_cs(fit_m2, expo_sd, "Model 2")
cat("\n===== Cause-specific Cox: Model 2 =====\n")
print(cs_model2, digits = 3, row.names = FALSE)

# ============================================================================
# Model 3：5 暴露 + 人口学 + 药物
# ============================================================================
f_m3   <- as.formula(paste0("Surv(time_yr, event_cs) ~ ",
                            paste(expo_sd, collapse = " + "), " + ", cov_form_m3))
fit_m3 <- coxph(f_m3, data = charls)
cs_model3 <- extract_cs(fit_m3, expo_sd, "Model 3")
cat("\n===== Cause-specific Cox: Model 3（+ 药物）=====\n")
print(cs_model3, digits = 3, row.names = FALSE)

cs_all <- bind_rows(cs_model1, cs_model2, cs_model3)

# ---- 【F2】分析人数记账 ----
ledger <- data.frame(
  model    = c("Model 1", "Model 2", "Model 3"),
  n        = c(cs_model1$n[1], fit_m2$n, fit_m3$n),
  n_event  = c(cs_model1$n_event[1], fit_m2$nevent, fit_m3$nevent),
  n_input  = nrow(charls),
  row.names = NULL
)
ledger$n_dropped <- ledger$n_input - ledger$n
cat("\n===== 分析人数记账 =====\n"); print(ledger, row.names = FALSE)

# ---- 药物系数（参考用）----
med_effect <- extract_cs(fit_m3, c("bpmed", "diamed", "lipidmed"), "Model 3")
cat("\n----- 药物协变量系数（参考）-----\n")
print(med_effect, digits = 3, row.names = FALSE)

# ---- 分类协变量系数（核查哑变量是否正确展开）----
cat("\n----- marriage / education 哑变量系数（Model 2）-----\n")
s2 <- summary(fit_m2)$coefficients
print(round(s2[grep("^marriage|^education", rownames(s2)), , drop = FALSE], 4))

# ============================================================================
# 比例风险假设检验（Schoenfeld 残差）
# ============================================================================
ph_test <- cox.zph(fit_m3)
cat("\n===== 比例风险假设检验（Model 3）=====\n")
print(ph_test)

# ============================================================================
# 导出
# ============================================================================
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
write.csv(cs_all,  paste0(OUT_DIR, "TableS20_cs_cox_models.csv"), row.names = FALSE)
write.csv(ledger,  paste0(OUT_DIR, "cs_cox_n_ledger.csv"),        row.names = FALSE)
write.csv(as.data.frame(ph_test$table),
          paste0(OUT_DIR, "cs_cox_ph_test.csv"))

cat("\n结果已保存至：", OUT_DIR, "\n")