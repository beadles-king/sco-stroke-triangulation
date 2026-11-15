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
# 敏感性分析：肌酐法 eGFR vs CKD-EPI 2021 肌酐+胱抑素 C 联合方程 eGFR_C
# —— 修订版（对应 Table S22）
#
# 相对原 30_charls_cystatin_c.R 的四处修改：
#
#  (1) 【A6】所有 n 与事件数改为报告"进模型人数"（成列删除之后）。
#      原脚本 run_fg_by_age() 用 nrow(sub) 计数，是删除前的行数，
#      导致 S22 的 ≥70 岁 1207/94 大于 S19 的 1178/92，看似不可能。
#      修订后两表口径一致，并额外输出删除前后的对照供表注引用。
#
#  (2) 【F1】marriage / education 经 build_design() 转哑变量，
#      不再作为单一线性项进入设计矩阵。
#
#  (3) 【F2】Models 1/2/3 各自报告分析人数——Model 3 因追加三个药物
#      变量而人数更少，此前未报告。
#
#  (4) 年龄分层同时给出两种设定：
#      single = 单暴露（原脚本的实际设定，产生已发表的 1·24 / 0·96）
#      mutual = 五暴露互相调整（与 Table S19 的年龄层设定一致，可比）
#      并可选用全样本 SD 标准化，使之与 S19 完全同尺度。
#
# 依赖：21_charls_model_utils.R（同目录）
# ============================================================================

setwd(RAWDAT)

library(readr)
library(dplyr)
library(tidyr)
library(survival)
library(cmprsk)

source(file.path(CODE_DIR, "21_21_charls_model_utils.R"))

OUT_DIR  <- file.path(RAWDAT, "CHARLS")
# "subsample" = 子样本内标准化（复现已发表结果）
# "whole"     = 全队列 SD 标准化（与 Table S19 同尺度）
SD_BASIS <- "subsample"

# ============================================================================
# 1. 数据与子样本
# ============================================================================
charls <- read.csv("charls_dat.csv")
charls$female <- as.integer(charls$gender == 2)

# 全队列 SD（备用基准）
whole_sd <- vapply(c(EXPO_RAW, "eGFR_C"),
                   function(v) sd(charls[[v]], na.rm = TRUE), numeric(1))
# 无胱抑素C测量
charls_sens <- charls %>%
  filter(!is.na(eGFR), !is.na(eGFR_C)) 

# 胱抑素C为不可能值
charls_sens <- charls_sens %>%
  filter(eGFR_C >= 15 & eGFR_C <= 200)

cat(sprintf("\n胱抑素 C 子样本：n = %d（全队列 %d 的 %.1f%%）\n",
            nrow(charls_sens), nrow(charls),
            nrow(charls_sens) / nrow(charls) * 100))

std <- function(x, v) {
  if (SD_BASIS == "whole") (x - mean(charls[[v]], na.rm = TRUE)) / whole_sd[[v]]
  else                     as.numeric(scale(x))
}
charls_sens <- charls_sens %>%
  mutate(SBP_sd    = std(SBP,    "SBP"),
         BMI_sd    = std(BMI,    "BMI"),
         FG_sd     = std(FG,     "FG"),
         LDL_sd    = std(LDL,    "LDL"),
         eGFR_sd   = std(eGFR,   "eGFR"),
         eGFR_C_sd = std(eGFR_C, "eGFR_C"),
         age_grp   = factor(if_else(age >= 70, "Age >= 70", "Age < 70"),
                            levels = c("Age < 70", "Age >= 70")))

cat(sprintf("标准化基准：%s\n", SD_BASIS))

# ============================================================================
# 2. 分布层面对比（不涉及协变量，仍用全子样本 n = 6109）
# ============================================================================
desc_compare <- charls_sens %>%
  summarise(n = n(),
            eGFR_mean = mean(eGFR),     eGFR_sd_raw   = sd(eGFR),
            eGFR_C_mean = mean(eGFR_C), eGFR_C_sd_raw = sd(eGFR_C),
            diff_mean = mean(eGFR - eGFR_C), diff_sd = sd(eGFR - eGFR_C),
            pearson_r = cor(eGFR, eGFR_C),
            spearman_r = cor(eGFR, eGFR_C, method = "spearman"))
cat("\n===== eGFR vs eGFR_C 分布对比 =====\n"); print(desc_compare)

ttest_res  <- t.test(charls_sens$eGFR, charls_sens$eGFR_C, paired = TRUE)
wilcox_res <- wilcox.test(charls_sens$eGFR, charls_sens$eGFR_C, paired = TRUE)
pooled_sd  <- sqrt((sd(charls_sens$eGFR)^2 + sd(charls_sens$eGFR_C)^2) / 2)

paired_test <- data.frame(
  comparison    = "eGFR vs eGFR_C (paired)",
  n             = nrow(charls_sens),
  mean_diff     = unname(ttest_res$estimate),
  CI_low        = ttest_res$conf.int[1],
  CI_high       = ttest_res$conf.int[2],
  t_pvalue      = ttest_res$p.value,
  wilcox_pvalue = wilcox_res$p.value,
  SMD           = (mean(charls_sens$eGFR) - mean(charls_sens$eGFR_C)) / pooled_sd,
  row.names = NULL
)
cat("\n===== 配对差异与 SMD =====\n"); print(paired_test)

ba_summary <- charls_sens %>%
  transmute(age_grp, diff_egfr = eGFR - eGFR_C) %>%
  group_by(age_grp) %>%
  summarise(n = n(), mean_diff = mean(diff_egfr), sd_diff = sd(diff_egfr),
            LoA_low  = mean_diff - 1.96 * sd_diff,
            LoA_high = mean_diff + 1.96 * sd_diff, .groups = "drop")
cat("\n===== Bland-Altman（按年龄层）=====\n"); print(ba_summary)

# ============================================================================
# 3. Models 1–3，两套 eGFR 定义，各自报告分析人数
# ============================================================================
run_one <- function(expo_focus, model, label) {
  covs <- if (model == "M3") COV_M3 else COV_M2
  expo <- if (model == "M1") expo_focus
  else c("SBP_sd", "FG_sd", "LDL_sd", expo_focus, "BMI_sd")
  des <- build_design(charls_sens, expo = expo, covs = covs, label = label)
  fit <- fit_crr(des)
  list(des = des,
       res = extract_fg(fit, expo_focus) %>%
         mutate(model = model, focus = expo_focus,
                n = des$n, n_stroke = des$n_stroke))
}

cat("\n===== 各模型分析人群 =====\n")
runs <- list(
  e_m1 = run_one("eGFR_sd",   "M1", "eGFR  M1"),
  e_m2 = run_one("eGFR_sd",   "M2", "eGFR  M2"),
  e_m3 = run_one("eGFR_sd",   "M3", "eGFR  M3"),
  c_m1 = run_one("eGFR_C_sd", "M1", "eGFR_C M1"),
  c_m2 = run_one("eGFR_C_sd", "M2", "eGFR_C M2"),
  c_m3 = run_one("eGFR_C_sd", "M3", "eGFR_C M3")
)

res_all <- bind_rows(lapply(runs, `[[`, "res")) %>%
  mutate(method = if_else(focus == "eGFR_sd",
                          "eGFR (Cr)", "eGFR_C (Cr+CysC)")) %>%
  select(model, method, n, n_stroke, SHR, CI_low, CI_high, P)
cat("\n===== Fine-Gray SHR（每 1 SD）=====\n"); print(res_all, digits = 3)

ledger_models <- do.call(n_ledger, lapply(runs, `[[`, "des"))
ledger_models$analysis <- names(runs)
cat("\n===== 分析人数记账（删除前 vs 进模型）=====\n")
print(ledger_models, row.names = FALSE)

# ============================================================================
# 4. 两套 eGFR 定义的模型间比较（Z 检验）
# ============================================================================
cor_egfr <- cor(charls_sens$eGFR_sd, charls_sens$eGFR_C_sd)
cat(sprintf("\neGFR_sd 与 eGFR_C_sd 的相关：r = %.4f\n", cor_egfr))

z_compare <- function(a, b, r = cor_egfr) {
  d  <- a$coef - b$coef
  se_n <- sqrt(a$se^2 + b$se^2)
  se_c <- sqrt(a$se^2 + b$se^2 - 2 * r * a$se * b$se)
  data.frame(diff_coef = d, SHR_ratio = exp(d),
             z_naive = d / se_n, p_naive = 2 * pnorm(-abs(d / se_n)),
             z_corrected = d / se_c, p_corrected = 2 * pnorm(-abs(d / se_c)),
             row.names = NULL)
}

compare_models <- bind_rows(
  z_compare(runs$e_m1$res, runs$c_m1$res) %>% mutate(model = "Model 1"),
  z_compare(runs$e_m2$res, runs$c_m2$res) %>% mutate(model = "Model 2"),
  z_compare(runs$e_m3$res, runs$c_m3$res) %>% mutate(model = "Model 3")
) %>% select(model, diff_coef, SHR_ratio, z_corrected, p_corrected,
             z_naive, p_naive)
cat("\n===== 两法比较（p_corrected 为推荐值）=====\n")
print(compare_models, digits = 3)

# ============================================================================
# 5. 年龄分层：两种设定并列，n 为进模型人数
#    single = 单暴露（原脚本设定，对应已发表的 1·24 / 0·96）
#    mutual = 五暴露互相调整（与 Table S19 年龄层同设定，可直接比较）
# ============================================================================
age_strat_one <- function(expo_focus, spec) {
  bind_rows(lapply(levels(charls_sens$age_grp), function(lv) {
    sub  <- charls_sens[charls_sens$age_grp == lv, ]
    expo <- if (spec == "single") expo_focus
    else c("SBP_sd", "FG_sd", "LDL_sd", expo_focus, "BMI_sd")
    des  <- build_design(sub, expo = expo, covs = COV_M2,
                         label = sprintf("%s %s %s", expo_focus, spec, lv))
    if (des$n_stroke < 20) {
      return(data.frame(spec = spec, focus = expo_focus, group = lv,
                        n = des$n, n_stroke = des$n_stroke,
                        SHR = NA, CI_low = NA, CI_high = NA, P = NA,
                        row.names = NULL))
    }
    extract_fg(fit_crr(des), expo_focus) %>%
      mutate(spec = spec, focus = expo_focus, group = lv,
             n = des$n, n_stroke = des$n_stroke, n_before = des$n_before,
             ev_before = des$ev_before) %>%
      select(spec, focus, group, n_before, ev_before, n, n_stroke,
             SHR, CI_low, CI_high, P)
  }))
}

cat("\n===== 年龄分层敏感性 =====\n")
age_strat <- bind_rows(
  age_strat_one("eGFR_sd",   "single"),
  age_strat_one("eGFR_C_sd", "single"),
  age_strat_one("eGFR_sd",   "mutual"),
  age_strat_one("eGFR_C_sd", "mutual")
) %>%
  mutate(method = if_else(focus == "eGFR_sd",
                          "eGFR (Cr)", "eGFR_C (Cr+CysC)")) %>%
  select(spec, method, group, n_before, ev_before, n, n_stroke,
         SHR, CI_low, CI_high, P)
print(age_strat, digits = 3, row.names = FALSE)

cat("\n提示：Table S22 的年龄层请填 n / n_stroke（进模型口径）；\n")
cat("      n_before / ev_before 仅用于说明与旧版数字的差异来源。\n")
cat("      若要与 Table S19 直接对照，取 spec = \"mutual\" 且\n")
cat("      将脚本顶部 SD_BASIS 设为 \"whole\"。\n")

# ============================================================================
# 6. 输出
# ============================================================================
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
write.csv(desc_compare,   paste0(OUT_DIR, "Sens_eGFR_desc.csv"),         row.names = FALSE)
write.csv(paired_test,    paste0(OUT_DIR, "Sens_eGFR_paired_SMD.csv"),   row.names = FALSE)
write.csv(ba_summary,     paste0(OUT_DIR, "Sens_eGFR_BlandAltman.csv"),  row.names = FALSE)
write.csv(res_all,        paste0(OUT_DIR, "Sens_eGFR_FineGray.csv"),     row.names = FALSE)
write.csv(compare_models, paste0(OUT_DIR, "Sens_eGFR_modelcompare.csv"), row.names = FALSE)
write.csv(age_strat,      paste0(OUT_DIR, "Sens_eGFR_agestrat.csv"),     row.names = FALSE)
write.csv(ledger_models,  paste0(OUT_DIR, "Sens_eGFR_n_ledger.csv"),     row.names = FALSE)

cat("\n结果已保存至：", OUT_DIR, "\n")