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
# CHARLS 亚组分析 —— 全样本 SD 标准化版（对应 Table S19）
#
# 合并来源：subgroup_wholeSD.R（分析口径）+ CHARLS_亚组.R（分离检测与导出结构）
#
# 依赖：21_charls_model_utils.R
#       【不依赖 24_charls_fine_gray.R】本脚本自行读取 charls_dat.csv 并生成 female 与
#       SBP_sd…BMI_sd 五列，生成方式与 24_charls_fine_gray.R 中的 scale() 完全一致
#       （全队列基准），因此与主分析同尺度、可复现。
#       若环境中已有 24_charls_fine_gray.R 建好的 charls，可设 REUSE_EXISTING <- TRUE。
#
# 五处实质修正：
#   (1) 【A6】n 与事件数在成列删除【之后】计数，与真正进模型的人群一致。
#       旧 CHARLS_亚组.R 在删除之前计数（≥70 岁报 1212/95，实际拟合 1178/92），
#       这也是 Table S22 的 1207/94 看似大于 S19 的 1178/92 的根源。
#   (2) 【A1】交互检验显式纳入分层变量主效应。旧脚本在年龄层只放交互项、
#       不放 age_group 主效应，属误设模型，其输出即正文那个 p=0·004；
#       正确设定给出 p=0·005，与 Table S19 一致。
#   (3) 协变量集固定为 Model 2（不含药物），与 Table S19 表注一致。
#       旧 CHARLS_亚组.R 实际用的是 Model 3 的 cov_vars_m3。
#   (4) 【F1】marriage / education 经 build_design() 转哑变量。
#   (5) 【F2】每个模型报告自身的分析人数；MODEL 可切 M2 / M3。
#
# 输出：
#   subgroup_wholeSD.csv           各层 SHR、层内原始 SD、分析 n
#   interaction_wholeSD.csv        交互检验
#   subgroup_n_reconciliation.csv  两套 n 口径对照（供 S19 / S22 表注引用）
# ============================================================================

library(cmprsk)
library(dplyr)
library(tidyr)

source(file.path(CODE_DIR, "21_21_charls_model_utils.R"))

# ============================================================================
# 0. 配置
# ============================================================================
DATA_PATH      <- file.path(RAWDAT, "charls_dat.csv")
OUT_DIR        <- file.path(RAWDAT, "CHARLS")
REUSE_EXISTING <- FALSE   # TRUE = 复用环境中已有的 charls（需含 *_sd 与 female）
MODEL          <- "M2"    # "M2" = Table S19 设定；"M3" = 追加药物调整

cov_vars <- if (MODEL == "M3") COV_M3 else COV_M2

# ============================================================================
# 1. 数据准备（与 24_charls_fine_gray.R 的 scale() 定义一致）
# ============================================================================
if (REUSE_EXISTING) {
  stopifnot(exists("charls"))
  miss <- setdiff(c(EXPO_SD, "female"), names(charls))
  if (length(miss) > 0)
    stop("环境中的 charls 缺少：", paste(miss, collapse = ", "),
         "；请设 REUSE_EXISTING <- FALSE")
} else {
  charls <- read.csv(DATA_PATH)
  charls$female <- as.integer(charls$gender == 2)
  for (i in seq_along(EXPO_RAW)) {
    charls[[EXPO_SD[i]]] <- as.numeric(scale(charls[[EXPO_RAW[i]]]))
  }
}

charls <- charls %>%
  mutate(
    age_group = factor(if_else(age >= 70, "Age >= 70", "Age < 70"),
                       levels = c("Age < 70", "Age >= 70")),
    sex_group = factor(if_else(female == 1, "Female", "Male"),
                       levels = c("Male", "Female")),
    res_group = factor(if_else(rural == 1, "Rural", "Urban"),
                       levels = c("Urban", "Rural"))
  )

subgroups <- list(age = "age_group", sex = "sex_group", rural = "res_group")
# 分层变量与协变量中的同源项（分层/交互时须剔除，避免共线）
drop_map  <- list(age_group = "age", sex_group = "female", res_group = "rural")

# ============================================================================
# 2. 统一分析人群：对完整协变量集做成列删除
#    所有层内拟合与交互检验共用同一分析集，保证 n 与交互 p 口径一致
# ============================================================================
analysis_ok <- stats::complete.cases(
  charls[, c("time_yr", "new_stroke", EXPO_SD, cov_vars,
             unlist(subgroups))]
)
charls_ana <- charls[analysis_ok, ]

cat(sprintf("\n===== 分析人群（%s 协变量集）=====\n", MODEL))
cat(sprintf("  原始行数              : n = %d, stroke = %d, death = %d\n",
            nrow(charls), sum(charls$new_stroke == 1),
            sum(charls$new_stroke == 2)))
cat(sprintf("  成列删除后（进模型）  : n = %d, stroke = %d, death = %d\n",
            nrow(charls_ana), sum(charls_ana$new_stroke == 1),
            sum(charls_ana$new_stroke == 2)))
cat("  各协变量缺失数：\n")
print(sapply(charls[, cov_vars], function(x) sum(is.na(x))))

# ============================================================================
# 3. 两套 n 口径对照（供 Table S19 / S22 表注引用）
# ============================================================================
recon <- bind_rows(lapply(names(subgroups), function(sg) {
  sg_col <- subgroups[[sg]]
  bind_rows(lapply(levels(charls[[sg_col]]), function(lv) {
    pre  <- charls[!is.na(charls[[sg_col]])         & charls[[sg_col]] == lv, ]
    post <- charls_ana[!is.na(charls_ana[[sg_col]]) & charls_ana[[sg_col]] == lv, ]
    data.frame(subgroup = sg_col, level = lv,
               n_before_del    = nrow(pre),
               events_before   = sum(pre$new_stroke == 1),
               n_analysis      = nrow(post),
               events_analysis = sum(post$new_stroke == 1),
               n_dropped       = nrow(pre) - nrow(post),
               row.names = NULL)
  }))
}))
cat("\n===== n 口径对照（未剔缺失 vs 进模型）=====\n")
print(recon, row.names = FALSE)

# ============================================================================
# 4. 层内拟合：全样本 SD，5 暴露互相调整
#    build_design() 负责哑变量、零方差协变量剔除、完全分离检测
# ============================================================================
fit_subgroup <- function(sg_col, level) {
  d   <- charls_ana[!is.na(charls_ana[[sg_col]]) & charls_ana[[sg_col]] == level, ]
  des <- build_design(d, expo = EXPO_SD, covs = cov_vars,
                      label = sprintf("%s = %s", sg_col, level))
  fit <- fit_crr(des)
  s   <- summary(fit)$coef
  
  b  <- s[EXPO_SD, "coef"]
  se <- s[EXPO_SD, "se(coef)"]
  p  <- z_to_p(b, se)                        # ← 取代 s[EXPO_SD, "p-value"]
  
  data.frame(
    subgroup = sg_col, level = level,
    n_before = des$n_before, events_before = des$ev_before,
    n_total  = des$n, n_stroke = des$n_stroke, n_death = des$n_death,
    variable = EXPO_SD,
    SHR      = exp(b),
    CI_low   = exp(b - 1.96 * se),
    CI_high  = exp(b + 1.96 * se),
    P        = p,
    P_exact  = sprintf("%.3e", p),
    P_disp   = fmt_p(p),
    raw_sd   = vapply(EXPO_RAW, function(v) sd(des$d[[v]], na.rm = TRUE), numeric(1)),
    dropped_cov = paste(des$dropped_cov, collapse = ";"),
    dropped_col = paste(des$dropped_col, collapse = ";"),
    row.names   = NULL
  )
}

cat("\n===== 层内拟合 =====\n")
subgroup_res <- bind_rows(lapply(names(subgroups), function(sg) {
  sg_col <- subgroups[[sg]]
  bind_rows(lapply(levels(charls[[sg_col]]),
                   function(lv) fit_subgroup(sg_col, lv)))
}))

# 标准化基准核对
sd_basis <- rbind(
  `全队列(标准化基准)` = vapply(EXPO_RAW, function(v) sd(charls[[v]],     na.rm = TRUE), numeric(1)),
  `分析人群`           = vapply(EXPO_RAW, function(v) sd(charls_ana[[v]], na.rm = TRUE), numeric(1))
)
cat("\n===== 标准化基准 SD =====\n"); print(round(sd_basis, 3))

cat("\n===== 各层原始暴露 SD =====\n")
print(as.data.frame(subgroup_res %>%
                      distinct(subgroup, level, variable, raw_sd) %>%
                      pivot_wider(names_from = variable, values_from = raw_sd)),
      digits = 3, row.names = FALSE)

# ============================================================================
# 5. 交互检验：全样本单一模型，与层内点估计共用同一 SD 与同一分析人群
#    关键：分层变量主效应 grp 显式纳入，并剔除与其同源的协变量
# ============================================================================
interaction_test <- function(sg_col) {
  lv2     <- levels(charls_ana[[sg_col]])[2]
  cov_use <- setdiff(cov_vars, drop_map[[sg_col]])
  
  des <- build_design(charls_ana, expo = EXPO_SD, covs = cov_use,
                      extra_cols = sg_col, verbose = FALSE,
                      label = paste("interaction", sg_col))
  ind <- as.integer(des$d[[sg_col]] == lv2)
  
  bind_rows(lapply(EXPO_SD, function(e) {
    extra <- cbind(grp = ind, inter = des$d[[e]] * ind)
    fit   <- fit_crr(des, extra_X = extra)
    s     <- summary(fit)$coef
    b  <- s["inter", "coef"]
    se <- s["inter", "se(coef)"]
    p  <- z_to_p(b, se)                      # ← 取代 s["inter", "p-value"]
    data.frame(subgroup = sg_col, exposure = e,
               n_model       = des$n, events_model = des$n_stroke,
               interact_coef = b,
               interact_se   = se,
               P_interaction = p,
               P_exact       = sprintf("%.3e", p),
               P_disp        = fmt_p(p),
               row.names = NULL)
  }))
}

cat("\n===== 交互检验 =====\n")
interaction_res <- bind_rows(lapply(unlist(subgroups, use.names = FALSE),
                                    interaction_test))

# ============================================================================
# 6. 输出
# ============================================================================
cat(sprintf("\n===== 亚组点估计（全队列 SD，%s）=====\n", MODEL))
print(subgroup_res[, c("subgroup", "level", "n_total", "n_stroke",
                       "variable", "SHR", "CI_low", "CI_high", "P")],
      digits = 3, row.names = FALSE)

cat("\n===== 交互检验（含分层变量主效应）=====\n")
print(interaction_res, digits = 3, row.names = FALSE)

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
write.csv(subgroup_res,    file.path(OUT_DIR, "subgroup_wholeSD.csv"),           row.names = FALSE)
write.csv(interaction_res, file.path(OUT_DIR, "interaction_wholeSD.csv"),        row.names = FALSE)
write.csv(recon,           file.path(OUT_DIR, "subgroup_n_reconciliation.csv"),  row.names = FALSE)

cat("\n已写出 subgroup_wholeSD.csv / interaction_wholeSD.csv /",
    "subgroup_n_reconciliation.csv\n")
cat("提示：Table S19 的 n 与事件数请取 n_total / n_stroke（进模型口径）；\n")
cat("      n_before / events_before 仅用于说明与旧数字的差异来源。\n")