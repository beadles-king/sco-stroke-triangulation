# =============================================================================
# Reproducibility repository - Wang et al., transportability-aware triangulation
# of metabolic stroke evidence across the Shanghai Cooperation Organisation.
# Paths are resolved by config.R; see README.md.
# =============================================================================
if (!exists("RAWDAT")) {
  if (file.exists("config.R")) source("config.R") else
    stop("Run this script from the code/ directory (the folder containing config.R).")
}

###############################################################################
# Figure S20：CHARLS 机器学习佐证分析（XGBoost + LASSO）
#
# 本脚本由 31_charls_shap.R 精简润色而来，只保留进入 Figure S20 的内容。
#
# 已删除的非发表内容：
#   混淆矩阵与 Youden 阈值、Brier score、准确度/敏感度/特异度/PPV/NPV/F1、
#   模型性能综合对比表、SHAP 蜂群图（附录第 275 段已写明 beeswarm 不展示）、
#   waterfall 个体解释图、SBP×BMI / SBP×FG 交互图、代谢变量排序对比图、
#   cv.glmnet 的 lambda 路径图、以及从未被使用的 shap_long_trimmed。
#
# 相对原脚本的修改：
#   【M1】配色统一。原脚本条形图 XGBoost=蓝、LASSO=红，而 ROC 图的
#         scale_color_manual 按因子字母序上色（LASSO 在前），实际把蓝给了
#         LASSO、红给了 XGBoost，与条形图相反。现固定 XGBoost=COL_XGB、
#         LASSO=COL_LAS，三个面板一致。
#   【M2】依赖图去掉着色映射。原脚本用 color_feature = "auto"，五个面板各自
#         挑选交互最强的变量上色，而共享色条只写 "Feature value"，读者无从
#         知道每个面板的颜色代表哪个变量。附录对本图的主张只涉及函数形状，
#         颜色无贡献，故改为灰点 + 单色平滑线。
#   【M3】两个模型共用同一套交叉验证折。原脚本调用了两次 createFolds，
#         XGBoost 与 LASSO 的 CV AUC 落在不同划分上，不可配对比较。
#   【M4】发表格式：Helvetica、base_size 8、中点小数号、变量名对齐
#         Table S14、坐标轴带单位、面板标题移除（面板字母在拼图阶段加）。
#   【M5】fmt_p 自检。本项目已两次因 fmt_p 版本不对导致图上 p 值错误，
#         故在出图前断言其行为，版本不对直接停。
#   【M6】marriage 与 education 改为哑变量，参考水平取编码最小值，与
#         v6 §3【F1】对 Fine-Gray 模型所做的修正一致。此前二者以单一
#         数值列进入 glmnet，对线性模型而言是把名义变量拟合成了线性项。
#   【M7】可复现性。此前"设了 set.seed 仍每次不同"的根源不在 R 的随机数：
#         (i) 未固定 nthread，OpenMP 线程数随运行时环境变化，直方图算法
#             多线程求和的浮点累加顺序随之改变，分裂点出现毫厘之差；
#         (ii) which.max 在平坦的 CV AUC 曲线上是刀刃判断，上述毫厘之差
#             足以让最优轮次跳变几十轮，模型、SHAP、AUC 随之全变。
#         现固定 NTHREAD、显式向 XGBoost 传 seed（subsample 与
#         colsample_bytree 用的是它自己的 RNG，不受 set.seed 管辖）、
#         每个随机步骤就近设种子、cv.glmnet 改传显式 foldid。
#         调参规则未改动：先修病因，不顺手换估计量。脚本会打印曲线平坦度
#         诊断，若修完仍不稳，可把 NROUNDS_RULE 切到 "1se"。
#
# 输出：FigS20_A_xgb_importance / FigS20_B_lasso_importance /
#       FigS20_C_roc / FigS20_D_dependence，各出 600 与 300 dpi 两档 TIF
###############################################################################

setwd(RAWDAT)

library(tidyverse)
library(xgboost)
library(SHAPforxgboost)
library(glmnet)
library(caret)
library(pROC)
library(patchwork)
library(ragg)

# ---- 格式化工具：与 21_charls_model_utils.R / 05_paf_attribution.R 中的定义必须一致 --------
if (file.exists(file.path(CODE_DIR, "21_21_charls_model_utils.R"))) source(file.path(CODE_DIR, "21_21_charls_model_utils.R"))
if (file.exists(file.path(CODE_DIR, "05_paf_attribution.R"))) source(file.path(CODE_DIR, "05_paf_attribution.R"))
stopifnot(exists("dec"), exists("fmt_p"))
setwd(RAWDAT)

# 【M5】fmt_p 自检：两档规则（<0·0001 / 0·0001–0·01 四位 / ≥0·01 两位）
stopifnot(identical(
  unname(sapply(c(0.00004, 0.00045, 0.00631, 0.27, 0.03394), fmt_p)),
  c("<0\u00b70001", "0\u00b70005", "0\u00b70063", "0\u00b727", "0\u00b703")
))
cat("fmt_p 自检通过\n")

# ---- 可复现性设定（见【M7】）-----------------------------------------------
SEED     <- 1111
NTHREAD  <- 1           # 固定线程数；不固定则浮点累加顺序随运行时环境变化
NROUNDS_RULE <- "max"   # "max" = which.max（原规则）；"1se" = 最优 1 个标准误内的最小轮次

# ---- 出图常数 --------------------------------------------------------------
COL_XGB  <- "#2C5985"   # XGBoost，全文主蓝
COL_LAS  <- "#A63603"   # LASSO，全文橙红
COL_PT   <- "#9E9E9E"   # 依赖图散点
FS_BASE  <- 8
OUT_DIR  <- file.path(RAWDAT, "CHARLS/ML")

save_panel <- function(p, stem, w_mm, h_mm) {
  for (r in c(600, 300)) {
    f <- file.path(OUT_DIR,
                   if (r == 600) paste0(stem, ".tif") else paste0(stem, "_300dpi.tif"))
    agg_tiff(f, width = w_mm, height = h_mm, units = "mm",
             res = r, compression = "lzw", background = "white")
    print(p); dev.off()
    cat(sprintf("已保存 %s (%.0f x %.0f mm @ %d dpi)\n", f, w_mm, h_mm, r))
  }
}

theme_pub <- function() {
  theme_minimal(base_size = FS_BASE, base_family = "Helvetica") +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "grey92", linewidth = 0.25),
      axis.text        = element_text(size = FS_BASE - 1, colour = "grey20"),
      axis.title       = element_text(size = FS_BASE),
      axis.line        = element_line(colour = "grey35", linewidth = 0.3),
      axis.ticks       = element_line(colour = "grey35", linewidth = 0.3),
      plot.title       = element_text(size = FS_BASE, face = "bold", hjust = 0),
      plot.margin      = margin(2, 6, 2, 2)
    )
}

# 坐标轴数字用中点小数号
lab_dec <- function(d) function(x) dec(sprintf(paste0("%.", d, "f"), x))

# ============================================================================
# 1. 数据准备
# ============================================================================
charls <- read.csv('charls_dat.csv')

feature_names <- c("SBP", "FG", "LDL", "eGFR", "BMI",
                   "age", "gender", "marriage", "education",
                   "rural", "smoke", "drink", "heart disease", "depression")

# 图上显示名，与 Table S14 一致
LABELS <- c(
  SBP = "SBP", FG = "FG", LDL = "LDL", eGFR = "eGFR", BMI = "BMI",
  age = "Age", gender = "Sex", marriage = "Marital status",
  education = "Education", rural = "Residence", smoke = "Smoking",
  drink = "Drinking", `heart disease` = "Heart disease",
  depression = "Depression"
)

ml_dat <- charls %>%
  rename(`heart disease` = heartd) %>%
  filter(new_stroke != 2) %>%                       # 竞争死亡不进入二分类
  mutate(stroke_binary = as.integer(new_stroke == 1)) %>%
  select(stroke_binary, all_of(feature_names)) %>%
  drop_na()

cat(sprintf("\nML 样本：n = %d，卒中 %d 例（%.2f%%）\n",
            nrow(ml_dat), sum(ml_dat$stroke_binary),
            mean(ml_dat$stroke_binary) * 100))
cat("注意：这一分析人群与 Models 1–2（7896）和 Model 3（7710）都不同，\n")
cat("      图注中必须写明本图的 n。\n")

# ---- 【M6】marriage 与 education 改为哑变量，参考水平 = 编码最小值 ----------
# 与 v6 §3【F1】对 Fine-Gray 模型所做的修正一致。
# marriage: 1 已婚/同居(参考)、2 分居/离异、3 丧偶、4 未婚
# education: 1 初中以下(参考)、2 高中/中专、3 大专及以上
# 标签按图面宽度取短形式；各水平相对的参考类别在图注中说明
DUMMY_SPEC <- list(
  marriage = c(`2` = "Separated/divorced",
               `3` = "Widowed",
               `4` = "Never married"),
  education = c(`2` = "Upper secondary",
                `3` = "Tertiary")
)

for (v in names(DUMMY_SPEC)) {
  obs <- sort(unique(ml_dat[[v]]))
  exp_lv <- as.numeric(c(min(obs), names(DUMMY_SPEC[[v]])))
  if (!identical(as.numeric(obs), exp_lv))
    stop(sprintf("%s 的水平是 {%s}，与预期的 {%s} 不符，请先核对编码。",
                 v, paste(obs, collapse = ","), paste(exp_lv, collapse = ",")))
  for (L in names(DUMMY_SPEC[[v]])) {
    nm <- paste0(v, "_", L)
    ml_dat[[nm]] <- as.integer(ml_dat[[v]] == as.numeric(L))
    LABELS[nm] <- DUMMY_SPEC[[v]][[L]]
    n_lv <- sum(ml_dat[[nm]])
    cat(sprintf("  哑变量 %-12s n = %5d (%.1f%%)  参考水平 = %s\n",
                nm, n_lv, n_lv / nrow(ml_dat) * 100, min(obs)))
    if (n_lv < 20)
      warning(sprintf("哑变量 %s 只有 %d 人，估计可能不稳。", nm, n_lv), call. = FALSE)
  }
  ml_dat[[v]] <- NULL
}

model_vars <- setdiff(names(ml_dat), "stroke_binary")
X <- as.matrix(ml_dat[, model_vars])
y <- ml_dat$stroke_binary
stopifnot(all(apply(X, 2, function(z) length(unique(z))) > 1))   # 无零方差列
cat(sprintf("进模型变量 %d 个（哑变量化前 %d 个）\n",
            ncol(X), length(feature_names)))

set.seed(SEED)
train_idx <- createDataPartition(y, p = 0.7, list = FALSE)
X_train <- X[train_idx, ];  y_train <- y[train_idx]
X_test  <- X[-train_idx, ]; y_test  <- y[-train_idx]
cat(sprintf("训练集 %d | 测试集 %d\n", nrow(X_train), nrow(X_test)))

# ============================================================================
# 2. XGBoost：5 折交叉验证调轮次 + 最终模型
# ============================================================================
# 【M7】nthread 与 seed 都写进参数：前者固定浮点累加顺序，后者控制
# subsample / colsample_bytree —— 它们用 XGBoost 自己的 RNG，不受 set.seed 管辖
XGB_PARAMS <- list(
  objective        = "binary:logistic",
  eval_metric      = "auc",
  max_depth        = 4,
  eta              = 0.05,
  subsample        = 0.8,
  colsample_bytree = 0.8,
  min_child_weight = 10,
  nthread          = NTHREAD,
  seed             = SEED
)

spw <- sum(y_train == 0) / sum(y_train == 1)
dtrain <- xgb.DMatrix(X_train, label = y_train, nthread = NTHREAD)
dtest  <- xgb.DMatrix(X_test,  label = y_test,  nthread = NTHREAD)

set.seed(SEED)
xgb_cv <- xgb.cv(params = c(XGB_PARAMS, list(scale_pos_weight = spw)),
                 data = dtrain, nrounds = 500, nfold = 5,
                 early_stopping_rounds = 30, stratified = TRUE, verbose = 0)

# ---- 轮次选择与曲线平坦度诊断 ----
elog  <- xgb_cv$evaluation_log
i_max <- which.max(elog$test_auc_mean)
se_at_max <- elog$test_auc_std[i_max] / sqrt(5)
within_1se <- which(elog$test_auc_mean >= elog$test_auc_mean[i_max] - se_at_max)
i_1se <- min(within_1se)

cat(sprintf("\nCV AUC 最优 %.5f 于第 %d 轮；1 个标准误内的轮次范围 %d-%d（共 %d 轮）\n",
            elog$test_auc_mean[i_max], i_max,
            min(within_1se), max(within_1se), length(within_1se)))
if (length(within_1se) > 50)
  cat("  曲线在最优点附近很平：which.max 对微小数值扰动敏感。若跨会话仍不稳，\n",
      "  可把 NROUNDS_RULE 切到 \"1se\"。\n", sep = "")

best_nrounds <- if (NROUNDS_RULE == "1se") i_1se else i_max
cat(sprintf("采用规则 \"%s\"，最优迭代轮次 = %d\n", NROUNDS_RULE, best_nrounds))

set.seed(SEED)
xgb_model <- xgb.train(params = c(XGB_PARAMS, list(scale_pos_weight = spw)),
                       data = dtrain, nrounds = best_nrounds, verbose = 0)

xgb_roc <- roc(y_test, predict(xgb_model, dtest), quiet = TRUE)
xgb_auc <- as.numeric(auc(xgb_roc))

# ============================================================================
# 3. LASSO
# ============================================================================
preproc   <- preProcess(X_train, method = c("center", "scale"))
X_train_s <- predict(preproc, X_train)
X_test_s  <- predict(preproc, X_test)

w_train <- ifelse(y_train == 1, sum(y_train == 0) / sum(y_train == 1), 1)

# 【M7】显式 foldid：cv.glmnet 内部抽折依赖调用时刻的随机数流位置，
# 传入固定的 foldid 后，结果与前面执行过多少随机操作无关
set.seed(SEED)
foldid_lasso <- sample(rep_len(1:10, length(y_train)))

lasso_cv <- cv.glmnet(X_train_s, y_train, family = "binomial", alpha = 1,
                      weights = w_train, foldid = foldid_lasso,
                      type.measure = "auc")
lasso_model <- glmnet(X_train_s, y_train, family = "binomial", alpha = 1,
                      weights = w_train, lambda = lasso_cv$lambda.min)

lasso_roc <- roc(y_test, predict(lasso_model, X_test_s, type = "response")[, 1],
                 quiet = TRUE)
lasso_auc <- as.numeric(auc(lasso_roc))

# ============================================================================
# 4. 【M3】两模型共用同一套 5 折，CV AUC 方可配对比较
# ============================================================================
set.seed(SEED)
cv_folds <- createFolds(y, k = 5, list = TRUE, returnTrain = FALSE)

cv_auc_pair <- vapply(seq_along(cv_folds), function(k) {
  test_i <- cv_folds[[k]]
  Xtr <- X[-test_i, ]; ytr <- y[-test_i]
  Xte <- X[ test_i, ]; yte <- y[ test_i]
  w   <- sum(ytr == 0) / sum(ytr == 1)
  
  set.seed(SEED + k)
  m_xgb <- xgb.train(params = c(XGB_PARAMS, list(scale_pos_weight = w)),
                     data = xgb.DMatrix(Xtr, label = ytr, nthread = NTHREAD),
                     nrounds = best_nrounds, verbose = 0)
  a_xgb <- as.numeric(auc(roc(yte,
                              predict(m_xgb, xgb.DMatrix(Xte, nthread = NTHREAD)), quiet = TRUE)))
  
  pp <- preProcess(Xtr, method = c("center", "scale"))
  set.seed(SEED + k)
  fid <- sample(rep_len(1:5, length(ytr)))
  inner <- cv.glmnet(predict(pp, Xtr), ytr, family = "binomial", alpha = 1,
                     weights = ifelse(ytr == 1, w, 1), foldid = fid,
                     type.measure = "auc")
  a_las <- as.numeric(auc(roc(yte,
                              predict(inner, predict(pp, Xte), s = inner$lambda.min,
                                      type = "response")[, 1], quiet = TRUE)))
  c(xgb = a_xgb, lasso = a_las)
}, numeric(2))

delong <- roc.test(xgb_roc, lasso_roc, method = "delong")

cat("\n===== 图注所需数值 =====\n")
cat(sprintf("XGBoost 测试集 AUC = %.4f  ->  %s\n", xgb_auc, dec(sprintf("%.2f", xgb_auc))))
cat(sprintf("LASSO   测试集 AUC = %.4f  ->  %s\n", lasso_auc, dec(sprintf("%.2f", lasso_auc))))
cat(sprintf("DeLong z = %.3f, p = %.4f  ->  p=%s\n",
            delong$statistic, delong$p.value, fmt_p(delong$p.value)))
cat(sprintf("5 折 CV AUC（同一套折）XGBoost %.3f ± %.3f | LASSO %.3f ± %.3f\n",
            mean(cv_auc_pair["xgb", ]), sd(cv_auc_pair["xgb", ]),
            mean(cv_auc_pair["lasso", ]), sd(cv_auc_pair["lasso", ])))
cat(sprintf("LASSO 非零系数 %d / %d 个变量\n",
            sum(coef(lasso_model)[-1] != 0), ncol(X)))

# ---- 跨会话可复现性指纹：换个会话重跑，这几行必须逐字相同 ----
cat("\n===== 可复现性指纹（跨会话比对用）=====\n")
cat(sprintf("nthread=%d  seed=%d  规则=%s  轮次=%d\n",
            NTHREAD, SEED, NROUNDS_RULE, best_nrounds))
cat(sprintf("AUC  xgb=%.10f  lasso=%.10f\n", xgb_auc, lasso_auc))
cat(sprintf("DeLong z=%.10f\n", as.numeric(delong$statistic)))
cat("R:", R.version.string, "| xgboost:",
    as.character(packageVersion("xgboost")), "\n")

# 与附录第 275 段核对
if (abs(round(lasso_auc, 2) - 0.66) > 1e-9 || abs(round(xgb_auc, 2) - 0.65) > 1e-9)
  warning("AUC 与附录第 275 段的「0·66 和 0·65」不符，请同步更新附录。", call. = FALSE)

# ============================================================================
# 5. SHAP
# ============================================================================
shap_values <- shap.values(xgb_model = xgb_model, X_train = X_train)
shap_long   <- shap.prep(shap_contrib = shap_values$shap_score, X_train = X_train)
cat(sprintf("指纹续：SHAP 绝对值总和 = %.10f\n", sum(abs(shap_values$shap_score))))

imp_xgb <- shap_values$shap_score %>% abs() %>% colMeans() %>%
  enframe("feature", "value") %>%
  mutate(label = LABELS[feature]) %>%
  arrange(value) %>%
  mutate(label = factor(label, levels = label))

imp_las <- coef(lasso_model) %>% as.matrix() %>% as.data.frame() %>%
  rownames_to_column("feature") %>% rename(coefficient = s0) %>%
  filter(feature != "(Intercept)") %>%
  mutate(value = abs(coefficient), label = LABELS[feature]) %>%
  arrange(value) %>%
  mutate(label = factor(label, levels = label))

cat("\nXGBoost SHAP 排序：", paste(rev(as.character(imp_xgb$label)), collapse = " > "), "\n")
cat("LASSO   系数排序：", paste(rev(as.character(imp_las$label)), collapse = " > "), "\n")

METAB <- c("SBP", "BMI", "FG", "LDL", "eGFR")
rank_xgb <- rev(as.character(imp_xgb$label)); rank_xgb <- rank_xgb[rank_xgb %in% METAB]
rank_las <- rev(as.character(imp_las$label)); rank_las <- rank_las[rank_las %in% METAB]
cat("五个暴露的排序 — XGBoost:", paste(rank_xgb, collapse = " > "), "\n")
cat("                LASSO  :", paste(rank_las, collapse = " > "), "\n")
if (!identical(rank_xgb, c("SBP", "BMI", "FG", "LDL", "eGFR")))
  warning("XGBoost 的暴露排序与附录第 275 段不符，请同步更新附录。", call. = FALSE)

# ============================================================================
# 6. 面板 A / B：变量重要性
# ============================================================================
bar_panel <- function(df, fill, xlab, digits) {
  ggplot(df, aes(x = value, y = label)) +
    geom_col(fill = fill, width = 0.65) +
    scale_x_continuous(labels = lab_dec(digits), n.breaks = 3,
                       expand = expansion(c(0, 0.04))) +
    labs(x = xlab, y = NULL) +
    theme_pub() +
    theme(panel.grid.major.y = element_blank())
}

p_A <- bar_panel(imp_xgb, COL_XGB, "Mean absolute SHAP value", 2)
p_B <- bar_panel(imp_las, COL_LAS, "Absolute standardised coefficient", 1)

save_panel(p_A, "FigS20_A_xgb_importance",   82, 78)
save_panel(p_B, "FigS20_B_lasso_importance", 82, 78)

# ============================================================================
# 7. 面板 C：ROC。【M1】颜色按模型固定，不随因子顺序变化
# ============================================================================
roc_df <- bind_rows(
  tibble(fpr = 1 - xgb_roc$specificities,   tpr = xgb_roc$sensitivities,   model = "XGBoost"),
  tibble(fpr = 1 - lasso_roc$specificities, tpr = lasso_roc$sensitivities, model = "LASSO")
)
auc_lab <- c(
  XGBoost = paste0("XGBoost (AUC ", dec(sprintf("%.2f", xgb_auc)), ")"),
  LASSO   = paste0("LASSO (AUC ",   dec(sprintf("%.2f", lasso_auc)), ")")
)

p_C <- ggplot(roc_df, aes(fpr, tpr, colour = model)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              colour = "grey60", linewidth = 0.3) +
  geom_line(linewidth = 0.5) +
  scale_colour_manual(values = c(XGBoost = COL_XGB, LASSO = COL_LAS),
                      labels = auc_lab, breaks = c("XGBoost", "LASSO"),
                      name = NULL) +
  scale_x_continuous(labels = lab_dec(2), breaks = c(0, 0.5, 1), expand = c(0, 0)) +
  scale_y_continuous(labels = lab_dec(2), breaks = c(0, 0.5, 1), expand = c(0, 0)) +
  annotate("text", x = 0.97, y = 0.06, hjust = 1, size = 2.5,
           colour = "grey20", family = "Helvetica",
           label = paste0("DeLong p=", fmt_p(delong$p.value))) +
  labs(x = "1 \u2212 specificity", y = "Sensitivity") +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  theme_pub() +
  theme(legend.position = c(0.98, 0.22), legend.justification = c(1, 0),
        legend.text = element_text(size = FS_BASE - 1),
        legend.key.height = unit(8, "pt"), legend.key.width = unit(12, "pt"),
        legend.background = element_blank())

p_C <- p_C + theme(plot.margin = margin(2, 9, 2, 2))   # 末端刻度 1·00 需要额外右边距
save_panel(p_C, "FigS20_C_roc", 54, 54)

# ============================================================================
# 8. 面板 D：五个暴露的 SHAP 依赖图
#    【M2】不作着色映射；x 轴按 0–99·5 百分位裁视野（平滑仍用全部数据）
# ============================================================================
UNITS <- c(SBP = "SBP, mm Hg", BMI = "BMI, kg/m\u00b2", FG = "FG, mg/dL",
           LDL = "LDL, mg/dL", eGFR = "eGFR, mL/min per 1\u00b773 m\u00b2")

dep_panel <- function(v) {
  d  <- shap_long %>% filter(variable == v)
  hi <- quantile(d$rfvalue, 0.995, na.rm = TRUE)
  ggplot(d, aes(rfvalue, value)) +
    geom_hline(yintercept = 0, colour = "grey80", linewidth = 0.25) +
    geom_point(colour = COL_PT, alpha = 0.25, size = 0.25, stroke = 0) +
    geom_smooth(method = "loess", formula = y ~ x, se = FALSE,
                colour = COL_XGB, linewidth = 0.5) +
    scale_x_continuous(labels = lab_dec(0)) +
    scale_y_continuous(labels = lab_dec(2)) +
    coord_cartesian(xlim = c(min(d$rfvalue, na.rm = TRUE), hi)) +
    labs(x = UNITS[[v]], y = "SHAP value", title = v) +
    theme_pub()
}

# 五个依赖图各自存盘，由拼图脚本排进 2-3-3 版面的后两行
DEP_STEMS <- c(SBP = "FigS20_D_dep_SBP", BMI = "FigS20_E_dep_BMI",
               FG  = "FigS20_F_dep_FG",  LDL = "FigS20_G_dep_LDL",
               eGFR = "FigS20_H_dep_eGFR")
for (v in METAB) save_panel(dep_panel(v), DEP_STEMS[[v]], 54, 50)

cat("\n八个面板已出齐（A B / C D E / F G H），拼图在 Python 端完成。\n")