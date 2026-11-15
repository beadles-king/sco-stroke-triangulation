# =============================================================================
# Reproducibility repository - Wang et al., transportability-aware triangulation
# of metabolic stroke evidence across the Shanghai Cooperation Organisation.
# Paths are resolved by config.R; see README.md.
# =============================================================================
if (!exists("RAWDAT")) {
  if (file.exists("config.R")) source("config.R") else
    stop("Run this script from the code/ directory (the folder containing config.R).")
}

# ============================================================
# Step 3b-i (REVISED): 非线性检验 —— 基于 Fine-Gray 亚分布风险的 RCS
#
# 替换原脚本中的整个 "Step 3b-i" 区块：
#   从  dd <- datadist(charls)
#   到  print(rcs_nonlinear)
# 全部删除，用本文件内容替换。
#
# 原实现用 rms::cph 拟合 cause-specific Cox（竞争死亡被删失），
# 与主分析 Fine-Gray 不同尺度。本实现改为在亚分布风险尺度上拟合样条：
#   survival::finegray() 构造 IPCW 加权风险集 → coxph(weights = fgwt, robust)
#   该拟合与 cmprsk::crr 在数学上等价，但支持样条基与 anova/预测。
#
# 依赖：survival (>= 3.0), Hmisc, dplyr, ggplot2, patchwork
# 注：可从原脚本移除 library(rms)；rcspline.eval 来自 Hmisc（rms 的依赖）。
# ============================================================

library(survival)
library(Hmisc)
library(dplyr)
library(ggplot2)
library(patchwork)
library(ragg)

stopifnot(exists("charls"),
          nrow(charls) == 8012,
          all(EXPO_SD %in% names(charls)))
cat(sprintf("RCS 输入校验通过：charls n=%d\n", nrow(charls)))

# ---- 输出位置：改这一行即可 ----
OUT_DIR   <- file.path(RAWDAT, "CHARLS/RCS")                   # 例如 file.path(DATA_ROOT, "figures")
out_path  <- function(f) {
  if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)
  file.path(OUT_DIR, f)
}


exposures_cont <- c("SBP", "FG", "LDL", "eGFR", "BMI")

cov_vars <- c("age", "female", "marriage", "education",
              "rural", "smoke", "drink", "heartd", "depression")
charls <- charls %>%
  dplyr::mutate(
    marriage  = factor(marriage,  levels = sort(unique(marriage))),
    education = factor(education, levels = sort(unique(education)))
  )
cat("\n===== marriage / education 水平分布 =====\n")
print(table(charls$marriage,  useNA = "ifany"))
print(table(charls$education, useNA = "ifany"))

NK   <- 4      # 节点数，与原分析一致（rms 默认分位 .05/.35/.65/.95）
NORM <- 2      # rcspline.eval 的标准化方式，网格预测时必须保持一致

# 参考值：显式指定，避免 log(HR)=0 的位置含糊
# （原稿附录"HR=1 在 26 kg/m2"之所以自相矛盾，就是因为参考点未定义）
ref_values <- c(SBP = 115, FG = 90, LDL = 100, eGFR = 100, BMI = 22)


# ------------------------------------------------------------
# 0. 数据准备：事件因子 + 行 id
# ------------------------------------------------------------
charls <- charls %>%
  mutate(
    id        = dplyr::row_number(),
    status_fg = factor(new_stroke, levels = c(0, 1, 2),
                       labels = c("censor", "stroke", "death"))
  )

# finegray() 不接受 time <= 0 与缺失值
stopifnot(all(charls$time_yr > 0, na.rm = TRUE))


# ------------------------------------------------------------
# 1. 构建 RCS 基（必须在 finegray 展开之前构建）
#    finegray 会按事件时间切分行；样条基是个体层面的常量，
#    先建好再随 "~ ." 一并带入展开后的数据框。
# ------------------------------------------------------------
KNOTS_ON <- "analysis"     # "analysis"（推荐）或 "full"（复现旧结果）

# 先算出分析集（与 fg_base 的成列删除口径一致）
.need_cols <- c("time_yr", "status_fg", "id", cov_vars,
                paste0(exposures_cont, "_sd"), exposures_cont)
.ana_idx   <- stats::complete.cases(charls[, .need_cols])
cat(sprintf("\nRCS 分析人群：%d -> %d（剔除 %d）；事件 %d -> %d\n",
            nrow(charls), sum(.ana_idx), nrow(charls) - sum(.ana_idx),
            sum(charls$new_stroke == 1),
            sum(charls$new_stroke[.ana_idx] == 1)))

basis_info <- list()
for (expo in exposures_cont) {
  x_knot <- if (KNOTS_ON == "analysis") charls[[expo]][.ana_idx] else charls[[expo]]
  kn <- Hmisc::rcspline.eval(x_knot, nk = NK, inclx = TRUE, knots.only = TRUE)
  B  <- Hmisc::rcspline.eval(charls[[expo]], knots = kn, inclx = TRUE, norm = NORM)
  colnames(B) <- paste0(expo, "_rcs", seq_len(ncol(B)))
  charls[, colnames(B)] <- as.data.frame(B)
  basis_info[[expo]] <- list(knots = kn, cols = colnames(B))
  cat(sprintf("%-5s knots (%s): %s\n", expo, KNOTS_ON,
              paste(round(kn, 2), collapse = ", ")))
}


# ------------------------------------------------------------
# 2. 构造 Fine-Gray 加权风险集
# ------------------------------------------------------------
all_basis_cols <- unlist(lapply(basis_info, `[[`, "cols"), use.names = FALSE)
expo_sd_cols   <- paste0(exposures_cont, "_sd")

fg_base <- charls[, c("time_yr", "status_fg", "id",
                      cov_vars, all_basis_cols, expo_sd_cols)]
fg_base <- fg_base[stats::complete.cases(fg_base), ]
stopifnot(nrow(fg_base) == 7896, sum(fg_base$status_fg == "stroke") == 565)

cat("\nfinegray 输入 n =", nrow(fg_base), "\n")
rcs_ledger <- data.frame(
  analysis        = "RCS Fine-Gray (Model 1 设定)",
  n_before        = nrow(charls),
  n_analysis      = nrow(fg_base),
  events_before   = sum(charls$new_stroke == 1),
  events_analysis = sum(fg_base$status_fg == "stroke"),
  knots_on        = KNOTS_ON,
  row.names = NULL
)
print(rcs_ledger)
write.csv(rcs_ledger, out_path("RCS_n_ledger.csv"), row.names = FALSE)

fg_dat <- survival::finegray(Surv(time_yr, status_fg) ~ .,
                             data = fg_base, etype = "stroke")

cat("finegray 展开后行数 =", nrow(fg_dat),
    "| 事件数 =", sum(fg_dat$fgstatus), "\n\n")


# ------------------------------------------------------------
# 3. 拟合函数
#    mutual = FALSE : 单暴露样条 + 协变量（对应 Model 1 的设定，
#                     与原 cause-specific 版本可比，作为主分析）
#    mutual = TRUE  : 目标暴露样条 + 其余四暴露线性项 + 协变量
#                     （对应 Model 2 的设定，作为敏感性分析）
#    ties = "breslow" 与 crr 的处理方式一致；本队列因中点插补
#    存在大量结点时间，tie 处理方式的选择会影响估计。
# ------------------------------------------------------------
fit_rcs_fg <- function(expo, mutual = FALSE) {
  cols  <- basis_info[[expo]]$cols
  other <- setdiff(exposures_cont, expo)
  rhs   <- c(cols, cov_vars, if (mutual) paste0(other, "_sd"))
  f <- as.formula(paste("Surv(fgstart, fgstop, fgstatus) ~",
                        paste(rhs, collapse = " + ")))
  survival::coxph(f, data = fg_dat, weights = fgwt,
                  id = id, robust = TRUE, ties = "breslow")
}


# 联合 Wald 检验（稳健方差）
wald_joint <- function(fit, terms) {
  idx <- match(terms, names(coef(fit)))
  b   <- coef(fit)[idx]
  V   <- vcov(fit)[idx, idx, drop = FALSE]
  st  <- as.numeric(t(b) %*% solve(V) %*% b)
  df  <- length(idx)
  c(chisq = st, df = df, p = pchisq(st, df, lower.tail = FALSE))
}


# ------------------------------------------------------------
# 4. 主循环：拟合 + 非线性检验
# ------------------------------------------------------------
rcs_fg <- list()

for (expo in exposures_cont) {
  cols <- basis_info[[expo]]$cols
  fit  <- fit_rcs_fg(expo, mutual = FALSE)
  
  w_overall   <- wald_joint(fit, cols)        # 总体关联（df = 3）
  w_nonlinear <- wald_joint(fit, cols[-1])    # 非线性（df = 2）
  
  rcs_fg[[expo]] <- list(fit = fit,
                         overall = w_overall,
                         nonlinear = w_nonlinear)
  
  cat(sprintf("%-5s  overall chi2=%6.2f (df=%d) P=%.4f | nonlinear chi2=%6.2f (df=%d) P=%.4f\n",
              expo, w_overall["chisq"], w_overall["df"], w_overall["p"],
              w_nonlinear["chisq"], w_nonlinear["df"], w_nonlinear["p"]))
}

# FDR 校正（5 个暴露为一个检验家族，与原设定一致）
rcs_nonlinear <- data.frame(
  variable    = exposures_cont,
  chisq       = round(sapply(exposures_cont, function(e) rcs_fg[[e]]$nonlinear["chisq"]), 2),
  df          = sapply(exposures_cont, function(e) rcs_fg[[e]]$nonlinear["df"]),
  P_nonlinear = round(sapply(exposures_cont, function(e) rcs_fg[[e]]$nonlinear["p"]), 4),
  P_overall   = round(sapply(exposures_cont, function(e) rcs_fg[[e]]$overall["p"]), 4),
  row.names   = NULL
)
rcs_nonlinear$P_fdr <- round(p.adjust(rcs_nonlinear$P_nonlinear, method = "fdr"), 4)

cat("\n===== Table S16 (Fine-Gray RCS 非线性检验) =====\n")
print(rcs_nonlinear)


# ------------------------------------------------------------
# 5. 预测曲线：相对于显式参考值的 log(SHR) 及 95% CI
# ------------------------------------------------------------
predict_rcs_fg <- function(expo, ref = ref_values[[expo]], n_grid = 200,
                           lo = 0.01, hi = 0.99) {
  info <- basis_info[[expo]]
  fit  <- rcs_fg[[expo]]$fit
  x    <- charls[[expo]]
  
  grid <- seq(quantile(x, lo, na.rm = TRUE),
              quantile(x, hi, na.rm = TRUE), length.out = n_grid)
  
  # 参考点与网格一起变换，避免单点调用 rcspline.eval 的边界问题
  Ball <- Hmisc::rcspline.eval(c(ref, grid), knots = info$knots,
                               inclx = TRUE, norm = NORM)
  Bref <- Ball[1, , drop = FALSE]
  Bg   <- Ball[-1, , drop = FALSE]
  D    <- sweep(Bg, 2, Bref, "-")          # 对比矩阵
  
  idx  <- match(info$cols, names(coef(fit)))
  beta <- coef(fit)[idx]
  V    <- vcov(fit)[idx, idx, drop = FALSE]
  
  lp <- as.numeric(D %*% beta)
  se <- sqrt(pmax(rowSums((D %*% V) * D), 0))
  
  data.frame(
    exposure = expo, x = grid, ref = ref,
    log_shr  = lp,
    lower    = lp - 1.96 * se,
    upper    = lp + 1.96 * se,
    SHR      = exp(lp),
    SHR_low  = exp(lp - 1.96 * se),
    SHR_high = exp(lp + 1.96 * se)
  )
}

pred_all <- bind_rows(lapply(exposures_cont, predict_rcs_fg))


# ------------------------------------------------------------
# 6. 曲线特征：最低风险点、过参考线的位置、指定点的 SHR
#    （正文与附录需要引用的拐点数值由此产生）
# ------------------------------------------------------------
curve_features <- pred_all %>%
  group_by(exposure) %>%
  summarise(
    ref            = first(ref),
    x_min_risk     = x[which.min(log_shr)],                       # 风险最低点
    SHR_at_min     = round(min(SHR), 3),
    x_cross_ref    = {                                            # 曲线回到 log(SHR)=0 的位置
      s <- sign(log_shr); j <- which(diff(s) != 0)
      if (length(j) == 0) NA_real_ else round(x[j[length(j)]], 1)
    },
    x_p10 = round(quantile(charls[[first(exposure)]][.ana_idx], .10, na.rm = TRUE), 1),
    x_p90 = round(quantile(charls[[first(exposure)]][.ana_idx], .90, na.rm = TRUE), 1),
    .groups = "drop"
  )

cat("\n===== 曲线特征（用于正文拐点表述） =====\n")
print(as.data.frame(curve_features))

# 指定值处的 SHR（按需修改 at_values）
at_values <- list(SBP = c(110, 120, 140, 160, 180),
                  FG  = c(80, 90, 97, 110, 126, 160),
                  LDL = c(80, 100, 130, 160, 190),
                  eGFR= c(45, 60, 80, 100, 120),
                  BMI = c(18, 20, 22, 25, 28, 32))

shr_at_points <- bind_rows(lapply(exposures_cont, function(e) {
  p <- pred_all %>% filter(exposure == e)
  do.call(rbind, lapply(at_values[[e]], function(v) {
    k <- which.min(abs(p$x - v))
    data.frame(exposure = e, x_target = v, x_used = round(p$x[k], 2),
               SHR = round(p$SHR[k], 3),
               CI = sprintf("%.2f-%.2f", p$SHR_low[k], p$SHR_high[k]))
  }))
}))
cat("\n===== 指定值处的 SHR（相对参考值） =====\n")
print(shr_at_points)


# ---- 成对对比：SHR(x2 vs x1) 及 delta 法 95% CI --------------
contrast_rcs <- function(expo, x1, x2) {
  info <- basis_info[[expo]]
  fit  <- rcs_fg[[expo]]$fit
  B    <- Hmisc::rcspline.eval(c(x1, x2), knots = info$knots,
                               inclx = TRUE, norm = NORM)
  d    <- B[2, ] - B[1, ]
  idx  <- match(info$cols, names(coef(fit)))
  b    <- coef(fit)[idx]; V <- vcov(fit)[idx, idx, drop = FALSE]
  est  <- sum(d * b)
  se   <- sqrt(as.numeric(t(d) %*% V %*% d))
  data.frame(exposure = expo, x1 = x1, x2 = x2,
             SHR   = round(exp(est), 3),
             low   = round(exp(est - 1.96 * se), 3),
             high  = round(exp(est + 1.96 * se), 3),
             P     = signif(2 * pnorm(-abs(est / se)), 3),
             per_unit_pct = round((exp(est / (x2 - x1)) - 1) * 100, 1))
}

# 正文与附录需要的全部对比
contrasts_needed <- rbind(
  contrast_rcs("BMI", 20, 25),   # GBD TMREL 窗口内部梯度 → 1.80
  contrast_rcs("BMI", 25, 28),   # 平坦段斜率
  contrast_rcs("BMI", 22, 28),
  contrast_rcs("FG",  90, 110),  # 陡峭段斜率
  contrast_rcs("FG",  110, 126), # 平坦段斜率
  contrast_rcs("FG",  90, 97),   # FG 的 TMREL 窗口内部梯度 → 1.19
  contrast_rcs("SBP", 115, 140),
  contrast_rcs("SBP", 115, 180)
)
print(contrasts_needed)

# 目标值处精确 SHR（取代 shr_at_points 的网格近似）
exact_shr <- function(expo, xs, ref = ref_values[[expo]])
  do.call(rbind, lapply(xs, function(x) contrast_rcs(expo, ref, x)))
print(exact_shr("BMI", c(18, 20, 25, 28, 32)))
print(exact_shr("FG",  c(80, 97, 110, 126, 160)))
print(exact_shr("SBP", c(110, 140, 160, 180)))
print(exact_shr("LDL", c(80, 130, 160, 190)))
print(exact_shr("eGFR", c(45, 60, 80, 120)))


# ------------------------------------------------------------
# 7.0 全文格式约定
#   小数号  中点 U+00B7
#   p 值    p < 0.0001         -> "<0·0001"
#           0.0001 <= p < 0.01 -> 保持四位小数
#           p >= 0.01          -> 两位小数
#   区间    en dash U+2013；负号 minus U+2212
#   字号    正文 8 pt，面板字母 11 pt 粗体
#   幅面    180 mm 双栏上限，高度 <= 240 mm
# ------------------------------------------------------------
DEC       <- "\u00B7"
EN_DASH   <- "\u2013"
MINUS     <- "\u2212"
BASE_SIZE <- 8
FONT      <- "Helvetica"     # Helvetica 等宽替代
COL_LINE  <- "#2C5985"             # 与 Figure 5 同色
FIG_W_MM  <- 180
FIG_H_MM  <- 118


mm2in <- function(x) x / 25.4

# 小数点 -> 中点；ASCII 连字符 -> 真减号（S19 的 y 轴有负值，必需）
dec <- function(x) {
  x <- gsub("\\.", DEC, as.character(x))
  gsub("(?<![0-9A-Za-z])-(?=[0-9])", MINUS, x, perl = TRUE)
}

fmt_p <- function(p) {
  vapply(p, function(pi) {
    if (is.na(pi))   return(NA_character_)
    if (pi < 0.0001) return(paste0("<0", DEC, "0001"))
    if (pi < 0.01)   return(dec(sprintf("%.4f", pi)))
    dec(sprintf("%.2f", pi))
  }, character(1))
}

# "p<0·0001" 不带等号，其余带等号
p_inline <- function(prefix, p) {
  paste0(prefix, if (p < 0.0001) "" else "=", fmt_p(p))
}

# 兜底：把成图 gtable 内所有 text grob 的小数号与负号统一
dotify <- function(g) {
  if (inherits(g, "gTree") || inherits(g, "gtable")) {
    g$children <- lapply(g$children, dotify)
  }
  if (inherits(g, "text") && !is.null(g$label)) g$label <- dec(g$label)
  g
}


# ------------------------------------------------------------
# 7.1 面板元数据
#   顺序按效应强度：SBP、BMI、FG、LDL、eGFR
#   （与 Figure 4、Figure S18、Table S15–S18 一致；
#    不改动 exposures_cont，故建模与 FDR 家族均不受影响）
# ------------------------------------------------------------
plot_order <- c("SBP", "BMI", "FG", "LDL", "eGFR")

axis_labs <- c(
  SBP  = "Systolic blood pressure (mmHg)",
  BMI  = "Body-mass index (kg/m\u00B2)",
  FG   = "Fasting plasma glucose (mg/dL)",
  LDL  = "LDL cholesterol (mg/dL)",
  eGFR = paste0("eGFR (mL/min per 1", DEC, "73 m\u00B2)")
)

theme_pub <- function() {
  theme_classic(base_size = BASE_SIZE, base_family = FONT) +
    theme(
      plot.title  = element_text(size = BASE_SIZE, face = "plain", hjust = 0,
                                 margin = margin(b = 3)),
      axis.title  = element_text(size = BASE_SIZE),
      axis.text   = element_text(size = BASE_SIZE, colour = "black"),
      axis.line   = element_line(linewidth = 0.3, colour = "black"),
      axis.ticks  = element_line(linewidth = 0.3, colour = "black"),
      plot.margin = margin(2, 4, 2, 2, "mm")
    )
}


# ------------------------------------------------------------
# 7.2 单面板
# ------------------------------------------------------------
plot_rcs_panel <- function(expo) {
  
  d  <- pred_all %>% filter(exposure == expo)
  kn <- basis_info[[expo]]$knots
  pn <- rcs_nonlinear$P_nonlinear[rcs_nonlinear$variable == expo]
  pf <- rcs_nonlinear$P_fdr[rcs_nonlinear$variable == expo]
  
  # 只画落在绘图区（1st–99th 百分位）内的节点
  kn <- kn[kn >= min(d$x) & kn <= max(d$x)]
  
  lab <- paste0(p_inline("Nonlinearity p", pn), "\n",
                p_inline("FDR p", pf))
  
  ggplot(d, aes(x = x)) +
    geom_hline(yintercept = 0, linetype = "dashed",
               colour = "grey45", linewidth = 0.3) +
    geom_vline(xintercept = d$ref[1], linetype = "dotted",
               colour = "grey45", linewidth = 0.3) +
    geom_ribbon(aes(ymin = lower, ymax = upper),
                fill = COL_LINE, alpha = 0.18) +
    geom_line(aes(y = log_shr), colour = COL_LINE, linewidth = 0.5) +
    geom_rug(data = data.frame(k = kn), aes(x = k), inherit.aes = FALSE,
             sides = "b", length = unit(2, "mm"),
             linewidth = 0.4, colour = "grey20") +
    scale_x_continuous(
      labels = function(v) dec(format(v, trim = TRUE, scientific = FALSE)),
      expand = expansion(mult = c(0.02, 0.02))
    ) +
    scale_y_continuous(
      breaks = scales::pretty_breaks(n = 5),
      labels = function(v) dec(sprintf("%.1f", v)),
      expand = expansion(mult = c(0.05, 0.22))
    ) +
    annotate("text", x = -Inf, y = Inf, label = lab,
             hjust = -0.04, vjust = 1.25, size = BASE_SIZE / .pt,
             lineheight = 1.05, colour = "grey25", family = FONT) +
    labs(x     = axis_labs[[expo]],
         y     = "log SHR") +
    theme_pub()
}

rcs_plots_fg <- lapply(plot_order, plot_rcs_panel)
names(rcs_plots_fg) <- plot_order


# ------------------------------------------------------------
# 7.3 拼接：3 列 × 2 行，第 6 格留空，面板字母 A–E
# ------------------------------------------------------------
figS19 <- patchwork::wrap_plots(rcs_plots_fg, ncol = 3) +
  patchwork::plot_spacer() +
  patchwork::plot_layout(ncol = 3) +
  patchwork::plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 11, face = "bold", family = FONT),
        plot.tag.position = c(0, 1))

print(figS19)


# ------------------------------------------------------------
# 7.4 导出：矢量母版（字体嵌入）+ 600/300 dpi LZW TIFF
# ------------------------------------------------------------
save_figure <- function(pw, stem, w_mm = FIG_W_MM, h_mm = FIG_H_MM) {
  draw <- function() {
    grid.newpage()
    grid.draw(dotify(patchwork::patchworkGrob(pw)))
  }
  # 原先此处另出一份 cairo_pdf 矢量母版, 已移除: cairo 把每个字形的步进四舍五入
  # 到整数 point (0.353 mm), 再用 TJ 微调补偿, 同一单词内字距不匀 (最坏 1 pt)。
  # 投稿件本就是下面两份 TIF, ragg 在目标分辨率排字, 600 dpi 误差 <= 0.042 mm。
  agg_tiff(out_path(paste0(stem, ".tif")), width = w_mm, height = h_mm, units = "mm",
           res = 600, compression = "lzw", background = "white")
  draw(); dev.off()
  agg_tiff(out_path(paste0(stem, "_300dpi.tif")), width = w_mm, height = h_mm,
           units = "mm", res = 300, compression = "lzw", background = "white")
  draw(); dev.off()
  invisible(NULL)
}

save_figure(figS19, "FigureS19")


# ------------------------------------------------------------
# 7.5 导出配套数据表（数值一律带格式化显示列，供直接贴表）
# ------------------------------------------------------------
# Table S16：总体关联 + 非线性 + FDR
tabS16_out <- rcs_nonlinear %>%
  mutate(
    reference        = ref_values[variable],
    P_overall_disp   = fmt_p(P_overall),
    P_nonlinear_disp = fmt_p(P_nonlinear),
    P_fdr_disp       = fmt_p(P_fdr)
  ) %>%
  select(variable, reference, chisq, df,
         P_overall, P_overall_disp,
         P_nonlinear, P_nonlinear_disp,
         P_fdr, P_fdr_disp)
write.csv(tabS16_out, out_path("TableS16_RCS_FineGray.csv"), row.names = FALSE)
cat("\n===== Table S16（含显示列） =====\n"); print(tabS16_out)

# Table S18：指定取值点的 SHR（用 contrast_rcs 的精确值，非网格近似）
tabS18_out <- bind_rows(lapply(plot_order, function(e) {
  xs <- setdiff(at_values[[e]], ref_values[[e]])
  do.call(rbind, lapply(xs, function(v) {
    r <- contrast_rcs(e, ref_values[[e]], v)
    data.frame(
      exposure  = e,
      reference = ref_values[[e]],
      value     = v,
      SHR       = r$SHR, CI_low = r$low, CI_high = r$high, P = r$P,
      display   = sprintf("%s (%s%s%s)",
                          dec(sprintf("%.2f", r$SHR)),
                          dec(sprintf("%.2f", r$low)), EN_DASH,
                          dec(sprintf("%.2f", r$high))),
      P_disp    = fmt_p(r$P)
    )
  }))
}))
write.csv(tabS18_out, out_path("TableS18_RCS_contrasts.csv"), row.names = FALSE)
cat("\n===== Table S18（含显示列） =====\n"); print(tabS18_out)

# 正文引用的成对对比（BMI 20→25 等）
contrasts_out <- contrasts_needed %>%
  mutate(display = sprintf("%s (%s%s%s)",
                           dec(sprintf("%.2f", SHR)),
                           dec(sprintf("%.2f", low)), EN_DASH,
                           dec(sprintf("%.2f", high))),
         P_disp  = fmt_p(P))
write.csv(contrasts_out, out_path("RCS_contrasts_for_text.csv"), row.names = FALSE)
cat("\n===== 正文成对对比（含显示列） =====\n"); print(contrasts_out)

write.csv(curve_features, out_path("RCS_curve_features.csv"), row.names = FALSE)


# ------------------------------------------------------------
# 8. 敏感性分析 A：互调整设定（对应 Model 2）
# ------------------------------------------------------------
rcs_mutual <- bind_rows(lapply(exposures_cont, function(expo) {
  cols <- basis_info[[expo]]$cols
  fit  <- fit_rcs_fg(expo, mutual = TRUE)
  w    <- wald_joint(fit, cols[-1])
  data.frame(variable = expo, chisq = round(w["chisq"], 2),
             P_nonlinear = round(w["p"], 4), row.names = NULL)
}))
rcs_mutual$P_fdr <- round(p.adjust(rcs_mutual$P_nonlinear, method = "fdr"), 4)
cat("\n===== 敏感性分析 A：互调整（Model 2 设定）下的非线性检验 =====\n")
print(rcs_mutual)


# ------------------------------------------------------------
# 9. 敏感性分析 B：保留原 cause-specific Cox 结果作为对照
#    （不进正文，仅在附录说明两种尺度结论一致/不一致）
# ------------------------------------------------------------
rcs_cs <- bind_rows(lapply(exposures_cont, function(expo) {
  cols <- basis_info[[expo]]$cols
  f <- as.formula(paste("Surv(time_yr, new_stroke == 1) ~",
                        paste(c(cols, cov_vars), collapse = " + ")))
  fit <- survival::coxph(f, data = charls, ties = "breslow")
  w   <- wald_joint(fit, cols[-1])
  data.frame(variable = expo, P_nonlinear = round(w["p"], 4), row.names = NULL)
}))
rcs_cs$P_fdr <- round(p.adjust(rcs_cs$P_nonlinear, method = "fdr"), 4)
cat("\n===== 敏感性分析 B：cause-specific Cox（原设定）对照 =====\n")
print(rcs_cs)


# ------------------------------------------------------------
# 10. 一致性核验：finegray+coxph 应重现 crr 的线性估计
#     若两者 SHR 差异 >0.02，说明 tie 处理或权重设定需要复查
# ------------------------------------------------------------
check_lin <- survival::coxph(
  as.formula(paste("Surv(fgstart, fgstop, fgstatus) ~",
                   paste(c(expo_sd_cols, cov_vars), collapse = " + "))),
  data = fg_dat, weights = fgwt, id = id, robust = TRUE, ties = "breslow")

cat("\n===== 核验：finegray+coxph vs crr（Model 2 线性设定）=====\n")
print(data.frame(
  variable     = expo_sd_cols,
  SHR_finegray = round(exp(coef(check_lin)[expo_sd_cols]), 3),
  SHR_crr      = round(exp(fg_m2$coef[expo_sd_cols]), 3),
  row.names    = NULL
))
# 两列差异 > 0·02 时需复查 tie 处理与权重设定。