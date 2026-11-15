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
# 敏感性分析 2（回归稀释）：数据整理 —— 子集化重构版 v2
#
# 数据来源分工（三个文件都要）：
#   CHARLS_敏感.csv   W3 测量（r3systo1–3、r3mbmi）、r3agey、波次状态、访谈年月
#                     —— 纳入链的主力，与旧版一致
#   W3_blood.csv      W3 空腹血检
#   charls_dat.csv    主分析合格队列（n = 8012）：作 ID 白名单，并提供 W1 侧的
#                     暴露与协变量。不需要重新导出，主分析代码无须改动。
#
# 相对旧版的三处结构性修改：
#
#  【R1】纳入链的起点改为"主队列成员"。第 1 步即 filter(ID %in% charls$ID)，
#        此后 W1 侧的 SBP/BMI/FG/LDL/eGFR 与全部协变量一律从 charls_dat.csv
#        join 进来，不在本脚本内重算。由此 charls_sens 在构造上必然是
#        charls 的子集，Figure S17 的分支画法才成立。连带效果：
#          - W1 的 SBP（第 2、3 次读数均值）、eGFR（含女性系数 1·012、用
#            r1agey）与主分析逐字节一致，旧版【S1】【S2】【S3】自动继承
#          - W1 血检的 join 不再需要（FG/LDL/eGFR 已在主数据中）
#          - W1 侧的生理学不可能值过滤已在主数据完成，此处只洗 W3 侧
#          - 协变量沿用主数据已编码、已按 >20% NA 筛过的同一套，两边模型
#            的协变量集合保证相同
#          - r3agey >= 45 的过滤删除：r1agey >= 45 已由主队列保证
#
#  【R2】ID 一律以定宽零填充字符串参与 join / 比对。旧版
#        mutate(ID = as.numeric(ID)) 会吃掉前导零，使省份编码 01–09 的
#        参与者在与主队列比对时全部失配。
#        注意 pad_id() 必须走字符串路径：formatC(format = "d") 会把值强制
#        转成 integer，12 位 ID 超出 .Machine$integer.max 会全部变成 NA。
#
#  【R3】删除 FIX_EGFR_W1_AGE 开关。eGFR_w1 继承主数据（已用 r1agey）；
#        eGFR_w3 用 r3agey（W3 抽血时年龄，无可选项）。协变量 age 由
#        顶部 AGE_COVAR 指定，默认 r3agey（本分析的基线是 W3）。
#
# 输出：charls_dat_敏感.csv           列名与主分析一致，供 29_charls_dilution_fine_gray.R 使用
#       sens2_flow.csv               逐步排除记录（Figure S17 敏感性分支②用）
#       sens2_exposure_compare.csv   新旧两套 SBP 定义的分布对照
# ============================================================================

setwd(RAWDAT)
library(tidyverse)

# ---- 开关 ------------------------------------------------------------------
FIX_SBP_READINGS <- TRUE      # S1：W3 收缩压取第 2、3 次读数（W1 侧已由主数据固定）
FIX_EGFR_FEMALE  <- TRUE      # S2：W3 eGFR 的女性系数 1·012（W1 侧已由主数据固定）
AGE_COVAR        <- "r3agey"  # 协变量 age 的来源："r3agey"（基线年龄）或 "r1agey"

OUT_DIR  <- file.path(RAWDAT, "CHARLS")
ID_WIDTH <- 12                # CHARLS 个人 ID 位数；pad_id() 会自检

stopifnot(AGE_COVAR %in% c("r3agey", "r1agey"))

# ---- ID 规范化：定宽零填充字符串（纯字符串路径，勿用 formatC format="d"）----
pad_id <- function(x, width = ID_WIDTH) {
  s <- if (is.numeric(x)) sprintf("%.0f", x) else as.character(x)
  s <- gsub("[^0-9]", "", s)
  s[!nzchar(s)] <- NA_character_
  out <- str_pad(s, width = width, side = "left", pad = "0")
  bad <- !is.na(out) & nchar(out) != width
  if (any(bad)) {
    warning(sprintf("pad_id(): %d 个 ID 规范化后长度不是 %d 位（实际长度 %s），请检查 ID_WIDTH。",
                    sum(bad), width,
                    paste(sort(unique(nchar(out[bad]))), collapse = "/")),
            call. = FALSE)
  }
  out
}

num_ <- function(x) as.numeric(ifelse(str_starts(as.character(x), "\\."), NA, x))

# ---- 逐步计数记录 ----------------------------------------------------------
flow <- list()
note <- function(step, n) {
  flow[[length(flow) + 1]] <<- data.frame(step = step, n = n)
  cat(sprintf("  %-46s n = %d\n", step, n))
}

# ============================================================================
# 0. 读入并断言（依赖全局对象的脚本一律先断言，见交接文档 v6 §4.1）
# ============================================================================
charls <- read.csv('charls_dat.csv', colClasses = c(ID = "character"))
charls$ID <- pad_id(charls$ID)

stopifnot(nrow(charls) == 8012)
stopifnot(!any(duplicated(charls$ID)))
cat(sprintf("主分析队列：n = %d，stroke = %d，competing death = %d\n",
            nrow(charls), sum(charls$new_stroke == 1), sum(charls$new_stroke == 2)))

har  <- read_csv('CHARLS_敏感.csv', show_col_types = FALSE)
w3_b <- read_csv('W3_blood.csv',    show_col_types = FALSE)

har$ID  <- pad_id(har$ID)
w3_b$ID <- pad_id(w3_b$ID)
stopifnot(!any(duplicated(har$ID)))

cat("\n===== 纳入流程（敏感性分析 2）=====\n")
note("主分析合格队列 charls_dat.csv", nrow(charls))

# ============================================================================
# 1. 限定在主队列内，再加 W3 / W4 条件
# ============================================================================
n_matched <- sum(charls$ID %in% har$ID)
if (n_matched < nrow(charls)) {
  warning(sprintf("有 %d 名主队列成员在 CHARLS_敏感.csv 中匹配不到，请核对 ID 来源。",
                  nrow(charls) - n_matched), call. = FALSE)
}

har <- har %>% filter(ID %in% charls$ID)
note("在 CHARLS_敏感.csv 中可匹配", nrow(har))

har <- har %>% filter(r3iwstat == "1.resp, alive" & r3stroke == "0.No")
note("W3 应答存活且无卒中", nrow(har))

har <- har %>% filter(r4iwstat %in% c("1.resp, alive", "5.nr, died this wv"))
note("W4 有随访信息", nrow(har))

# ---- W3 收缩压：三次读数 ----
har <- har %>%
  mutate(across(c(r1systo1, r1systo2, r1systo3,
                  r3systo1, r3systo2, r3systo3), num_),
         systo_w3 = if (FIX_SBP_READINGS) rowMeans(cbind(r3systo2, r3systo3), na.rm = TRUE)
         else                  rowMeans(cbind(r3systo1, r3systo2), na.rm = TRUE),
         # 旧定义，仅供分布对照
         systo_w1_old = rowMeans(cbind(r1systo1, r1systo2), na.rm = TRUE),
         systo_w3_old = rowMeans(cbind(r3systo1, r3systo2), na.rm = TRUE))

har <- har %>% filter(!is.na(systo_w3)); note("有 W3 收缩压", nrow(har))

# ---- W3 BMI ----
har <- har %>% mutate(r3mbmi = num_(r3mbmi)) %>% filter(!is.na(r3mbmi))
note("有 W3 BMI", nrow(har))

# ---- 年龄与访谈时间 ----
har <- har %>%
  mutate(r1agey = num_(r1agey),
         r3agey = num_(r3agey),
         across(c(r3iwy, r3iwm, r4iwy, r4iwm), num_),
         r3iwt = r3iwy + r3iwm / 12,
         r4iwt = r4iwy + r4iwm / 12,
         age_w3 = r3agey)
stopifnot(all(!is.na(har$age_w3)))

# ============================================================================
# 2. W3 空腹血检
# ============================================================================
w3_b <- w3_b %>%
  mutate(bl_fasting_num = as.numeric(str_extract(bl_fasting, "^\\d+"))) %>%
  filter(bl_fasting_num == 1) %>%
  drop_na(bl_glu, bl_crea, bl_ldl)

dat <- inner_join(har, w3_b, by = "ID")
note("合并 W3 空腹血检", nrow(dat))

# ============================================================================
# 3. 并入主队列的 W1 侧变量（暴露、协变量，全部沿用主分析定义）
# ============================================================================
dat <- inner_join(
  dat %>% select(ID, r1agey, r3agey, age_w3, r3iwt, r4iwt,
                 systo_w3, systo_w1_old, systo_w3_old,
                 r3mbmi, bl_glu, bl_crea, bl_ldl,
                 r3iwstat, r3stroke, r4iwstat, r4stroke),
  charls %>% select(ID, gender, marriage, education, rural, smoke, drink,
                    heartd, depression, bpmed, diamed, lipidmed,
                    SBP_w1 = SBP, BMI_w1 = BMI, FG_w1 = FG,
                    LDL_w1 = LDL, eGFR_w1 = eGFR),
  by = "ID"
)
stopifnot(nrow(dat) == n_distinct(dat$ID))

dat$age <- dat[[AGE_COVAR]]
cat(sprintf("  协变量 age 取自：%s\n", AGE_COVAR))

# ============================================================================
# 4. W3 的 eGFR（CKD-EPI 2021 肌酐方程，与主分析同式，年龄用 r3agey）
# ============================================================================
fem_coef <- if (FIX_EGFR_FEMALE) 1.012 else 1.000

egfr_ckdepi <- function(scr, sex, age) {
  kappa <- if_else(sex == 2, 0.7, 0.9)
  alpha <- if_else(sex == 2, -0.241, -0.302)
  sexc  <- if_else(sex == 2, fem_coef, 1)
  r <- scr / kappa
  142 * pmin(r, 1)^alpha * pmax(r, 1)^(-1.200) * (0.9938)^age * sexc
}

dat <- dat %>% mutate(eGFR_w3 = egfr_ckdepi(bl_crea, gender, age_w3))

# ============================================================================
# 5. W3 侧生理学范围清洗（阈值与主分析一致），再取两波均值
#    W1 侧不必再洗：charls_dat.csv 已在该步之后导出
# ============================================================================
rng <- function(x, lo, hi) ifelse(x >= lo & x <= hi, x, NA_real_)

dat <- dat %>%
  mutate(
    systo_w3 = rng(systo_w3, 70, 250),
    r3mbmi   = rng(r3mbmi,   12,  60),
    bl_glu   = rng(bl_glu,   40, 500),
    bl_ldl   = rng(bl_ldl,   20, 400),
    eGFR_w3  = rng(eGFR_w3,  15, 200),
    SBP  = (SBP_w1  + systo_w3) / 2,
    BMI  = (BMI_w1  + r3mbmi)   / 2,
    FG   = (FG_w1   + bl_glu)   / 2,
    LDL  = (LDL_w1  + bl_ldl)   / 2,
    eGFR = (eGFR_w1 + eGFR_w3)  / 2
  )

n_before <- nrow(dat)
dat <- dat %>% drop_na(SBP, BMI, FG, LDL, eGFR)
note(sprintf("W3 暴露在生理范围内（排除 %d）", n_before - nrow(dat)), nrow(dat))

# ============================================================================
# 6. 结局与时间（基线 = W3，随访至 W4）
# ============================================================================
r4iwt_mean <- mean(dat$r4iwt, na.rm = TRUE)

dat <- dat %>%
  mutate(
    r4iwt_f = coalesce(r4iwt, r4iwt_mean),
    new_stroke = case_when(
      r4stroke == "1.Yes"              ~ 1L,
      r4iwstat == "5.nr, died this wv" ~ 2L,
      TRUE                             ~ 0L
    ),
    time_yr = case_when(
      new_stroke == 2 ~ (r4iwt_f + r3iwt) / 2 - r3iwt,   # 竞争死亡取中点
      TRUE            ~ r4iwt_f - r3iwt
    )
  )

dat <- dat %>% filter(time_yr > 0)
note("随访时间 > 0", nrow(dat))

# ============================================================================
# 7. 精简与导出
# ============================================================================
charls_sens <- dat %>%
  select(ID, age, gender, marriage, education, rural, smoke, drink,
         heartd, depression, bpmed, diamed, lipidmed,
         SBP, BMI, FG, eGFR, LDL, new_stroke, time_yr,
         r1agey, r3agey, SBP_w1, systo_w3, systo_w1_old, systo_w3_old,
         eGFR_w1, eGFR_w3)

# ---- 子集关系断言：Figure S17 的分支画法依赖这一条 ----
stopifnot(all(charls_sens$ID %in% charls$ID))
cat(sprintf("\n子集校验通过：%d 名敏感性分析对象全部来自 %d 人的主队列\n",
            nrow(charls_sens), nrow(charls)))

# ---- 协变量缺失（Figure S17 分支图注用）----
confounders <- c("marriage", "education", "rural", "smoke", "drink",
                 "heartd", "depression", "bpmed", "diamed", "lipidmed")
na_report <- charls_sens %>%
  summarise(across(all_of(confounders), ~ sum(is.na(.)))) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "n_missing") %>%
  mutate(pct = round(n_missing / nrow(charls_sens) * 100, 2))
cat("\n===== 协变量缺失 =====\n"); print(as.data.frame(na_report))

cat(sprintf("\n敏感性分析 2 样本：N = %d | 卒中 %d | 竞争死亡 %d | 无事件 %d\n",
            nrow(charls_sens),
            sum(charls_sens$new_stroke == 1),
            sum(charls_sens$new_stroke == 2),
            sum(charls_sens$new_stroke == 0)))
cat("注意：以上为暴露完整者的人数；进模型人数还要再扣除协变量缺失，\n")
cat("      由 29_charls_dilution_fine_gray.R 的 build_design() 报告。\n")

# ---- 新旧 SBP 定义对照 ----
cmp <- data.frame(
  quantity = c("SBP W1 (主数据)", "SBP W3"),
  new_mean = c(mean(charls_sens$SBP_w1,   na.rm = TRUE),
               mean(charls_sens$systo_w3, na.rm = TRUE)),
  old_mean = c(mean(charls_sens$systo_w1_old, na.rm = TRUE),
               mean(charls_sens$systo_w3_old, na.rm = TRUE)),
  new_sd   = c(sd(charls_sens$SBP_w1,   na.rm = TRUE),
               sd(charls_sens$systo_w3, na.rm = TRUE)),
  old_sd   = c(sd(charls_sens$systo_w1_old, na.rm = TRUE),
               sd(charls_sens$systo_w3_old, na.rm = TRUE))
)
cmp$mean_diff <- cmp$new_mean - cmp$old_mean
cat("\n===== 新旧 SBP 定义对照（新 = 修正后）=====\n")
print(cmp, digits = 4, row.names = FALSE)

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
write.csv(charls_sens, 'charls_dat_敏感.csv', row.names = FALSE)
write.csv(charls_sens[1:10, ], 'charls_dat_敏感_struc.csv', row.names = FALSE)
write.csv(bind_rows(flow), file.path(OUT_DIR, "sens2_flow.csv"), row.names = FALSE)
write.csv(cmp, file.path(OUT_DIR, "sens2_exposure_compare.csv"), row.names = FALSE)

cat("\n已写出 charls_dat_敏感.csv / sens2_flow.csv / sens2_exposure_compare.csv\n")
