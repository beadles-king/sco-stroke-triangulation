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
library(readr)
library(tidyverse)

har <- read_csv('Har_charls.csv')
w1_b <- read_csv('W1_blood.csv')


##har处理
#w1-参与随访且无中风
har <- har %>%
  filter(
    r1iwstat == "1.resp, alive" &
    r1stroke == '0.No'
  )

#w2-此波失访且此后未参与
har <- har %>%
  filter(
    !(r2iwstat %in% c("4.nr, alive", "9.nr, dk if alive or died") & 
        !r3iwstat %in% c("1.resp, alive", "5.nr, died this wv", "6.nr, died prev wv") & 
        !r4iwstat %in% c("1.resp, alive", "5.nr, died this wv", "6.nr, died prev wv"))
  )

#有收缩压记录
har <- har %>%
  mutate(
    r1systo2 = as.numeric(ifelse(str_starts(r1systo2, "\\."), NA, r1systo2)),
    r1systo3 = as.numeric(ifelse(str_starts(r1systo3, "\\."), NA, r1systo3)),  # 修正：右侧改回 r1systo3
    systo = rowMeans(cbind(r1systo2, r1systo3), na.rm = TRUE)                  # 改：取第2、3次均值
  ) %>%
  filter(!is.na(systo))   # is.na() 同时捕获两次皆缺失产生的 NaN

#有bmi记录
har <- har %>%
  mutate(
    r1mbmi = as.numeric(ifelse(str_starts(r1mbmi, "\\."), NA, r1mbmi))
  ) %>%
  filter(!is.na(r1mbmi))

#有年龄记录且大于等于45岁
har <- har %>%
  mutate(
    r1agey = as.numeric(ifelse(str_starts(r1agey, "\\."), NA, r1agey))
  ) %>%
  filter(!is.na(r1agey)) 
har <- har %>% filter(r1agey >= 45)

##合并har和w1_b
dat <- inner_join(
  har %>% mutate(ID_w1 = as.numeric(ID_w1)),
  w1_b %>% mutate(ID = as.numeric(ID)),
  by = c("ID_w1" = "ID")
)

##FG、LDL、Crea缺失
dat <- dat %>%
  drop_na(newglu, newcrea, newldl)

##非空腹
dat <- dat %>%
  filter(qc1_va003 == 1) 

#确定访谈时间
dat <- dat %>%
  mutate(
    r1iwt = r1iwy + r1iwm/12,
    r2iwt = r2iwy + r2iwm/12,
    r3iwt = r3iwy + r3iwm/12,
    r4iwt = r4iwy + r4iwm/12)


##数据整理
#性别整理
dat <- dat %>%
  mutate(
    ragender = as.numeric(str_extract(ragender, "^\\d+"))
  )

#学历整理
dat <- dat %>%
  mutate(
    raeducl = as.numeric(str_extract(raeducl, "^\\d+"))
  )

#婚姻状态整理
dat <- dat %>%
  mutate(
    r1mstat = case_when(
      r1mstat %in% c("1.married", "3.partnered") ~ 1,
      r1mstat %in% c("4.separated", "5.divorced") ~ 2,
      r1mstat == "7.widowed" ~ 3,
      r1mstat == "8.never married" ~ 4
    )
  )

#城乡居住地整理
dat <- dat %>%
  mutate(
    h1rural = as.numeric(str_extract(h1rural, "^\\d+"))
  )

#体力运动整理
dat <- dat %>%
  mutate(
    r1mdact_c = ifelse(str_starts(r1mdact_c, "\\."), NA, 
                      as.numeric(str_extract(r1mdact_c, "^\\d+")))
  )

#饮酒状态整理
dat <- dat %>%
  mutate(
    r1drinkev = ifelse(str_starts(r1drinkev, "\\."), NA, 
                       as.numeric(str_extract(r1drinkev, "^\\d+")))
  )

#吸烟状态整理
dat <- dat %>%
  mutate(
    r1smokev = ifelse(str_starts(r1smokev, "\\."), NA, 
                       as.numeric(str_extract(r1smokev, "^\\d+")))
  )

#收入状态整理
dat <- dat %>%
  mutate(
    hh1itot = ifelse(str_starts(hh1itot, "\\."), NA, as.numeric(hh1itot)),
    hh1itot = ntile(hh1itot, 4) - 1
  )

#心脏病史整理
dat <- dat %>%
  mutate(
    r1hearte = ifelse(str_starts(r1hearte, "\\."), NA, 
                      as.numeric(str_extract(r1hearte, "^\\d+")))
  )

#抑郁状态整理
dat <- dat %>%
  mutate(
    r1cesd10 = ifelse(str_starts(r1cesd10, "\\."), NA, as.numeric(r1cesd10)),
    r1cesd10 = ifelse(is.na(r1cesd10), NA,
                        ifelse(r1cesd10 >= 10, 1, 0))
  )

#高血压药整理
dat <- dat %>%
  mutate(
    r1rxhibp_c = ifelse(str_starts(r1rxhibp_c, "\\."), NA, 
                      as.numeric(str_extract(r1rxhibp_c, "^\\d+")))
  )

#糖尿病药整理
dat <- dat %>%
  mutate(
    r1rxdiab_c = ifelse(str_starts(r1rxdiab_c, "\\."), NA, 
                        as.numeric(str_extract(r1rxdiab_c, "^\\d+")))
  )

#血脂药整理
dat <- dat %>%
  mutate(
    r1rxdyslip_c = ifelse(str_starts(r1rxdyslip_c, "\\."), NA, 
                        as.numeric(str_extract(r1rxdyslip_c, "^\\d+")))
  )


#计算eGFR
dat <- dat %>%
  mutate(
    eGFR = case_when(
      ragender == 2 & newcrea <= 0.7 ~
        142 * (newcrea / 0.7)^(-0.241) * (0.9938)^r1agey * 1.012,   # ← 加 1.012
      ragender == 2 & newcrea > 0.7 ~
        142 * (newcrea / 0.7)^(-1.200) * (0.9938)^r1agey * 1.012,   # ← 加 1.012
      ragender == 1 & newcrea <= 0.9 ~
        142 * (newcrea / 0.9)^(-0.302) * (0.9938)^r1agey,
      ragender == 1 & newcrea > 0.9 ~
        142 * (newcrea / 0.9)^(-1.200) * (0.9938)^r1agey
    )
  )

#计算eGFR_C（CKD-EPI 2021 Crea+CysC联合方程，敏感性分析用）
#注意：请确认 CHARLS 中 CysC 变量名（单位 mg/L）
dat <- dat %>%
  mutate(
    cystatinc = as.numeric(ifelse(str_starts(as.character(cystatinc), "\\."), 
                               NA, cystatinc)),
    # 性别参数
    .kappa    = if_else(ragender == 2, 0.7, 0.9),
    .alpha    = if_else(ragender == 2, -0.219, -0.144),
    .sex_coef = if_else(ragender == 2, 0.963, 1),
    .scr_r    = newcrea / .kappa,
    .scys_r   = cystatinc / 0.8,
    
    eGFR_C = 135 *
      pmin(.scr_r,  1)^.alpha *
      pmax(.scr_r,  1)^(-0.544) *
      pmin(.scys_r, 1)^(-0.323) *
      pmax(.scys_r, 1)^(-0.778) *
      (0.9961)^r1agey *
      .sex_coef
  ) %>%
  select(-.kappa, -.alpha, -.sex_coef, -.scr_r, -.scys_r)

# CysC缺失率报告（预期约25-30%）
cat(sprintf("CysC missing: %.1f%%; eGFR_C missing: %.1f%%\n",
            mean(is.na(dat$cystatinc)) * 100,
            mean(is.na(dat$eGFR_C)) * 100))


##缺失处理，NA>20%删去
confounders <- c("r1mstat", "raeducl", "h1rural", "hh1itot", 
                 "r1smokev", "r1drinkev", "r1mdact_c", 
                 "r1hearte", "r1cesd10", "r1rxhibp_c", "r1rxdiab_c", "r1rxdyslip_c")
# 报告缺失率
na_report <- dat %>%
  summarise(across(all_of(confounders), 
                   ~ mean(is.na(.)) * 100)) %>%
  pivot_longer(everything(), 
               names_to = "variable", 
               values_to = "na_pct") %>%
  mutate(na_pct = round(na_pct, 2),
         action = ifelse(na_pct > 20, "建议剔除", "保留"))
print(na_report)
vars_to_remove <- na_report %>% filter(action == "建议剔除") %>% pull(variable)
dat <- dat %>% select(-all_of(vars_to_remove))


##时间、结局变量定义
# Step 1: 各波次iwt均值填补NA
w_mean <- c(
  mean(dat$r1iwt, na.rm = TRUE),
  mean(dat$r2iwt, na.rm = TRUE),
  mean(dat$r3iwt, na.rm = TRUE),
  mean(dat$r4iwt, na.rm = TRUE)
)

dat <- dat %>%
  mutate(
    r1iwt_f = coalesce(r1iwt, w_mean[1]),
    r2iwt_f = coalesce(r2iwt, w_mean[2]),
    r3iwt_f = coalesce(r3iwt, w_mean[3]),
    r4iwt_f = coalesce(r4iwt, w_mean[4])
  )

# Step 2: 逐行定义结局变量
dat <- dat %>%
  rowwise() %>%
  mutate(
    .iwt  = list(c(r1iwt_f, r2iwt_f, r3iwt_f, r4iwt_f)),
    .stat = list(c(r1iwstat, r2iwstat, r3iwstat, r4iwstat)),
    
    stroke_wave = {
      s <- which(c(r2stroke, r3stroke, r4stroke) == "1.Yes")
      if (length(s) > 0) s[1] + 1L else NA_integer_
    },
    death_wave = {
      d <- which(c(r2iwstat, r3iwstat, r4iwstat) == "5.nr, died this wv")
      if (length(d) > 0) d[1] + 1L else NA_integer_
    },
    prev_alive_stroke = if (!is.na(stroke_wave) & stroke_wave > 1) {
      a <- which(.stat[1:(stroke_wave - 1)] == "1.resp, alive")
      if (length(a) > 0) max(a) else NA_integer_
    } else NA_integer_,
    prev_alive_death = if (!is.na(death_wave) & death_wave > 1) {
      a <- which(.stat[1:(death_wave - 1)] == "1.resp, alive")
      if (length(a) > 0) max(a) else NA_integer_
    } else NA_integer_,
    last_alive = {
      a <- which(.stat == "1.resp, alive")
      if (length(a) > 0) max(a) else NA_integer_
    },
    
    stroke_first = !is.na(stroke_wave) & (is.na(death_wave) | stroke_wave < death_wave),
    
    # 先安全提取时间点，NA索引返回NA
    t_stroke    = if (!is.na(stroke_wave)) .iwt[[stroke_wave]] else NA_real_,
    t_prev_s    = if (!is.na(prev_alive_stroke)) .iwt[[prev_alive_stroke]] else NA_real_,
    t_death     = if (!is.na(death_wave)) .iwt[[death_wave]] else NA_real_,
    t_prev_d    = if (!is.na(prev_alive_death)) .iwt[[prev_alive_death]] else NA_real_,
    t_last      = if (!is.na(last_alive)) .iwt[[last_alive]] else NA_real_,
    
    new_stroke = case_when(
      stroke_first                        ~ 1L,
      !is.na(death_wave) & !stroke_first  ~ 2L,
      TRUE                                ~ 0L
    ),
    time = case_when(
      stroke_first                        ~ (t_stroke + t_prev_s) / 2,
      !is.na(death_wave) & !stroke_first  ~ (t_death + t_prev_d) / 2,
      TRUE                                ~ t_last
    )
  ) %>%
  select(-.iwt, -.stat, -stroke_wave, -death_wave,
         -prev_alive_stroke, -prev_alive_death, -last_alive, -stroke_first,
         -t_stroke, -t_prev_s, -t_death, -t_prev_d, -t_last) %>%
  ungroup()

#确定时间变量
dat <- dat %>%
  mutate(time_yr = time - r1iwt)


# ---- 统一清洗：剔除所有暴露的生理学不可能值 ----
n_before <- nrow(dat)

dat <- dat %>%
  filter(
    systo  >= 70  & systo  <= 250,     # mmHg
    newglu   >= 40  & newglu   <= 500,     # mg/dL
    newldl  >= 20  & newldl  <= 400,     # mg/dL
    eGFR >= 15  & eGFR <= 200,     # mL/min/1.73m²
    r1mbmi  >= 12  & r1mbmi  <= 60       # kg/m²
  )

cat(sprintf("Excluded %d participants with implausible metabolic values; n: %d → %d\n",
            n_before - nrow(dat), n_before, nrow(dat)))

# 验证各暴露清洗后范围
sapply(dat[, c("systo", "newglu", "newldl", "eGFR", "r1mbmi")], summary)



##总体简化
charls <- dat %>%
  select(ID, r1agey, ragender, r1mstat, raeducl, h1rural, r1smokev, r1drinkev, 
         r1hearte, r1cesd10,
         r1rxhibp_c, r1rxdiab_c, r1rxdyslip_c,    # ← 新增：药物变量
         systo, r1mbmi, newglu, eGFR, eGFR_C, newldl,   # ← 新增 eGFR_C
         new_stroke, time_yr) %>%
  rename(
    age        = r1agey,
    gender     = ragender,
    marriage   = r1mstat,
    education  = raeducl,
    rural      = h1rural,
    smoke      = r1smokev,
    drink      = r1drinkev,
    bpmed      = r1rxhibp_c,        # ← 新增
    diamed     = r1rxdiab_c,        # ← 新增
    lipidmed   = r1rxdyslip_c,      # ← 新增
    heartd     = r1hearte,
    depression = r1cesd10,
    SBP        = systo,
    BMI        = r1mbmi,
    FG         = newglu,
    LDL        = newldl
  )




##存储
write.csv(charls,'charls_dat.csv')




