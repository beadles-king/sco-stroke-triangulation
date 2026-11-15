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
# Figure S10: 亚型代谢归因聚类地图（IS / ICH / SAH）
#
# 与 62_fig2b_map.R 同源，逐项对齐其发表标准写法：
#   - 出图一律 ragg::agg_tiff，600 dpi 与 300 dpi 各渲染一次
#     （cairo_pdf / pdf() / png() / ggsave 全部禁用，见 62_fig2b_map.R 顶部说明）
#   - 海洋 #DCEBF5，非 SCO 国家 grey88，国界 grey35 linewidth 0.15
#   - 国名用 st_point_on_surface + geom_text_repel，白色描边，size 2.5
#   - Mantel 标注置于面板左下角，size 2.6，走 dec() / fmt_p()
#   - theme_void(base_size = 8)，面板边框 grey35 linewidth 0.3
#   - 面板字母（A/B/C）不在此处加，由 Python 拼图阶段统一加
#
# 依赖：05_paf_attribution.R 已运行，全局有 is_res / ich_res / sah_res，以及 dec() / fmt_p()
#
# Mantel 数值不在本脚本内重算，一律读 07_paf_spatial_tests.R 写出的
# spatial_PAF_cluster_test_results.csv。原因：mantel() 是置换检验，
# 图与表若各算各的必然分叉（2026-08 已发生一次：表里 IS p=0·0011，
# 重跑得 0·0002）。单一数据源是唯一可靠的做法。
# 前提：07_paf_spatial_tests.R 已加 set.seed()，否则那份 CSV 本身也不可复现。
#
# 三个面板各自存盘，拼图在 Python 端完成。
###############################################################################
setwd(file.path(RAWDAT, "PAF"))
library(sf)
library(rnaturalearth)
library(ggrepel)
# vegan 不再需要：Mantel 检验不在本脚本内运行
library(geosphere)  # 大圆距离
library(dplyr)
library(ggplot2)
library(cowplot)    # 仅用于抽取共享图例

if (!requireNamespace("ragg", quietly = TRUE)) install.packages("ragg")

# ---- 版面常数 --------------------------------------------------------------
PANEL_W <- 128     # 单面板宽 mm（三面板竖排，右侧留给共享图例）
PANEL_H <-  75     # 单面板高 mm；3*75 + 2*3 = 231 mm，在 240 上限内
XLIM    <- c(20, 150)
YLIM    <- c(0, 78)
OUT_DIR <- "."

save_panel <- function(p, stem, w_mm, h_mm) {
  for (r in c(600, 300)) {
    f <- file.path(OUT_DIR,
                   if (r == 600) paste0(stem, ".tif") else paste0(stem, "_300dpi.tif"))
    ragg::agg_tiff(f, width = w_mm, height = h_mm, units = "mm",
                   res = r, compression = "lzw", background = "white",
                   scaling = 1)
    print(p); dev.off()
    cat(sprintf("已保存: %s (%.1f x %.1f mm @ %d dpi)\n", f, w_mm, h_mm, r))
  }
}

# ---------------------------------------------------------------------------
# 0. 十国首都坐标（Mantel 检验用，与 62_fig2b_map.R、附录一致）
# ---------------------------------------------------------------------------
capitals <- data.frame(
  country = c("Belarus", "China", "India", "Iran", "Kazakhstan",
              "Kyrgyzstan", "Pakistan", "Russia", "Tajikistan", "Uzbekistan"),
  iso3    = c("BLR", "CHN", "IND", "IRN", "KAZ",
              "KGZ", "PAK", "RUS", "TJK", "UZB"),
  # 坐标与 07_paf_spatial_tests.R 逐位一致（BLR、KAZ 两处此前有约 1 km 的分叉）
  lon     = c(27.5667, 116.4074, 77.2090, 51.3890, 71.4491,
              74.5698, 73.0479, 37.6173, 68.7870, 69.2401),
  lat     = c(53.9000, 39.9042, 28.6139, 35.6892, 51.1605,
              42.8746, 33.6844, 55.7558, 38.5598, 41.2995),
  stringsAsFactors = FALSE
)

cluster_colors <- c("C1" = "#E6550D", "C2" = "#7FBFA0")

# ---------------------------------------------------------------------------
# 1. 读入空间检验结果 —— 图上标注与 Table S12 的唯一数据源
# ---------------------------------------------------------------------------
SPATIAL_CSV <- "spatial_PAF_cluster_test_results.csv"

spatial <- read.csv(SPATIAL_CSV, stringsAsFactors = FALSE)
row_of <- c(IS  = "Ischemic stroke",
            ICH = "Intracerebral hemorrhage",
            SAH = "Subarachnoid hemorrhage")
stopifnot(all(row_of %in% spatial$Outcome))

cat("===== 本图标注所用的空间检验结果 =====\n")
print(spatial[, c("Outcome", "Mantel_r", "Mantel_P")], row.names = FALSE)
cat("提醒：这几个 p 必须与 Table S12 当前版本一致。若 07_paf_spatial_tests.R\n")
cat("      刚重跑过，请先把 Table S12 更新到同一次运行的结果。\n\n")

# ---------------------------------------------------------------------------
# 1b. 偏 Mantel（控制 GBD region/super-region 层级），同样只读不算
#     C2 表用的结局名与空间检验表不同，此处显式映射，勿改。
# ---------------------------------------------------------------------------
C2_CSV <- file.path(RAWDAT, "大改/C2_hierarchy_table.csv")
h2_all <- read.csv(C2_CSV, stringsAsFactors = FALSE)
c2_row_of <- c(IS = "Ischemic stroke", ICH = "ICH", SAH = "SAH")
stopifnot(all(c2_row_of %in% h2_all$outcome))

cat("===== 偏 Mantel（控制 GBD 层级）=====\n")
print(h2_all[h2_all$test %in% c("geo~PAF", "hier~PAF", "geo~PAF | hier"),
             c("outcome", "test", "statistic", "p")], row.names = FALSE)
cat("提醒：这些值必须与新增的层级检验附录表一致。\n\n")

# ---------------------------------------------------------------------------
# 2. 底图（一次读入，三个面板共用）
# ---------------------------------------------------------------------------
world <- ne_countries(scale = "medium", returnclass = "sf")

# ---------------------------------------------------------------------------
# 3. 单面板绘图
# ---------------------------------------------------------------------------
# 含 China 的簇一律重贴为 C1，与 draw_dual_heatmap() 的逻辑一致；
# 62_fig2b_map.R 的注释已指出 res_paf$clusters 里存的是重贴前的标签。
align_c1 <- function(clusters) {
  if ("China" %in% names(clusters) && clusters["China"] != 1) {
    out <- ifelse(clusters == 1, 2, 1)
    names(out) <- names(clusters)
    return(out)
  }
  clusters
}

# Mantel 只读不算。r 是确定性的，因此仍就地复算一次作为核对：
# 若 r 对不上，说明距离矩阵或 PAF 矩阵的构造与空间检验脚本不同，
# 那是真错误，必须停下。p 是置换得出的，只取 CSV 的值，不重算。
mantel_of <- function(paf_m, key) {
  paf_z  <- t(scale(t(paf_m)))                 # 按 RF 行标准化，与主聚类同法
  d_feat <- dist(t(paf_z), method = "euclidean")
  ord    <- match(labels(d_feat), capitals$country)
  d_geo  <- as.dist(distm(as.matrix(capitals[ord, c("lon", "lat")]),
                          fun = distHaversine) / 1000)
  r_local <- suppressWarnings(
    cor(as.vector(d_feat), as.vector(d_geo), method = "pearson"))
  
  ref <- spatial[spatial$Outcome == row_of[[key]], ]
  if (abs(r_local - ref$Mantel_r) > 1e-3)
    stop(sprintf(paste0("%s 的 Mantel r 与空间检验脚本不一致：",
                        "本脚本算出 %.4f，CSV 里是 %.3f。",
                        "r 是确定性的，对不上说明距离矩阵的构造不同，先查清再出图。"),
                 key, r_local, ref$Mantel_r))
  # 同一统计量在两份 CSV 中的 r 必须一致
  h2 <- h2_all[h2_all$outcome == c2_row_of[[key]], ]
  r_c2 <- h2$statistic[h2$test == "geo~PAF"]
  if (abs(r_c2 - ref$Mantel_r) > 1e-3)
    stop(sprintf(paste0("%s 的 geo~PAF 在两份 CSV 中不一致：%.4f vs %.3f。",
                        "两脚本读入的 PAF 矩阵或距离定义不同，先查清再出图。"),
                 key, r_c2, ref$Mantel_r))
  pr <- h2[h2$test == "geo~PAF | hier", ]
  cat(sprintf("%-4s r=%.3f p=%.4f（取自 CSV）| 控制层级后 r=%.3f p=%.4f\n",
              key, r_local, ref$Mantel_P, pr$statistic, pr$p))
  list(statistic = ref$Mantel_r, signif = ref$Mantel_P,
       p_stat = pr$statistic, p_signif = pr$p)
}

make_map <- function(res_obj, key, stem, show_legend = FALSE) {
  clusters <- align_c1(res_obj$res_paf$clusters)
  mt       <- mantel_of(res_obj$paf_mat, key)
  
  mantel_lab <- paste(
    paste0("Mantel r=", dec(sprintf("%.2f", mt$statistic)),
           ", p=", fmt_p(mt$signif)),
    paste0("adjusted for estimation hierarchy: r=", dec(sprintf("%.2f", mt$p_stat)),
           ", p=", fmt_p(mt$p_signif)),
    sep = "\n")
  
  cl_df <- data.frame(iso3    = capitals$iso3,
                      country = capitals$country,
                      stringsAsFactors = FALSE)
  cl_df$cluster <- paste0("C", clusters[cl_df$country])
  
  w <- world %>% left_join(cl_df, by = c("iso_a3" = "iso3"))
  
  lab_pts <- w %>%
    filter(!is.na(cluster)) %>%
    st_point_on_surface() %>%
    mutate(lon = st_coordinates(.)[, 1],
           lat = st_coordinates(.)[, 2]) %>%
    st_drop_geometry()
  
  p <- ggplot(w) +
    geom_sf(aes(fill = cluster), colour = "grey35", linewidth = 0.15) +
    scale_fill_manual(values = cluster_colors, na.value = "grey88",
                      breaks = c("C1", "C2"),
                      name = "Risk-structure\ncluster") +
    geom_text_repel(data = lab_pts,
                    aes(x = lon, y = lat, label = iso_a3),
                    size = 2.5, fontface = "bold", colour = "grey10",
                    bg.color = "white", bg.r = 0.12,
                    min.segment.length = 0.2, segment.size = 0.25,
                    segment.colour = "grey30", max.overlaps = Inf,
                    seed = 42) +
    annotate("text", x = XLIM[1] + 2, y = YLIM[1] + 3, hjust = 0, vjust = 0,
             label = mantel_lab, size = 2.6, colour = "grey10",
             lineheight = 1.05) +
    coord_sf(xlim = XLIM, ylim = YLIM, expand = FALSE) +
    theme_void(base_size = 8) +
    theme(
      panel.background = element_rect(fill = "#DCEBF5", colour = NA),
      panel.border     = element_rect(fill = NA, colour = "grey35", linewidth = 0.3),
      legend.title     = element_text(size = 7.5, face = "bold"),
      legend.text      = element_text(size = 7),
      legend.key.size  = unit(9, "pt"),
      legend.position  = if (show_legend) "right" else "none",
      plot.margin      = margin(2, 2, 2, 2)
    )
  
  if (!show_legend) {
    save_panel(p + guides(fill = "none"), stem, PANEL_W, PANEL_H)
  }
  invisible(p)
}

# ---------------------------------------------------------------------------
# 4. 三个面板
#    面板顺序与 Figure S9 一致：A 缺血性，B 脑出血，C 蛛网膜下腔出血
#    面板字母与亚型标题由 Python 拼图阶段统一加，此处不画标题
# ---------------------------------------------------------------------------
p_is  <- make_map(is_res,  "IS",  "FigS10_A_ischaemic")
p_ich <- make_map(ich_res, "ICH", "FigS10_B_ich")
p_sah <- make_map(sah_res, "SAH", "FigS10_C_sah")

# ---------------------------------------------------------------------------
# 5. 独立共享图例（与 Figure S15 的做法一致）
# ---------------------------------------------------------------------------
leg <- cowplot::get_legend(
  p_is + theme(legend.position = "right",
               legend.box.margin = margin(0, 0, 0, 0))
)
save_panel(cowplot::ggdraw(leg), "FigS10_legend", 34, 26)

# ---------------------------------------------------------------------------
# 6. 成员名单核对：应与 Table S11 一致
#    IS  C1 = China, India, Pakistan, Tajikistan, Uzbekistan
#    ICH C1 = China, India, Pakistan
#    SAH C1 = China, India, Pakistan
# ---------------------------------------------------------------------------
cat("\n===== 聚类成员核对（对照 Table S11）=====\n")
expected <- list(
  IS  = c("China", "India", "Pakistan", "Tajikistan", "Uzbekistan"),
  ICH = c("China", "India", "Pakistan"),
  SAH = c("China", "India", "Pakistan")
)
for (nm in names(expected)) {
  cl <- align_c1(switch(nm, IS = is_res, ICH = ich_res, SAH = sah_res)$res_paf$clusters)
  c1 <- sort(names(cl)[cl == 1])
  ok <- identical(c1, sort(expected[[nm]]))
  cat(sprintf("%-4s C1 = {%s}  %s\n", nm, paste(c1, collapse = ", "),
              if (ok) "与 Table S11 一致" else "*** 与 Table S11 不一致 ***"))
  if (!ok) stop(sprintf("%s 的聚类成员与 Table S11 对不上，先查清再出图。", nm))
}

