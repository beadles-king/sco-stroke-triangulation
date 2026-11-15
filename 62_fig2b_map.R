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
# Figure 2B: 全中风代谢归因聚类地图（单幅）
# 追加到 05_paf_attribution.R 末尾运行，依赖 main_res（第 296 行 run_cause 的返回值）
#
# 注意：draw_dual_heatmap() 内部把"中国所在簇"重贴为 C1，但没有把改过的
#       标签写回 res_paf$clusters。因此这里必须重做一次同样的重贴，
#       否则地图配色会与热图注释条相反。
###############################################################################

library(sf)
library(rnaturalearth)
library(ggrepel)
library(geosphere)  # 大圆距离
# vegan 不再需要：Mantel 检验不在本脚本内运行，见下方【M1】

# ---- 出图: 一律走 ragg -------------------------------------------------------
# cairo_pdf / pdf() 会把每个字形的步进四舍五入到整数 point (0.353 mm) 或不嵌
# 字体, 前者导致同一单词内字距不匀 (最坏 1 pt), 后者栅格化时字面被替换。
# ragg 在目标分辨率的像素网格上排字, 600 dpi 下取整误差 <= 0.042 mm, 且全部
# 图形出自同一引擎、同一份 Helvetica, 与其余图观感一致。
if (!requireNamespace("ragg", quietly = TRUE)) install.packages("ragg")

save_panel <- function(p, stem, w_mm, h_mm) {
  for (r in c(600, 300)) {
    f <- if (r == 600) paste0(stem, ".tif") else paste0(stem, "_300dpi.tif")
    ragg::agg_tiff(f, width = w_mm, height = h_mm, units = "mm",
                   res = r, compression = "lzw", background = "white",
                   scaling = 1)
    print(p); dev.off()
    cat(sprintf("已保存: %s (%.1f x %.1f mm @ %d dpi)\n", f, w_mm, h_mm, r))
  }
}


# ---------------------------------------------------------------------------
# 0. 十国首都坐标（Mantel 检验用，与附录一致）
# ---------------------------------------------------------------------------
capitals <- data.frame(
  country = c("Belarus", "China", "India", "Iran", "Kazakhstan",
              "Kyrgyzstan", "Pakistan", "Russia", "Tajikistan", "Uzbekistan"),
  iso3    = c("BLR", "CHN", "IND", "IRN", "KAZ",
              "KGZ", "PAK", "RUS", "TJK", "UZB"),
  # 坐标与 07_paf_spatial_tests.R 逐位一致。此前 BLR 与 KAZ 两处与统计脚本
  # 相差约 1 km，虽在 r 的 1e-3 容差内，但同一坐标表出现两个版本本身是隐患。
  lon     = c(27.5667, 116.4074, 77.2090, 51.3890, 71.4491,
              74.5698, 73.0479, 37.6173, 68.7870, 69.2401),
  lat     = c(53.9000, 39.9042, 28.6139, 35.6892, 51.1605,
              42.8746, 33.6844, 55.7558, 38.5598, 41.2995),
  stringsAsFactors = FALSE
)

# ---------------------------------------------------------------------------
# 1. 取回聚类结果并重贴标签（与 draw_dual_heatmap 的逻辑一致）
# ---------------------------------------------------------------------------
clusters <- main_res$res_paf$clusters
k        <- main_res$res_paf$k

if (k == 2 && "China" %in% names(clusters) && clusters["China"] != 1) {
  clusters <- ifelse(clusters == 1, 2, 1)
}

cluster_colors <- setNames(
  c("#E6550D", "#7FBFA0", "#9970AB", "#F4A261", "#8AB17D")[1:k],
  paste0("C", 1:k)
)

# ---------------------------------------------------------------------------
# 2. Mantel 检验结果：只读不算
#
# 【M1】此前本脚本自行调用 mantel(permutations = 9999)，与
#       07_paf_spatial_tests.R（现为 99999 次、set.seed(2023)）是两次独立运行，
#       图上标注与 Table S12 必然分叉。现改为读同一份 CSV。
#       r 是确定性的，仍就地复算一次作为核对：r 对不上说明距离矩阵的构造
#       与空间检验脚本不同，那是真错误，必须停。
# ---------------------------------------------------------------------------
SPATIAL_CSV <- "spatial_PAF_cluster_test_results.csv"

spatial <- read.csv(SPATIAL_CSV, stringsAsFactors = FALSE)
ref <- spatial[spatial$Outcome == "Stroke (Total)", ]
stopifnot(nrow(ref) == 1)

paf_m  <- main_res$paf_mat
paf_z  <- t(scale(t(paf_m)))                 # 按 RF 行标准化，与主聚类同法
d_feat <- dist(t(paf_z), method = "euclidean")

ord    <- match(labels(d_feat), capitals$country)
d_geo  <- as.dist(distm(as.matrix(capitals[ord, c("lon", "lat")]),
                        fun = distHaversine) / 1000)   # km

r_local <- cor(as.vector(d_feat), as.vector(d_geo), method = "pearson")
if (abs(r_local - ref$Mantel_r) > 1e-3)
  stop(sprintf(paste0("Mantel r 与空间检验脚本不一致：本脚本算出 %.4f，",
                      "CSV 里是 %.3f。r 是确定性的，对不上说明距离矩阵的",
                      "构造不同，先查清再出图。"), r_local, ref$Mantel_r))

cat(sprintf("\nMantel r = %.3f（与 CSV 一致），p = %.5f（取自 CSV）\n",
            r_local, ref$Mantel_P))

mantel_lab <- paste0("Mantel r=", dec(sprintf("%.2f", ref$Mantel_r)),
                     ", p=", fmt_p(ref$Mantel_P))

# ---------------------------------------------------------------------------
# 2b. 偏 Mantel（控制 GBD region/super-region 层级）：同样只读不算
#
# 【M2】仅标注无条件的 Mantel 已构成选择性呈现：层级对 PAF 剖面的解释力
#       强于地理，且总卒中在控制层级后关联降至不显著。图上必须两值并列。
#       数据源为 C2 脚本写出的 C2_hierarchy_table.csv，与新增附录表同源。
# ---------------------------------------------------------------------------
C2_CSV <- file.path(RAWDAT, "大改/C2_hierarchy_table.csv")
h2 <- read.csv(C2_CSV, stringsAsFactors = FALSE)
h2 <- h2[h2$outcome == "Stroke (Total)", ]
stopifnot(all(c("geo~PAF", "geo~PAF | hier") %in% h2$test))

# 同一统计量在两份 CSV 中的 r 必须一致；不一致说明两脚本的输入矩阵不同。
r_c2 <- h2$statistic[h2$test == "geo~PAF"]
if (abs(r_c2 - ref$Mantel_r) > 1e-3)
  stop(sprintf(paste0("C2 表与空间检验表的 geo~PAF 不一致：%.4f vs %.3f。",
                      "两脚本读入的 PAF 矩阵或距离定义不同，先查清再出图。"),
               r_c2, ref$Mantel_r))

pr <- h2[h2$test == "geo~PAF | hier", ]
partial_lab <- paste0("adjusted for estimation hierarchy: r=", dec(sprintf("%.2f", pr$statistic)),
                      ", p=", fmt_p(pr$p))

annot_lab <- paste(mantel_lab, partial_lab, sep = "\n")
cat("图上标注：\n", annot_lab, "\n")
cat("请确认第一行与 Table S12、第二行与新增的层级检验表完全一致。\n")

# ---------------------------------------------------------------------------
# 3. 底图
# ---------------------------------------------------------------------------
world <- ne_countries(scale = "medium", returnclass = "sf")

cl_df <- data.frame(iso3    = capitals$iso3,
                    country = capitals$country,
                    stringsAsFactors = FALSE)
cl_df$cluster <- paste0("C", clusters[cl_df$country])

world <- world %>%
  left_join(cl_df, by = c("iso_a3" = "iso3"))

lab_pts <- world %>%
  filter(!is.na(cluster)) %>%
  st_point_on_surface() %>%
  mutate(lon = st_coordinates(.)[, 1],
         lat = st_coordinates(.)[, 2]) %>%
  st_drop_geometry()

# ---------------------------------------------------------------------------
# 4. 绘图
# ---------------------------------------------------------------------------
p_map <- ggplot(world) +
  geom_sf(aes(fill = cluster), colour = "grey35", linewidth = 0.15) +
  scale_fill_manual(values = cluster_colors, na.value = "grey88",
                    breaks = paste0("C", 1:k),
                    name = "Risk-structure\ncluster") +
  geom_text_repel(data = lab_pts,
                  aes(x = lon, y = lat, label = iso_a3),
                  size = 2.5, fontface = "bold", colour = "grey10",
                  bg.color = "white", bg.r = 0.12,
                  min.segment.length = 0.2, segment.size = 0.25,
                  segment.colour = "grey30", max.overlaps = Inf) +
  annotate("text", x = 22, y = 3, hjust = 0, vjust = 0,
           label = annot_lab, size = 2.6, colour = "grey10",
           lineheight = 1.05) +
  coord_sf(xlim = c(20, 150), ylim = c(0, 78), expand = FALSE) +
  theme_void(base_size = 8) +
  theme(
    panel.background = element_rect(fill = "#DCEBF5", colour = NA),
    panel.border     = element_rect(fill = NA, colour = "grey35", linewidth = 0.3),
    legend.title     = element_text(size = 7.5, face = "bold"),
    legend.text      = element_text(size = 7),
    legend.key.size  = unit(9, "pt"),
    legend.position  = "right",
    plot.margin      = margin(2, 2, 2, 2)
  ) + 
  guides(fill = "none")

# Figure 2 的面板 B; 面板 A 由 05_paf_attribution.R 出, 两者在拼图阶段合成 Figure 2
save_panel(p_map, "Fig2B_stroke_cluster_map", w_mm = 180, h_mm = 105)