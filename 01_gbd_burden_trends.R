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
# 方案B改进版: 时间序列特征提取（含整体负担水平）+ 层次聚类
# 特征: slope + pct_change + volatility + acceleration + mean_level × 3 measures
# 运行前请确保 general_year 已按原代码生成
#
# 本版改动仅限绘图部分（分析逻辑未动），已删去 PCA 图：
#   - 出图尺寸改为 180 mm 宽，字号 7-8 pt，投稿时无需缩放
#   - 分面标题只保留首行，x 轴标题只保留末行
#   - 簇内国名按字母序排列
#   - 树状图去标题、枝条改中性色、虚线框标注簇编号
#   - 配色改为色觉友好组合
###############################################################################
setwd(RAWDAT)
library(tidyverse)
library(cluster)
library(dendextend)
library(cowplot)

# ---- 出图: 一律走 ragg ------------------------------------------------------
# 原先主图/亚型图用 ggsave 出 png + pdf, 树状图用 png()。三者在 macOS 上分别
# 走 quartz / grDevices pdf / quartz 三套排版, 且 pdf() 不嵌字体, 栅格化时会被
# 替换成 Nimbus Sans, 与其余图 (真 Helvetica) 不是同一副字面。统一用 ragg 后,
# 全部图形出自同一引擎、同一份字体, 与 S4-S9 / S18-S19 观感一致。
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

# ============================================================================
# 0. 数据加载与预处理（沿用你的原始逻辑）
# ============================================================================
general <- read_csv('general.csv')
general_globalrate <- read_csv('general_globalrate.csv')
general <- rbind(general_globalrate, general)

general_year <- general %>%
  filter(age_name == 'Age-standardized',
         year %in% 2001:2023,
         sex_name == 'Both',
         metric_name == 'Rate') %>%
  mutate(location_name = case_when(
    location_name == "Russian Federation" ~ "Russia",
    location_name == "Iran (Islamic Republic of)" ~ "Iran",
    TRUE ~ location_name
  ))

# ============================================================================
# 通用工具
# ============================================================================
build_wide_matrix <- function(data, cause, measure) {
  data %>%
    filter(cause_name == cause, measure_name == measure,
           location_name != "Global") %>%
    select(location_name, year, val) %>%
    pivot_wider(names_from = year, values_from = val) %>%
    column_to_rownames("location_name")
}

# ============================================================================
# 特征提取 + 聚类（含 mean_level）
# ============================================================================
feature_based_cluster <- function(data, cause, k_range = 2:4, k_force = NULL) {
  measures <- c("ASIR" = "Incidence",
                "ASMR" = "Deaths",
                "ASDR" = "DALYs (Disability-Adjusted Life Years)")
  
  feat_all <- list()
  
  for (mname in names(measures)) {
    wide <- build_wide_matrix(data, cause, measures[mname])
    mat <- as.matrix(wide)
    years <- as.numeric(colnames(mat))
    
    feats <- t(apply(mat, 1, function(ts) {
      # 1. 线性斜率（标准化：每年变化占基线%）
      slope <- coef(lm(ts ~ years))[2] / ts[1] * 100
      # 2. 总变化率
      pct_change <- (ts[length(ts)] - ts[1]) / ts[1] * 100
      # 3. 波动性
      diffs <- diff(ts)
      volatility <- sd(diffs) / (abs(mean(diffs)) + 1e-10)
      # 4. 趋势加速度（后半 vs 前半斜率变化）
      mid <- length(ts) %/% 2
      slope1 <- coef(lm(ts[1:mid] ~ years[1:mid]))[2]
      slope2 <- coef(lm(ts[(mid+1):length(ts)] ~ years[(mid+1):length(ts)]))[2]
      accel <- (slope2 - slope1) / (abs(slope1) + 1e-10)
      # 5. 整体负担水平（23年均值）
      mean_level <- mean(ts)
      
      c(slope = slope, 
        volatility = volatility, 
        pct_change = pct_change, 
        acceleration = accel, 
        mean_level = mean_level * 2)
      
    }))
    
    colnames(feats) <- paste0(mname, "_", colnames(feats))
    feat_all[[mname]] <- feats
  }
  
  # 合并: 10国 * 15特征
  feat_matrix <- do.call(cbind, feat_all)
  feat_scaled <- scale(feat_matrix)
  
  cat("\n特征矩阵 (标准化前, 保留2位小数):\n")
  print(round(feat_matrix, 2))
  
  d <- dist(feat_scaled, method = "euclidean")
  hc <- hclust(d, method = "ward.D2")
  
  # 轮廓系数
  sil <- sapply(k_range, function(ki) {
    cl <- cutree(hc, k = ki)
    mean(silhouette(cl, d)[, "sil_width"])
  })
  names(sil) <- paste0("k=", k_range)
  cat("\n轮廓系数:", paste(names(sil), round(sil, 3), sep = ": ", collapse = "  "), "\n")
  
  best_k <- if (!is.null(k_force)) k_force else k_range[which.max(sil)]
  #best_k <- 3
  cat("使用 k =", best_k, "\n\n")
  
  # 打印各k结果
  for (ki in k_range) {
    clusters <- cutree(hc, k = ki)
    sizes <- table(clusters)
    cat("k =", ki, " (组大小:", paste(sizes, collapse = "-"), ")\n")
    for (ci in sort(unique(clusters))) {
      cat("  Cluster", ci, ":", paste(names(clusters[clusters == ci]), collapse = ", "), "\n")
    }
    cat("\n")
  }
  
  clusters_final <- cutree(hc, k = best_k)
  list(hc = hc, clusters = clusters_final, k = best_k,
       dist = d, feat_matrix = feat_matrix, feat_scaled = feat_scaled)
}

# ============================================================================
# 树状图
# ============================================================================
plot_dend <- function(hc_obj, k, filename,
                      width_mm = 180, height_mm = 92) {   # 85 → 92，补偿加高的下边距
  dend <- as.dendrogram(hc_obj) %>%
    set("branches_lwd", 1.1) %>%
    set("labels_cex", 1)
  
  # 虚线框按树状图左→右顺序绘制，簇编号需据此排序，不能直接写 1:k
  cl        <- cutree(hc_obj, k = k)
  ord       <- order.dendrogram(dend)
  grp       <- cl[ord]                                        # 左→右每片叶子的簇编号
  lab_order <- unique(grp)
  centres   <- sapply(lab_order, function(g) mean(which(grp == g)))  # 各簇中心的叶片坐标
  
  draw <- function() {
    par(mar = c(8, 4, 0.5, 9.3))                             # 下边距 6 → 8
    plot(dend, main = NULL, ylab = "Height")
    rect.dendrogram(dend, k = k, border = "grey40", lty = 2)
    # 簇编号改用 mtext 置于虚线框下方，按边距行定位，不受国名长度影响
    mtext(paste0("Cluster ", lab_order), side = 1, at = centres,
          line = 6.3, cex = 0.95, col = "grey20")
  }
  stem <- sub("\\.(png|pdf|tif|tiff)$", "", filename)
  for (r in c(600, 300)) {
    f <- if (r == 600) paste0(stem, ".png") else paste0(stem, "_300dpi.png")
    ragg::agg_png(f, width = width_mm, height = height_mm, units = "mm",
                  res = r, pointsize = 8, background = "white", scaling = 1)
    draw(); dev.off()
    cat("已保存:", f, "\n")
  }
}

# ============================================================================
# 分面趋势图
# ============================================================================
plot_faceted_trends <- function(data, cause, clusters,
                                measure_labels = c("ASIR" = "Incidence",
                                                   "ASMR" = "Deaths",
                                                   "ASDR" = "DALYs (Disability-Adjusted Life Years)")) {
  
  # 色觉友好配色：避免同一簇内出现红/绿配对，且不使用低对比度的纯黄
  country_colors <- c(
    "Global" = "black",
    "Belarus" = "#661100", "China" = "#D55E00", "India" = "#117733",
    "Iran" = "#AA4499", "Kazakhstan" = "#E69F00", "Kyrgyzstan" = "#CC6677",
    "Pakistan" = "#999933", "Russia" = "#0072B2",
    "Tajikistan" = "#888888", "Uzbekistan" = "#44AA99"
  )
  
  country_linetypes <- c(
    "Global" = "dashed",
    "Belarus" = "solid", "China" = "solid", "India" = "solid",
    "Iran" = "solid", "Kazakhstan" = "solid", "Kyrgyzstan" = "solid",
    "Pakistan" = "solid", "Russia" = "solid",
    "Tajikistan" = "solid", "Uzbekistan" = "solid"
  )
  
  # 图例顺序: Global最上, 其余按字母序
  country_order <- c("Global", "Belarus", "China", "India", "Iran",
                     "Kazakhstan", "Kyrgyzstan", "Pakistan", "Russia",
                     "Tajikistan", "Uzbekistan")
  
  cluster_df <- data.frame(location_name = names(clusters),
                           cluster = clusters, stringsAsFactors = FALSE)
  
  cluster_labels <- cluster_df %>%
    group_by(cluster) %>%
    summarise(label = paste0("Cluster ", cluster[1], ": ",
                             paste(sort(location_name), collapse = ", ")),
              .groups = "drop") %>%
    mutate(label = str_wrap(label, width = 30))
  
  cluster_df <- cluster_df %>% left_join(cluster_labels, by = "cluster")
  
  plot_list <- list()
  
  for (i in seq_along(measure_labels)) {
    mname <- names(measure_labels)[i]
    mval  <- measure_labels[i]
    
    country_data <- data %>%
      filter(cause_name == cause, measure_name == mval,
             location_name != "Global") %>%
      select(location_name, year, val) %>%
      left_join(cluster_df, by = "location_name")
    
    global_base <- data %>%
      filter(cause_name == cause, measure_name == mval,
             location_name == "Global") %>%
      select(year, val)
    
    cluster_levs <- cluster_labels$label[order(cluster_labels$cluster)]
    country_data$label <- factor(country_data$label, levels = cluster_levs)
    
    global_data <- do.call(rbind, lapply(cluster_levs, function(lev) {
      global_base %>% mutate(label = lev)
    }))
    global_data$label <- factor(global_data$label, levels = cluster_levs)
    
    ylab_text <- paste0(mname, ", per 100 000")
    
    p <- ggplot() +
      geom_line(data = global_data %>% mutate(location_name = "Global"),
                aes(x = year, y = val, color = location_name,
                    linetype = location_name),
                linewidth = 0.5, alpha = 0.9) +
      geom_line(data = country_data,
                aes(x = year, y = val, color = location_name,
                    linetype = location_name),
                linewidth = 0.5) +
      facet_wrap(~ label, nrow = 1) +
      scale_color_manual(values = country_colors, breaks = country_order,
                         limits = country_order) +
      scale_linetype_manual(values = country_linetypes, breaks = country_order,
                            limits = country_order) +
      labs(y = ylab_text, x = "Year", color = "", linetype = "") +
      theme_minimal(base_size = 8) +
      theme(
        text = element_text(size = 8),
        strip.text = element_text(size = 7.5, face = "bold"),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.key.width = unit(10, "pt"),
        panel.grid.minor = element_blank(),
        plot.margin = margin(2, 2, 2, 2)
      ) +
      scale_y_continuous(n.breaks = 4) +
      guides(color = guide_legend(ncol = 1, override.aes = list(linewidth = 0.7)),
             linetype = guide_legend(ncol = 1))
    
    # x 轴文字与标题只在末行显示
    if (i < length(measure_labels)) {
      p <- p + theme(axis.text.x = element_blank(),
                     axis.ticks.x = element_blank(),
                     axis.title.x = element_blank())
    }
    # 分面标题只在首行显示
    if (i > 1) {
      p <- p + theme(strip.text = element_blank())
    }
    
    plot_list[[mname]] <- p
  }
  
  # 首行多出分面标题、末行多出 x 轴文字与标题，用 rel_heights 补偿使绘图区等高
  combined <- plot_grid(
    plot_list[["ASIR"]] + theme(legend.position = "none"),
    plot_list[["ASMR"]] + theme(legend.position = "none"),
    plot_list[["ASDR"]] + theme(legend.position = "none"),
    ncol = 1, align = "v", axis = "lr",
    rel_heights = c(1.28, 1, 1.19)
  )
  
  legend <- get_legend(
    plot_list[["ASIR"]] +
      theme(legend.position = "right",
            legend.spacing.y = unit(2, "pt"),
            legend.key.height = unit(10, "pt"),
            legend.title = element_blank())
  )
  
  plot_grid(combined, legend, ncol = 2, rel_widths = c(1, 0.17))
}

setwd(file.path(RAWDAT, "总体、年龄、性别负担/总体"))
# ============================================================================
# 运行：Stroke 主分析
# ============================================================================
cat("\n========== Stroke 特征聚类（含 mean_level）==========\n")
res_stroke <- feature_based_cluster(general_year, "Stroke")

# 树状图
plot_dend(res_stroke$hc, res_stroke$k, "stroke_feature_dendrogram.png")

# 分面趋势图
fig_stroke <- plot_faceted_trends(general_year, "Stroke", res_stroke$clusters)
save_panel(fig_stroke, "Figure1_stroke_clustered", w_mm = 180, h_mm = 130)


# ============================================================================
# 运行：三个亚型
# ============================================================================
subtypes <- c("Ischemic stroke", "Intracerebral hemorrhage", "Subarachnoid hemorrhage")
subtype_short <- c("IS", "ICH", "SAH")

res_list <- list()
for (si in seq_along(subtypes)) {
  stype <- subtypes[si]
  sshort <- subtype_short[si]
  
  cat("\n\n==========", stype, "==========\n")
  res <- feature_based_cluster(general_year, stype)
  res_list[[sshort]] <- res
  
  plot_dend(res$hc, res$k, paste0("Supp_", sshort, "_dendrogram.png"))
  
  fig <- plot_faceted_trends(general_year, stype, res$clusters)
  save_panel(fig, paste0("SuppFig_", sshort, "_clustered"),
             w_mm = 180, h_mm = 130)
}

cat("\n========== 全部完成 ==========\n")