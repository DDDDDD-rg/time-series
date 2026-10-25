# PM11 水质时间序列：五个指标、三个水层、干湿季和 TS 标记
# 用法：在 RStudio 打开本文件，点击 Source，一次生成五张图。
# 唯一外部依赖：ggplot2；若未安装，先运行 install.packages("ggplot2")。
# 原始 CSV 和已有 Excel 均只读，不会修改。

# 1. 可调整设置 ---------------------------------------------------------------
project_dir <- "C:/Users/xzy/Desktop/CUHK_Study/GNBF6010_Research_Projects/Research/Msc main project-time series"
input_file <- file.path(project_dir, "data", "water_quality", "Bluff Island_PM11_marine_water_quality.csv")
batch_file <- file.path(project_dir, "results", "all_hosts", "repeat_groups",
                        "season_analysis_input.csv")

number_of_years <- 5L
end_year <- 2025L              # 当前原始数据截至 2025 年；起始年自动向前推四年。
wet_months <- 4:9              # 湿季：4–9 月；干季：10–3 月。
save_pdf <- TRUE               # 同时导出可编辑的矢量 PDF；PNG 始终导出。
figure_width <- 12.5           # 单位：英寸。
figure_height <- 4.8
png_dpi <- 320
base_font <- "sans"

# 默认每次新建一个带时间戳的 figures 子目录，保留以往出图。
# 此环境变量仅供指定另一个输出位置；普通使用无需设置。
output_dir <- Sys.getenv("PM11_PLOT_OUTPUT_DIR", unset = "")
if (!nzchar(output_dir)) {
  output_dir <- file.path(project_dir, "figures", paste0(
    "PM11_timeseries_", end_year - number_of_years + 1L, "_", end_year,
    "_", format(Sys.time(), "%Y%m%d_%H%M%S")
  ))
}

depth_levels <- c("Surface Water", "Middle Water", "Bottom Water")
depth_colours <- c("Surface Water" = "#D64B42", "Middle Water" = "#E69F00",
                   "Bottom Water" = "#2866A0")
depth_shapes <- c("Surface Water" = 16, "Middle Water" = 15, "Bottom Water" = 17)

# DIN 沿用此前作图表的简称；数据取官方 Total Inorganic Nitrogen 列。
# 不把处理后的硝酸盐、亚硝酸盐和氨氮相加来重新计算总无机氮。
metric_columns <- c(
  DIN = "Total Inorganic Nitrogen (mg/L)",
  Nitrate = "Nitrate Nitrogen (mg/L)",
  Nitrite = "Nitrite Nitrogen (mg/L)",
  Ammonia = "Ammonia Nitrogen (mg/L)",
  DO = "Dissolved Oxygen (mg/L)"
)
metric_titles <- c(DIN = "Total inorganic nitrogen", Nitrate = "Nitrate nitrogen",
                   Nitrite = "Nitrite nitrogen", Ammonia = "Ammonia nitrogen",
                   DO = "Dissolved oxygen")

# 2. 检查依赖和输入 -----------------------------------------------------------
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop('请先安装 ggplot2：install.packages("ggplot2")', call. = FALSE)
}
if (utils::packageVersion("ggplot2") < "3.4.0") {
  stop('请更新 ggplot2 至 3.4.0 或以上：install.packages("ggplot2")', call. = FALSE)
}
suppressPackageStartupMessages(library(ggplot2))
stopifnot(number_of_years >= 1L, end_year >= 1900L,
          all(wet_months %in% 1:12))
if (!file.exists(input_file)) stop("找不到原始数据：", input_file, call. = FALSE)

raw_data <- utils::read.csv(input_file, check.names = FALSE, colClasses = "character",
                            fileEncoding = "UTF-8-BOM", stringsAsFactors = FALSE)
needed_columns <- c("Station", "Dates", "Sample No", "Depth", unname(metric_columns))
if (!all(needed_columns %in% names(raw_data))) {
  stop("原始 CSV 缺少列：", paste(setdiff(needed_columns, names(raw_data)), collapse = ", "))
}
raw_data <- raw_data[raw_data$Station == "PM11", , drop = FALSE]
raw_data$Date <- as.Date(raw_data$Dates, format = "%Y/%m/%d")
if (anyNA(raw_data$Date)) stop("原表存在不能解析的日期，请检查 Dates 列。")
if (!all(raw_data$Depth %in% depth_levels)) stop("原表存在未识别的水层名称。")
record_key <- paste(raw_data$Date, raw_data$Depth, raw_data[["Sample No"]], sep = "_")
if (anyDuplicated(record_key)) stop("原表存在重复采样记录，请先核对，避免重复计数。")

start_year <- end_year - number_of_years + 1L
start_date <- as.Date(sprintf("%04d-01-01", start_year))
end_date <- as.Date(sprintf("%04d-12-31", end_year))
if (max(raw_data$Date) < as.Date(sprintf("%04d-12-01", end_year))) {
  stop("end_year 超出了原始数据覆盖的完整日历年，请修改设置。")
}
selected_data <- raw_data[raw_data$Date >= start_date & raw_data$Date <= end_date,
                          , drop = FALSE]
if (!nrow(selected_data)) stop("所选时间范围内没有 PM11 数据。")
selected_data$Month <- as.Date(format(selected_data$Date, "%Y-%m-01"))

# 3. 检出限处理和逐月数据 -----------------------------------------------------
# 如 <0.002，作图数值用 0.001；同时保留原值和检出限标记供核对。
# 空白/NA 保留缺失；遇到不能识别的文字时停止，不静默改成零。
parse_measurement <- function(x, column_name) {
  text <- trimws(x)
  missing <- is.na(text) | toupper(text) %in% c("", "NA", "N/A", "NAN", "-", "--")
  censored <- !missing & grepl("^(<=|<|\u2264)", text)
  number_text <- sub("^(<=|<|\u2264)\\s*", "", text)
  values <- suppressWarnings(as.numeric(number_text))
  bad <- !missing & !is.finite(values)
  if (any(bad)) stop(column_name, " 有无法解析的值：", paste(unique(text[bad]), collapse = ", "))
  values[missing] <- NA_real_
  if (any(values < 0, na.rm = TRUE)) stop(column_name, " 存在负浓度。")
  reporting_limit <- ifelse(censored, values, NA_real_)
  values[censored] <- values[censored] / 2
  data.frame(Value = values, Below_reporting_limit = censored,
             Reporting_limit = reporting_limit)
}

long_data <- do.call(rbind, lapply(names(metric_columns), function(metric) {
  original_values <- selected_data[[metric_columns[[metric]]]]
  parsed <- parse_measurement(original_values, metric_columns[[metric]])
  data.frame(Date = selected_data$Date, Month = selected_data$Month,
             Depth = selected_data$Depth, Metric = metric, Raw_value = original_values,
             parsed, stringsAsFactors = FALSE)
}))
rownames(long_data) <- NULL

# 同一年同一月同一水层若有多次采样，先做月内算术平均。
# 2021–2025 窗口内实际上每个有记录的年月水层都仅有一次采样。
# 不跨年份平均，也不把三个水层合并为重复来计算 SEM。
group_key <- interaction(long_data$Metric, long_data$Depth, long_data$Month, drop = TRUE)
monthly_data <- do.call(rbind, lapply(split(long_data, group_key), function(z) {
  data.frame(Month = z$Month[1], Depth = z$Depth[1], Metric = z$Metric[1],
             Value = if (all(is.na(z$Value))) NA_real_ else mean(z$Value, na.rm = TRUE),
             N_sampling = nrow(z), N_valid = sum(!is.na(z$Value)),
             N_below_reporting_limit = sum(z$Below_reporting_limit),
             Source_dates = paste(sort(unique(as.character(z$Date))), collapse = ";"),
             stringsAsFactors = FALSE)
}))

# 补齐每一个日历月，缺测月保留 NA；绘图时连接前后已有观测，不填补缺测值。
calendar <- seq(start_date, as.Date(sprintf("%04d-12-01", end_year)), by = "month")
calendar_grid <- expand.grid(Month = calendar, Depth = depth_levels,
                             Metric = names(metric_columns), stringsAsFactors = FALSE)
plot_data <- merge(calendar_grid, monthly_data, by = c("Month", "Depth", "Metric"),
                   all.x = TRUE, sort = FALSE)
for (column in c("N_sampling", "N_valid", "N_below_reporting_limit")) {
  plot_data[[column]][is.na(plot_data[[column]])] <- 0L
}
plot_data$Source_dates[is.na(plot_data$Source_dates)] <- ""
plot_data$Season <- ifelse(as.integer(format(plot_data$Month, "%m")) %in% wet_months,
                           "Wet", "Dry")
plot_data$Depth <- factor(plot_data$Depth, levels = depth_levels)
plot_data$Plot_date <- plot_data$Month + 14L  # 月中定位，不代表实际采样日。
plot_data <- plot_data[order(match(plot_data$Metric, names(metric_columns)),
                            plot_data$Depth, plot_data$Month), ]
rownames(plot_data) <- NULL

# 4. 读取 TS 年月，并生成干湿季区间 ------------------------------------------
# 优先从已有 PERMANOVA 输入结果读取；无此文件时用已核对的年月表。
if (file.exists(batch_file)) {
  batch_input <- utils::read.csv(batch_file, check.names = FALSE,
                                 fileEncoding = "UTF-8-BOM", stringsAsFactors = FALSE)
  if (!all(c("TS", "Sampling_month") %in% names(batch_input))) {
    stop("TS 来源文件缺少 TS / Sampling_month 列。")
  }
  ts_dates <- unique(batch_input[, c("TS", "Sampling_month")])
} else {
  ts_dates <- data.frame(TS = paste0("TS", 1:8),
                         Sampling_month = c("2024-04", "2024-07", "2024-10", "2025-02",
                                            "2025-04", "2025-08", "2025-11", "2026-02"))
  message("未找到 TS 来源表，使用脚本中已核对的 TS 年月。")
}
if (anyNA(ts_dates) || anyDuplicated(ts_dates$TS)) stop("同一个 TS 对应多个或缺失年月。")
ts_dates$Plot_date <- as.Date(paste0(ts_dates$Sampling_month, "-15"))
if (anyNA(ts_dates$Plot_date)) stop("TS 年月应为 YYYY-MM 格式。")
ts_dates$Label <- paste(ts_dates$TS,
                        month.abb[as.integer(format(ts_dates$Plot_date, "%m"))], sep = "\n")
ts_dates <- ts_dates[order(ts_dates$Plot_date), ]
ts_in_range <- ts_dates[ts_dates$Plot_date >= start_date & ts_dates$Plot_date <= end_date, ]
ts_outside <- ts_dates[ts_dates$Plot_date < start_date | ts_dates$Plot_date > end_date, ]

# 按每月干湿季分类生成连续区间，干季自然跨越相邻两年。
season_grid <- data.frame(Start = calendar,
                          Season = ifelse(as.integer(format(calendar, "%m")) %in% wet_months,
                                          "Wet", "Dry"))
season_grid$End <- c(calendar[-1], as.Date(sprintf("%04d-01-01", end_year + 1L)))
season_grid$Run <- cumsum(c(TRUE, tail(season_grid$Season, -1) != head(season_grid$Season, -1)))
season_intervals <- do.call(rbind, lapply(split(season_grid, season_grid$Run), function(z) {
  data.frame(Start = min(z$Start), End = max(z$End), Season = z$Season[1])
}))
season_intervals$End[season_intervals$End > end_date] <- end_date
season_intervals$Midpoint <- season_intervals$Start +
  floor(as.numeric(season_intervals$End - season_intervals$Start) / 2)
wet_intervals <- season_intervals[season_intervals$Season == "Wet", ]

# 5. 单指标绘图函数 -----------------------------------------------------------
make_metric_plot <- function(metric) {
  z <- plot_data[plot_data$Metric == metric, ]
  if (all(is.na(z$Value))) stop(metric, " 在所选范围内全部缺失，无法出图。")
  # 仅将已有观测用于画线和画点。缺测月仍保留在 plot_data 中，不新增观测点。
  observed_data <- z[!is.na(z$Value), , drop = FALSE]
  # 中层和底层先画，表层最后画，保持原有颜色、线宽和点大小。
  other_depth_data <- observed_data[observed_data$Depth != "Surface Water", , drop = FALSE]
  surface_data <- observed_data[observed_data$Depth == "Surface Water", , drop = FALSE]
  max_value <- max(z$Value, na.rm = TRUE)
  y_top <- if (max_value > 0) max_value * 1.18 else 1
  y_breaks <- pretty(c(0, max_value * 1.05), n = if (metric == "DO") 7 else 5)
  y_breaks <- y_breaks[y_breaks >= 0 & y_breaks <= y_top * 0.93]
  y_label <- if (metric == "DO") expression(Concentration~(mg~L^{-1})) else
    expression(Concentration~(mg~N~L^{-1}))

  caption <- paste0(
    "Shading: wet season (April-September). White background: dry season (October-March).\n",
    "Lines connect available monthly observations across missing months. Values below reporting limits are shown at half the limit."
  )
  if (!identical(sort(wet_months), 4:9)) {
    caption <- paste0("Shading: wet season. White background: dry season.\n",
                       "Lines connect available monthly observations across missing months. Values below reporting limits are shown at half the limit.")
  }

  p <- ggplot(z, aes(x = Plot_date, y = Value, colour = Depth, shape = Depth, group = Depth)) +
    geom_rect(data = wet_intervals, aes(xmin = Start, xmax = End, ymin = 0, ymax = Inf),
               inherit.aes = FALSE, fill = "#D9DEE5", alpha = 0.75, colour = NA) +
    geom_line(data = other_depth_data, linewidth = 0.55, na.rm = TRUE, show.legend = TRUE) +
    geom_point(data = other_depth_data, size = 1.85, stroke = 0.25, na.rm = TRUE, show.legend = TRUE) +
    geom_line(data = surface_data, linewidth = 0.55, na.rm = TRUE, show.legend = FALSE) +
    geom_point(data = surface_data, size = 1.85, stroke = 0.25, na.rm = TRUE, show.legend = FALSE) +
    geom_text(data = season_intervals, aes(x = Midpoint, y = y_top * 0.98, label = Season),
               inherit.aes = FALSE, family = base_font, size = 3.25, colour = "#515863", vjust = 1) +
    scale_colour_manual(values = depth_colours, breaks = depth_levels, drop = FALSE) +
    scale_shape_manual(values = depth_shapes, breaks = depth_levels, drop = FALSE) +
    scale_x_date(breaks = as.Date(sprintf("%04d-01-01", start_year:end_year)),
                  date_labels = "%Y", minor_breaks = calendar, expand = expansion(mult = 0)) +
    scale_y_continuous(breaks = y_breaks, expand = expansion(mult = 0)) +
    coord_cartesian(xlim = c(start_date, end_date), ylim = c(0, y_top), clip = "off") +
    labs(title = metric_titles[[metric]],
          subtitle = sprintf("PM11, Bluff Island, %d-%d", start_year, end_year),
          x = "Year", y = y_label, colour = NULL, shape = NULL, caption = caption) +
    theme_classic(base_size = 12, base_family = base_font) +
    theme(
      plot.title = element_text(size = 17, face = "bold", margin = margin(b = 5)),
      plot.subtitle = element_text(size = 11, colour = "#515863", margin = margin(b = 13)),
      axis.title.y = element_text(size = 12, margin = margin(r = 9)),
      axis.title.x = element_text(size = 12, margin = margin(t = if (nrow(ts_in_range)) 44 else 10)),
      axis.text = element_text(colour = "#20252B", size = 11),
      axis.line = element_line(colour = "#20252B", linewidth = 0.45),
      axis.ticks = element_line(colour = "#20252B", linewidth = 0.4),
      axis.ticks.length = grid::unit(2, "mm"),
      legend.position = "top", legend.justification = "left",
      legend.text = element_text(size = 10.5),
      legend.key.width = grid::unit(1.2, "cm"),
      legend.margin = margin(0, 0, 7, 0),
      plot.caption = element_text(size = 8.5, colour = "#515863", hjust = 0,
                                   lineheight = 1.15, margin = margin(t = 12)),
      plot.margin = margin(14, 18, 12, 14)
    )

  # 月中箭头表示研究采样批次所在月份，不声称具体采样日。
  if (nrow(ts_in_range)) {
    p <- p +
      geom_segment(data = ts_in_range,
                     aes(x = Plot_date, xend = Plot_date, y = -0.12 * y_top, yend = -0.022 * y_top),
                     inherit.aes = FALSE, colour = "#20252B", linewidth = 0.35,
                     arrow = grid::arrow(length = grid::unit(1.8, "mm"), type = "closed")) +
      geom_text(data = ts_in_range, aes(x = Plot_date, y = -0.16 * y_top, label = Label),
                  inherit.aes = FALSE, family = base_font, colour = "#20252B", size = 3,
                  lineheight = 0.95, vjust = 1)
  }
  p
}

# 6. 导出五张图及对应作图数据 -------------------------------------------------
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
plots <- setNames(lapply(names(metric_columns), make_metric_plot), names(metric_columns))
for (i in seq_along(plots)) {
  stem <- sprintf("%02d_%s_PM11_%d_%d", i, names(plots)[i], start_year, end_year)
  png_args <- list(filename = file.path(output_dir, paste0(stem, ".png")),
                    plot = plots[[i]], width = figure_width, height = figure_height,
                    units = "in", dpi = png_dpi, bg = "white")
  do.call(ggsave, png_args)
  if (save_pdf) {
    pdf_device <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
    ggsave(file.path(output_dir, paste0(stem, ".pdf")), plot = plots[[i]],
             width = figure_width, height = figure_height, units = "in", device = pdf_device,
             bg = "white")
  }
}
utils::write.csv(plot_data, file.path(output_dir, "PM11_monthly_plot_data.csv"),
                  row.names = FALSE, na = "", fileEncoding = "UTF-8")
utils::write.csv(ts_dates, file.path(output_dir, "TS_month_mapping.csv"),
                  row.names = FALSE, fileEncoding = "UTF-8")

run_notes <- c(
  sprintf("PM11 time series: %d-%d", start_year, end_year),
  paste("Source:", input_file),
  paste("TS source:", if (file.exists(batch_file)) batch_file else "verified mapping in script"),
  paste("Wet months:", paste(sort(wet_months), collapse = ", "),
        "; dry months:", paste(setdiff(1:12, wet_months), collapse = ", ")),
  "Points represent year-specific monthly values, not a five-year seasonal average.",
  "Multiple samples within one year/month/depth are averaged when present. No SEM is computed.",
  "Missing months remain NA in the exported data. Lines connect existing observations across missing months; no values are imputed and no points are added.",
  "Censored values are shown at half their reporting limits. DIN uses the original official TIN column.",
  "Nitrogen concentrations: mg N/L. Dissolved oxygen: mg/L, not percent saturation.",
  "TS arrows indicate calendar months, positioned at the 15th for display.",
  paste("TS shown:", paste(ts_in_range$TS, collapse = ", ")),
  paste("TS outside the plotting range:", paste(ts_outside$TS, collapse = ", ")),
  paste("R:", R.version.string),
  paste("ggplot2:", utils::packageVersion("ggplot2"))
)
writeLines(run_notes, file.path(output_dir, "plot_notes.txt"), useBytes = TRUE)
capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo.txt"))
message("完成：五张图已导出至 ", normalizePath(output_dir, winslash = "/"))
if (nrow(ts_outside)) message("时间范围之外的 TS：", paste(ts_outside$TS, collapse = ", "))
if (interactive()) print(plots[[1]])
