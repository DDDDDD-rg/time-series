# 直观展示：按DO由低到高排列各批次，分别显示每个宿主的所选MC相对丰度。
# 图示仅为描述，不运行新检验，不把整体PERMANOVA的p值赋给单个MC。
project_dir <- "C:/Users/xzy/Desktop/CUHK_Study/GNBF6010_Research_Projects/Research/Msc main project-time series"
source(file.path(project_dir, "scripts", "permanova_plot_style.R"), local=TRUE)
analysis_scope <- getOption("mc.analysis_scope", "all_hosts")
latest_result_pointer <- file.path(project_dir, "results", "00_latest_PERMANOVA_MC8.txt")
latest_figure_pointer <- file.path(project_dir, "figures", "00_latest_PERMANOVA_MC8.txt")
default_results_root <- if (file.exists(latest_result_pointer))
  readLines(latest_result_pointer, n = 1L, warn = FALSE, encoding = "UTF-8") else
  file.path(project_dir, "results", "PERMANOVA_MC8")
default_figures_root <- if (file.exists(latest_figure_pointer))
  readLines(latest_figure_pointer, n = 1L, warn = FALSE, encoding = "UTF-8") else
  file.path(project_dir, "figures", "PERMANOVA_MC8")
preview_dir <- getOption("mc.figure_dir", file.path(default_figures_root, analysis_scope))
suppressPackageStartupMessages(library(ggplot2))
saved <- readRDS(file.path(getOption("mc.results_root", default_results_root), analysis_scope,
                           "repeat_groups", "season_environment_analysis.rds"))
data <- as.data.frame(if(!is.null(saved$environment_inputs)) saved$environment_inputs[["DO_mg_L"]] else saved$environment_input)
if(analysis_scope=="all_hosts" && !all(c("Seawater","Sediment") %in% data$Host_type))
  stop("Old results exclude environmental samples; run the updated analysis first.")
result <- saved$environment[saved$environment$Factor == "DO_mg_L", , drop = FALSE]
mc_cols <- grep("^MC", names(data), value = TRUE)
stopifnot(nrow(result) == 1, result$Status == "ok", result$Scope == analysis_scope,
          !anyNA(data[, c("TS", "Host_type", "DO_mg_L", mc_cols)]),
          all(as.matrix(data[, mc_cols]) >= 0), all(rowSums(data[, mc_cols]) <= 1 + 1e-8))

batch <- unique(data[, c("TS", "DO_mg_L")])
stopifnot(anyDuplicated(batch$TS) == 0)
batch <- batch[order(batch$DO_mg_L), ]
# 每根柱子对应一个Repeat_group，避免同宿主同TS的多组丰度叠加在同一根柱上。
stopifnot("Repeat_group_ID" %in% names(data), anyDuplicated(data$Repeat_group_ID) == 0,
          result$Analysis_unit == "Repeat_group均值")
data <- data[order(data$DO_mg_L, data$TS, data$Host_type, data$Repeat_group_ID), ]
data$Group_number <- ave(seq_len(nrow(data)), data$Host_type, data$TS, FUN = seq_along)
labels <- setNames(sprintf("%.1f\n%s\nR%s", data$DO_mg_L, data$TS, data$Group_number),
                   data$Repeat_group_ID)
long <- do.call(rbind, lapply(mc_cols, function(mc) data.frame(
  Host_type = data$Host_type, TS = data$TS, Repeat_group_ID = data$Repeat_group_ID,
  MC = mc, Abundance = data[[mc]])))
long$Repeat_group_ID <- factor(long$Repeat_group_ID, levels = data$Repeat_group_ID)
long$MC <- factor(long$MC, levels = mc_cols)
colours <- pv_subset(pv_mc, mc_cols)
caption <- paste0("Overall host-adjusted PERMANOVA: R2 = ", sprintf("%.3f", result$R2_raw),
  "; BH p = ", formatC(result$P_BH_raw, format = "g", digits = 3),
  ". This p-value is for the ", length(mc_cols), "-MC vector, not individual MCs.\n",
  "One bar = one Repeat_group mean. R labels identify different groups within a host and batch.\n",
  paste0("Samples share ",length(unique(data$TS))," batch-level DO records. Columns are ordered by DO, not time.\n"),
  "Only the ", length(mc_cols), " selected MCs are shown; relative abundances retain their original denominator and are not rescaled to 100%.")
plot <- ggplot(long, aes(Repeat_group_ID, Abundance, fill = MC)) +
  geom_col(width = 0.75, position = position_stack(reverse = TRUE)) +
  facet_wrap(~ Host_type, ncol = 1, scales = "free_x") +
  scale_fill_manual(values = colours, drop = FALSE) +
  scale_x_discrete(labels = labels) +
  scale_y_continuous(labels = function(x) paste0(round(100 * x), "%"),
                     breaks = seq(0, 1, by = 0.1), expand = expansion(mult = c(0, 0.05))) +
  labs(title = NULL,
       subtitle = "Within each host: Repeat_groups ordered from lower to higher batch DO",
       x = "Dissolved oxygen (mg/L), sampling batch and Repeat_group", y = "Relative abundance in the original MC table",
       fill = "Selected MC", caption = caption) +
  pv_theme(base_size = 14) +
  theme(strip.background = element_rect(fill = "grey95", colour = "grey70"),
        strip.text = element_text(size = 12), panel.spacing = grid::unit(0.6, "lines"),
        legend.key.height = grid::unit(0.65, "cm"), plot.title = element_blank(),
        plot.subtitle = element_blank(), plot.caption = element_text(hjust = 0, size = 9),
        plot.margin = margin(12, 14, 12, 12))
dir.create(preview_dir, recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(preview_dir, "DO_MC_composition_preview.png"), plot, width = 11.5, height = max(6, 1.6 * length(unique(data$Host_type)) + 2), dpi = 250, bg = "white")
saveRDS(list(plot = plot, data = data, batch = batch, result = result), file.path(preview_dir, "preview_data.rds"))
# 核对柱子的总高度与所选MC原始丰度之和一致，且批次排序采用真实DO。
built <- ggplot_build(plot)
totals <- rowSums(data[, mc_cols])
drawn <- unlist(lapply(split(built$data[[1]], built$data[[1]]$PANEL), function(panel) {
  tapply(panel$ymax, panel$x, max)
}), use.names = FALSE)
stopifnot(length(drawn) == nrow(data),
          isTRUE(all.equal(sort(as.numeric(drawn)), sort(totals), tolerance = 1e-8, check.attributes = FALSE)))
cat("PASS:", nrow(data), "Repeat_group bars; original abundances and DO sorting retained.\n")
print(batch)
