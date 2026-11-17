# 当前MC19 coral方案停用本历史入口；原代码留作查阅，原文件另有备份。
stop("本轮只做coral环境关联，不运行旧Host_type/TS调整模型。请使用run_permanova.R或04_plot_permanova.R。",call.=FALSE)

# 宿主和时间：可以单独打开本文件，从头运行。
# 设置在01_data_settings.R中；本脚本不依赖旧脚本留下的工作区对象。
project_dir <- "C:/Users/xzy/Desktop/CUHK_Study/GNBF6010_Research_Projects/Research/Msc main project-time series"
source(file.path(project_dir, "scripts", "PERMANOVA", "01_data_settings.R"), encoding = "UTF-8", local = TRUE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ===== 1. 宿主：比较同批次里的不同宿主 =====
# 模型：MC组成 ~ Host_type + TS；TS内置换；只解释Host_type行。
if (nlevels(analysis_data$Host_type) < 2 || nlevels(analysis_data$TS) < 2) {
  stop("宿主/时间检验至少需要两个宿主和两个采样批次。")
}
set.seed(random_seed)
perm_host <- permute::shuffleSet(nrow(analysis_data),
  control = permute::how(nperm = n_permutations, blocks = analysis_data$TS))
host_result <- test_factor(
  data = analysis_data, rhs = c("Host_type", "TS"), effect = "Host_type",
  permutations = perm_host, dispersion_group = "Host_type"
)
host_result$Permutation_design <- "在TS内置换重复组分析单位"

# ===== 2. 时间：比较同类宿主在不同批次的组成 =====
# 模型相同，但在Host内置换；只解释TS行。
# 沿用用户的独立重复组假设；不等于已确认的同一个体纵向分析。
set.seed(random_seed)
perm_time <- permute::shuffleSet(nrow(analysis_data),
  control = permute::how(nperm = n_permutations, blocks = analysis_data$Host_type))
time_result <- test_factor(
  data = analysis_data, rhs = c("Host_type", "TS"), effect = "TS",
  permutations = perm_time, dispersion_group = "TS"
)
time_result$Permutation_design <- "在Host_type内置换重复组分析单位"

# ===== 3. 保存两项总体检验 =====
factor_results <- bind_rows(host_result, time_result)
factor_results$Analysis_unit <- "Repeat_group均值"
factor_results$Scope <- analysis_scope
factor_results$P_BH_raw <- add_bh(factor_results$P_raw)
factor_results$P_BH_lingoes <- add_bh(factor_results$P_lingoes)
factor_results$Significance <- add_stars(factor_results$P_BH_raw)
# 为这两个总体问题提供BH值；未校正p也保留，不与旧两两比较的BH混用。
print(factor_results[, c("Factor", "N_units", "R2_raw", "P_raw", "P_BH_raw",
                        "P_BH_lingoes", "Significance", "Dispersion_note", "Status")])
write_csv_excel(factor_results, file.path(output_dir, "host_time_results.csv"), row.names = FALSE)
write_csv_excel(as.data.frame.matrix(table(analysis_data$Host_type, analysis_data$TS)),
          file.path(output_dir, "host_by_TS_coverage.csv"))
write_csv_excel(excluded_hosts, file.path(output_dir, "excluded_small_hosts.csv"), row.names = FALSE)
saveRDS(list(data = analysis_data, results = factor_results), file.path(output_dir, "host_time_analysis.rds"))
writeLines(capture.output(sessionInfo()), file.path(output_dir, "host_time_session.txt"))

# ===== 4. 可选展示：不是检验必需步骤 =====
# 想画图时把FALSE改成TRUE。检验始终使用Bray距离，不使用NMDS坐标。
make_nmds_plot <- getOption("mc.nmds_diagnostics", TRUE) # Full runner fits 2D once in 04.
if (make_nmds_plot) {
  set.seed(random_seed)
  nmds <- vegan::metaMDS(make_distance(analysis_data), k = 2, trymax = 100, trace = FALSE)
  grDevices::pdf(file.path(output_dir, "nmds_diagnostics.pdf"), width = 7, height = 7)
  tryCatch({
    plot(nmds, display = "sites", type = "p", main = paste("MC8 NMDS; stress =", round(nmds$stress, 3)))
    vegan::stressplot(nmds)
  }, finally = grDevices::dev.off())
  saveRDS(nmds, file.path(output_dir, "nmds.rds"))
}
message("宿主和时间分析完成。输出目录：", output_dir)
