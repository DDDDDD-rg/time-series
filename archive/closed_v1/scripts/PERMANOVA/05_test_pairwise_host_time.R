# 当前MC19 coral方案停用本历史入口；原代码留作查阅，原文件另有备份。
stop("本轮未安排宿主/TS两两检验，不运行旧MC8两两模型。请使用run_permanova.R或04_plot_permanova.R。",call.=FALSE)

# 宿主和采样批次的两两 PERMANOVA：独立运行，不重跑02的整体检验。
# 自动读取01的范围、Repeat_group均值、平方根Bray距离和置换次数。
project_dir <- "C:/Users/xzy/Desktop/CUHK_Study/GNBF6010_Research_Projects/Research/Msc main project-time series"
source(file.path(project_dir, "scripts", "PERMANOVA", "01_data_settings.R"), encoding = "UTF-8", local = TRUE)
source(file.path(project_dir, "scripts", "PERMANOVA", "pairwise_letters.R"), encoding = "UTF-8", local = TRUE)

# ===== 1. 比较设置 =====
# 两类分别做BH：Host_type比较宿主；TS比较采样批次。
# 如只需要其中一种，可写成 "Host_type" 或 "TS"。
pairwise_targets <- c("Host_type", "TS")
pairwise_dir <- file.path(output_dir, "pairwise")
stopifnot(length(pairwise_targets) > 0, !anyDuplicated(pairwise_targets),
          all(pairwise_targets %in% c("Host_type", "TS")))

# ===== 2. 准备一对：只保留两个组共同覆盖的另一因素 =====
# 宿主对 -> 共同TS；批次对 -> 共同宿主。
# 缺采组合不补0；每个Repeat_group保留一行。
source(file.path(project_dir,"scripts", "PERMANOVA","permanova_pairwise_helpers.R"), encoding="UTF-8", local=TRUE)

# ===== 7. 运行所选两类并保存；输出与02/03分别存放 =====
dir.create(pairwise_dir, recursive = TRUE, showWarnings = FALSE)
saved_families <- list()
for (effect in pairwise_targets) {
  family <- run_pairwise_family(
    analysis_data, effect,
    seed_offset = if (effect == "TS") 100000L else 0L
  )
  results <- family$results
  family$letters <- derive_pairwise_letters(results, family$levels)
  prefix <- if (effect == "Host_type") "pairwise_host" else "pairwise_TS"
  write_csv_excel(family$letters$table,
                  file.path(pairwise_dir, paste0(prefix, "_letters.csv")),
                  row.names = FALSE)
  write_csv_excel(results, file.path(pairwise_dir, paste0(prefix, "_results.csv")),
                  row.names = FALSE)
  write_csv_excel(family$coverage,
                  file.path(pairwise_dir, paste0(prefix, "_coverage.csv")),
                  row.names = FALSE)
  for (column in c("P_BH_raw", "P_BH_lingoes", "R2_raw")) {
    write_csv_excel(make_pair_matrix(results, family$levels, column),
                    file.path(pairwise_dir, paste0(prefix, "_", column, "_matrix.csv")))
  }
  display_order <- order(results$P_BH_raw, na.last = TRUE)
  print(results[display_order, c("Group1", "Group2", "N1", "N2",
                                "R2_raw", "P_raw", "P_BH_raw",
                                "Significance", "Status")], row.names = FALSE)
  saved_families[[effect]] <- family
}
saveRDS(list(data = analysis_data, families = saved_families,
             scope = analysis_scope, seed = random_seed,
             n_permutations_requested = n_permutations),
         file.path(pairwise_dir, "pairwise_host_time_analysis.rds"))
writeLines(capture.output(sessionInfo()),
            file.path(pairwise_dir, "pairwise_host_time_session.txt"))
message("两两比较完成。输出目录：", pairwise_dir)
