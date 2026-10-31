# 当前主分析：coral/ecosystem原始独立样本，19个MC，7项当期环境；不做季节或宿主检验。
# 手动source本文件也会建立独立结果目录；不会读取旧01或按Repeat_group平均。
project_dir <- "C:/Users/xzy/Desktop/CUHK_Study/GNBF6010_Research_Projects/Research/Msc main project-time series"
source(file.path(project_dir, "scripts", "PERMANOVA", "permanova_MC19_data_settings.R"), encoding="UTF-8", local=TRUE)
input <- read_mc19_samples()
output_dir <- mc19_new_run(project_dir, "PERMANOVA", analysis_scope, analysis_version, input_file,
  requested=getOption("mc19.run_id", NULL))
run_name <- basename(output_dir)
support_dir <- mc19_support_dir(output_dir)
dir.create(support_dir)
run_success <- FALSE
status_file <- file.path(support_dir, "00_run_status.txt")
writeLines(c("running", analysis_version), status_file)

tryCatch({
  scripts_used <- c("permanova_MC19_data_settings.R", "03_test_season_environment_batch_labels.R",
                    "03_test_season_environment.R", "run_permanova.R", "04_plot_permanova.R",
                    "permanova_MC19_plot_helpers.R", "07_plot_host_highlight_panels.R",
                    "plot_MC_environment_composition.R", "plot_DO_composition_preview.R", "../mc19_run_io.R")
  snapshot_dir <- file.path(support_dir, "scripts_used")
  dir.create(snapshot_dir)
  if (!all(file.copy(file.path(project_dir, "scripts", "PERMANOVA", scripts_used), snapshot_dir))) stop("脚本快照失败。")
  if (!file.copy(input_file, file.path(support_dir, "00_input_used_MC19_env7.xlsx"))) stop("输入快照失败。")
  writeLines(analysis_version, file.path(support_dir, "00_analysis_version.txt"))
  message(analysis_scope, "原始样本：", nrow(input$community), "；MC：", ncol(input$community),
          "；有效TS：", nrow(input$batches), "。无样本平均，无Host_type调整。")
  distance <- make_coral_distance(input$community)
  orders <- enumerate_ts_orders(nrow(input$batches))
  results <- vector("list", length(env_cols))
  null_f <- matrix(NA_real_, nrow=nrow(orders), ncol=length(env_cols), dimnames=list(NULL, env_cols))
  for (i in seq_along(env_cols)) {
    message("环境 ", i, "/", length(env_cols), "：", env_cols[i])
    fit <- fit_one_environment(env_cols[i], input, distance, orders)
    results[[i]] <- fit$result
    null_f[, i] <- fit$F_per_order
  }
  environment_results <- do.call(rbind, results)
  environment_results$P_BH_raw <- p.adjust(environment_results$P_raw, method="BH", n=length(env_cols))
  environment_results$Significance <- ifelse(environment_results$P_BH_raw <= 0.001, "***",
    ifelse(environment_results$P_BH_raw <= 0.01, "**",
      ifelse(environment_results$P_BH_raw <= 0.05, "*", "ns")))
  environment_results$Adjustment <- paste("BH across seven", analysis_scope, "environment tests")
  stopifnot(nrow(environment_results)==7L, all(environment_results$Status=="ok"))
  main_columns <- c("Scope", "Indicator", "Factor", "N_samples", "N_MC", "N_TS",
                    "R2_raw", "F_raw", "P_raw", "P_BH_raw", "Significance")
  write_csv_excel(environment_results[, main_columns, drop=FALSE],
                  file.path(output_dir, "environment_results.csv"))
  write_csv_excel(environment_results, file.path(support_dir, "environment_results_full.csv"))
  write_csv_excel(cbind(input$metadata, as.data.frame(input$community)),
                  file.path(support_dir, "environment_analysis_input.csv"))
  write_csv_excel(input$batches, file.path(support_dir, "environment_TS_values.csv"))
  write_csv_excel(input$coverage, file.path(support_dir, paste0(analysis_scope, "_sample_coverage.csv")))
  labels <- matrix(input$batches$TS[as.vector(orders)], nrow=nrow(orders), ncol=ncol(orders))
  colnames(labels) <- paste0("Target_", input$batches$TS)
  identity <- seq_len(nrow(orders)) == 1L
  write_csv_excel(data.frame(Order_ID=seq_len(nrow(orders)), Is_observed=identity, labels, check.names=FALSE),
                  file.path(support_dir, "TS_label_orders.csv"))
  write_csv_excel(data.frame(Order_ID=seq_len(nrow(orders)), Is_observed=identity, null_f, check.names=FALSE),
                  file.path(support_dir, "permutation_F_values.csv"))
  # 相关性按每TS一行计算；不重复共享值来增加相关信息。
  correlations <- cor(as.matrix(input$batches[env_cols]))
  write_csv_excel(data.frame(Factor=rownames(correlations), correlations, check.names=FALSE),
                  file.path(support_dir, "environment_correlations_TS.csv"))
  saveRDS(list(results=environment_results, input=input, distance=distance,
    TS_orders=orders, F_per_order=null_f, version=analysis_version, scope=analysis_scope),
    file.path(support_dir, paste0(analysis_scope, "_environment_analysis.rds")))
  sources <- c(input_file, file.path(snapshot_dir, basename(scripts_used)))
  write_csv_excel(data.frame(File=sources, MD5=unname(tools::md5sum(sources))),
                  file.path(support_dir, "00_input_and_script_manifest.csv"))
  writeLines(capture.output(sessionInfo()), file.path(support_dir, "00_session_info.txt"))
  summary_text <- c(paste("MC19 /", analysis_scope, "/ 7项当期环境"), paste("输入：", input_file),
    paste("原始样本：", nrow(input$community), "；有效TS：", nrow(input$batches)),
    "每个样本保留一行；无Repeat_group或任何样本平均；无Host_type模型。",
    "原始Ruegeria比例逐元素平方根后计算Bray–Curtis；19个MC不重新归一化。",
    "每种TS分配的R²、F、平方和与自由度均直接读取vegan::adonis2(permutations=0)。",
    paste("所有环境共用", nrow(orders), "种TS行排列，含原始顺序。"),
    "p值=所有TS排列中F不小于实际F的比例；7项环境p值统一BH校正。",
    "低于检出限采用输入表中的L/2；保留原始报告值与删失标记。",
    "此置换推断依赖TS环境标签在零假设下可交换；结果解释为探索性关联。",
    "TS8当期环境缺测，不补值；不运行季节、宿主或两两比较；作图由04读取保存结果。",
    "顶层environment_results.csv为主要结果；完整结果、输入、置换记录、日志和快照均在supporting_files。", "",
    capture.output(print(environment_results[, c("Factor", "N_samples", "N_TS", "R2_raw",
      "F_raw", "P_raw", "P_BH_raw", "Significance")], row.names=FALSE)))
  writeLines(enc2utf8(summary_text), file.path(support_dir, "00_summary.txt"), useBytes=TRUE)
  writeLines(c("complete", analysis_version, output_dir), status_file)
  # 使用MC19专属指针；不更改旧MC8结果指针。
  mc19_publish_result(project_dir, "PERMANOVA", analysis_scope, output_dir)
  run_success <- TRUE
  print(environment_results[, c("Factor", "N_samples", "N_TS", "R2_raw", "F_raw", "P_raw", "P_BH_raw")], row.names=FALSE)
  message("完成。主要结果：", file.path(output_dir, "environment_results.csv"))
}, error=function(e) {
  writeLines(c("failed", analysis_version, conditionMessage(e)), status_file)
  stop(conditionMessage(e), call.=FALSE)
}, interrupt=function(e) {
  writeLines(c("interrupted", analysis_version), status_file)
  stop("运行已中断，保留部分输出；最新完成指针未推进。", call.=FALSE)
})
