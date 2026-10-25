# MC8完整入口：含环境全样本 + 仅珊瑚 + HOST/Seawater/Sediment，检验、两两比较和作图。
# 在本地R/RStudio中source本文件即可。每次建立独立目录；不自动打开图窗/侧边栏。
(function() {
  project_dir <- "C:/Users/xzy/Desktop/CUHK_Study/GNBF6010_Research_Projects/Research/Msc main project-time series"
  input_file <- file.path(project_dir, "data", "PERMANOVA_input_MC8_20261011.xlsx")
  selected_mc <- c("MC6", "MC507", "MC21", "MC520", "MC10", "MC15", "MC0", "MC25")
  scopes <- c("all_hosts", "coral")
  run_id <- Sys.getenv("MC_PERMANOVA_RUN_ID", unset = format(Sys.time(), "%Y%m%d_%H%M%S"))
  if (!grepl("^[0-9]{8}_[0-9]{6}$", run_id)) stop("无效运行编号。")
  run_name <- paste0("PERMANOVA_MC8_", run_id)
  result_root <- file.path(project_dir, "results", run_name)
  figure_root <- file.path(project_dir, "figures", run_name)
  if (dir.exists(result_root) || dir.exists(figure_root)) stop("本次结果目录已存在。")
  required <- c("readxl", "dplyr", "vegan", "permute", "ggplot2", "viridisLite", "plotly", "htmlwidgets", "htmltools")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("缺少R包：", paste(missing, collapse = ", "))
  stopifnot(file.exists(input_file), identical(names(readxl::read_excel(input_file,
    sheet = "Community", n_max = 0))[-1L], selected_mc))
  input_md5 <- unname(tools::md5sum(input_file))
  dir.create(result_root, recursive = TRUE)
  dir.create(figure_root, recursive = TRUE)
  writeLines("MC8_env_inclusive_v1",file.path(result_root,"00_analysis_version.txt"))
  stopifnot(file.copy(input_file, file.path(result_root, "00_input_used_MC8.xlsx")))
  writeLines(c("running", result_root, figure_root), file.path(result_root, "00_run_status.txt"))
  start_time <- Sys.time()
  completed <- FALSE
  on.exit(if(!completed) writeLines(c("failed_or_interrupted",result_root,
    "See per-script logs / partial outputs; latest completed pointer was not advanced."),
    file.path(result_root,"00_run_status.txt")),add=TRUE)
  option_keys <- c("mc.analysis_scope", "mc.results_root", "mc.figures_root", "mc.figure_dir", "mc.batch_output_dir", "mc.pairwise_figure_dir", "mc.nmds_diagnostics")
  old_options <- setNames(lapply(option_keys, getOption), option_keys)
  on.exit(options(old_options), add = TRUE)
  script_names <- c("01_data_settings.R","02_test_host_time.R","03_test_season_environment.R",
    "03_test_season_environment_batch_labels.R","04_plot_permanova.R","05_test_pairwise_host_time.R",
    "06_plot_pairwise_permanova.R","07_plot_host_highlight_panels.R","08_test_sample_type.R",
    "permanova_pairwise_helpers.R","permanova_plot_style.R","pairwise_letters.R",
    "plot_DO_composition_preview.R","run_permanova.R")
  snapshot_dir <- file.path(result_root,"scripts_used")
  dir.create(snapshot_dir)
  stopifnot(all(file.copy(file.path(project_dir,"scripts",script_names),snapshot_dir)))
  all_results <- all_pairs <- pair_summaries <- all_plots <- list()
  write_csv <- function(data, file) {
    utils::write.csv(data, file, row.names = FALSE, fileEncoding = "UTF-8")
    content <- readBin(file, "raw", n = file.info(file)$size)
    writeBin(c(as.raw(c(239, 187, 191)), content), file)
  }
  esc <- function(x) {
    x <- gsub("&", "&amp;", as.character(x), fixed = TRUE)
    x <- gsub("<", "&lt;", x, fixed = TRUE)
    x <- gsub(">", "&gt;", x, fixed = TRUE)
    gsub('"', "&quot;", x, fixed = TRUE)
  }
  html_table <- function(data) {
    head <- paste0("<tr>", paste0("<th>", esc(names(data)), "</th>", collapse = ""), "</tr>")
    rows <- vapply(seq_len(nrow(data)), function(i)
      paste0("<tr>", paste0("<td>", esc(unlist(data[i, ], use.names = FALSE)), "</td>", collapse = ""), "</tr>"), character(1))
    paste0("<div class='table-wrap'><table><thead>", head, "</thead><tbody>", paste(rows, collapse = "\n"), "</tbody></table></div>")
  }
  p_text <- function(p) ifelse(is.finite(p), formatC(p, format = "g", digits = 4), "—")
  labels <- c(Host_type="宿主（调整TS）", TS="TS（调整宿主）", Season="干湿季（调整宿主）",
    Salinity_psu="盐度", SST_C="水温", pH="pH", DO_mg_L="溶解氧", Chla_ug_L="叶绿素a",
    DIN_mgN_L="DIN", PSM_mg_L="悬浮固体", Turbidity_NTU="浊度")
  for (scope in scopes) {
    scope_result <- file.path(result_root, scope)
    scope_figure <- file.path(figure_root, scope)
    dir.create(scope_result)
    options(mc.analysis_scope = scope, mc.results_root = result_root,
      mc.figures_root = figure_root, mc.figure_dir = scope_figure,
      mc.batch_output_dir = file.path(scope_result, "repeat_groups"), mc.nmds_diagnostics=FALSE)
    run_env <- new.env(parent = globalenv())
    for (script in c("02_test_host_time.R", "03_test_season_environment.R",
      "05_test_pairwise_host_time.R", "04_plot_permanova.R", "plot_DO_composition_preview.R")) {
      message("[", scope, "] 开始 ", script, " · ", format(Sys.time(), "%H:%M:%S"))
      log_file <- file.path(scope_result, paste0(sub("\\.R$", "", script), ".log"))
      connection <- file(log_file, open = "wt", encoding = "UTF-8")
      output_depth <- sink.number()
      sink(connection, type = "output")
      sink(connection, type = "message")
      failure <- tryCatch({
        source(file.path(project_dir, "scripts", script), encoding = "UTF-8", local = run_env)
        NULL
      }, error = function(e) e)
      sink(type = "message")
      while (sink.number() > output_depth) sink(type = "output")
      close(connection)
      if (!is.null(failure)) {
        writeLines(paste("failed:", scope, script, conditionMessage(failure)), file.path(result_root, "00_run_status.txt"))
        stop(scope, " / ", script, "：", conditionMessage(failure), "；日志：", log_file)
      }
      message("[", scope, "] 完成 ", script, " · ", format(Sys.time(), "%H:%M:%S"))
    }
    host <- readRDS(file.path(scope_result, "host_time_analysis.rds"))
    environmental <- readRDS(file.path(scope_result, "repeat_groups", "season_environment_analysis.rds"))
    pairwise <- readRDS(file.path(scope_result, "pairwise", "pairwise_host_time_analysis.rds"))
    combined <- dplyr::bind_rows(host$results, environmental$season, environmental$environment)
    stopifnot(nrow(combined) == 11L, all(combined$Status == "ok"), all(combined$N_MC == 8L),
      all(combined$Input_MD5 == input_md5), identical(grep("^MC", names(host$data), value = TRUE), selected_mc))
    combined$P_report <- ifelse(combined$Factor == "Season", combined$P_raw, combined$P_BH_raw)
    combined$Adjustment <- ifelse(combined$Factor == "Season", "单项检验，未校正",
      ifelse(combined$Factor %in% c("Host_type", "TS"), "宿主/TS两项BH", "环境八项BH"))
    all_results[[scope]] <- combined
    for (effect in c("Host_type", "TS")) {
      family <- pairwise$families[[effect]]
      results <- family$results
      stopifnot(nrow(results) == choose(length(family$levels), 2L), all(results$N_MC == 8L),
        all(results$Input_MD5 == input_md5), !any(grepl("^error:", results$Status)))
      all_pairs[[paste(scope, effect)]] <- results
      pair_summaries[[paste(scope, effect)]] <- data.frame(Scope = scope, Factor = effect,
        Planned_pairs = nrow(results), Completed_pairs = sum(is.finite(results$P_BH_raw)),
        Skipped_pairs = sum(!is.finite(results$P_BH_raw)),
        Significant_pairs_BH = sum(results$P_BH_raw <= 0.05, na.rm = TRUE),
        Groups_without_letters = paste(family$letters$excluded, collapse = ";"))
    }
    all_plots <- c(all_plots, setNames(run_env$plots, paste(scope, names(run_env$plots), sep = ":")))
  }
  options(mc.analysis_scope="sample_type",mc.results_root=result_root)
  sample_env <- new.env(parent=globalenv())
  source(file.path(project_dir,"scripts","08_test_sample_type.R"),encoding="UTF-8",local=sample_env)
  sample_saved <- readRDS(file.path(result_root,"sample_type","sample_type_analysis.rds"))
  all_pairs[["sample_type Sample_type"]] <- sample_saved$family$results
  pair_summaries[["sample_type Sample_type"]] <- data.frame(Scope="sample_type",Factor="Sample_type",
    Planned_pairs=nrow(sample_saved$family$results),
    Completed_pairs=sum(is.finite(sample_saved$family$results$P_BH_raw)),
    Skipped_pairs=sum(!is.finite(sample_saved$family$results$P_BH_raw)),
    Significant_pairs_BH=sum(sample_saved$family$results$P_BH_raw<=0.05,na.rm=TRUE),
    Groups_without_letters="")
  if(any(grepl("^error:",sample_saved$family$results$Status))) stop("Three-category pairwise test failed; inspect sample_type outputs.")
  # Four detailed host/TS heatmaps plus the new three-category heatmap.
  options(mc.pairwise_figure_dir=file.path(figure_root,"pairwise_heatmaps"))
  heatmap_env <- new.env(parent=globalenv())
  source(file.path(project_dir,"scripts","06_plot_pairwise_permanova.R"),encoding="UTF-8",local=heatmap_env)
  overall <- dplyr::bind_rows(all_results,sample_saved$result)
  pair_summary <- dplyr::bind_rows(pair_summaries)
  pairs <- dplyr::bind_rows(all_pairs)
  first <- c("Scope", "Factor", "N_MC", "N_units", "N_hosts", "N_TS", "F_raw", "R2_raw", "P_raw",
    "P_BH_raw", "P_report", "Adjustment", "Significance", "P_BH_lingoes", "P_disp_raw", "P_disp_lingoes", "Status", "Warnings")
  overall <- overall[, c(first, setdiff(names(overall), first))]
  write_csv(overall, file.path(result_root, "00_all_tests_summary.csv"))
  write_csv(pair_summary, file.path(result_root, "00_pairwise_summary.csv"))
  write_csv(pairs, file.path(result_root, "01_all_pairwise_results.csv"))
  write_csv(pairs[is.finite(pairs$P_BH_raw) & pairs$P_BH_raw <= 0.05, ], file.path(result_root, "01_significant_pairwise_results.csv"))
  pdf_device <- if(isTRUE(capabilities("cairo"))) grDevices::cairo_pdf else grDevices::pdf
  pdf_device(file.path(figure_root, "00_all_final_results.pdf"), width = 12, height = 9.6,
    onefile=TRUE,family=if(isTRUE(capabilities("cairo"))) "Arial" else "sans")
  tryCatch({ for (p in all_plots) print(if(isTRUE(capabilities("cairo"))) p else p + ggplot2::theme(text=ggplot2::element_text(family="sans"))) }, finally = grDevices::dev.off())
  figure_files <- list.files(figure_root, recursive = TRUE, full.names = TRUE)
  stopifnot(length(all_plots) >= 6L, all(file.info(figure_files)$size > 0))
  key_results <- c("MC8 PERMANOVA：本次完整结果", paste("MC：", paste(selected_mc, collapse = ", ")),
    paste("输入：", input_file), paste("Input MD5：", input_md5),
    "所有MC零值保留；原始相对丰度不按所选8个MC重新归一化。",
    "分析单位为Repeat_group均值；平方根转换后计算Bray–Curtis。",
    "宿主/TS两项BH；每个范围的环境八项BH；季节为单项原始p。", "")
  html_sections <- character()
  for (scope in scopes) {
    x <- all_results[[scope]]
    title <- if (scope == "all_hosts") "全样本（含Seawater、Sediment）" else "仅珊瑚"
    env_significant <- x$Factor[!x$Factor %in% c("Host_type", "TS", "Season") & x$P_report <= 0.05]
    key_results <- c(key_results, paste0("【", title, "】"),
      vapply(seq_len(nrow(x)), function(i) sprintf("%s：R2=%.4f；p=%s（%s）；N=%d，TS=%d。",
        labels[x$Factor[i]], x$R2_raw[i], p_text(x$P_report[i]), x$Adjustment[i], x$N_units[i], x$N_TS[i]), character(1)),
      paste("BH显著环境指标：", if (length(env_significant)) paste(labels[env_significant], collapse = "、") else "无"), "")
    display <- data.frame(检验 = unname(labels[x$Factor]), 分析单位数 = x$N_units, TS数 = x$N_TS,
      R2 = sprintf("%.4f", x$R2_raw), 原始p = p_text(x$P_raw), 报告p = p_text(x$P_report),
      校正范围 = x$Adjustment, 显著性 = x$Significance, PERMDISP原始p = p_text(x$P_disp_raw), check.names = FALSE)
    relative_figures <- paste0("../../figures/", run_name, "/", scope, "/")
    pngs <- list.files(file.path(figure_root, scope), pattern = "_final\\.png$")
    gallery <- paste(vapply(pngs, function(name) paste0("<a href='", relative_figures, name,
      "'><img loading='lazy' src='", relative_figures, name, "' alt='", name, "'></a>"), character(1)), collapse = "\n")
    html_sections <- c(html_sections, paste0("<section><h2>", title, "</h2>", html_table(display),
      "<p><a href='", relative_figures, "final_results.pdf'>本范围展示图PDF</a> · ",
      "<a href='", scope, "/pairwise/pairwise_host_results.csv'>宿主两两比较</a> · ",
      "<a href='", scope, "/pairwise/pairwise_TS_results.csv'>TS两两比较</a> · ",
      "<a href='", relative_figures, "DO_MC_composition_preview.png'>DO梯度组成图</a></p><div class='gallery'>", gallery, "</div></section>"))
  }
  broad <- sample_saved$result
  html_sections <- c(html_sections,paste0("<section><h2>HOST / Seawater / Sediment</h2>",
    html_table(broad[,c("Factor","N_units","R2_raw","P_raw","P_disp_raw","Status")]),
    html_table(sample_saved$family$results[,c("Group1","Group2","N1","N2","R2_raw","P_raw","P_BH_raw","Status")]),
    "<p>HOST合并全部生物宿主；保留所有117个已观测Repeat_group。模型调整TS；在TS内置换。三对比较单独BH。</p></section>"))
  key_results <- c(key_results,sprintf("【HOST / Seawater / Sediment】N=%d；R2=%.4f；原始p=%s。三对比较另做BH。",
    broad$N_units,broad$R2_raw,p_text(broad$P_raw)), "【两两比较汇总】", capture.output(print(pair_summary, row.names = FALSE)), "",
    "分类因素的PERMDISP提示需与组成检验一并阅读；连续环境指标不人为分组。",
    "宿主/TS图不显示两两字母；两两结果由五张独立热图展示。",
    "环境使用7个批次；TS8缺测仍保留于宿主/TS/季节分析。",
    paste("结果目录：", result_root), paste("图目录：", figure_root),
    paste("总耗时（分钟）：", round(as.numeric(difftime(Sys.time(), start_time, units = "mins")), 2)))
  writeLines(enc2utf8(key_results), file.path(result_root, "00_key_results.txt"), useBytes = TRUE)
  html <- paste0("<!doctype html><html lang='zh-CN'><meta charset='UTF-8'><meta name='viewport' content='width=device-width,initial-scale=1'>",
    "<title>MC8 PERMANOVA汇总</title><style>body{font:15px/1.65 system-ui,sans-serif;color:#233443;background:#f4f6f8;margin:0}main{max-width:1400px;margin:auto;padding:28px}section{background:white;padding:24px;margin:22px 0;border-radius:10px}h1,h2{color:#173954}.table-wrap{overflow:auto}table{border-collapse:collapse;width:100%;font-variant-numeric:tabular-nums;font-size:14px}th,td{padding:9px 11px;border-bottom:1px solid #e3e7ec;text-align:left;white-space:nowrap}th{background:#eaf0f5}tr:nth-child(even){background:#f8fafb}a{color:#1261a0}.gallery{display:grid;grid-template-columns:repeat(auto-fit,minmax(360px,1fr));gap:18px}.gallery img{width:100%;border:1px solid #e3e7ec;border-radius:6px}.meta{color:#536474}</style><main>",
    "<h1>MC8 PERMANOVA · 完整结果</h1><p class='meta'>", run_name, " · ", paste(selected_mc, collapse = " · "),
    "</p><p>400个输入样本；按Repeat_group取均值后分析。全样本含环境：115组、14类（当前环境指标101组）；仅珊瑚39组（环境33组）；三类样本比较117组。</p>",
    "<p><a href='00_all_tests_summary.csv'>总体与环境检验完整数值</a> · <a href='00_key_results.txt'>中文结论摘要</a> · ",
    "<a href='01_significant_pairwise_results.csv'>显著两两比较</a> · ",
    "<a href='../../figures/", run_name, "/00_all_final_results.pdf'>全部展示图PDF</a></p>",
    "<section><h2>两两比较完成情况</h2>", html_table(pair_summary),
    "<p>宿主对与TS对分别BH校正；跳过的比较保留在完整表。图不显示两两字母；另保存五张两两比较热图。</p></section>",
    paste(html_sections, collapse = "\n"), "<section><h2>方法与数据版本</h2><p>Repeat_group均值 → 平方根相对丰度 → Bray–Curtis。",
    "宿主检验调整TS并在TS内置换；TS检验调整宿主并在宿主内置换（各9999次）。干湿季和环境调整宿主，完整枚举批次标签分配。",
    "宿主/TS两项BH，环境每个范围八项BH；季节是单项检验。原始Bray为主，Lingoes为敏感性检查。</p>",
    "<p>相对丰度保留原总Ruegeria分母，零值保留。连续环境指标不人为分组做PERMDISP。分类因素的PERMDISP结果见表格及原始输出。",
    "环境来自7个批次；TS8环境缺测。以前不含环境样本的历史结果保留；本轮另存新目录。</p><p>输入MD5：", input_md5, "</p></section></main></html>")
  writeLines(enc2utf8(html), file.path(result_root, "00_PERMANOVA_summary.html"), useBytes = TRUE)
  sources <- c(input_file,file.path(snapshot_dir,script_names))
  write_csv(data.frame(file = sources, MD5 = unname(tools::md5sum(sources))), file.path(result_root, "00_input_and_script_manifest.csv"))
  writeLines(enc2utf8(result_root), file.path(project_dir, "results", "00_latest_PERMANOVA_MC8.txt"), useBytes = TRUE)
  writeLines(enc2utf8(figure_root), file.path(project_dir, "figures", "00_latest_PERMANOVA_MC8.txt"), useBytes = TRUE)
  writeLines(c("complete", enc2utf8(result_root), enc2utf8(figure_root)), file.path(result_root, "00_run_status.txt"), useBytes = TRUE)
  completed <- TRUE
  cat(paste(key_results, collapse = "\n"), "\n\n优先查看：00_PERMANOVA_summary.html；完整数值：00_all_tests_summary.csv。\n")
})()
