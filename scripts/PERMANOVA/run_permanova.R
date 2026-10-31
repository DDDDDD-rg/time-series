# MC19 coral/ecosystem完整入口：用户手动source后，先做7项环境PERMANOVA，再读取保存结果作图。
(function() {
  # 通常只改这两项：FALSE可只分析；可选3D需要plotly和htmlwidgets。
  analysis_scope <- getOption("mc19.analysis_scope", "coral")  # coral / ecosystem
  old_scope <- options(mc19.analysis_scope=analysis_scope)
  on.exit(options(old_scope), add=TRUE)
  make_plots <- TRUE
  make_interactive_3d <- FALSE
  project_dir <- "C:/Users/xzy/Desktop/CUHK_Study/GNBF6010_Research_Projects/Research/Msc main project-time series"
  packages <- c("readxl","vegan",if(make_plots) "ggplot2",
                if(make_plots && make_interactive_3d) c("plotly","htmlwidgets"))
  missing <- packages[!vapply(packages,requireNamespace,logical(1),quietly=TRUE)]
  if (length(missing)) stop("请先安装R包：",paste(missing,collapse=", "))
  worker <- new.env(parent=globalenv())
  source(file.path(project_dir,"scripts", "PERMANOVA","03_test_season_environment.R"),encoding="UTF-8",local=worker)
  if (!isTRUE(worker$run_success)) stop("分析尚未完整完成，未开始作图。")
  pipeline_status <- file.path(worker$support_dir,"00_pipeline_status.txt")
  writeLines("analysis_complete",pipeline_status)
  if (make_plots) {
    old <- options(mc19.result_dir=worker$output_dir,mc19.interactive_3d=make_interactive_3d)
    on.exit(options(old),add=TRUE)
    tryCatch({
      source(file.path(project_dir,"scripts", "PERMANOVA","04_plot_permanova.R"),encoding="UTF-8",local=worker)
      writeLines("analysis_and_plots_complete",pipeline_status)
    },error=function(e) {
      writeLines(c("analysis_complete_plot_failed",conditionMessage(e)),pipeline_status)
      stop("分析结果已保存；作图未完成。修正后只运行04即可，无需重跑分析。原因：",conditionMessage(e),call.=FALSE)
    })
  } else {
    writeLines("analysis_only_complete",pipeline_status)
  }
})()
