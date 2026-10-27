# MC19硝酸盐先导：逐TS留一敏感性；由用户运行。
# 读取既有8个完整模型，各留出TS1–TS7，共新增56次拟合。
# 不重拟合完整模型，不修改input或原结果，不模拟，不进行P/q校准。
# 历史56次已完成，无需重跑；当前GLMM.R运行19MC/env4的Wald P/BH，不自动运行本入口。
# LOTO翻转次数只作描述，不筛除模型，不决定TS随机项的去留。

local({
  project_dir <- "C:/Users/xzy/Desktop/CUHK_Study/GNBF6010_Research_Projects/Research/Msc main project-time series"
  source_results_dir <- file.path(project_dir, "results",
    "GLMM_MC19_nitrate_pilot_20261027_120255")
  make_plots <- TRUE

  engine <- new.env(parent = globalenv())
  sys.source(file.path(project_dir, "scripts", "glmm_MC19_pilot_helpers.R"), envir = engine)
  sys.source(file.path(project_dir, "scripts", "glmm_MC19_loto_helpers.R"), envir = engine)
  cfg <- engine$mc19_pilot_config(project_dir)
  run_dir <- engine$mc19_loto_run(cfg, source_results_dir)
  if (make_plots) tryCatch(engine$mc19_loto_plots(cfg, run_dir), error = function(e) {
    engine$mc19_write_text(conditionMessage(e),
      file.path(run_dir, "supporting_files", "plot_error.txt"))
    engine$mc19_write_text("analysis_complete_plot_failed",
      file.path(run_dir, "supporting_files", "run_status.txt"))
    warning("56次留一拟合已结束，作图需检查：", conditionMessage(e), call. = FALSE)
  })
  cat("\n逐TS留一流程结束。主表：", file.path(run_dir, "00_result.csv"),
    "\n完整估计来自原8模型；新增56次敏感性拟合，P/q未校准。\n", sep = "")
  invisible(run_dir)
})
