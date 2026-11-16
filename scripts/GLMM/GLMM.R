# Single-indicator GLMM entry; profile selects TS and environmental indicators.
# Default: full eight-TS/five-nutrient preparation. TS7_temperature.R selects seven TS/six indicators.
# Fresh Beta-binomial/logit fits and BH across the selected full family; diagnostics/plots follow.

local({
  analysis_scope <- getOption("mc19.analysis_scope", "coral")  # coral / ecosystem
  project_dir <- getOption("mc19.project_dir","C:/Users/xzy/Desktop/CUHK_Study/GNBF6010_Research_Projects/Research/Msc main project-time series")
  run_diagnostics <- isTRUE(getOption("mc19.run_diagnostics",TRUE))        # 读取保存模型做诊断，不重拟合。
  run_plots <- isTRUE(getOption("mc19.run_plots",TRUE))
  n_sim <- 2000L                 # 仅诊断模拟，绝不用于beta P值。
  simulation_seed <- 20261027L

  engine <- new.env(parent = globalenv())
  for (helper in c("glmm_MC19_pilot_helpers.R", "glmm_MC19_wald_helpers.R",
      "glmm_MC19_wald_plot_helpers.R")) sys.source(file.path(project_dir, "scripts", "GLMM", helper), envir = engine)
  cfg <- engine$mc19_wald_config(project_dir, analysis_scope)
  cfg$n_sim <- n_sim; cfg$seed <- simulation_seed
  engine$mc19_require(c("readxl", "glmmTMB", if (run_diagnostics) "DHARMa"))
  run_dir <- engine$mc19_wald_run(cfg)
  notices <- character()
  if (run_diagnostics) tryCatch(engine$mc19_wald_diagnostics(cfg, run_dir),
    error = function(e) { notices <<- c(notices, paste("diagnostics:", conditionMessage(e))) })
  if (run_plots) tryCatch(engine$mc19_wald_plots(cfg, run_dir),
    error = function(e) { notices <<- c(notices, paste("plots:", conditionMessage(e))) })
  if (length(notices)) {
    engine$mc19_write_text(notices, file.path(run_dir, "supporting_files", "pipeline_notices.txt"))
    engine$mc19_write_text("Wald_results_saved_pipeline_needs_review", file.path(run_dir, "supporting_files", "run_status.txt"))
    warning(paste(notices, collapse = "\n"), call. = FALSE)
  }
  cat("\nWald流程结束。主表：", file.path(run_dir, "00_result.csv"),
    sprintf("\n%s: %d single-indicator models; %d TS; BH n=%d. No bootstrap/permutation/LOTO.\n",analysis_scope,cfg$n_models,cfg$expected_TS,cfg$q_family_size), sep = "")
  invisible(run_dir)
})
