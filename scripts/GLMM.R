# 当前入口：19MC x 4营养盐，仅TS随机截距 = 76个模型。
# 固定公式：cbind(MC_reads, Other_reads) ~ E_z + (1 | TS)。
# 由用户运行。glmmTMB自带Wald P；全部76项统一BH；不bootstrap/置换/LOTO。
# 98个coral原始样本含TS4/TS6，不平均，无Host/Repeat/水温项。
# 默认核对后仅复用120255的4个有TS硝酸盐拟合，新增72个。

local({
  project_dir <- "C:/Users/xzy/Desktop/CUHK_Study/GNBF6010_Research_Projects/Research/Msc main project-time series"
  reuse_pilot_fits <- TRUE       # FALSE时从头拟合76个；不影响模型公式。
  run_diagnostics <- TRUE        # 读取保存模型做诊断，不重拟合。
  run_plots <- TRUE
  n_sim <- 2000L                 # 仅诊断模拟，绝不用于beta P值。
  simulation_seed <- 20261027L

  engine <- new.env(parent = globalenv())
  for (helper in c("glmm_MC19_pilot_helpers.R", "glmm_MC19_wald_helpers.R",
      "glmm_MC19_wald_plot_helpers.R")) sys.source(file.path(project_dir, "scripts", helper), envir = engine)
  cfg <- engine$mc19_wald_config(project_dir)
  cfg$reuse_pilot <- reuse_pilot_fits; cfg$n_sim <- n_sim; cfg$seed <- simulation_seed
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
    "\n98coral样本，76个有TS模型，Wald beta检验统一BH；未运行bootstrap、置换或LOTO。\n", sep = "")
  invisible(run_dir)
})
