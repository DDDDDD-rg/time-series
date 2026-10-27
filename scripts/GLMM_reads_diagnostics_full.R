# 当前MC19 x 4营养盐/TS-only Wald诊断。只读取保存的76个有TS模型，不重新拟合。
# 模拟来自已保存的Beta-binomial模型；TS随机截距重新抽样。
# 不重新拟合，不自动补零膨胀项，不置换，不bootstrap。
# GLMM.R默认调用本诊断功能；单独重做诊断时运行本文件。
# 原8先导/153110 bootstrap/LOTO结果不由本入口读写。

local({
  project_dir <- "C:/Users/xzy/Desktop/CUHK_Study/GNBF6010_Research_Projects/Research/Msc main project-time series"
  results_dir <- ""              # 空字符串读env4/TS-only Wald专用latest。
  n_sim <- 2000L
  simulation_seed <- 20261027L
  refresh_plots <- TRUE

  engine <- new.env(parent = globalenv())
  for (helper in c("glmm_MC19_pilot_helpers.R", "glmm_MC19_wald_helpers.R",
      "glmm_MC19_wald_plot_helpers.R")) sys.source(file.path(project_dir, "scripts", helper), envir = engine)
  cfg <- engine$mc19_wald_config(project_dir)
  cfg$n_sim <- n_sim
  cfg$seed <- simulation_seed
  run <- engine$mc19_wald_read_run(cfg, results_dir)
  engine$mc19_wald_diagnostics(cfg, run$run_dir)
  if (refresh_plots) engine$mc19_wald_plots(cfg, run$run_dir)
  cat("诊断更新完成；未重新拟合。\n")
})
