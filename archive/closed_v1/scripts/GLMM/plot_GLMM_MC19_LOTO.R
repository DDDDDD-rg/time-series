# 只重画已经完成的逐TS留一结果；不拟合，不模拟，不读取旧MC8指针。

local({
  project_dir <- "C:/Users/xzy/Desktop/CUHK_Study/GNBF6010_Research_Projects/Research/Msc main project-time series"
  loto_results_dir <- ""         # 空字符串读取LOTO专用latest；也可填LOTO运行目录。

  engine <- new.env(parent = globalenv())
  sys.source(file.path(project_dir, "scripts", "GLMM", "glmm_MC19_pilot_helpers.R"), envir = engine)
  sys.source(file.path(project_dir, "scripts", "GLMM", "glmm_MC19_loto_helpers.R"), envir = engine)
  cfg <- engine$mc19_pilot_config(project_dir)
  engine$mc19_loto_plots(cfg, loto_results_dir)
})
