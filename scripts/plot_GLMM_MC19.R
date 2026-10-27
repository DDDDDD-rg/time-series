# 当前MC19 x 4营养盐/TS-only Wald作图：只读保存的76个有TS模型及P/q，不拟合、不模拟。
# 5份主要合并PDF；诊断PDF位于supporting_files/diagnostics。
# 历史bootstrap图用plot_GLMM_MC19_full_data.R；历史LOTO图用plot_GLMM_MC19_LOTO.R。

local({
  project_dir <- "C:/Users/xzy/Desktop/CUHK_Study/GNBF6010_Research_Projects/Research/Msc main project-time series"
  results_dir <- ""              # 空字符串读env4/TS-only Wald独立latest。

  engine <- new.env(parent = globalenv())
  for (helper in c("glmm_MC19_pilot_helpers.R", "glmm_MC19_wald_helpers.R",
      "glmm_MC19_wald_plot_helpers.R")) sys.source(file.path(project_dir, "scripts", helper), envir = engine)
  cfg <- engine$mc19_wald_config(project_dir)
  engine$mc19_wald_plots(cfg, results_dir)
})
