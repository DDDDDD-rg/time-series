# HISTORICAL bootstrap 153110 replot only. Current Wald plots use plot_GLMM_MC19.R.
# No fitting, simulation or LOTO execution; existing historical results are preserved.
local({
  project_dir <- "C:/Users/xzy/Desktop/CUHK_Study/GNBF6010_Research_Projects/Research/Msc main project-time series"
  results_dir <- "" # Full-data inference's dedicated latest, or an explicit full-data directory.
  engine <- new.env(parent = globalenv())
  for (helper in c("glmm_MC19_pilot_helpers.R", "glmm_MC19_loto_helpers.R",
      "glmm_MC19_inference_helpers.R")) sys.source(file.path(project_dir, "scripts", helper), envir = engine)
  cfg <- engine$mc19_pilot_config(project_dir)
  engine$mc19_inference_plots(cfg, results_dir)
})
