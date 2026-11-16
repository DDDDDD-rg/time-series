# Selected workflow: audited reuse for seven-TS temperature; legacy full eight-TS remains separate.
project_dir <- getOption("mc19.project_dir","C:/Users/xzy/Desktop/CUHK_Study/GNBF6010_Research_Projects/Research/Msc main project-time series")
if(identical(getOption("mc19.pipeline","TS8_env5"),"TS7_temperature")) {
  source(file.path(project_dir,"scripts","PERMANOVA","permanova_MC19_data_settings.R"),encoding="UTF-8",local=TRUE)
  source(file.path(project_dir,"scripts","mc19_temperature_incremental.R"),encoding="UTF-8",local=TRUE)
  input <- read_mc19_samples()
  output_dir <- mc19_incremental_permanova_run(project_dir,analysis_scope,input,input_file)
  support_dir <- file.path(output_dir,"supporting_files")
  run_success <- TRUE
} else source(file.path(project_dir,"scripts","PERMANOVA","03_test_season_environment_batch_labels.R"),encoding="UTF-8",local=TRUE)
