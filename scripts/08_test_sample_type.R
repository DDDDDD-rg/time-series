# Broad comparison: HOST / Seawater / Sediment, adjusted for TS.
# HOST pools all observed biological hosts, including singleton host types.
# Sourcing this script performs tests ONLY when the user runs it.
project_dir <- "C:/Users/xzy/Desktop/CUHK_Study/GNBF6010_Research_Projects/Research/Msc main project-time series"
local({
  old_scope <- getOption("mc.analysis_scope")
  on.exit(options(mc.analysis_scope=old_scope),add=TRUE)
  options(mc.analysis_scope="sample_type")
  source(file.path(project_dir,"scripts","01_data_settings.R"),encoding="UTF-8",local=TRUE)
  source(file.path(project_dir,"scripts","permanova_pairwise_helpers.R"),encoding="UTF-8",local=TRUE)
  source(file.path(project_dir,"scripts","permanova_plot_style.R"),encoding="UTF-8",local=TRUE)
  dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
  stopifnot(nlevels(analysis_data$Sample_type)==3L)
  set.seed(random_seed)
  permutations <- permute::shuffleSet(nrow(analysis_data),control=permute::how(
    nperm=n_permutations,blocks=analysis_data$TS))
  result <- test_factor(analysis_data,c("Sample_type","TS"),"Sample_type",
    permutations,dispersion_group="Sample_type")
  result$Scope <- "sample_type"
  result$Analysis_unit <- "Repeat_group均值"
  result$Permutation_design <- "Sample_type adjusted for TS; profile permutations within TS"
  result$P_BH_raw <- result$P_BH_lingoes <- NA_real_
  result$P_report <- result$P_raw
  result$Adjustment <- "三类样本总体：单项原始p；三对比较另做BH"
  result$Significance <- add_stars(result$P_raw)
  write_csv_excel(result,file.path(output_dir,"sample_type_results.csv"),row.names=FALSE)
  write_csv_excel(analysis_data,file.path(output_dir,"sample_type_analysis_input.csv"),row.names=FALSE)
  write_csv_excel(as.data.frame.matrix(table(analysis_data$Sample_type,analysis_data$TS)),
    file.path(output_dir,"sample_type_by_TS_coverage.csv"))
  if(result$Status!="ok") stop("Three-category omnibus test failed: ",result$Status)
  family <- run_pairwise_family(analysis_data,"Sample_type",seed_offset=200000L)
  pairwise_dir <- file.path(output_dir,"pairwise")
  dir.create(pairwise_dir,showWarnings=FALSE)
  write_csv_excel(family$results,file.path(pairwise_dir,"pairwise_sample_type_results.csv"),row.names=FALSE)
  write_csv_excel(family$coverage,file.path(pairwise_dir,"pairwise_sample_type_coverage.csv"),row.names=FALSE)
  for(column in c("P_BH_raw","P_BH_lingoes","R2_raw")) write_csv_excel(
    make_pair_matrix(family$results,family$levels,column),
    file.path(pairwise_dir,paste0("pairwise_sample_type_",column,"_matrix.csv")))
  saveRDS(list(data=analysis_data,result=result,family=family,analysis_version=analysis_version),
    file.path(output_dir,"sample_type_analysis.rds"))
  writeLines(capture.output(sessionInfo()),file.path(output_dir,"sample_type_session.txt"))
  print(result[,c("Factor","N_units","R2_raw","P_raw","Status")])
  print(family$results[,c("Group1","Group2","N1","N2","R2_raw","P_BH_raw","Status")])
})
