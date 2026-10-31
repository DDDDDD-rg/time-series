# Select scope before sourcing: options(mc19.analysis_scope="ecosystem")
# Default: coral. Explicit result directories must match this scope.
# 7项环境的逐样本MC组成图；每柱一个原始样本，19个MC不重新归一化。
# 其余Ruegeria补为Other_Ruegeria；不是样本均值。
# 结果路径由共用helper解析，兼容supporting_files和旧版平铺目录。
project_dir <- "C:/Users/xzy/Desktop/CUHK_Study/GNBF6010_Research_Projects/Research/Msc main project-time series"
source(file.path(project_dir,"scripts", "PERMANOVA","permanova_MC19_data_settings.R"),encoding="UTF-8",local=TRUE)
source(file.path(project_dir,"scripts", "PERMANOVA","permanova_MC19_plot_helpers.R"),encoding="UTF-8",local=TRUE)
plot_mc19_suite("composition")
