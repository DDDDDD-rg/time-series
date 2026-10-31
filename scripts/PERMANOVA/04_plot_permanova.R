# Select scope before sourcing: options(mc19.analysis_scope="ecosystem")
# Default: coral. Explicit result directories must match this scope.
# 当前全部配套图：读取MC19 coral/ecosystem已完成结果；不重跑PERMANOVA或BH。
# 二维坐标统一缓存；输出原始样本点、7项环境图、组成图和描述性宿主分面。
# 兼容新supporting_files布局及旧版平铺结果；不需要为路径变化重跑分析。
project_dir <- "C:/Users/xzy/Desktop/CUHK_Study/GNBF6010_Research_Projects/Research/Msc main project-time series"
source(file.path(project_dir,"scripts", "PERMANOVA","permanova_MC19_data_settings.R"),encoding="UTF-8",local=TRUE)
source(file.path(project_dir,"scripts", "PERMANOVA","permanova_MC19_plot_helpers.R"),encoding="UTF-8",local=TRUE)
plot_mc19_suite(c("main","host","composition"))
