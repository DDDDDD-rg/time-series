# 描述性宿主高亮图：只读04保存的同一套二维坐标；不做宿主检验。
# 没有坐标缓存时停止，请先运行04。
# 结果/坐标路径由共用helper解析，兼容supporting_files和旧版平铺目录。
project_dir <- "C:/Users/xzy/Desktop/CUHK_Study/GNBF6010_Research_Projects/Research/Msc main project-time series"
source(file.path(project_dir,"scripts","permanova_MC19_data_settings.R"),encoding="UTF-8",local=TRUE)
source(file.path(project_dir,"scripts","permanova_MC19_plot_helpers.R"),encoding="UTF-8",local=TRUE)
plot_mc19_suite("host")
