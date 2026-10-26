# 当前入口只运行coral的7项环境PERMANOVA；保留文件名便于原有使用习惯。
# 每个原始样本一行；不执行季节检验、Host_type调整或Repeat_group平均。
# R²/F由vegan拟合，脚本汇总TS标签置换P；主表在顶层，其余放supporting_files。
project_dir <- "C:/Users/xzy/Desktop/CUHK_Study/GNBF6010_Research_Projects/Research/Msc main project-time series"
source(file.path(project_dir, "scripts", "03_test_season_environment_batch_labels.R"), encoding="UTF-8", local=TRUE)
