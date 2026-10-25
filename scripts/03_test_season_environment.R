# 干湿季和环境的统一入口：保留每个Repeat_group，不做第二次平均。
# all_hosts或coral在01_data_settings.R中设置。
# 完整分析代码在下方引用的脚本；两个03入口运行的是同一个版本。
project_dir <- "C:/Users/xzy/Desktop/CUHK_Study/GNBF6010_Research_Projects/Research/Msc main project-time series"
source(file.path(project_dir, "scripts", "03_test_season_environment_batch_labels.R"), encoding = "UTF-8", local = TRUE)
