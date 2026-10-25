# R script entry points

当前分析脚本保留本目录，模型内容和计算设置没有因整理而修改。

| 任务 | 入口 |
|---|---|
| 完整PERMANOVA分析 | [run_permanova.R](run_permanova.R) |
| 当前Z_TS_SD_Q结构的GLMM拟合 | [GLMM.R](GLMM.R) |
| 读取已有GLMM结果并复核 | [GLMM_reads_diagnostics_full.R](GLMM_reads_diagnostics_full.R) |
| PERMANOVA总体图与交互3D HTML | [04_plot_permanova.R](04_plot_permanova.R) |
| 两两比较热图 | [06_plot_pairwise_permanova.R](06_plot_pairwise_permanova.R) |
| 宿主高亮分面 | [07_plot_host_highlight_panels.R](07_plot_host_highlight_panels.R) |
| GLMM宿主系数热图 | [08_plot_GLMM_host_beta_heatmap.R](08_plot_GLMM_host_beta_heatmap.R) |
| HOST/Seawater/Sediment图 | [09_plot_sample_type.R](09_plot_sample_type.R) |
| DO组成预览 | [plot_DO_composition_preview.R](plot_DO_composition_preview.R) |
| PM11环境时间序列 | [plot_PM11_timeseries_5years.R](plot_PM11_timeseries_5years.R) |

其余编号脚本和 `permanova_*`、`pairwise_letters.R` 为上述流程使用的单步分析或辅助函数。

RStudio项目 `scripts.Rproj`、`.RData`、`.Rhistory` 和 `.Rproj.user/` 保留原位置。
历史 `.bak`、before版本及备份目录集中在 [archive/script_backups/](../archive/script_backups/)。
