# 当前MC19 coral分析与作图

2026-10-27核对：v2分析与作图已于2026-10-26完成，当前运行编号20261026_224333。统计量来源为vegan；本轮整理目录与发布脚本，不重跑完整置换。

在RStudio中打开run_permanova.R，点击Source。不要依次运行所有编号脚本。

首次缺包时，在R控制台运行：
```r
install.packages(c("readxl", "vegan", "ggplot2"))
```

- 默认run_permanova.R先分析、再生成PNG和PDF；make_plots=FALSE时只分析。
- 只重新作图：运行04_plot_permanova.R，读取MC19专属最新完成结果，不重跑PERMANOVA或BH。
- 只重画宿主描述性高亮：运行07_plot_host_highlight_panels.R，复用04保存的二维坐标。
- 只重画7项环境组成图：运行plot_MC_environment_composition.R；旧plot_DO_composition_preview.R转到同一新入口。
- 如需交互3D，将run_permanova.R中的make_interactive_3d改为TRUE，另装plotly、htmlwidgets；默认关闭。
- 分析完成但图报错时，分析结果保留；排查后只运行04即可。

## 统计设置

- 输入：data/PERMANOVA_input_MC19_env7_20261026.xlsx；当前完整环境coral为98个样本、19个MC、7个有效TS。
- 每个原始样本保留一行；不使用Repeat_group、不平均、不加入Host_type，不执行宿主、TS、季节或两两检验。
- MC比例逐元素平方根后由vegan::vegdist计算Bray–Curtis，19个MC不重新归一化。
- 7项环境分别建立单因子模型。7个TS的环境行共用同一套5040种完整排列，包含实际分配；每次把TS环境值扩展到其原始样本，MC样本行固定。
- 每一种分配调用vegan::adonis2(distance ~ Env_value, permutations=0, by=NULL)，直接读取R²、F、平方和及自由度。脚本不再手写这些统计量的计算公式，也不使用adonis2的样本行置换P值。
- P值保持为全部TS排列中F不小于实际F的比例，含实际分配，保留重复环境值排列的重数；比较容差沿用v1。7项P值仍统一由stats::p.adjust(method="BH")校正。
- 盐度、水温、pH、硝酸盐氮、亚硝酸盐氮、氨氮、DIP清单及检测限L/2规则不变；TS8当期环境缺测仍排除。
- 每项环境会拟合5040次vegan模型，7项共35280次，因此会比v1矩阵实现慢；每完成500次显示进度。

## 结果目录

新分析运行后，每个results运行目录的顶层只有主表和一个配套文件夹：

```text
results/PERMANOVA_MC19_env7_coral_samples_时间戳/
├─ environment_results.csv          # 主要结果：7行，含R²、F、P、BH q和样本数
└─ supporting_files/
   ├─ environment_results_full.csv  # 完整统计量、模型与输入指纹等
   ├─ 00_summary.txt
   ├─ 00_run_status.txt
   ├─ 00_pipeline_status.txt        # 完整入口运行时写入
   ├─ 00_input_used_MC19_env7.xlsx
   ├─ coral_environment_analysis.rds
   ├─ environment_analysis_input.csv
   ├─ environment_TS_values.csv
   ├─ coral_sample_coverage.csv
   ├─ TS_label_orders.csv
   ├─ permutation_F_values.csv
   ├─ environment_correlations_TS.csv
   ├─ 00_analysis_version.txt
   ├─ 00_input_and_script_manifest.csv
   ├─ 00_session_info.txt
   ├─ scripts_used/
   └─ plot_NMDS2D.rds               # 作图时产生；可选3D另有缓存
```

主表保留Indicator、Factor、N_samples、N_MC、N_TS、R2_raw、F_raw、P_raw、P_BH_raw、Significance；平方和、自由度、排列数量等完整字段放在配套表与RDS。作图直接读完整RDS，不依赖主表删减的列。

新计算版本标记为MC19_env7_coral_individual_TSlabels_v2_vegan。作图入口同时兼容v1平铺目录及v2的supporting_files目录，缓存按其实际来源版本检查；不会把旧结果标记成新计算结果。当前结果保留原路径；MC8和MC19 v1历史结果/图已归档，原路径保留兼容联接；MC19专属latest指针继续指向同一运行根目录。

## 配套图与入口

图仍输出到figures/同一运行名/，不自动打开窗口。包括7项R²/q值、TS着色NMDS、7项环境着色NMDS、环境时间线、19MC原始样本点、MC×样本热图、宿主描述性高亮和7项环境排序组成柱。显著性读取保存的BH结果。

- 完整入口：run_permanova.R → 03_test_season_environment.R → 03_test_season_environment_batch_labels.R。
- 共用统计与路径：permanova_MC19_data_settings.R；共用作图读取器：permanova_MC19_plot_helpers.R。
- 旧01仅作MC8历史设置；02/05/06/08/09历史分支的停止保护保留。
- 本次修改前的脚本与导航备份：archive/script_backups/MC19_before_vegan_statistics_20261026/。
- 完整TS标签置换仍依赖时间标签可交换；本次只替换统计量计算实现与结果目录布局。

## 图目录与自动合并

```text
figures/PERMANOVA_MC19_env7_coral_samples_时间戳/
├─ 00_1_总体结果与NMDS.pdf         # 10页
├─ 00_2_环境与MC时间变化.pdf     # 3页
├─ 00_3_环境梯度下MC组成.pdf     # 7页
└─ supporting_files/            # 单图PDF/PNG、坐标、排序、来源清单、日志和脚本快照
```

04、07和环境组成图仍使用同一套读取器，单图统一写入supporting_files。每次作图完成后调用bundle_MC19_figures.py刷新合并PDF；只重画某个分支时，刷新具有完整来源的分组。原页拼接保留矢量内容，来源缺失时不覆盖已有合并文件。图latest指针始终保存运行根目录。

合并需要Python和pypdf。本机自动检测可用Python（包含已安装的Codex附带解释器）；也可设置`options(mc19.python="解释器完整路径")`或环境变量`MC19_PYTHON`。缺少工具会明确停止作图；修复后只运行04，不需重跑分析。可单独执行`python scripts/bundle_MC19_figures.py --figure-dir "图运行目录"`，兼容旧版平铺和新supporting_files来源。

当前已完成结果与本轮验证详见[更新日志](../CHANGELOG.md)。

以下为原有脚本导航，包含尚未迁移的MC8历史模块：

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
