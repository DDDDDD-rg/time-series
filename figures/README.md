# 图目录导航

2026-10-28核对：当前76模型TS-only Wald图已完成，位于[GLMM_MC19_env4_TS_Wald_coral_20261027_233529/](GLMM_MC19_env4_TS_Wald_coral_20261027_233529/)。顶层5份合并PDF：总览2页，NO3/NO2/NH4/DIP各21页；supporting_files/diagnostics内4份诊断PDF各39页。合计9份242页，与同批76项BH结果对应。原始点、固定部分/边际曲线与条件TS拟合值保留。[运行及重画说明](../scripts/README.md)。

以下为历史流程及结果记录，当前运行按上述TS-only入口执行。

2026-10-27当前Wald作图代码已改、未运行：GLMM_MC19_env4_Wald_coral_时间戳/顶层一份3页总览及四份各22页营养盐效应/曲线PDF；4份诊断补充在supporting_files/diagnostics。当前plot_GLMM_MC19.R只读152模型Wald结果，不拟合或模拟。历史153110 bootstrap图另用plot_GLMM_MC19_full_data.R。[重要代码说明](../docs/GLMM_MC19_env4_Wald_code_guide_20261027.md)。

2026-10-27：完整数据推断作图源码已修改、尚未运行。新推断结束后，输出`GLMM_MC19_nitrate_full_data_时间戳/00_full_data_inference.pdf`七页：全98样本、完整beta/近似区间、4MC曲线及beta P/q；注明高硝酸盐集中于TS4/TS6。重画入口`plot_GLMM_MC19_full_data.R`或默认模式的`plot_GLMM_MC19.R`，只读保存结果。重画原先导图则将后者plot_mode设为"pilot"。[运行说明](../scripts/README.md)。

逐TS留一已由用户完成：[140106五页对照图](GLMM_MC19_nitrate_LOTO_20261027_140106/00_LOTO_stability.pdf)。已核对4个MC页面；图注与横轴重叠的源码已修正，原PDF未重画。单独刷新用plot_GLMM_MC19_LOTO.R，不拟合或模拟；原120255两份PDF不覆盖。[说明](../scripts/README.md)。

2026-10-27：MC19 GLMM先导120255已由用户运行，图在[同名目录](GLMM_MC19_nitrate_pilot_20261027_120255/)，分别为7页效应/曲线、17页诊断PDF。QQ页重复标题的排版源码已修正，原PDF尚未重画；将plot_GLMM_MC19.R的plot_mode设为"pilot"后可刷新，不拟合、不模拟。[运行说明](../scripts/README.md)。

2026-10-27核对。当前PERMANOVA为MC19 coral / env7 / v2，运行编号20261026_224333。

| 当前图 | 入口 |
|---|---|
| 总体结果与NMDS（10页） | [00_1_总体结果与NMDS.pdf](PERMANOVA_MC19_env7_coral_samples_20261026_224333/00_1_总体结果与NMDS.pdf) |
| 环境与MC时间变化（3页） | [00_2_环境与MC时间变化.pdf](PERMANOVA_MC19_env7_coral_samples_20261026_224333/00_2_环境与MC时间变化.pdf) |
| 环境梯度下MC组成（7页） | [00_3_环境梯度下MC组成.pdf](PERMANOVA_MC19_env7_coral_samples_20261026_224333/00_3_环境梯度下MC组成.pdf) |
| 原始单图、PNG、坐标和来源清单 | [supporting_files/](PERMANOVA_MC19_env7_coral_samples_20261026_224333/supporting_files/) |
| 当前主要结果 | [environment_results.csv](../results/PERMANOVA_MC19_env7_coral_samples_20261026_224333/environment_results.csv) |

## 历史PERMANOVA

- [MC19 v1（2026-10-26 21:51）](../archive/PERMANOVA/figures/PERMANOVA_MC19_env7_coral_samples_20261026_215151/)。
- [MC8（2026-10-14）](../archive/PERMANOVA/figures/PERMANOVA_MC8_20261014_183535/)。
- [其他历史PERMANOVA图、两两热图与高亮图](../archive/PERMANOVA/figures/)。

## 其他已有图与可编辑源文件

| 内容 | 位置 |
|---|---|
| 既有GLMM宿主热图 | [GLMM_host_beta_heatmap_20261014_185246/](GLMM_host_beta_heatmap_20261014_185246/) |
| DO系数森林图 | [GLMM_summary/](GLMM_summary/) |
| 历史HOST/Seawater/Sediment图 | [sample_type_20261014_191318/](sample_type_20261014_191318/) |
| PM11时间序列 | [PM11_timeseries_2021_2025_surface_top_20261009_180827/](PM11_timeseries_2021_2025_surface_top_20261009_180827/) |
| 环境图及AI/Prism/EMF源文件 | [environment_assets/](environment_assets/) |
| 本目录新增Prism源文件 | [TS-ENV折线图.prism](TS-ENV折线图.prism) |

机器读取的latest指针保持原路径与内容并设为隐藏；历史目录旧位置保留隐藏兼容联接。HTML和同名_files依赖目录整组保存。[更新日志](../CHANGELOG.md)记录本次整理与验证。
