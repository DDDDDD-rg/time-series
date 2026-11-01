# 当前图导航

更新：2026-10-31。主要展示新版talk图、原GLMM显著关系曲线及两部分分析图。

已集中整理：[MC19_阶段汇报结果图_v01](MC19_阶段汇报结果图_v01/00_阅读说明.md)。含39个PDF、28个PNG、2个交互3D HTML及其依赖，按八类内容分文件夹、以Coral/Ecosystem清楚命名；包内文件清单记录原路径和SHA-256。原输出保留原位。

| 内容 | Coral | Ecosystem |
|---|---|---|
| 19MC热图＋森林图 | [查看](MC19_talk_coral_20261031_153730/MC19_heatmap_forest.pdf) | [查看](MC19_talk_ecosystem_20261031_153730/MC19_heatmap_forest.pdf) |
| NO₂二维NMDS | [查看](MC19_talk_coral_20261031_153730/NMDS2D_Nitrite_mgN_L.pdf) | [查看](MC19_talk_ecosystem_20261031_153730/NMDS2D_Nitrite_mgN_L.pdf) |
| 氨氮二维NMDS | [查看](MC19_talk_coral_20261031_153730/NMDS2D_Ammonia_mgN_L.pdf) | [查看](MC19_talk_ecosystem_20261031_153730/NMDS2D_Ammonia_mgN_L.pdf) |
| NO₂／氨氮三维交互图（右侧切换） | [查看](MC19_talk_coral_20261031_153730/NMDS3D_nutrients.html) | [查看](MC19_talk_ecosystem_20261031_153730/NMDS3D_nutrients.html) |
| NMDS诊断 | [查看](MC19_talk_coral_20261031_153730/NMDS_diagnostics.csv) | [查看](MC19_talk_ecosystem_20261031_153730/NMDS_diagnostics.csv) |
| 原GLMM所有显著曲线 | [查看](GLMM_MC19_env4_TS_Wald_coral_20261027_233529/02_all_significant_nutrient_curves.pdf) | [查看](GLMM_MC19_env4_TS_Wald_ecosystem_20261029_123804/02_all_significant_nutrient_curves.pdf) |
| 两部分效应 | [查看](MC19_two_part_coral_legacy_plot001/two_part_effects.pdf) | [查看](MC19_two_part_ecosystem_legacy_plot001/two_part_effects.pdf) |
| 两部分曲线 | [查看](GLMM_MC19_two_part_coral_20261031_153730/two_part_curves.pdf) | [查看](GLMM_MC19_two_part_ecosystem_20261031_153730/two_part_curves.pdf) |
| 两部分残差审查 | [查看](GLMM_MC19_two_part_coral_20261031_153730/residual_review.pdf) | [查看](GLMM_MC19_two_part_ecosystem_20261031_153730/residual_review.pdf) |

- [N/P四指标随TS折线](MC19_talk_coral_20261031_153730/NP_by_TS.pdf)：两个范围共用水质表。
- Coral氨氮是对照面板。Ecosystem二维stress约.220且最佳解未重复；三维约.157已重复确认。三维轴对PDF/PNG是投影。HTML与依赖目录一起保存。
- [旧图归档与清单](../archive/obsolete_outputs_20261031_155512/README.md)：旧独立热图、旧主图备份、旧PERMANOVA完整展示套图。旧套图其它环境/组成等补充完整保留，原路径用隐藏联接兼容。
- [Coral原诊断](GLMM_MC19_env4_TS_Wald_coral_20261027_233529/supporting_files/diagnostics/) / [ecosystem原诊断](GLMM_MC19_env4_TS_Wald_ecosystem_20261029_123804/supporting_files/diagnostics/)。
- [PM11背景图](PM11_timeseries_2021_2025_surface_top_20261009_180827/)、[可编辑环境源文件](environment_assets/)、[Prism](TS-ENV折线图.prism)。
- [历史PERMANOVA图](../archive/PERMANOVA/figures/) / [历史GLMM图](../archive/GLMM_figures/)。
- [当前结果](../results/README.md)。


## 2026-10-31：字体/布局统一更新

当前MC19主图、两部分曲线/残差、NMDS Shepard与原GLMM诊断已用本地R重画并覆盖原路径；字体放大加粗，热图＋森林图压缩空白。旧版备份和发布清单位于[style_refresh_run001](../archive/style_refresh_run001/README.md)。

| 范围 | 2D stress | 3D stress |
|---|---:|---:|
| Coral | 0.178246 | 0.133440 |
| Ecosystem | 0.220116 | 0.156899 |

同范围NO2/氨氮使用同一套各维度坐标。3D轴名及统计说明固定在图外，原始P和BH q分别标记；本次只改样式，没有重新计算NMDS或统计结果。


3D更新：轴线和刻度已固定到同一组空间棱上，随立方体整体旋转，不再自动换棱。每个范围使用一张整合图，右侧选择NO2-N或Ammonia-N；切换保持当前视角与全部坐标，只更新颜色/色标/原PERMANOVA统计说明。原单指标HTML入口也已覆盖为可切换版本。

## 当前PPT横版布局

热图＋森林图保持16:9：热图加宽、森林面板收窄，MC编号、英文筛选依据和完整置信区间保留。NMDS静态图（含3D轴对投影及Shepard）、N/P折线、two-part森林图已统一为16:9。原GLMM显著曲线每页最多2个MC，同MC同行同页，Coral共3页、ecosystem共4页；主图完整范围、右侧Zoom保留。

已经清楚的横版two-part曲线/残差和GLMM补充诊断保留原布局。全部仅基于保存结果重画，P/q、预测和NMDS坐标不变，交互3D保持现有版本。覆盖前备份、源码及核对记录见[本次归档](../archive/ppt_landscape_run001/README.md)。


## 2026-11-01图形更新

已用本地R从保存模型和缓存优化并覆盖GLMM曲线、two-part森林/曲线、NP和交互3D，并同步MC19_阶段汇报结果图_v01。Two-part森林按det/pos分面，正值曲线显示百分比并在需要时附独立Zoom；原始点和完整CI保留。NO2-N/DIP检出限分别0.002 mg N/L、0.002 mg P/L（原分析L/2=.001），时间折线有水平虚线、环境响应曲线有垂直虚线。3D默认视角放大、固定轴与统计说明及营养盐切换保留，HTML应与旁边lib目录一起传送。既有refresh_MC19_figure_style.R继续读取保存缓存；本次执行和备份见研究根archive/figure_polish_run001/README.md。未重跑分析。
