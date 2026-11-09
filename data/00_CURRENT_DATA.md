# 当前输入与数据来源

更新：2026-11-02。当前GLMM、PERMANOVA及后续图统一使用8TS、19MC、五项营养盐。

| 文件/目录 | 用途 |
|---|---|
| [MC19_TS8_env5.xlsx](MC19_TS8_env5.xlsx) | 当前共用输入；Community、Metadata、Environment_TS、ReadMe四表 |
| [mc/MC_abundance_absolute.tsv](mc/MC_abundance_absolute.tsv) | 原始MC计数；分母加总全部Ruegeria MC含Unassigned |
| [water_quality/](water_quality/) | EPD PM11水质原始资料 |
| [source_tables/](source_tables/) | 原始样本元数据、DNA与网络等资料 |
| [PERMANOVA_input_MC19_env7_20261026.xlsx](PERMANOVA_input_MC19_env7_20261026.xlsx) | 原输入保留为历史记录，当前程序不再读取 |
| [Environment_noDO_Nsplit_20261026.xlsx](Environment_noDO_Nsplit_20261026.xlsx)、[TS_environment_noDO_Nsplit_20261026.csv](TS_environment_noDO_Nsplit_20261026.csv) | 原环境整理及历史参考，保持不变 |

TS8 MC和宿主元数据已保留：coral共116个样本（TS8为18），ecosystem共356个（TS8为46）。TS8环境对应2026-02，待填。

填写Environment_TS第9行：C9为实际水样日期，E9:H9为NO3-N、NO2-N、Ammonia-N、DIP原始报告值。DIN自动由三项无机氮相加；低于报告限的成分先取L/2。未齐全时8TS程序停止，不排除TS8后运行。

PM11表层按生物样本采样月匹配，原始报告值与删失标记保留。未使用DO、Chl-a；盐度、水温、pH历史字段仅追溯，不进入五指标模型。

现有7TS结果表和图保留，本次未重新运行分析。完整入口见[运行说明](../scripts/README.md)。历史MC8、Cooccurence.xlsx和其他原数据保持原位。
