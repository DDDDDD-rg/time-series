# 当前输入与数据来源

更新：2026-10-28。当前GLMM与PERMANOVA共用MC19工作簿；原始输入内容未修改。

| 文件/目录 | 用途 |
|---|---|
| [PERMANOVA_input_MC19_env7_20261026.xlsx](PERMANOVA_input_MC19_env7_20261026.xlsx) | 当前样本元数据、19MC比例与TS环境；Community/Metadata/Environment_TS/ReadMe四表 |
| [mc/MC_abundance_absolute.tsv](mc/MC_abundance_absolute.tsv) | GLMM原始MC计数；分母为全部Ruegeria MC reads |
| [mc/](mc/) | 原始MC绝对和相对丰度表 |
| [Environment_noDO_Nsplit_20261026.xlsx](Environment_noDO_Nsplit_20261026.xlsx)、[TS_environment_noDO_Nsplit_20261026.csv](TS_environment_noDO_Nsplit_20261026.csv) | 当期/历史环境整理，7项分析指标及DIN参考 |
| [source_tables/](source_tables/) | 原始DNA、样本元数据、指示种/网络等来源 |
| [asv/](asv/) | 原始ASV资料 |
| [water_quality/](water_quality/) | PM11水质CSV和营养盐制图资料 |
| [Cooccurence.xlsx](Cooccurence.xlsx) | 原网络资料工作簿，保留 |

当前有效环境覆盖TS1–TS7；coral逐样本分析共98个样本。GLMM用4种营养盐，PERMANOVA用7项环境；具体筛选、标准化和检测限L/2规则见[运行说明](../scripts/README.md)。

## 历史输入与归类

- [GLMM_input_MC8_20261011.xlsx](GLMM_input_MC8_20261011.xlsx)、[PERMANOVA_input_MC8_20261011.xlsx](PERMANOVA_input_MC8_20261011.xlsx)为历史MC8输入。保留原位置，供历史源码读取。
- [archive/](archive/)存放更早MC11输入及备份。
- `20260928.pptx`属于汇报文件，已移到[presentations/20260928.pptx](../presentations/20260928.pptx)。

旧版输入导航保存在[本次文档备份](../archive/organization_project_20261028_214616/documentation_before/data/00_CURRENT_DATA.md)。
