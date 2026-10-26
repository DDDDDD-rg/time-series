# Time-series research project

2026-10-26当前PERMANOVA：coral原始样本、19个MC、7项当期环境；入口为[scripts/run_permanova.R](scripts/run_permanova.R)，输入为[data/PERMANOVA_input_MC19_env7_20261026.xlsx](data/PERMANOVA_input_MC19_env7_20261026.xlsx)。v2已于2026-10-26 22:43完成分析与作图；统计量由vegan计算，保留TS标签置换与P值汇总规则。最新运行和整理记录见[更新日志](CHANGELOG.md)。

当前results运行目录顶层只放environment_results.csv和supporting_files/；figures对应目录顶层为3份合并PDF和supporting_files/，单图与配套文件集中保存。配套作图兼容旧目录。运行及字段说明见[scripts/README.md](scripts/README.md)。

以下保留2026-10-14的MC8与GLMM入口记录；其中MC8 PERMANOVA结果和图为历史运行。

## 常用入口

| 用途 | 文件或目录 |
|---|---|
| 当前MC19 PERMANOVA | [主结果表](results/PERMANOVA_MC19_env7_coral_samples_20261026_224333/environment_results.csv)、[三份合并图入口](figures/00_LATEST_PERMANOVA.md) |
| 当前GLMM结果 | [结果总览](results/GLMM_Z_TS_SD_Q_MC8_20261013_185720/00_key_results.html)、[中文摘要](results/GLMM_Z_TS_SD_Q_MC8_20261013_185720/00_key_results.txt) |
| GLMM数值检查 | [各MC检查表](results/GLMM_Z_TS_SD_Q_MC8_20261013_185720/00_model_checks.csv) |
| GLMM效应结果 | [DO系数](results/GLMM_Z_TS_SD_Q_MC8_20261013_185720/00_DO_significance.csv)、[全部系数](results/GLMM_Z_TS_SD_Q_MC8_20261013_185720/00_beta_significance.csv) |
| 历史MC8 PERMANOVA | [结果总览](results/PERMANOVA_MC8_20261014_183535/00_PERMANOVA_summary.html)、[关键结果](results/PERMANOVA_MC8_20261014_183535/00_key_results.txt) |
| 历史MC8两两比较 | [比较族汇总](results/PERMANOVA_MC8_20261014_183535/00_pairwise_summary.csv)、[完整两两结果](results/PERMANOVA_MC8_20261014_183535/01_all_pairwise_results.csv) |
| 历史MC8 PERMANOVA图 | [完整PDF](figures/PERMANOVA_MC8_20261014_183535/00_all_final_results.pdf)、[图片目录](figures/PERMANOVA_MC8_20261014_183535/) |
| 其他当前图和源文件 | [图目录说明](figures/README.md) |
| 本次整理记录 | [移动清单](archive/organization_20261014/move_plan.csv)、[整理说明](archive/organization_20261014/README.md) |

## 目录用途

| 目录 | 内容 |
|---|---|
| [data/](data/) | 当前输入、原始MC/ASV表、样本元数据和水质资料 |
| [scripts/](scripts/) | 当前R脚本、辅助函数和RStudio项目；[运行入口说明](scripts/README.md) |
| [results/](results/) | 最新GLMM/PERMANOVA结果、latest/current指针和当前结果工作簿 |
| [figures/](figures/) | 当前图、环境作图源文件和GLMM人工汇总图 |
| [docs/](docs/) | PERMANOVA、批次置换和两两比较的方法说明 |
| [presentations/](presentations/) | 此项目目录原有的PPT |
| [archive/](archive/) | 历史运行、旧图、脚本备份和整理记录；[归档目录说明](archive/README.md) |

## 输入与运行

- 历史MC8 PERMANOVA输入：[PERMANOVA_input_MC8_20261011.xlsx](data/PERMANOVA_input_MC8_20261011.xlsx)。
- GLMM输入：[GLMM_input_MC8_20261011.xlsx](data/GLMM_input_MC8_20261011.xlsx)。
- 完整PERMANOVA入口：[scripts/run_permanova.R](scripts/run_permanova.R)。
- 当前GLMM完整拟合入口：[scripts/GLMM.R](scripts/GLMM.R)。
- 已有GLMM结果复核入口：[scripts/GLMM_reads_diagnostics_full.R](scripts/GLMM_reads_diagnostics_full.R)。

当前输入、当前分析结果、当前脚本和RStudio项目保留原路径。历史结果和旧图已归档，旧目录位置留有隐藏的兼容联接，供既有脚本、缓存和旧链接使用。联接与归档目录指向同一份文件，不产生第二份模型或结果。

HTML若带有同名 `_files` 目录，保存和分享时应整组保留。MC、候选结构、诊断和脚本快照继续保存在各运行目录内部。
