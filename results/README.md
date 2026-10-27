# 结果目录导航

2026-10-28核对：最新已完成GLMM为[GLMM_MC19_env4_TS_Wald_coral_20261027_233529](GLMM_MC19_env4_TS_Wald_coral_20261027_233529/00_result.csv)。98coral/7TS、76个TS随机截距模型，Wald P统一BH n=76；8项q<0.05。[全精度明细](GLMM_MC19_env4_TS_Wald_coral_20261027_233529/supporting_files/all_results.csv)、[诊断总表](GLMM_MC19_env4_TS_Wald_coral_20261027_233529/supporting_files/diagnostics/diagnostic_overview.csv)。主表重要参数已前置。历史152项两结构结果、先导、bootstrap和LOTO保留本地；不与本批q混用。

以下为历史流程及结果记录，当前运行按上述TS-only入口执行。

2026-10-27当前GLMM代码已切到19MC×4营养盐×两结构152项Wald/BH，尚未运行。运行后新建GLMM_MC19_env4_Wald_coral_时间戳/，顶层152行00_result.csv、简述和supporting_files；复用120255原8拟合，新增144个。153110仍是已完成的bootstrap历史结果，不与新Wald P/q混用。[当前运行及代码解读](../docs/GLMM_MC19_env4_Wald_code_guide_20261027.md)。

2026-10-27：完整数据beta推断源码已修改、尚未运行。GLMM.R默认续接原8模型，全部98样本（含TS4/TS6）、两结构并列，不以LOTO翻转筛除。运行后新建`GLMM_MC19_nitrate_full_data_时间戳/`，主表含beta_P_value/beta_q_value；配套目录保存全部bootstrap及LOTO补充。999次零模型参数bootstrap，BH当前8项；不改旧批结果。[运行说明](../scripts/README.md)。

逐TS留一已由用户完成：[140106八行汇总](GLMM_MC19_nitrate_LOTO_20261027_140106/00_result.csv)，64行估计和56个新模型在supporting_files。56均数值可用，1个TS方差边界标记；关系对TS4/TS6敏感，P/q未校准。只读复用120255完整模型，LOTO latest与pilot分离。[说明](../scripts/README.md)。

2026-10-27：MC19 GLMM硝酸盐先导已由用户完成[120255主表](GLMM_MC19_nitrate_pilot_20261027_120255/00_result.csv)。8模型均收敛/Hessian通过，但均有诊断标记；P/q未校准。每批顶层仅主表、简述和supporting_files；新先导不读取旧MC8指针。[运行说明](../scripts/README.md)。

整理核对：2026-10-27。日常优先查看下列当前结果；历史批次统一从归档入口查找。

## 当前结果

| 内容 | 入口 |
|---|---|
| 当前MC19 GLMM：TS随机截距、Wald/BH n=76 | [GLMM_MC19_env4_TS_Wald_coral_20261027_233529/](GLMM_MC19_env4_TS_Wald_coral_20261027_233529/) |
| 当前MC19 PERMANOVA：coral、7项当期环境、v2 | [PERMANOVA_MC19_env7_coral_samples_20261026_224333/](PERMANOVA_MC19_env7_coral_samples_20261026_224333/) |
| MC19主要结果表 | [environment_results.csv](PERMANOVA_MC19_env7_coral_samples_20261026_224333/environment_results.csv) |
| MC19详细输出、运行说明、输入与脚本快照 | [supporting_files/](PERMANOVA_MC19_env7_coral_samples_20261026_224333/supporting_files/) |
| 最近已完成的历史MC8 GLMM结果 | [GLMM_Z_TS_SD_Q_MC8_20261013_185720/](GLMM_Z_TS_SD_Q_MC8_20261013_185720/) |
| GLMM导航 | [00_LATEST_GLMM.md](00_LATEST_GLMM.md) |
| PERMANOVA导航 | [00_LATEST_PERMANOVA.md](00_LATEST_PERMANOVA.md) |

MC19 v2已核对运行状态为complete、流程状态为analysis_and_plots_complete。主结果表和supporting_files保持本次运行原有结构。

## 汇总工作簿

| 内容 | 文件 |
|---|---|
| TS及历史环境数据 | [TS1-TS7_环境指标_SPSS.xlsx](TS1-TS7_环境指标_SPSS.xlsx) |
| GLMM人工汇总 | [Z_TS_SD_Q_模型结果与主要结论_20261013_185720.xlsx](Z_TS_SD_Q_模型结果与主要结论_20261013_185720.xlsx) |

## 历史结果

- [MC19 v1（2026-10-26 21:51）](../archive/PERMANOVA/results/PERMANOVA_MC19_env7_coral_samples_20261026_215151/)。
- [MC8 PERMANOVA（2026-10-14）](../archive/PERMANOVA/results/PERMANOVA_MC8_20261014_183535/)。
- [全部PERMANOVA历史批次](../archive/PERMANOVA/results/)、[GLMM MC8历史筛查](../archive/GLMM_MC8/)、[早期GLMM](../archive/GLMM_early/)。

## 程序指针与旧路径

`00_latest_*.txt`、`00_current_*.txt`保留原位置和原内容，设为隐藏，由现有脚本继续读取。历史运行的旧目录位置保留隐藏的兼容联接，实际数据位于archive，同一批结果只有一份物理文件。

各当前运行和两个Excel保留原路径。已有模型、结果表、诊断、HTML配套文件和脚本快照整批保留；后续新运行仍可由原脚本写入results。

本次[整理记录](../archive/organization_results_20261026_225450/README.md)包含移动清单、整理前文件校验和整理后验证。

本次图目录整理、脚本与验证更新见[更新日志](../CHANGELOG.md)。
