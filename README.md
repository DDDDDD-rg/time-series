# Time-series research project

更新：2026-10-31。Coral（98样本）和ecosystem（310样本，排除seawater/sediment）均有7个有效TS，已完成PERMANOVA、76项GLMM与后续两部分分析，共6批当前结果。

| 内容 | 入口 |
|---|---|
| 全部当前结果与数值明细 | [results/README.md](results/README.md) |
| 当前图：横版热图＋森林图、NMDS、显著关系和两部分结果 | [figures/README.md](figures/README.md) |
| GLMM、PERMANOVA与两部分分析入口 | [scripts/README.md](scripts/README.md) |
| talk图与两部分的check/fit/plot说明 | [README_MC19_next_steps.md](scripts/GLMM/README_MC19_next_steps.md) |
| run/plot编号及最新完成批次发现 | [README_run_names.md](scripts/README_run_names.md) |
| 输入与来源 | [data/00_CURRENT_DATA.md](data/00_CURRENT_DATA.md) |
| 方法与交接 | [docs/README.md](docs/README.md) |
| 本次及历史更新 | [CHANGELOG.md](CHANGELOG.md) |

GLMM按范围分别使用`cbind(MC_reads, Other_reads) ~ E_z + (1 | TS)`，Beta-binomial/logit，每次一个营养盐，76项统一BH；coral有8个、ecosystem有14个q<0.05组合。两部分分析在这8/14组合中拆分检出和正值占比，分别固定16/28个分模型的探索性BH族；诊断和科学解释继续保留审查状态。

PERMANOVA使用19MC原始比例平方根和Bray–Curtis，按7TS枚举5040种环境标签排列，7项统一BH；两个范围均无q<0.05。保留原始样本，不平均，不含Host或Repeat项。

当前talk主图为大号粗体、PPT横版；热图保留19MC×4指标，森林图仅展示存在显著结果的指标，完整参数仍在CSV。3D NMDS整合NO2/氨氮切换，固定棱轴与统计说明，沿用保存坐标和stress。

## 保存与复现

GitHub保留当前源码、文献/问题笔记、完整数值明细、诊断摘要和当前图（HTML连同本地依赖）。原始Excel/计数输入、模型RDS、NMDS缓存、诊断模拟、运行锁与历史大输出继续留本地。重新拟合或按保存模型重画需准备这些输入/缓存并调整入口project_dir。

scripts按[GLMM](scripts/GLMM/README.md)和[PERMANOVA](scripts/PERMANOVA/README.md)分组；独立PM11作图和RStudio项目留外层。新结果采用run编号，重画采用plot修订，方法版本与来源单独记录；读取器按分析类型/范围选择最新完成批次并回退失效指针。

[archive](archive/README.md)保存旧模型/图、测试、修改前备份及整理清单。旧结果目录的兼容联接服务于本地历史路径；GitHub通过提交历史保留先前发布图。原始输入和已保存统计结果不因展示样式或目录调整而重算。
