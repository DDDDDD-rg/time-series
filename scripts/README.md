# MC19时间序列：运行说明

2026-11-16联合GLMM更新：当前[七TS普通联合入口](GLMM/GLMM_MC19_joint_TS7.R)同时考虑NO3、NO2、总氨氮、DIP和温度，不加DIN；每范围19模型、95个条件系数，默认check。温度逐TS核对PM11表层原CSV；只读共线性与隔离检查通过，正式fit由用户运行。使用独立`joint_glmm_TS7_env5_{scope}_runNNN`结果目录，旧四营养盐联合结果保留。详见[运行说明](GLMM/GLMM_MC19_joint_TS7_README.md)。下方四营养盐联合入口说明为历史准备记录；六指标单环境流程仍独立保留DIN。

2026-11-16最新：按用户新要求，七TS入口恢复DIN并保留DIP，共NO3-N、NO2-N、总氨氮、DIN、DIP、温度六项；使用独立TS7_env6版本，全量重新拟合每范围114项GLMM、重新计算六项PERMANOVA，two-part仅跟进新GLMM的BH显著组合。结果/作图仍分开，默认只读check；本轮只改脚本，未运行正式分析。见[当前运行说明](TS7_temperature_README.md)。以下2026-11-15五指标增量记录为历史方案，已被本项替代。

2026-11-15最新：七TS温度入口现为“四项营养盐＋温度”的增量版，排除DIN。核对并复用旧76个营养盐GLMM，只补缺少缓存的温度模型；95项统一BH。PERMANOVA复用五项原始结果、重算五项BH。结果与作图分别使用results/plots，见[最新运行说明](TS7_temperature_README.md)。下方六指标/114项为此前准备记录，已由本决定替代。

2026-11-15新增：用户指定的当前七TS全流程补入温度，使用[TS7_temperature.R](TS7_temperature.R)，详见[七TS温度全流程运行说明](TS7_temperature_README.md)。覆盖单指标GLMM、PERMANOVA、two-part、全套静态图与交互NMDS；六指标为五项营养盐加水温。默认只读check，正式计算由用户运行；联合GLMM脚本和结果不改动。下方八TS主流程保留其完整输入要求。

2026-11-15新增：用户明确要求准备当前七TS的普通联合GLMM，独立入口为[GLMM_MC19_joint_TS7.R](GLMM/GLMM_MC19_joint_TS7.R)，详见[运行与代码解读](GLMM/GLMM_MC19_joint_TS7_README.md)。该入口同时考虑四项营养盐，默认check模式，不执行拟合；正式fit时为每范围19模型、76个条件系数。它不改动下述八TS主流程和完整性要求，已有模型及结果不覆盖。

更新：2026-11-02。正式流程统一为TS1–TS8、19MC和五项营养盐，直接替代旧7TS输入。现有7TS结果表和图保留。当前TS8环境数据待填，本次仅准备输入与源码，没有重新运行分析。

## 填写输入

共用工作簿：[MC19_TS8_env5.xlsx](../data/MC19_TS8_env5.xlsx)。Community与Metadata已含TS8样本。Environment_TS每TS一行，TS8在第9行。

只填写Environment_TS中的实际水样日期（C9）和四项原始报告值（E9:H9）。列顺序为NO3-N、NO2-N、Ammonia-N、DIP。普通值如0.008，低于报告限如<0.002。日期需属于生物样本的采样月份；TS8目前为2026-02。

数值列、DIN、删失标记和Metadata环境列由公式计算。R读取时也直接从Environment_TS原始报告值计算并按TS匹配，不依赖Excel公式缓存。缺失不是0，未填齐时停止，不回退为7TS。保留的盐度、水温、pH历史字段不进入模型。

DIN = NO3-N + NO2-N + Ammonia-N，单位mg N/L；成分低于报告限先取L/2。DIN的删失标记表示至少一个组成成分使用L/2，没有另设DIN直接测量报告限。

## 分析设置

| 项目 | 当前设置 |
|---|---|
| 范围 | coral和ecosystem分别分析；ecosystem包含coral并排除海水/沉积物 |
| 样本 | 所有环境完整时，coral 116、ecosystem 356；不按宿主平均 |
| 环境 | NO3-N、NO2-N、Ammonia-N、DIN、DIP；每次一个环境指标 |
| 原GLMM | 19×5=95个Beta-binomial/logit模型；`cbind(MC_reads, Other_reads) ~ E_z + (1 | TS)` |
| 标准化与BH | 8TS等权中心/SD；每范围95项BH，失败项不缩小检验族 |
| PERMANOVA | sqrt原始MC比例、Bray–Curtis；8! = 40320种整TS环境标签排列；每范围5项BH |
| two-part | 新8TS原GLMM中q<0.05的组合；检出Bernoulli/logit、正值Beta/logit；合并2×组合数校正BH |
| MC分母 | 全部原Ruegeria MC reads含Unassigned；19MC不重新归一化 |

## 数据齐全后的运行顺序

在RStudio打开scripts.Rproj。每个范围依次运行下列入口，再切换另一个范围。不要运行全部历史脚本。

```r
options(mc19.analysis_scope="coral") # 另一个范围用"ecosystem"
```

1. `GLMM/GLMM.R`：原GLMM、诊断和配套图。完整8TS重新拟合，不复用旧先导。
2. `PERMANOVA/run_permanova.R`：整体检验和配套图。
3. `GLMM/plot_MC19_talk_figures.R`：热图+森林图、五指标时间线、五指标NMDS 2D/3D。
4. `GLMM/GLMM_MC19_two_part.R`：默认check核对任务；设置`options(mc19.two_part.stage="fit")`后运行两部分拟合和图。

只重画用`GLMM/plot_GLMM_MC19.R`、`PERMANOVA/04_plot_permanova.R`，或two-part的plot阶段。来源必须是已完成的8TS新批次。

新结果使用递增run/plot编号，见[命名规则](README_run_names.md)。最新结果定位会检查8TS及五指标版本，旧7TS结果继续留在原目录供查阅。

## 输入与文件读取检查

`test_mc19_input.R`只核对配置、DIN和缺失/日期规则；`test_mc19_run_io.R`使用临时文件验证编号、范围和版本。两者均不拟合模型或运行置换检验。预备阶段另用独立合成TS8环境的工作簿验证两范围完整数据读取；该副本仅在archive，不是研究数据。

历史先导、LOTO、bootstrap和MC8入口保留为历史源码，不属于此次8TS运行顺序。主流程未新增Host、Repeat、水温或覆盖协变量。
