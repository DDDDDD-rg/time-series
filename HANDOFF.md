# Time series收束交接：可保留资产与新方向接口

日期：2026-11-17。状态：用户已终止此研究方向，整体留在原time series文件夹内归档，不删除数据。不再默认推进TS8、联合GLMM、收缩后验或dbRDA。候选新方向为“珊瑚相关细菌的昼夜厌氧呼吸”，具体对象、机制和实验设计尚未确定。

## 1. 以后先读什么

当前文件是面向后续研究的交接入口。完整冻结材料在[archive/closed_v1](archive/closed_v1)，保留原data、scripts、results、figures、presentations、docs、archive和运行注册表结构。

| 需要找的内容 | 冻结位置 |
|---|---|
| 最后PPT与中文talk要点 | [presentations](archive/closed_v1/presentations) |
| 全部模型、统计表与检查结果 | [results](archive/closed_v1/results) |
| 原始MC/ASV、宿主、环境输入 | [data](archive/closed_v1/data) |
| 可复用分析代码及运行说明 | [scripts](archive/closed_v1/scripts) |
| 作图缓存、六指标交互图 | [figures](archive/closed_v1/figures) |
| 历史模型、布局备份及此前审查 | [archive](archive/closed_v1/archive) |
| 本次文件清单、校验与路径对照 | [archive/closeout_records_v1](archive/closeout_records_v1) |

共享研究知识库继续保留在项目根目录Research_Knowledge，供其他方向使用。相关入口：[模型审查](../../Research_Knowledge/questions/MC19_TS7_PPT_data_figure_audit_v1.md)、[联合分析记录](../../Research_Knowledge/questions/MC19_joint_TS7_results_review_v1.md)、[原研究状态](../../Research_Knowledge/questions/Time_series_current.md)。

## 2. 最值得留下的判断

1. 旧项目留下了MC层级的环境关联线索和可追溯数据，适合用于挑选比较对象、查找样本及审查功能假说；现有结果没有测定昼夜厌氧呼吸。
2. MC6有最直接的外部反硝化研究背景；MC15与MC0/MC10代表株的功能差异可用于构建比较框架。MC98/MC80在旧现场模型中与MC15方向相反，值得保留，但其具体呼吸功能仍需独立鉴定。
3. 现有分析的独立环境条件只有七套TS值。增加拟合复杂度没有增加环境重复；新设计应优先取得实际的昼夜、氧微环境和过程速率证据。
4. 新方向可从Ruegeria切入，也可以扩展到其他珊瑚相关细菌；旧19MC候选名单不应自动限定新方向。

## 3. 实际有哪些数据，还缺什么

### 可直接重用的原始资产

| 资产 | 价值与注意事项 |
|---|---|
| [MC原始计数](archive/closed_v1/data/mc/MC_abundance_absolute.tsv)及[原比例表](archive/closed_v1/data/mc/MC_abundance_relative.tsv) | 可回溯MC检出与组成；当前模型分母为全部Ruegeria MC reads，含Unassigned；19个候选不另归一至100% |
| [Ruegeria ASV表](archive/closed_v1/data/asv) | 可回到标记层级重新审查群体分辨率；表格本身不等于参考基因组或功能基因数据 |
| [TS1–TS8宿主元数据](archive/closed_v1/data/source_tables/metadata_TS1_TS8.xlsx) | 本轮逐表确认字段为New_ID、sample_ID、type、host niche；可追踪样本/宿主；这些表没有采样时刻和昼夜标签 |
| [共用工作簿](archive/closed_v1/data/MC19_TS8_env5.xlsx) | 已保留TS8生物信息；文件名含TS8，不代表八TS环境分析完成 |
| [PM11原始水质](archive/closed_v1/data/water_quality/Bluff%20Island_PM11_marine_water_quality.csv) | 含营养盐、温度、DO等原始字段；有站点、日期、水层。DO虽然存在，不能替代珊瑚表面/胃腔/组织内的氧剖面 |
| [DNA、网络和ASV来源资料](archive/closed_v1/data/source_tables) | 包含data_DNA_time-series_2026Aug19.xlsx、ISA_asv.xlsx、Network.xlsx，保留筛选与原始追溯价值；本轮未重新定义这些文件的全部实验来源 |

当前正式范围：TS1–TS7；coral 98、ecosystem 310个样本。Ecosystem包含coral，排除海水/沉积物；二者不是独立重复。TS8生物样本已进入准备输入，完整环境数据仍未补齐。

环境采用邻近Bluff Island的PM11表层水，按生物采样月份匹配；不是现场同一时刻测得的珊瑚微环境。七TS温度为25.8、31.4、28.6、16.7、24.6、26.6、25.2°C。

TAN是总氨氮，字段Ammonia_mgN_L；非离子氨U-NH3是另一测量项。DIN=NO3-N+NO2-N+TAN。营养盐低于报告限记录保留原“<L”标记，数值用L/2，温度不采用L/2。

### 不能从当前表直接获得的证据

- 已进入模型的元数据没有采样时刻、光照周期或昼夜配对；不能把季节TS直接重新标成昼夜。
- 未见本项目已完成的珊瑚微环境连续氧测量、细菌呼吸速率、15N过程示踪或昼夜细菌功能表达分析。
- MC/ASV丰度不能直接给出反硝化活性。若另有未纳入本目录的样本日志、菌株库、基因组、转录组或传感器资料，需单独核对后接入。

## 4. 终止时到底做到了哪一步

本轮重新核对保存的单环境与联合系数表：以下为已运行输出；诊断/敏感性/先验完成状态结合此前实际文件审查记录保留，没有重跑模型。

| 分析 | 冻结状态与主要信息 |
|---|---|
| 单环境Beta-binomial GLMM | 每范围19MC×6环境=114项；coral 10项、ecosystem 18项通过BH，涉及6/10个MC |
| PERMANOVA | 19MC原比例平方根变换后Bray–Curtis距离；按整个TS分配环境标签；每范围六项BH，均无BH显著结果 |
| two-part | 只分解父GLMM显著组合；coral 20、ecosystem 36个组件；同数据筛选后的探索性分解 |
| 五指标联合GLMM | NO3、NO2、TAN、DIP、温度同时入模，不加DIN；每范围19模型、95条件系数；本轮CSV确认coral 1项、ecosystem 14项通过BH |
| 联合残差诊断 | 38模型完成，每模型1000次unconditional模拟；均存在某类review标记，需要按零值、尾部等具体问题解释 |
| 逐TS及极端样本敏感性 | 42个逐TS任务＋18个样本影响任务，60/60完成 |
| 先验预测 | 两范围×19MC×tau 0.5/1/2，114/114项完成；tau2偏宽，后续先验调整未执行 |
| NMDS稳定性 | Coral二维稳定，stress≈0.178246；ecosystem二维局部点不稳定，stress≈0.220116；ecosystem三维稳定，stress≈0.156899，五次fresh均重复找到布局 |
| 未完成 | 完整八TS、dbRDA、收缩后验与后验先验敏感性、控制实验/实际昼夜过程验证 |

这些未完成项随旧方向一同关闭，不是新方向的默认待办清单。

### 最后版本的文件入口

- 单环境：[coral](archive/closed_v1/results/MC19_GLMM_TS7_env6_coral_run001)、[ecosystem](archive/closed_v1/results/MC19_GLMM_TS7_env6_ecosystem_run001)。不要用旧四营养盐或温度增量批次代替。
- 联合：[coral](archive/closed_v1/results/joint_glmm_TS7_env5_coral_run001)、[ecosystem](archive/closed_v1/results/joint_glmm_TS7_env5_ecosystem_run001)；文件00_coefficients.csv和01_model_checks.csv。
- 检查：[残差](archive/closed_v1/results/joint_glmm_TS7_env5_diagnostics_run001)、[敏感性](archive/closed_v1/results/joint_TS7_sensitivity_run001)、[先验预测](archive/closed_v1/results/joint_TS7_prior_predictive_run001)。
- 排序：[二维稳定性](archive/closed_v1/results/nmds_2d_stability_v1)、[ecosystem三维稳定性](archive/closed_v1/results/nmds_3d_ecosystem_stability_v1)、[三维交互图](archive/closed_v1/figures/TS7_env6_PPT_v1/ecosystem/auxiliary/NMDS3D_six_indicators.html)。

## 5. 重点MC：哪部分可转移，哪部分仍是线索

| 对象 | 旧项目实际线索 | 对新方向的潜在用途 |
|---|---|---|
| MC6 | 五指标联合中温度两范围负向；ecosystem方向/幅度较稳定，coral逐TS敏感 | 最值得先核查菌株与参考基因组对应；文献MC6-AB11有直接反硝化背景；当前温度关联不是昼夜响应证据 |
| MC15 | 单环境中NO2/TAN负向，两范围已检出后比例均下降；ecosystem检出也下降；未检出集中TS6–7 | 适合与具有反硝化能力的代表株比较；联合MC15—NO2留出TS5后变号，不能当作稳定的独立氮效应 |
| MC98、MC80 | 单环境中与MC15形成相反方向；部分two-part结果涉及检出和正值比例 | 是现场对照候选；不能仅据正相关把它们称作反硝化菌。MC98—TAN联合方向一致但精度敏感 |
| MC0-A5、MC10-B4 | 功能线索来自外部代表株文献；没有在旧现场模型中建立其昼夜作用 | 可用来讨论菌株功能比较、适应策略与MC15-BG7的对照，须确认菌株是否实际可获得 |
| MC25、MC2s1、MC95、MC43 | 温度关联与残差/零值/尾部问题需合并解释；MC95/43极端样本未主导温度方向 | 保留为次级探索线索；没有现成的昼夜厌氧呼吸机制证据 |

MC编号跨项目、分群版本及标记库未必对应。进入新实验前，核对参考菌株/基因组与标记树；不要把某代表株的功能直接扩展到整个MC。旧筛选理由含宿主关联、网络位置和文献功能，不等于全体19MC都经过过程实验验证。

## 6. 需要继承的方法经验

- 区分环境条件重复与同条件下的生物样本数；新研究若有同一珊瑚的昼夜多次取样，需明确珊瑚个体、实验缸和时间的层级。
- 单环境关联、联合条件关联、预测重要性是不同问题；不按P值给环境因素排“最重要”顺序。
- 当前NO3与DIP、NO2与TAN强共变；五指标营养盐VIF约8.55–99.91，温度约1.619。加入多个变量并未解决七套环境条件的信息限制。
- 数值收敛、残差匹配与逐TS稳定性分别检查。旧显著单环境组合28项中27项有某类诊断标记；review不等于全部作废，但需保留具体限制。
- 组成丰度、检出、功能基因、表达和过程速率各是不同证据层。新方向最好从一开始就定义需要的是哪层变化。
- MC丰度用全部Ruegeria reads作分母；如果新研究要比较细菌总体绝对活性，需考虑菌量和归一化，旧相对比例口径不能自动沿用。

## 7. 与候选新方向最直接相关的文献

本节是承接线索，不是系统文献综述。本轮只针对新方向入口做少量原文/摘要核查。

| 文献与读取范围 | 支持什么 | 使用边界 |
|---|---|---|
| [Xiang et al., 2026, Nature Communications](https://doi.org/10.1038/s41467-026-72571-w)；既有定向精读笔记 | Ruegeria MC1/6/13/18 specialist背景及低氧15NO3功能比较；MC0亦有non-specialist过程证据 | 不能等同于本项目已测昼夜；先核对MC/菌株对应。见[笔记](../../Research_Knowledge/papers/Xiang_2026_RuegeriaDenitrifiers.md) |
| [Xie et al., 2026, The ISME Journal](https://doi.org/10.1093/ismejo/wrag106)；既有定向笔记 | MC10-B4、MC0-A5与MC15-BG7的功能和宿主相关差异，提供代表株比较背景 | 菌株结果与现场MC相对丰度分别保留。见[笔记](../../Research_Knowledge/papers/Xie_2026_MC10Traits.md) |
| [Yang et al., 2024, Chemosphere](https://pubmed.ncbi.nlm.nih.gov/38777199/)；本轮核对摘要/题录 | Pocillopora damicornis的15N培养实验报告反硝化及anammox产生N2，并报告昼夜节律 | 对新方向最直接的起点；实验条件、具体时段及重复设计仍需阅读全文/补充材料；不能替代本地过程测定 |
| [Linsmayer et al., 2020, Scientific Reports](https://pmc.ncbi.nlm.nih.gov/articles/PMC7669893/)；本轮核对摘要及氧测量结果段 | 珊瑚扩散边界层氧可随昼夜显著变化，支持直接测微环境而非只用海水DO | 论文也讨论珊瑚发酵；宿主发酵不等同于细菌厌氧呼吸。暂不将其视为细菌功能的直接证明 |

Yang与Linsmayer为新增入口核查，未升级为完整精读笔记、未下载批量资料、未复现。既有[氮与共生生理主题](../../Research_Knowledge/topics/nitrogen_probiotics_physiology.md)及[宿主/环境层级](../../Research_Knowledge/topics/MC_time_host_environment.md)可继续查阅。

## 8. 新方向开始前需要确定的四件事

1. 研究对象与空间位置：活体珊瑚整体、表面黏液/组织、胃腔，还是分离菌株？优先Ruegeria还是更广的细菌群体？
2. 要测的过程：以硝酸盐/亚硝酸盐为电子受体的呼吸，完整反硝化、DNRA或其他过程？发酵与anammox另作过程分类，不笼统当作同一种“厌氧呼吸”。
3. 要解释的昼夜变化：菌量/组成、功能表达、单位生物量活性，还是实际转化速率？需要对应测量，不能只用一个丰度结果代替。
4. 可获得的材料：现成菌株、参考基因组、留存珊瑚样本、采样日志、光照与氧微电极/光学传感器、过程产物/同位素测量条件。先盘点，再选可执行的研究问题。

可讨论的起始问题：光暗转换引起的珊瑚局部氧变化，是否伴随特定细菌氮呼吸活性的变化？此句是候选问题，不是已经确定的研究假设或授权实验。

## 9. 归档与将来恢复

- 所有材料留在当前time series目录内部，无模型、原始表、图、PPT或未知用途材料被删除。
- closed_v1为本次冻结内容。外层README/HANDOFF是阅读入口；历史绝对路径字段与源码参数保留原状，用于追溯。
- 现有历史目录联接随整体移动修正目标，仍指向同一批归档材料；不是额外复制一套数据。路径清单与修正记录随档保存。
- 本轮不保证冻结脚本可直接从新路径续跑：它们可能含原绝对路径、旧RDS路径和源文件定位。以后若恢复，应先复制归档到独立工作目录，按清单处理路径并进行只读输入检查，避免在唯一归档上直接拟合或覆盖。
- 未提交或推送Git。共享方法笔记继续保留；旧time series状态明确为终止，未来不要根据历史“下一步”自动恢复计算。

本轮核查：重新读取当前PPT/讲述稿位置、原始目录及元数据字段、单环境与联合CSV；结合此前已完成的诊断、敏感性和NMDS实际检查记录整理。没有重新拟合、模拟、置换或开展新方向实验。
