方法记录。当前MC8结果：[总览](../results/PERMANOVA_MC8_20261014_183535/00_PERMANOVA_summary.html)。旧MC11结果见[历史归档](../archive/PERMANOVA/results/PERMANOVA_MC11_history/)。
当前分析范围与流程见[新版PERMANOVA说明](README_permanova_readable.md)。下文保留已有方法记录。

# 宿主和采样批次两两多重比较

2026-10-06 新增 scripts/PERMANOVA/05_test_pairwise_host_time.R。它自动读取01，在选定范围内独立进行两两比较；不要求先运行02。

2026-10-11已对新版8个MC完成全宿主与仅珊瑚两套两两检验，结果随本轮整体分析保存。

04_plot_permanova.R也已增加宿主和TS图例的字母显示，默认范围为all_hosts，对应用户提供的两张图。

## 怎么运行

1. 在 scripts/PERMANOVA/01_data_settings.R 选择 analysis_scope：all_hosts 或 coral，并设置置换次数。
2. 打开 scripts/PERMANOVA/05_test_pairwise_host_time.R，从头运行。
3. 默认 pairwise_targets = c("Host_type", "TS")，两类都做；需要单类时改成 "Host_type" 或 "TS"。
4. 在 <本次MC8结果目录>/<scope>/pairwise 查看结果。旧整体检验结果在原目录。
5. 04选择与01/05相同范围，再从头运行；它读取05结果自动给图例添加字母。04默认show_pairwise_letters=TRUE；FALSE保留无字母图例。

与其他脚本一样，project_dir 当前指向 Desktop 下的 Msc main project-time series 主目录。

## 每一对怎么比较

| 因素 | 取用的数据 | 模型与读取项 | 置换限制 |
|---|---|---|---|
| 两种宿主 | 只用两种宿主均出现的TS，保留其中全部Repeat_group | MC组成 ~ Host_type + TS，读Host_type | 同一TS内置换组 |
| 两个采样批次 | 只用两个TS均出现的宿主，保留其中全部Repeat_group | MC组成 ~ Host_type + TS，读TS | 同一Host_type内置换组 |

每组仍是01算好的Repeat_group均值，不做Host_type×TS第二次平均。8个MC相对丰度开平方根后计算Bray–Curtis；raw为主要距离版本，Lingoes作敏感性检查。

沿用现有02对批次比较的独立重复组假设。脚本没有建立跨TS个体配对。

为了保留当前跨批次／跨宿主调整设计，每对至少需要2个共同TS或共同宿主、两组各至少2个Repeat_group、满秩模型与正的残差自由度。不足的组合记录Status和NA；coverage文件记录共同范围以及被排除的缺采组合。

## 多重校正的范围

- all_hosts：12类宿主产生66对，宿主比较的BH按66项；8个TS产生28对，批次比较的BH按28项。
- coral：3类宿主产生3对，宿主比较的BH按3项；8个TS仍产生28对。
- 两类分别校正。raw和Lingoes分别校正，同一类比较中的跳过项仍计入计划比较数。
- 独立coral分析的3项BH与从all_hosts结果中抽出3项不同：从all_hosts表筛珊瑚时，应保留那张表原有的66项BH值。
- PERMDISP另作一组诊断并单独BH，保留未校正p供核查。

BH用于控制这组比较的错误发现率。[R的p.adjust说明](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/p.adjust.html)

## 结果怎么读

| 列 | 含义 |
|---|---|
| Group1 / Group2 | 本次比较的两个宿主或两个TS |
| N1 / N2 | 共同覆盖范围内，两组各有多少个Repeat_group |
| TS_used / Hosts_used | 本对实际使用的TS和宿主 |
| N_common_blocks / Common_blocks | 共同TS或共同宿主的数量及名单 |
| Residual_df | 残差自由度 |
| N_label_assignments | 在分层设计中两组标签可能的分配数 |
| N_perm | 实际使用的行置换数 |
| R2_raw | 在本对子集中，这一比较项占总距离变异的比例 |
| P_raw | 未经多重校正的PERMANOVA p |
| P_BH_raw | 主要使用的BH校正p |
| P_BH_lingoes | Lingoes敏感性版本的BH校正p |
| BH_family_size | 本次BH使用的计划比较数 |
| Significance | 根据P_BH_raw生成的星号或ns |
| P_disp_raw / P_disp_lingoes | 该对的未校正PERMDISP p |
| P_disp_BH_raw / P_disp_BH_lingoes | 离散度诊断单独BH后的p |
| Status / Warnings | 跳过、错误和计算警告 |
| Dispersion_note_BH | 结合离散度诊断的阅读提示 |

整体检验回答多个组是否存在总体差异；这里逐对回答“哪些组之间有差异”。判断某一对使用P_BH_raw ≤ 0.05，并结合R²、样本数和离散度结果阅读。

adonis2的by="margin"检验在其他模型项调整后的边际效应；分层置换沿用当前设计。[vegan的adonis2说明](https://vegandevs.github.io/vegan/reference/adonis.html)
PERMDISP检验组内离散度；显著时要结合组内变异解释组成差异。[vegan的betadisper说明](https://vegandevs.github.io/vegan/reference/betadisper.html)

## 输出文件

两类分别以 pairwise_host 和 pairwise_TS 开头：

- *_results.csv：每对一行，包括未完成的组合。
- *_coverage.csv：每对在另一因素上的覆盖数量，Used_in_pair表示是否进入本对。
- *_P_BH_raw_matrix.csv：对称BH校正p矩阵，适合查看成对关系。
- *_P_BH_lingoes_matrix.csv：敏感性版本p矩阵。
- *_R2_raw_matrix.csv：对称R²矩阵。
- *_letters.csv：各组对应的a、ab、b等字母、校正范围与标注状态。
- pairwise_host_time_analysis.rds：数据、结果及设置。
- pairwise_host_time_session.txt：运行环境。

矩阵对角线为NA；未完成的对也为NA，原因见results的Status。

干湿季只有Dry和Wet，现有季节检验已是两组比较。DO等连续环境指标检验的是组成与数值变化的关联，继续使用03的8项环境BH结果。

## NMDS上的字母怎么读

字母放在图例中的组名后，例如TS1 [a]、TS2 [ab]、TS3 [b]。这些是说明格式的例子，不是当前数据的结果。

- 两组共享至少一个字母：BH校正后未检出显著差异。
- 两组没有共同字母：BH校正p≤0.05。
- a、b、c不表示丰度或效应大小排序；一个组可能需要多个字母。

以a、ab、b为例：第一组和第二组共享a，第二组和第三组共享b，但第一组和第三组没有共同字母。图例比在椭圆中心堆放字母更容易辨认。

字母由pairwise_letters.R的insert-and-absorb逻辑产生，输入是P_BH_raw；它不拟合新模型。生成后逐对验证共享字母关系与显著性矩阵一致。

若部分比较因为共同TS／宿主不足而缺失，脚本只给最大的完整可比子集生成字母，并列按原图例顺序选择，其余组显示n/a。BH仍沿用全部计划比较数，不对子集重新校正。具体缺失组合见结果Status和coverage；n/a不参与字母比较。例如Sea anemone与部分宿主只有一个共同TS，会触发现行的最低共同批次规则。

新图保留整体检验注释，并增加字母含义；图目录同时保存pairwise_legend_letters.csv。02和05的Repeat_group、宿主/TS和MC丰度必须一致，04会先核对，以免将不同数据版本的字母标到一起。

CLD的共同字母表示未检出差异；本脚本把这一含义写在图注中。[CLD方法说明](https://rvlenth.github.io/emmeans/reference/CLD.emmGrid.html)
