方法记录。当前MC8结果：[总览](../results/PERMANOVA_MC8_20261014_183535/00_PERMANOVA_summary.html)。旧MC11结果见[历史归档](../archive/PERMANOVA/results/PERMANOVA_MC11_history/)。
当前分析范围与流程见[新版PERMANOVA说明](README_permanova_readable.md)。下文保留已有方法记录。

# 季节和环境分析：保留Repeat_group，不做第二次平均

01默认设置为 `analysis_scope <- "all_hosts"`。改为 `"coral"` 即分析三种珊瑚；两种范围都使用同一套处理规则。

## 运行顺序

1. 在 `scripts/01_data_settings.R` 设置分析范围。01由02、03自动读取。
2. 运行 `scripts/02_test_host_time.R` 检验宿主、采样批次。02原本就按Repeat_group分析，算法未改。
3. 运行 `scripts/03_test_season_environment.R`，或直接运行 `scripts/03_test_season_environment_batch_labels.R`，检验干湿季和8项环境指标。两个入口运行同一版代码。
4. 在04中选择同一分析范围后，运行 `scripts/04_plot_permanova.R`。它读取02结果及 `<本次MC8结果目录>/<scope>/repeat_groups` 的新03结果。

## 组成数据如何处理

01仍将同一个Repeat_group_ID内的原始相对丰度取均值。

03的 `make_repeat_data()` 只筛选目标TS并排序，每个Repeat_group保留一行。同一宿主、同一TS有多个重复组时，保留全部组，不再按Host_type与TS取均值，不要求完整宿主覆盖，不补缺采格子。

随后对各组的8个MC丰度开平方根，并计算Bray–Curtis距离。原始相对丰度不对所选8个MC重新归一化。

| 范围 | 干湿季输入 | 环境输入 | 宿主数 |
|---|---:|---:|---:|
| all_hosts | 101个Repeat_group，TS1–TS8 | 89个Repeat_group，TS1–TS7 | 12 |
| coral | 39个Repeat_group，TS1–TS8 | 33个Repeat_group，TS1–TS7 | 3 |

这些数量按当前输入核对。01沿用全局只有一个Repeat_group的宿主排除规则。

## 检验和置换

每项指标单独检验，模型为 `distance ~ 指标 + Condition(Host_type)`。统计量对应adonis2模型 `distance ~ Host_type + 指标` 的指标边际效应。

代码使用标准0/1类别编码表示宿主，调整效果与Condition(Host_type)相同。这样可以避开当前vegan在某些负约束轴置换中计算分类质心的报错，保留全部置换分配。

MC行和宿主保持原位，只在TS之间重新分配季节或环境标签。某个TS的所有Repeat_group行共享同一个置换值，各TS行数可以不同。两种分析范围都仍有7个环境批次和8个季节批次。

季节完整枚举70种不同分配。环境完整枚举各指标全部不同的批次值分配；7个互异值有5040种。8项环境检验统一做BH校正，原始Bray为主要结果，Lingoes为敏感性检查。

TS8环境缺测，因此只进入宿主、时间和季节分析。DIP部分为检测限记录、精确批次仅4个，未纳入8项环境主检验。

批次标签检验采用TS之间可交换的假设，宿主作为调整项。重复组之间的组成差异保留，各宿主、各批次拥有的组数可以不同。

## 输出和绘图

03结果在 `<本次MC8结果目录>/<scope>/repeat_groups`，其中scope为all_hosts或coral。

- `season_results.csv` 与 `environment_results.csv`：按Repeat_group重新计算的结果。
- `season_analysis_input.csv` 与 `environment_analysis_input.csv`：每行一个Repeat_group，保留组号和n_records。
- `host_batch_coverage.csv`：每个宿主每个TS的Repeat_group数量，0表示未采到。
- `season_environment_analysis.rds`：结果和对应的作图输入。

04使用新的结果和输入。DO组成预览中，每根柱子对应一个Repeat_group，同宿主同批次的不同组分别显示。

旧均值版脚本和被更新的结果保存在archive的本次备份目录。旧版p、R2和图对应旧输入，应与新结果分别使用。
