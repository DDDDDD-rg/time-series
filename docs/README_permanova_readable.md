# 当前入口：MC19 coral原始样本环境分析

2026-10-27核对：v2已于2026-10-26完成，当前运行编号20261026_224333。vegan负责每次拟合的R²、F、平方和与自由度，脚本负责TS环境标签置换及P值汇总；现场状态为analysis_and_plots_complete。

在RStudio打开scripts/PERMANOVA/run_permanova.R并点击Source，默认完成7项环境分析及配套图；只重画图运行scripts/PERMANOVA/04_plot_permanova.R。完整说明见[scripts/README.md](../scripts/README.md)。

仅coral、19个MC、7项当期环境；原始样本逐行保留，无Repeat_group平均或Host_type模型。7个TS仍完整枚举5040种环境标签分配，MC样本行固定；每次调用vegan::adonis2(permutations=0)，P值仍按F的排列分布汇总，7项统一BH。其他筛选、转换、检测限和作图设置不变。

新的results运行目录顶层只保留environment_results.csv主表和supporting_files/配套文件夹；完整结果、输入、快照、置换记录、状态、日志、RDS和NMDS缓存都在配套文件夹。04、07及环境组成图通过共用读取器兼容新布局和旧v1平铺结果；当前MC19结果保持原路径；历史MC8/v1结果及图已成组归档，旧路径有兼容联接。图运行目录保留3份合并PDF，单图和配套资料归入supporting_files；作图后自动刷新合并PDF。详见[更新日志](../CHANGELOG.md)。

下文是2026-10-14的MC8历史方法，不适用于当前一键入口；原代码已备份至archive/script_backups/。

# MC8 PERMANOVA：加入 Seawater 和 Sediment 的版本

本版本已于2026-10-14完成运行。当前结果：[PERMANOVA_MC8_20261014_183535](../results/PERMANOVA_MC8_20261014_183535/)。历史结果与图已成组归档，旧路径留有兼容联接。

## 运行方式

在 R/RStudio 中运行 `scripts/PERMANOVA/run_permanova.R`。它建立新的时间戳结果和图目录，完整执行以下内容，成功后更新 latest 指针。

| 范围 | 当前可用重复组 | 分析 |
|---|---:|---|
| all_hosts | 115：14类，包括 Seawater、Sediment | 宿主、TS、干湿季、8项环境指标、宿主/TS两两比较 |
| coral | 39：3类珊瑚 | 相同检验，仅用珊瑚 |
| sample_type | 117：HOST=103、Seawater=7、Sediment=7 | 三类总体比较和三对BH校正比较，均调整TS |

逐宿主分析沿用至少两个Repeat_group的要求，因此 Cuttlefish、Sargassum 这两个各只有一个组的宿主未进入14类检验。三类分析将其并入HOST，保留全部117组。

## 数据处理和缺采

输入为 `data/PERMANOVA_input_MC8_20261011.xlsx`，8个MC为 MC6、MC507、MC21、MC520、MC10、MC15、MC0、MC25。Repeat_group内原始比例先取均值，再平方根转换，计算Bray–Curtis；不把8个MC重新归一化为100%。

Seawater、Sediment仅在TS2–TS8采样。TS1没有这两类，不补零、不创建虚拟行；其他已有TS1样本保留。宿主/TS/季节检验使用所有已观测且符合原筛选规则的行。

每个环境指标独立筛选有效批次，不因另一指标缺失而删掉该行。目前8项主环境指标均为TS1–TS7有效、TS8缺测：all_hosts的环境分析101组，coral的环境分析33组。TS8仍参与宿主/TS/季节分析。DIP的检测限记录保留，但不纳入8项环境主检验。

## 检验及BH范围

| 问题 | 模型与置换 | p值 |
|---|---|---|
| 宿主/样本类型 | Host_type + TS，读取Host_type；TS内置换Repeat_group | 与TS总体检验合为两项BH |
| TS | Host_type + TS，读取TS；Host_type内置换Repeat_group | 同上 |
| 干湿季 | 调整Host_type；枚举TS层面的干湿标签 | 单项原始p |
| 环境 | 调整Host_type；逐个变量枚举其有效TS的环境标签 | 每个范围的8项BH |
| 三类样本总体 | Sample_type + TS，读取Sample_type；TS内置换 | 单项原始p |
| 三类两两比较 | 每对只用共同TS，调整TS；TS内置换 | 3对单独BH |

三类比较不同时调整细分Host_type，因为HOST/Seawater/Sediment是由细分类别确定的。HOST代表此次采集的生物宿主相关样本的混合组成。

原始Bray为主，Lingoes为敏感性检查。分类因素保留PERMDISP；两两比较的PERMDISP单独BH。每个两两比较族按全部计划组合数校正，无法估计的组合保留NA。14类的宿主对为91对，仅珊瑚为3对，各范围TS为28对，三类样本为3对。

## 图和HTML

所有图为英文、无图题。保留高对比色、白底、统一圆点及现有细线椭圆风格。

`04_plot_permanova.R`每个范围只导出两种三维HTML：`05_host_NMDS_3D.html`和`06_time_NMDS_3D.html`，不新增季节、环境或三类样本的三维图。

- 鼠标悬停图例名称：该组点放大，全部组保留原色、不灰化。
- 单击名称：选中/取消该组，可同时选中多组；未选组灰色淡化。
- 有选中组时再悬停其他组：暂时恢复全部原色并突出悬停组；移出后恢复已点击的选择。
- `Show all`清除选择。默认点径3、选中4、悬停5。未选背景点径2、透明度0.18。
- 坐标轴的位置、范围和点的XYZ坐标固定；可旋转和缩放，使用正交投影。

新全样本包含更多点，旧101组的NMDS坐标缓存会被拒绝，运行时需要重算；之后数据完全一致的展示修改可以复用。三维随机起点保留4个PSOCK工作进程。完整入口不再在02重复拟合一遍二维NMDS；04计算并保存坐标。

`06_plot_pairwise_permanova.R`读取最新完成结果，输出五张热图（原四张加三类样本热图），统一线性R²色标，必要时自动扩大上限。`07_plot_host_highlight_panels.R`也读取新结果，仍每行3列，14类时自动为5行。DO组成预览读取DO检验的实际数据。

HTML若不能嵌入全部依赖，需与同名 `_files` 目录一起保存和分享。

## 文件入口

- `run_permanova.R`：完整运行、独立结果目录、汇总与脚本快照。
- `01_data_settings.R`：读取、分组、筛选、共用检验函数。
- `02_test_host_time.R`：宿主和TS总体检验。
- `03_test_season_environment.R`及`03_test_season_environment_batch_labels.R`：季节和环境批次标签检验。
- `05_test_pairwise_host_time.R`及`permanova_pairwise_helpers.R`：两两比较。
- `08_test_sample_type.R`：HOST/Seawater/Sediment总体及两两比较。
- `04_plot_permanova.R`、`06_plot_pairwise_permanova.R`、`07_plot_host_highlight_panels.R`、`plot_DO_composition_preview.R`：作图。

根目录的`00_all_tests_summary.csv`有23项总体/季节/环境结果（含三类总体）。`00_PERMANOVA_summary.html`为总览；`00_pairwise_summary.csv`及`01_all_pairwise_results.csv`含5个比较族。新运行会保存输入、脚本副本、MD5和会话信息。
