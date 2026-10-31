# PERMANOVA脚本

- [run_permanova.R](run_permanova.R)：当前MC19 coral/env7完整分析与作图。
- [04_plot_permanova.R](04_plot_permanova.R)：只重画保存结果。
- [07_plot_host_highlight_panels.R](07_plot_host_highlight_panels.R)：描述性宿主高亮。
- [plot_MC_environment_composition.R](plot_MC_environment_composition.R)：环境梯度组成图；旧名plot_DO_composition_preview.R转至相同入口。
- `permanova_MC19_data_settings.R`、`permanova_MC19_plot_helpers.R`：当前共享定义。
- `03_test_season_environment.R`、`03_test_season_environment_batch_labels.R`：当前分析模块，由完整入口调用。
- `bundle_MC19_figures.py`：拼接保存的PDF，参数和结果目录不变。
- 其余编号/两两比较/样式脚本保留MC8历史方法；无需逐个运行。

输入和输出仍位于项目data/results/figures。[完整运行说明](../README.md)、[统计方法](../../docs/README_permanova_readable.md)。

## Coral / ecosystem 范围切换

在R控制台先设置一次，然后source本目录的主入口或配套入口：

```r
options(mc19.analysis_scope = "ecosystem")  # 全部生物宿主；排除seawater/sediment
# options(mc19.analysis_scope = "coral")  # 切回珊瑚；默认值
```

保留原始独立样本，不平均、不使用Repeat_group或Host项；缺当期环境的样本不填补。相同MC和指标，两范围分别校正BH。当前输入完整环境为coral 98样本、ecosystem 310样本，均7TS；范围标签写入结果，结果与图目录均带范围。

主入口run_permanova.R：19MC×7环境，平方根比例的Bray–Curtis；vegan计算F/R²，完整TS环境标签置换并在本范围7项BH。

结果在`results/MC19_PERMANOVA_<scope>_run001/`等递增批次，配套图在同批次的`_plot001/`等绘图修订；主结果environment_results.csv，其余在supporting_files。两范围各自使用`00_latest_PERMANOVA_MC19_env7_<scope>.txt`。04/07/组成图入口遵循同一选项，并检查保存结果范围。ecosystem宿主高亮图自动扩展颜色及面板高度。旧01_data_settings.R中的all_hosts属于历史MC8平均流程，不用于本分析。

本轮只修改代码；未运行模型、置换、模拟诊断或正式作图。


## 运行命名与最新结果

新批次采用run编号；只重画配套图采用plot修订，保留上一版。新旧latest指针与时间戳结果均兼容，分别筛选各范围已完成结果。[完整命名与读取规则](../README_run_names.md)。
