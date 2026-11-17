# MC19 PERMANOVA：8TS与五项营养盐

更新：2026-11-02。共用输入为[data/MC19_TS8_env5.xlsx](../../data/MC19_TS8_env5.xlsx)。Environment_TS已留TS8整行，四项原始报告值与实际水样日期待填。缺失时停止，不自动排除TS8。

```r
options(mc19.analysis_scope="coral") # 或"ecosystem"
```

主入口[run_permanova.R](run_permanova.R)调用03系列完成统计，再由04出图。19MC原始比例逐元素开平方，计算Bray–Curtis。每项环境单独检验；五指标为NO3-N、NO2-N、Ammonia-N、DIN、DIP。DIN在成分按L/2处理后相加。

完整8TS共有40320种TS标签排列（含实际顺序）；各指标共用同一排列，扩展回逐样本环境值，由vegan::adonis2(permutations=0)读取统计量。五项P在范围内统一BH。标签可交换性仍是该推断设计的条件；8TS样本仍不平均，不加Host或Repeat项。

| 操作 | 入口 |
|---|---|
| 完整统计和配套图 | run_permanova.R |
| 读取新结果重画 | 04_plot_permanova.R |
| 宿主描述图 | 07_plot_host_highlight_panels.R |
| 五指标环境梯度组成图 | plot_MC_environment_composition.R |
| 五指标talk NMDS 2D/3D | ../GLMM/plot_MC19_talk_figures.R |

新版本为`MC19_env5_TS8_<scope>_individual_TSlabels_v1_vegan`。图注、合并PDF列表及快照名称均同步到五指标/8TS，最新读取不会选择旧7TS结果。

现有7TS结果表、模型记录、NMDS坐标和图保留原目录。此轮没有运行置换检验或新NMDS。

[完整运行顺序](../README.md) · [命名规则](../README_run_names.md)
