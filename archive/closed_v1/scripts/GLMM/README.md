# MC19 GLMM：8TS与五项营养盐

更新：2026-11-02。当前输入为[data/MC19_TS8_env5.xlsx](../../data/MC19_TS8_env5.xlsx)，TS8环境值待填。已有7TS结果与图保留，本次未重新拟合。

主入口[GLMM.R](GLMM.R)：19MC×NO3-N、NO2-N、Ammonia-N、DIN、DIP，每范围95项。每次一个环境指标，8TS等权标准化；Beta-binomial/logit，固定公式`cbind(MC_reads, Other_reads) ~ E_z + (1 | TS)`。恒定离散度，无额外零膨胀，不平均、不加Host/Repeat/水温；每范围95项Wald P统一BH。

```r
options(mc19.analysis_scope="coral") # 或"ecosystem"
```

完整8TS从头拟合，不复用旧7TS先导。未填齐TS8环境时在创建分析批次前停止。原始reads分母沿用全部Ruegeria MC含Unassigned；19MC不重新归一化。

| 操作 | 入口 |
|---|---|
| 拟合、诊断和配套图 | GLMM.R |
| 已存8TS模型的模拟诊断 | GLMM_reads_diagnostics_full.R |
| 已存8TS结果重画 | plot_GLMM_MC19.R |
| 热图+森林图、NP时间线、NMDS | plot_MC19_talk_figures.R |
| 两部分分析 | GLMM_MC19_two_part.R；默认check，正式运行选择fit |

森林图为逐项95% Wald区间，符号按本批BH q；添加DIN不改变这一含义。五指标标签、DIN的mg N/L单位和8TS色标均已接入；完整原始点/CI和独立Zoom规则沿用。

输出按run/plot递增，当前版本为`MC19_env5_TS8_<scope>_Wald_v1`。最新读取只接受完成的8TS新批次。现有7TS结果不迁移、不删除，历史源码保留在原处。

[完整运行顺序](../README.md) · [两部分与talk入口](README_MC19_next_steps.md) · [命名规则](../README_run_names.md)
