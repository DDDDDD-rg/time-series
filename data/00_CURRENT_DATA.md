当前结果入口见[项目首页](../README.md)（2026-10-14更新）。以下保留已有数据转换记录。

# 当前分析输入 · 2026-10-11

现在使用这两份文件（均在本文件同一层）：

| 用途 | 最新文件 | 数据规模 |
|---|---|---|
| GLMM | [GLMM_input_MC8_20261011.xlsx](GLMM_input_MC8_20261011.xlsx) | Input_long：400个样本 × 8个MC = 3,200行 |
| PERMANOVA | [PERMANOVA_input_MC8_20261011.xlsx](PERMANOVA_input_MC8_20261011.xlsx) | Community：400行、8个MC；Metadata：400行 |

MC顺序：MC6、MC507、MC21、MC520、MC10、MC15、MC0、MC25。重复提及的MC507仅保留一次。

MC0按New_ID从原始MC绝对/相对丰度表匹配补入。其他7个保留MC的丰度及所有样本元数据已逐项核对。全部零值保留；相对丰度分母仍为全部Ruegeria reads，不按8个MC重新归一化。

每个MC保留400条记录。DO完整的348条记录属于TS1–TS7、103个Repeat_group；TS8的52条DO缺测记录仍在输入表中，保持空白。DO_c仍为DO_mg_L−6.3。

## 文件夹怎么找

| 位置 | 内容 |
|---|---|
| 本层 | 两份最新MC8输入、本说明 |
| source_tables/ | 原始DNA表、样本元数据、指示种/网络资料 |
| mc/ | 原始MC绝对与相对丰度TSV；未修改 |
| asv/ | 原始ASV资料；未修改 |
| water_quality/ | PM11原水质CSV及营养盐/Prism制图表 |
| archive/ | 旧MC11输入及之前的GLMM备份；保留原文件内容 |

Cooccurence.xlsx当前被其他程序占用，暂留data顶层；本次未关闭该程序或修改该文件。

## 脚本入口

- scripts/GLMM.R、scripts/GLMM_reads_diagnostics_full.R 已指向本层新版GLMM输入，默认读取全部8个MC。
- scripts/01_data_settings.R、scripts/run_permanova.R 已指向本层新版PERMANOVA输入；02、03、05继续通过01读取设置。
- scripts/plot_PM11_timeseries_5years.R 已指向water_quality中的原水质CSV。

MC8 PERMANOVA已重新运行：[最新结果总览](../results/PERMANOVA_MC8_20261011_225709/00_PERMANOVA_summary.html)。旧MC11结果归入results/PERMANOVA_MC11_history，旧图归入figures/PERMANOVA_MC11_history。当前MC8的GLMM已完成reads零值和波动四候选比较：[最新GLMM总览](../results/GLMM_joint_reads_MC8_20261012_001527/00_all_MC_joint_summary.html)。
