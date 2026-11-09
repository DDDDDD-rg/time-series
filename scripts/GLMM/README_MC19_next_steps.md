# MC19两部分分析与talk作图

更新：2026-11-02。只接收新的8TS/五营养盐原GLMM和PERMANOVA结果。旧7TS结果表与图保留。本次准备未运行研究拟合或新NMDS。

## Talk图

运行`plot_MC19_talk_figures.R`。先完成同范围的原GLMM与PERMANOVA，脚本读取保存统计量。

```r
options(mc19.analysis_scope="coral")
options(mc19.talk.parts=c("forest","NP","NMDS2D","NMDS3D"))
```

NP图包含NO3-N、NO2-N、Ammonia-N、DIN、DIP和全部8TS，NO2/DIP报告限虚线沿用。DIN无独立的直接测量报告限；空心点表示其至少一个成分使用L/2。

同范围同维度的五个NMDS环境面板共用坐标。交互3D右侧可在五指标间切换，轴线/刻度按固定空间棱绘制，轴名与统计说明在图外；每项显示保存的PERMANOVA R²、原始P和BH q。DIP单位mg P/L，其余mg N/L。

热图五指标全列出，旁边森林图继续只列存在显著MC的指标。所有点估计与完整95%CI保留，显著符号读取新批次95项BH结果。

## 两部分分析

```r
options(mc19.analysis_scope="coral",mc19.two_part.stage="check")
# Source GLMM_MC19_two_part.R，核对新的显著MC×环境组合
options(mc19.two_part.stage="fit")
# 再Source同一入口，正式拟合和作图
```

检出部分为Bernoulli/logit，使用全部样本；正值部分为Beta/logit，仅MC_reads>0的样本。两部分均`E_z+(1|TS)`，沿用新原GLMM的8TS中心/SD、原分母，不新增Host/Repeat/深度协变量。正值子集若缺某TS，如实记录有效TS数。

任务数由新原GLMM的q<0.05组合决定，每范围合并2×组合数校正BH，失败项保留原因和NA，不缩小检验族。两部分仍是同一数据经原GLMM筛选后的探索性分解。

森林图重点MC15、MC98、MC80、MC0、MC10排列在上方；`*`为det q<0.05，`#`为pos q<0.05。95%CI为逐项Wald区间，q不重算CI。曲线保留全部原始点、完整区间和独立Zoom。

只重画：设置`options(mc19.two_part.stage="plot",mc19.two_part.run_dir="")`，读取最新完成的8TS两部分结果。继续拟合需显式指定同版本目录并通过来源校验。

[运行顺序](../README.md) · [命名规则](../README_run_names.md)
