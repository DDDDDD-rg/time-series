# 七TS：六项环境指标，全量重跑与分开作图

入口：[TS7_temperature.R](TS7_temperature.R)。2026-11-16按用户要求加回DIN、保留DIP，使用新的 `TS7_env6` profile。本版全量重跑，替代此前“只补温度、复用旧结果”的默认入口。脚本已修改；本轮不执行正式分析。联合GLMM仍保持原状。

## 范围和模型

- TS1–TS7；coral 98、ecosystem 310个原始样本，分范围运行。ecosystem包含coral；分母仍是原175行Ruegeria MC（含Unassigned），不平均、不加入Host/Repeat调整。
- 六项：NO3-N、NO2-N、总氨氮TAN（图上简写NH4）、DIN、DIP、温度SST_C。
- DIN由三项氮组分相加，单位mg N/L。低于报告限的组分先按原规则取L/2，再求和；保留组分删失标记。DIN是计算值，不单设直接测量检出限。无需重新填写DIN或修改原始工作簿。
- 温度继续使用七TS的PM11表层原CSV逐日期核对值：25.8、31.4、28.6、16.7、24.6、26.6、25.2°C。
- GLMM：19MC×6指标＝114个模型/范围，全部重新拟合；每次只用一个环境变量，公式 `cbind(MC_reads, Other_reads) ~ E_z + (1 | TS)`，Beta-binomial/logit、恒定离散度、无零膨胀。E_z按七个TS等权中心化/标准化。114项Wald P统一BH，数值失败行保留，不缩小计划检验族。
- GLMM诊断：沿用保存模型的DHARMa及原有分布摘要检查；默认每模型2000次模拟。诊断模拟不用于计算环境系数P值。此版重新生成诊断，不从旧五指标批次迁入缓存。
- PERMANOVA：六个单环境检验重新计算，原19MC比例逐元素sqrt后用Bray–Curtis；7!=5040种整TS环境标签分配（含原顺序），统计量仍来自vegan::adonis2；六项原始P统一BH。不复用旧P/F/R²。
- two-part：只从本次114项GLMM中选数值有效且BH q<0.05的MC×环境组合，重新拟合检出Bernoulli/logit、正值Beta/logit，均为 `E_z + (1 | TS)`。当前所有组件P在 `2×入选组合数` 的族内统一BH。没有显著组合则明确跳过。已有组件仅支持同一新批次中断后的续跑，不自动跨旧批次复用。
- 本轮沿用既有方法，不新增联合GLMM、dbRDA、bootstrap、LOTO或其它敏感性模型。

## 运行方式

打开 `scripts.Rproj`。先清除手动固定的旧结果/作图路径和run编号：

```r
options(
  mc19.analysis_scope = "coral", # 另一套改为 "ecosystem"
  mc19.glmm.result_dir = "",
  mc19.result_dir = "",
  mc19.two_part.run_dir = "",
  mc19.talk.output_dir = "",
  mc19.figure_dir = NULL,
  mc19.run_id = NULL,
  mc19.run_diagnostics = TRUE
)
```

只读核对输入和配置：

```r
options(mc19.temperature.stage = "check")
source("TS7_temperature.R", encoding = "UTF-8")
```

默认就是check。核对七TS、样本/MC分母、六项环境、DIN求和、两套分析的环境和样本匹配、温度来源；不依赖旧模型结果，不创建研究结果目录，不拟合、不模拟、不置换、不做NMDS或图。

准备好正式计算时，只运行结果阶段：

```r
options(mc19.temperature.stage = "results")
source("TS7_temperature.R", encoding = "UTF-8")
```

按顺序执行GLMM及诊断→PERMANOVA→筛选后的two-part，关闭自动作图。每次results都会创建新批次并全量计算；只想重画时运行plots。

结果完成后单独作图：

```r
options(mc19.temperature.stage = "plots")
source("TS7_temperature.R", encoding = "UTF-8")
```

plots读取新的六指标结果，不重拟合GLMM/two-part，不重算PERMANOVA P或BH；NMDS坐标在作图阶段计算/缓存。图中的DIN、温度和DIP均按字段名加入，不按固定列号匹配检出限。

两范围都计算时，可在清除上述旧路径后使用：

```r
for (scope in c("coral", "ecosystem")) {
  options(mc19.analysis_scope = scope, mc19.temperature.stage = "results")
  source("TS7_temperature.R", encoding = "UTF-8")
}
```

两范围都作图时将循环中的 `"results"` 改成 `"plots"`。不要在同一范围旧结果尚未完成时启动第二个同类批次。

分项阶段保留：`GLMM`、`PERMANOVA`、`two_part`只生成相应结果；`plot_GLMM`、`plot_PERMANOVA`、`plot_two_part`分别作图；`talk`生成汇报图；`all`连续结果和图。推荐使用分开的results/plots。

## 输出和兼容性

- 方法版本带 `env6_TS7_full`，目录例如 `MC19_GLMM_TS7_env6_coral_run001`、`MC19_PERMANOVA_TS7_env6_coral_run001`、`MC19_two_part_TS7_env6_coral_run001`，并有独立latest指针。ecosystem同理。
- 旧四营养盐七TS结果、五指标 `TS7_temp` 增量结果、早期 `TS7_temperature` 结果及其图均保留；新流程不会把它们读作本次结果。历史五指标profile `TS7_temperature` 的定义保留供旧图读取。
- 默认八TS五营养盐准备流程 `TS8_env5` 不变。当前入口显式选七TS，不是因TS8缺失而自动回退。
- 既有绘图函数按六指标生成热图/森林图、显著响应曲线、诊断图、环境时间线、组成梯度和NMDS2D/3D、交互按钮及two-part。当前仍不加前四/后四TS虚线。
- 不自动修改PPT；主图合并与辅助图分开的整理要求继续沿用。此前 `temperature_review_v1` 是旧五指标检查包，内容不自动替换为本次结果。

## 验证和依赖

依赖沿用readxl、glmmTMB、DHARMa、vegan、ggplot2、plotly、htmlwidgets；合并PDF需Python/pypdf。脚本不自动安装包。

`test_TS7_temperature.R`已更新为六指标/114项，并覆盖DIN求和、输出版本隔离、六列图形、六指标按钮、结果/作图阶段分离。该回归脚本使用人工夹具与分析拦截，不应产生真实拟合或置换；本轮只做源码解析、只读输入检查和阶段调度检查，不执行完整渲染回归或正式统计分析。历史增量模块及其测试保留。

本次全量版本计算量明显高于之前的复用版本：两范围共228个GLMM及其诊断、12个完整标签检验，two-part数量由新GLMM结果决定。
