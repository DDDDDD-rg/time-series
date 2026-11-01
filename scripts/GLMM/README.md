## 当前状态与工程修复（2026-10-31）

Coral/ecosystem主分析均已运行完成，two-part为8/14组合、16/28分模型。当前主曲线为16:9、每页最多2MC（Coral3页、ecosystem4页），主图full range、右侧独立Zoom。历史段落中的“未运行”、16pt、每页3MC及旧目录示例仅描述当时状态；当前字号及版式以GLMM/README.md最新规则为准。

Two-part作图先执行只读深度校验，核对原始输入、任务键、模型系数、主/明细表、BH及预测/观测缓存；旧批次可直接验证，不修改文件。新拟合批次另写artifact_fingerprints.csv。只改字体不改变分析签名，分析函数改变会阻止混用；旧签名通过保存源码验证后仅在用户明确恢复fit时迁移。

批量样式刷新按kind/scope/source和完成状态定位plot修订及talk缓存，不再硬编码日期；要求预测快照与结果一致、GLMM与PERMANOVA输入一致。HTML先在短暂存目录导出到短lib依赖路径，核对资源后发布；过深路径明确报错，不发布残缺HTML。

本轮只修代码和说明，未拟合、未模拟、未重算NMDS，现有图/结果保留。新图注在下次用户重画时生效：空心点表示营养盐低于LOD、占比分母为总Ruegeria reads，two-part正值为0–1比例。现有图没有覆盖；新增Zoom和字体嵌入不是本轮执行内容。

---

以下为已有说明与历史更新记录。

# GLMM脚本

- [GLMM.R](GLMM.R)：当前76项TS随机截距Wald/BH完整入口。
- [GLMM_reads_diagnostics_full.R](GLMM_reads_diagnostics_full.R)：只重做保存模型的诊断。
- [plot_GLMM_MC19.R](plot_GLMM_MC19.R)：只重画当前主图与诊断图。
- `glmm_MC19_pilot_helpers.R`、`glmm_MC19_wald_helpers.R`、`glmm_MC19_wald_plot_helpers.R`：当前共享定义。
- `GLMM_leave_one_TS_out.R`、`glmm_MC19_loto_helpers.R`、`plot_GLMM_MC19_LOTO.R`：历史LOTO。
- `glmm_MC19_inference_helpers.R`、`plot_GLMM_MC19_full_data.R`：历史bootstrap流程辅助。
- `08_plot_GLMM_host_beta_heatmap.R`：旧MC8宿主图，样式从PERMANOVA/permanova_plot_style.R读取。

输入、结果和图均继续写入项目data/results/figures，不写入本文件夹。[完整运行说明](../README.md)、[代码指南](../../docs/GLMM_MC19_env4_Wald_code_guide_20261027.md)。

## Coral / ecosystem 范围切换

在R控制台先设置一次，然后source本目录的主入口或配套入口：

```r
options(mc19.analysis_scope = "ecosystem")  # 全部生物宿主；排除seawater/sediment
# options(mc19.analysis_scope = "coral")  # 切回珊瑚；默认值
```

保留原始独立样本，不平均、不使用Repeat_group或Host项；缺当期环境的样本不填补。相同MC和指标，两范围分别校正BH。当前输入完整环境为coral 98样本、ecosystem 310样本，均7TS；范围标签写入结果，结果与图目录均带范围。

主入口GLMM.R：19MC×4营养盐，公式`cbind(MC_reads, Other_reads) ~ E_z + (1 | TS)`，每范围76项Wald P单独BH。ecosystem从头拟合76项，不复用coral先导。

结果在`results/MC19_GLMM_<scope>_run001/`等递增批次，配套图在`figures/MC19_GLMM_<scope>_run001_plot001/`等绘图修订；旧时间戳结果继续可读；主结果00_result.csv，模型/诊断/输入在supporting_files。新latest为`00_latest_GLMM_MC19_env4_TS_Wald_<scope>.txt`，摘要为`00_LATEST_GLMM_<scope>.md`。coral无新指针时兼容读取旧无范围指针；ecosystem不会回退读coral。诊断与作图入口也使用同一选项；手动指定results_dir必须属于该范围。

本轮只修改代码；未运行模型、置换、模拟诊断或正式作图。

## 2026-10-29：全部显著关系曲线

Coral和ecosystem共用同一作图规则：从当前范围保存的76模型结果中选择数值有效且`beta_q_value < 0.05`的每一对MC×营养盐，不要求同一MC在多个指标显著。按MC分组及环境配置顺序，每页最多3个MC，同MC的显著指标同行（2–4列，列宽固定），输出`02_all_significant_nutrient_curves.pdf`及逐页PNG；面板/页码索引在supporting_files/significant_curve_index.csv。没有显著关系时输出一页英文说明，热图仍照常生成。

保持16pt、英文、无总图题、TS点色、原始全部样本、b_TS=0预测线及95%CI、检出限空心点。同MC的显著面板统一纵轴；少数极端值压缩主体时放大主图，完整范围小图保留全部观测值。热图及76项BH原值不变。重画时把旧三MC曲线和被替换的分页曲线移入supporting_files/previous_main_figures，避免旧页混入。只运行plot_GLMM_MC19.R即可，无需重跑GLMM。此次仅修改与无绘图检查。

2026-10-29版式修订：同一MC的所有显著指标放在同一行、同一页（最多4列，列宽固定，页面宽度随列数扩展）；每页最多3个MC，单指标MC独占一行。底部预留1.4英寸、行间1英寸，轴标签和三行简短图注分开。相同MC主图纵轴一致；极端值压缩曲线时主图改为放大范围，标记Full range的小图显示全部原始点与完整纵轴，自动选择较少遮挡曲线的上方角落。移除原1.8倍纵轴留白。只修改源码、检查分组及版式坐标，不运行正式绘图。

2026-10-29当前版：用户确认大图必须full range，小图显示局部区域。大图覆盖全部原始点及95%CI，同MC统一纵轴，不加1.8倍留白。需要放大时，在大图右侧独立留出Zoom小图区（包含边框、刻度、文字的区域也不与大图相交），相同x范围、缩小y范围；不再向大图绘图区内叠放小图。列距自动增至7.8英寸，主绘图区宽度和16pt基准字号保持。显著组合/同行同页/模型/Pq/热图均不变。supporting_files/curve_axis_ranges.csv记录主图和放大图范围。此次只修改与只读/坐标检查，现有图由用户重画。


## 运行命名与最新结果

当前主流程采用递增run编号，重新作图采用plot修订号，时间/方法版本/来源保存在运行记录中。新旧latest指针均兼容，自动跳过未完成批次，分别读取coral/ecosystem。[完整命名与读取规则](../README_run_names.md)。

## 2026-10-31：当前PPT横版规则

当前规则覆盖上面的历史版式说明：热图＋森林图为32×18in，热图色块区加宽至10in，森林图按保留的2/3个环境面板分配剩余宽度；保留19MC、四列热图与完整CI。原GLMM显著曲线每页最多2个MC，同MC同行同页，页面按16:9输出（当前19.2×10.8in），主图full range、右侧独立Zoom；Coral为3页、ecosystem为4页，页码索引同步更新。

NMDS二维与三维轴对投影为16×9in，Shepard为12.8×7.2in，N/P折线为18×10.125in，two-part森林图为20×11.25in。原本清楚的two-part曲线、残差及GLMM补充诊断保留横版。字体继续放大加粗。3D交互图保持既有切换、固定轴和点size3.3/opacity0.65。

`refresh_MC19_figure_style.R`复用保存预测、诊断及NMDS缓存。可用`options(mc19.style.diagnostics=FALSE, mc19.style.two_part.parts="forest")`跳过未改版的诊断/两部分曲线与残差；默认仍重画全部。此次已用本地R重画并覆盖现有主图，不重新拟合或计算P/q。
