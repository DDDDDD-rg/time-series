## 当前状态与工程修复（2026-10-31）

Coral/ecosystem主分析均已运行完成，two-part为8/14组合、16/28分模型。当前主曲线为16:9、每页最多2MC（Coral3页、ecosystem4页），主图full range、右侧独立Zoom。历史段落中的“未运行”、16pt、每页3MC及旧目录示例仅描述当时状态；当前字号及版式以GLMM/README.md最新规则为准。

Two-part作图先执行只读深度校验，核对原始输入、任务键、模型系数、主/明细表、BH及预测/观测缓存；旧批次可直接验证，不修改文件。新拟合批次另写artifact_fingerprints.csv。只改字体不改变分析签名，分析函数改变会阻止混用；旧签名通过保存源码验证后仅在用户明确恢复fit时迁移。

批量样式刷新按kind/scope/source和完成状态定位plot修订及talk缓存，不再硬编码日期；要求预测快照与结果一致、GLMM与PERMANOVA输入一致。HTML先在短暂存目录导出到短lib依赖路径，核对资源后发布；过深路径明确报错，不发布残缺HTML。

本轮只修代码和说明，未拟合、未模拟、未重算NMDS，现有图/结果保留。新图注在下次用户重画时生效：空心点表示营养盐低于LOD、占比分母为总Ruegeria reads，two-part正值为0–1比例。现有图没有覆盖；新增Zoom和字体嵌入不是本轮执行内容。

---

以下为已有说明与历史更新记录。

# 时间序列分析运行入口

更新：2026-10-29。当前保留MC19 GLMM和PERMANOVA两条流程。历史源码随所属流程放入对应文件夹供复核，历史运行产物已归档。不要按编号顺序执行全部脚本。

## 文件夹与主入口

| 目录 | 主入口 | 内容 |
|---|---|---|
| [GLMM/](GLMM/README.md) | [GLMM.R](GLMM/GLMM.R) | GLMM分析、诊断、绘图及先导/LOTO/bootstrap辅助源码 |
| [PERMANOVA/](PERMANOVA/README.md) | [run_permanova.R](PERMANOVA/run_permanova.R) | PERMANOVA分析、绘图、历史两两比较和PDF合并 |

`scripts.Rproj`和`plot_PM11_timeseries_5years.R`留在本层。当前源码、快照复制路径以及跨目录样式依赖均已更新；历史运行的scripts_used、RDS和指纹记录保留其运行时内容。

## 当前GLMM

在RStudio打开`GLMM/GLMM.R`后Source。输入为`data/PERMANOVA_input_MC19_env7_20261026.xlsx`及`data/mc/MC_abundance_absolute.tsv`，按New_ID对应所选范围原始样本（当前coral 98个，ecosystem 310个，均7TS），分母包含全部Ruegeria MC，不对19MC重新归一化。

```r
cbind(MC_reads, Other_reads) ~ E_z + (1 | TS)
```

19MC×硝酸盐、亚硝酸盐、氨氮、DIP，共76个Beta-binomial/logit模型；恒定离散度，无额外零膨胀。每模型一个营养盐，按7TS等权标准化；不平均，无Host/Repeat/水温项。直接提取Wald P，76项统一BH，失败项保留NA但不缩小校正族。

Coral默认`reuse_pilot_fits=TRUE`，核对120255先导的4个有TS硝酸盐拟合后新增72个；先导现位于archive/GLMM_MC19/results，原results路径保留兼容联接。无先导时设为FALSE从头拟合76个。Ecosystem强制关闭coral先导复用，从头拟合76项。主流程不执行bootstrap、置换或LOTO。

| 操作 | 入口 |
|---|---|
| 拟合、诊断、作图 | `GLMM/GLMM.R`；`run_diagnostics`、`run_plots`可分别关闭 |
| 只做诊断 | `GLMM/GLMM_reads_diagnostics_full.R`，读取保存模型，每模型2000次响应模拟，不重拟合 |
| 只重画 | `GLMM/plot_GLMM_MC19.R`，读取保存模型、P/q和诊断，不拟合、不模拟、不重算P/q |
| 当前结果 | [233529主表](../results/GLMM_MC19_env4_TS_Wald_coral_20261027_233529/00_result.csv)、[全精度明细](../results/GLMM_MC19_env4_TS_Wald_coral_20261027_233529/supporting_files/all_results.csv) |
| 方法与代码 | [代码指南](../docs/GLMM_MC19_env4_Wald_code_guide_20261027.md) |

分析、诊断和作图共用`GLMM/glmm_MC19_pilot_helpers.R`、`GLMM/glmm_MC19_wald_helpers.R`、`GLMM/glmm_MC19_wald_plot_helpers.R`。pilot helper仍是当前输入验证和诊断的依赖。

当前输出：results同名目录顶层`00_result.csv`、`00_readme.txt`、supporting_files；CSV重要参数在前，全精度明细用于后续读取。figures同名目录为效应热图和全部显著关系分页图，各有PDF/PNG；诊断、预测与被替换的五份主PDF在supporting_files。运行后使用MC19 TS-only专属latest指针。

在其它电脑运行需准备原输入、设置入口中的project_dir及安装readxl/glmmTMB/DHARMa。GitHub审阅材料不包含原始输入和模型缓存；单独诊断/重画仍需要本地保存模型。

## 当前PERMANOVA

完整入口为`PERMANOVA/run_permanova.R`；当前批次[PERMANOVA_MC19_env7_coral_samples_20261026_224333](../results/PERMANOVA_MC19_env7_coral_samples_20261026_224333/)已完成。98个coral原始样本、19MC、7个有效TS、7项环境；无Repeat_group平均或Host_type调整。

MC原始比例平方根变换、Bray–Curtis；5040种TS环境标签分配含实际分配。每项每次由vegan::adonis2读取统计量，P按F置换分布汇总，7项统一BH。不运行MC8的宿主/两两/sample_type分支。

| 操作 | 入口 |
|---|---|
| 完整分析与作图 | `PERMANOVA/run_permanova.R` → `PERMANOVA/03_test_season_environment.R` → `PERMANOVA/03_test_season_environment_batch_labels.R` |
| 只重画完整图 | `PERMANOVA/04_plot_permanova.R` |
| 描述性宿主高亮 | `PERMANOVA/07_plot_host_highlight_panels.R` |
| 环境梯度组成图 | `PERMANOVA/plot_MC_environment_composition.R`；旧名`PERMANOVA/plot_DO_composition_preview.R`转至相同入口 |
| 共享设置与作图 | `PERMANOVA/permanova_MC19_data_settings.R`、`PERMANOVA/permanova_MC19_plot_helpers.R` |
| PDF合并 | `PERMANOVA/bundle_MC19_figures.py`，依赖Python/pypdf，支持旧平铺和supporting_files布局 |
| 方法说明 | [README_permanova_readable.md](../docs/README_permanova_readable.md) |

默认输出三份合并PDF（10/3/7页）和supporting_files；只重画不重跑35280次完整置换。R依赖readxl/vegan/ggplot2；交互3D另需plotly/htmlwidgets。

## 保留的历史及辅助源码

- `GLMM/GLMM_leave_one_TS_out.R`、`GLMM/plot_GLMM_MC19_LOTO.R`、`GLMM/glmm_MC19_loto_helpers.R`：历史硝酸盐LOTO。
- `GLMM/plot_GLMM_MC19_full_data.R`、`GLMM/glmm_MC19_inference_helpers.R`：历史bootstrap结果读取和辅助定义。
- `PERMANOVA/01_data_settings.R`、`PERMANOVA/02_test_host_time.R`、`PERMANOVA/05_test_pairwise_host_time.R`、`PERMANOVA/06_plot_pairwise_permanova.R`、`PERMANOVA/08_test_sample_type.R`、`PERMANOVA/09_plot_sample_type.R`及`PERMANOVA/pairwise_letters.R`、`PERMANOVA/permanova_pairwise_helpers.R`、`PERMANOVA/permanova_plot_style.R`：MC8历史方法。`test`在这些名称中表示统计检验。
- `GLMM/08_plot_GLMM_host_beta_heatmap.R`：旧MC8宿主热图。
- `plot_PM11_timeseries_5years.R`：环境时间线作图；`scripts.Rproj`继续作为RStudio入口。

原`.RData`、`.Rhistory`移至archive/r_sessions；Python字节码缓存移至archive/tests。RStudio项目状态`.Rproj.user`保留原路径并隐藏。旧的完整脚本README快照在[archive/navigation/scripts_README_before_organization_20261028_214616.md](../archive/navigation/scripts_README_before_organization_20261028_214616.md)，早期记录中的默认参数以当时版本为准。

## 两个范围共用的切换选项

在控制台先设置：

```r
options(mc19.analysis_scope = "ecosystem")
# options(mc19.analysis_scope = "coral")  # 默认值
```

随后Source GLMM/GLMM.R或PERMANOVA/run_permanova.R；独立诊断/作图也读取同一选项。每次只处理所选范围。Ecosystem涵盖全部已采集生物宿主，排除seawater/sediment；不按宿主平均或合并reads，不使用Host项。结果目录、图目录和latest指针按范围隔离；GLMM在每范围76项BH，PERMANOVA在每范围7项BH。详细路径及兼容读取见两个子目录README。新增范围尚未运行。

## 2026-10-29：全部显著关系曲线

Coral和ecosystem共用同一作图规则：从当前范围保存的76模型结果中选择数值有效且`beta_q_value < 0.05`的每一对MC×营养盐，不要求同一MC在多个指标显著。按MC分组及环境配置顺序，每页最多3个MC，同MC的显著指标同行（2–4列，列宽固定），输出`02_all_significant_nutrient_curves.pdf`及逐页PNG；面板/页码索引在supporting_files/significant_curve_index.csv。没有显著关系时输出一页英文说明，热图仍照常生成。

保持16pt、英文、无总图题、TS点色、原始全部样本、b_TS=0预测线及95%CI、检出限空心点。同MC的显著面板统一纵轴；少数极端值压缩主体时放大主图，完整范围小图保留全部观测值。热图及76项BH原值不变。重画时把旧三MC曲线和被替换的分页曲线移入supporting_files/previous_main_figures，避免旧页混入。只运行plot_GLMM_MC19.R即可，无需重跑GLMM。此次仅修改与无绘图检查。

2026-10-29版式修订：同一MC的所有显著指标放在同一行、同一页（最多4列，列宽固定，页面宽度随列数扩展）；每页最多3个MC，单指标MC独占一行。底部预留1.4英寸、行间1英寸，轴标签和三行简短图注分开。相同MC主图纵轴一致；极端值压缩曲线时主图改为放大范围，标记Full range的小图显示全部原始点与完整纵轴，自动选择较少遮挡曲线的上方角落。移除原1.8倍纵轴留白。只修改源码、检查分组及版式坐标，不运行正式绘图。

2026-10-29当前版：用户确认大图必须full range，小图显示局部区域。大图覆盖全部原始点及95%CI，同MC统一纵轴，不加1.8倍留白。需要放大时，在大图右侧独立留出Zoom小图区（包含边框、刻度、文字的区域也不与大图相交），相同x范围、缩小y范围；不再向大图绘图区内叠放小图。列距自动增至7.8英寸，主绘图区宽度和16pt基准字号保持。显著组合/同行同页/模型/Pq/热图均不变。supporting_files/curve_axis_ranges.csv记录主图和放大图范围。此次只修改与只读/坐标检查，现有图由用户重画。


## MC19输出命名

当前原GLMM、PERMANOVA、两部分和talk入口共用[运行编号/绘图修订与最新结果规则](README_run_names.md)。共享定义为`mc19_run_io.R`，文件I/O回归检查为`test_mc19_run_io.R`。现有时间戳输出保持原样。
