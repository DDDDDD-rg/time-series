## 当前状态与工程修复（2026-10-31）

Coral/ecosystem主分析均已运行完成，two-part为8/14组合、16/28分模型。当前主曲线为16:9、每页最多2MC（Coral3页、ecosystem4页），主图full range、右侧独立Zoom。历史段落中的“未运行”、16pt、每页3MC及旧目录示例仅描述当时状态；当前字号及版式以GLMM/README.md最新规则为准。

Two-part作图先执行只读深度校验，核对原始输入、任务键、模型系数、主/明细表、BH及预测/观测缓存；旧批次可直接验证，不修改文件。新拟合批次另写artifact_fingerprints.csv。只改字体不改变分析签名，分析函数改变会阻止混用；旧签名通过保存源码验证后仅在用户明确恢复fit时迁移。

批量样式刷新按kind/scope/source和完成状态定位plot修订及talk缓存，不再硬编码日期；要求预测快照与结果一致、GLMM与PERMANOVA输入一致。HTML先在短暂存目录导出到短lib依赖路径，核对资源后发布；过深路径明确报错，不发布残缺HTML。

本轮只修代码和说明，未拟合、未模拟、未重算NMDS，现有图/结果保留。新图注在下次用户重画时生效：空心点表示营养盐低于LOD、占比分母为总Ruegeria reads，two-part正值为0–1比例。现有图没有覆盖；新增Zoom和字体嵌入不是本轮执行内容。

---

以下为已有说明与历史更新记录。

# MC19：talk图与两部分分析（2026-10-31）

用户已于2026-10-31完成正式两部分与talk运行，现有结果保留。下列为后续操作说明，更新命名规则不要求重新分析。原测试预览在archive/talk_code_QA_20261031。

## 先画森林图和NP折线

在RStudio设置并运行：

```r
project_dir <- "C:/Users/xzy/Desktop/CUHK_Study/GNBF6010_Research_Projects/Research/Msc main project-time series"
options(mc19.analysis_scope="coral", mc19.talk.parts=c("forest","NP"))
source(file.path(project_dir,"scripts/GLMM/plot_MC19_talk_figures.R"),encoding="UTF-8")
```

再把scope改成ecosystem可画该范围森林图。NP使用相同环境表，只需一份。原plot_GLMM_MC19.R也会新增01_MC19_heatmap_forest，但若只需新图，优先用上面独立入口，避免重画原曲线/诊断。

## 二维/三维NMDS

```r
options(mc19.analysis_scope="coral",mc19.talk.parts=c("NMDS2D","NMDS3D"),mc19.nmds.trymax=500L)
source(file.path(project_dir,"scripts/GLMM/plot_MC19_talk_figures.R"),encoding="UTF-8")
```

两个范围都输出NO₂与氨氮，coral氨氮标为comparison panel，P/q均读原PERMANOVA。默认读对应latest；可用mc19.glmm.result_dir指定GLMM源目录、mc19.result_dir指定PERMANOVA源目录。

二维兼容缓存且最佳解已重复时复用；ecosystem现有二维未重复最佳解，将继续多起点确认。三维需计算，输出HTML及三组轴对投影（明确标3D projection），保存stress/起点重复数/Shepard图。HTML与同名依赖目录一起保存。NMDS不重做PERMANOVA或BH。

默认图输出到figures/MC19_talk_<scope>_run001；若要复用本入口新坐标，可设置mc19.talk.output_dir为上次目录。不同输入指纹不能混用。依赖ggplot2、vegan；三维另需plotly、htmlwidgets，当前本机均已可用。

## 独立两部分分析

```r
options(mc19.analysis_scope="coral", mc19.two_part.stage="check")
source(file.path(project_dir,"scripts/GLMM/GLMM_MC19_two_part.R"),encoding="UTF-8")
```

check是默认，只读核对并打印8个coral或14个ecosystem组合。全部模型用E_z+(1|TS)；检出Bernoulli/logit，正值Beta/logit，原分母及原7TS标准化沿用。无新增Host/Repeat/SST或深度协变量。

正式拟合时：

```r
options(mc19.two_part.stage="fit")
source(file.path(project_dir,"scripts/GLMM/GLMM_MC19_two_part.R"),encoding="UTF-8")
```

每范围16/28任务，结果独立保存results/MC19_two_part_<scope>_run001，拟合后作图。每个分模型保存检查点。恢复时指定同一mc19.two_part.run_dir；源数据、代码、R和glmmTMB版本不匹配会拒绝混用。若想重新拟合，应使用新目录。

只改图（run_dir留空读取该范围最新完成批次，也可指定历史完整路径）：

```r
options(mc19.two_part.stage="plot",mc19.two_part.run_dir="")
source(file.path(project_dir,"scripts/GLMM/GLMM_MC19_two_part.R"),encoding="UTF-8")
```

主表00_result.csv并列原总体与两个分量；supporting_files保留输入、来源、模型、P/q/CI、数值标记、Pearson残差和测序深度；figures对应批次的_plot001等修订目录有全部组合效应/曲线/残差图。曲线设TS随机截距为零，非对随机效应积分的总体均值。单个失败项保留原因，P为NA；固定BH族不因失败缩小。

新BH在scope内合并两分量（16或28项）。同数据原显著筛选的影响保留，所选组合内BH不是整条筛选流程的FDR控制。观测检出同时受比例和读数深度影响，诊断需结合输入、TS覆盖和残差检查。

## 原版脚本准备阶段的验证记录

全部相关R语法通过；真实输入只读核对8/14组合；合成数据检出/Beta模型、不可估失败路径、完整汇总作图及检查点恢复通过。现有coral二维缓存重画分支通过；合成距离的3D HTML和六张投影通过。实际查看森林/NP测试图并修正标题越界。没有正式44模型拟合、PERMANOVA重跑或新正式3D计算。


## 当前命名与读取规则

本轮已接入共享run/plot编号与新旧latest兼容读取。talk可用`mc19.talk.output_dir="latest"`向同一完成批次补充图；两部分plot默认读最新完成结果。恢复拟合仍明确指定原目录并校验原代码。[完整规则与本轮只读验证](../README_run_names.md)。


## Two-part森林图样式与显著性标记

森林图基准字号22，轴文字20、横轴标题22、图例18，均加粗。MC编号旁`*`表示det的保存BH q<0.05，`#`表示pos的保存BH q<0.05，`*#`表示两部分均显著；不使用原始P或原总体GLMM q来标记。保留所有点和完整95%CI，同一MC的环境组合相邻。导出PDF、PNG与逐行标签核对表。

只更新森林图时设置`options(mc19.two_part.stage="plot",mc19.two_part.run_dir="",mc19.two_part.plot_parts="forest")`，再Source当前入口。默认plot_parts为`c("forest","curves","residual")`。新图保存在独立plot修订目录；forest-only修订不重画曲线或残差，原批次曲线和残差图继续保留。


## 2026-10-31：当前全部MC19图的统一字体与紧凑布局

当前8个物理输出目录已通过本地R脚本覆盖更新：41个PDF（含全部诊断页）、28个PNG、4个3D HTML。原GLMM曲线基准20pt且全粗体，保持全部原始点/完整CI与独立右侧Zoom。two-part曲线/残差与NMDS/NP采用约21pt基准及粗体轴/图例，增加边距以容纳大字体。热图＋森林图按筛选依据实际文字宽度安排左列，压缩空白，19行/多行依据保持对齐。

3D图保留原坐标，隐藏随视角变化的场景轴标题，改成图外固定颜色对应：X=NMDS1、Y=NMDS2、Z=NMDS3；颜色对应场景轴线及刻度。固定图外注释显示stress、PERMANOVA原始P的显著符号、BH q及ns。该批NO2/氨氮BH均未显著；coral氨氮是对照面板。

新增refresh_MC19_figure_style.R，仅重画当前保存结果：原GLMM直接用保存预测网格；NMDS采用cache_only=TRUE且不写回坐标/RDS/诊断表；诊断只读已有残差及模拟摘要。所有1083个数据/模型/NMDS缓存文件哈希校验未变。旧版图及QA在archive/style_refresh_run001，发布清单为publish_manifest.json。默认NMDS绘图复用已有缓存；只有显式mc19.nmds.refit_cached=TRUE才尝试继续确认已有解。


## 3D固定棱轴线与营养盐切换

`mc19_talk_3d_widget()`关闭原生自动选择显示棱的场景轴，用原坐标单位下固定的包围盒、三条轴线及刻度trace显示坐标；固定公共角为x最大/y最大/z最小。保留aspectmode=data。每个范围新增NMDS3D_nutrients.html，右侧NO2-N/Ammonia-N按钮只切换两条具有完全相同坐标和样本ID的点trace，并更新色标/统计说明，uirevision保留当前视角。原两个单指标文件同样包含切换按钮并覆盖更新，初始指标保持原文件含义。

浏览器核验408个样本在两指标下的坐标、ID、环境数值与保存数据逐项相等；三个相机方向下轴线/刻度数据不变，切换前后视角不变，统计说明对应原表。只重新导出HTML，未改静态PDF、NMDS缓存或P/q。旧版和检查记录在archive/3d_controls_run001。


3D散点微调：marker size由4改为3.3，opacity由0.8改为0.65；仅改变点的显示样式，两个范围的整合与兼容HTML均已覆盖。数据和布局逐项比较一致。


热图＋森林图现为16:9（32×18in）横向画布，24pt基准粗体；19MC及热图四指标完整保留。森林图仅保留至少一个MC的保存BH q<0.05的环境整列：coral为NO2/氨氮，ecosystem为NO2/氨氮/DIP；被保留的每列仍显示全部19MC及完整CI。其余森林面板的宽度重新分配；所有入口共用mc19_export_heatmap_forest，当前PDF/PNG已覆盖。


## 2026-11-01图形更新

已用本地R从保存模型和缓存优化并覆盖GLMM曲线、two-part森林/曲线、NP和交互3D，并同步MC19_阶段汇报结果图_v01。Two-part森林按det/pos分面，正值曲线显示百分比并在需要时附独立Zoom；原始点和完整CI保留。NO2-N/DIP检出限分别0.002 mg N/L、0.002 mg P/L（原分析L/2=.001），时间折线有水平虚线、环境响应曲线有垂直虚线。3D默认视角放大、固定轴与统计说明及营养盐切换保留，HTML应与旁边lib目录一起传送。既有refresh_MC19_figure_style.R继续读取保存缓存；本次执行和备份见研究根archive/figure_polish_run001/README.md。未重跑分析。
