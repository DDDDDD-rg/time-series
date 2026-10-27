# 当前MC19 coral分析与作图

## 当前GLMM：19MC × 4营养盐，仅TS随机截距，Wald P/BH（233529已完成，2026-10-28核对）

当前主入口仅拟合 `cbind(MC_reads, Other_reads) ~ E_z + (1 | TS)`，共76个模型。仍为98个coral原始样本/7TS，保留TS4/TS6；每次仅一个营养盐，Beta-binomial/logit，恒定离散度，无额外零膨胀。无Host/Repeat/水温，不平均。

```r
source("C:/Users/xzy/Desktop/CUHK_Study/GNBF6010_Research_Projects/Research/Msc main project-time series/scripts/GLMM.R", encoding = "UTF-8")
```

默认仅复用120255先导中的4个有TS硝酸盐模型，经样本/计数/环境/公式严格核对后新增72次拟合。`reuse_pilot_fits=FALSE`则从头拟合76个。原先导共享模块保留历史两结构定义；当前拟合由新的TS专用`mc19_wald_fit_one()`完成。

直接提取glmmTMB的Wald beta P，全部76项一次性BH；失败项留NA，仍保持n=76。TS-only是用户查看172226两结构结果后确定的后续方案，历史152项q不改写、不直接复制进新结果。

`run_diagnostics`和`run_plots`可分别关闭。默认诊断每模型模拟2000次响应、不重新拟合；诊断P与beta关联P分开。主流程不执行bootstrap、置换或LOTO，也不按显著性、少数极端值或LOTO变号筛样本/MC。

新输出：`results/GLMM_MC19_env4_TS_Wald_coral_时间戳/`。顶层是76行`00_result.csv`、`00_readme.txt`和`supporting_files/`。主表列序：MC、环境、β、95%区间、P、q、边际预测占比变化、OR、SE、TS标准差、数值/诊断状态；样本数、模型标识、z值、删失数、来源等放后面。完整明细同样将这些重点列前置。

figures同名目录顶层5份合并PDF：总览2页，NO3/NO2/NH4/DIP各21页。诊断完整时，配套4份诊断PDF各39页。曲线每MC一整页，保留原始98个点、固定部分/边际预测、条件TS拟合值及固定部分Wald区间。

只重做诊断运行`GLMM_reads_diagnostics_full.R`；只重画图运行`plot_GLMM_MC19.R`。两者仅接收新TS-only批次，默认读取`00_latest_GLMM_MC19_env4_TS_Wald.txt`，当前指针已指向用户完成的233529批次。读取按MC/环境/model键对齐，主表允许显示舍入；全精度明细用于作图和诊断。

172226两结构Wald、120255先导、140106 LOTO、153110 bootstrap及其原模型/图/快照保留。历史bootstrap/LOTO仍使用各自命名入口。代码指南见[重要代码解读](../docs/GLMM_MC19_env4_Wald_code_guide_20261027.md)。用户已完成233529批次；本次发布核对31个R文件语法、保存的Wald/BH及90个来源/模型指纹，未重拟合或模拟。11份修改前源码/文档备份在`archive/script_backups/GLMM_before_TS_only_Wald_20261027/`。

GitHub保存当前源码、主表、全精度明细、环境TS/标准化摘要、诊断总表及9份合并PDF。原始工作簿/计数表、模型RDS、诊断模拟对象和历史批次留本地；GitHub材料用于审阅当前结果。在其它电脑重新拟合需准备相同输入、设置各入口project_dir；若没有本地120255先导，将GLMM.R的reuse_pilot_fits设为FALSE。单独重做诊断或重画还需本地保存模型及latest指针。

以下保留PERMANOVA及历史GLMM记录。历史GLMM章节中的analysis_mode、plot_mode和bootstrap参数属于旧入口版本，当前运行以上述Wald流程为准。

2026-10-27核对：v2分析与作图已于2026-10-26完成，当前运行编号20261026_224333。统计量来源为vegan；本轮整理目录与发布脚本，不重跑完整置换。

在RStudio中打开run_permanova.R，点击Source。不要依次运行所有编号脚本。

首次缺包时，在R控制台运行：
```r
install.packages(c("readxl", "vegan", "ggplot2"))
```

- 默认run_permanova.R先分析、再生成PNG和PDF；make_plots=FALSE时只分析。
- 只重新作图：运行04_plot_permanova.R，读取MC19专属最新完成结果，不重跑PERMANOVA或BH。
- 只重画宿主描述性高亮：运行07_plot_host_highlight_panels.R，复用04保存的二维坐标。
- 只重画7项环境组成图：运行plot_MC_environment_composition.R；旧plot_DO_composition_preview.R转到同一新入口。
- 如需交互3D，将run_permanova.R中的make_interactive_3d改为TRUE，另装plotly、htmlwidgets；默认关闭。
- 分析完成但图报错时，分析结果保留；排查后只运行04即可。

## 统计设置

- 输入：data/PERMANOVA_input_MC19_env7_20261026.xlsx；当前完整环境coral为98个样本、19个MC、7个有效TS。
- 每个原始样本保留一行；不使用Repeat_group、不平均、不加入Host_type，不执行宿主、TS、季节或两两检验。
- MC比例逐元素平方根后由vegan::vegdist计算Bray–Curtis，19个MC不重新归一化。
- 7项环境分别建立单因子模型。7个TS的环境行共用同一套5040种完整排列，包含实际分配；每次把TS环境值扩展到其原始样本，MC样本行固定。
- 每一种分配调用vegan::adonis2(distance ~ Env_value, permutations=0, by=NULL)，直接读取R²、F、平方和及自由度。脚本不再手写这些统计量的计算公式，也不使用adonis2的样本行置换P值。
- P值保持为全部TS排列中F不小于实际F的比例，含实际分配，保留重复环境值排列的重数；比较容差沿用v1。7项P值仍统一由stats::p.adjust(method="BH")校正。
- 盐度、水温、pH、硝酸盐氮、亚硝酸盐氮、氨氮、DIP清单及检测限L/2规则不变；TS8当期环境缺测仍排除。
- 每项环境会拟合5040次vegan模型，7项共35280次，因此会比v1矩阵实现慢；每完成500次显示进度。

## 结果目录

新分析运行后，每个results运行目录的顶层只有主表和一个配套文件夹：

```text
results/PERMANOVA_MC19_env7_coral_samples_时间戳/
├─ environment_results.csv          # 主要结果：7行，含R²、F、P、BH q和样本数
└─ supporting_files/
   ├─ environment_results_full.csv  # 完整统计量、模型与输入指纹等
   ├─ 00_summary.txt
   ├─ 00_run_status.txt
   ├─ 00_pipeline_status.txt        # 完整入口运行时写入
   ├─ 00_input_used_MC19_env7.xlsx
   ├─ coral_environment_analysis.rds
   ├─ environment_analysis_input.csv
   ├─ environment_TS_values.csv
   ├─ coral_sample_coverage.csv
   ├─ TS_label_orders.csv
   ├─ permutation_F_values.csv
   ├─ environment_correlations_TS.csv
   ├─ 00_analysis_version.txt
   ├─ 00_input_and_script_manifest.csv
   ├─ 00_session_info.txt
   ├─ scripts_used/
   └─ plot_NMDS2D.rds               # 作图时产生；可选3D另有缓存
```

主表保留Indicator、Factor、N_samples、N_MC、N_TS、R2_raw、F_raw、P_raw、P_BH_raw、Significance；平方和、自由度、排列数量等完整字段放在配套表与RDS。作图直接读完整RDS，不依赖主表删减的列。

新计算版本标记为MC19_env7_coral_individual_TSlabels_v2_vegan。作图入口同时兼容v1平铺目录及v2的supporting_files目录，缓存按其实际来源版本检查；不会把旧结果标记成新计算结果。当前结果保留原路径；MC8和MC19 v1历史结果/图已归档，原路径保留兼容联接；MC19专属latest指针继续指向同一运行根目录。

## 配套图与入口

图仍输出到figures/同一运行名/，不自动打开窗口。包括7项R²/q值、TS着色NMDS、7项环境着色NMDS、环境时间线、19MC原始样本点、MC×样本热图、宿主描述性高亮和7项环境排序组成柱。显著性读取保存的BH结果。

- 完整入口：run_permanova.R → 03_test_season_environment.R → 03_test_season_environment_batch_labels.R。
- 共用统计与路径：permanova_MC19_data_settings.R；共用作图读取器：permanova_MC19_plot_helpers.R。
- 旧01仅作MC8历史设置；02/05/06/08/09历史分支的停止保护保留。
- 本次修改前的脚本与导航备份：archive/script_backups/MC19_before_vegan_statistics_20261026/。
- 完整TS标签置换仍依赖时间标签可交换；本次只替换统计量计算实现与结果目录布局。

## 图目录与自动合并

```text
figures/PERMANOVA_MC19_env7_coral_samples_时间戳/
├─ 00_1_总体结果与NMDS.pdf         # 10页
├─ 00_2_环境与MC时间变化.pdf     # 3页
├─ 00_3_环境梯度下MC组成.pdf     # 7页
└─ supporting_files/            # 单图PDF/PNG、坐标、排序、来源清单、日志和脚本快照
```

04、07和环境组成图仍使用同一套读取器，单图统一写入supporting_files。每次作图完成后调用bundle_MC19_figures.py刷新合并PDF；只重画某个分支时，刷新具有完整来源的分组。原页拼接保留矢量内容，来源缺失时不覆盖已有合并文件。图latest指针始终保存运行根目录。

合并需要Python和pypdf。本机自动检测可用Python（包含已安装的Codex附带解释器）；也可设置`options(mc19.python="解释器完整路径")`或环境变量`MC19_PYTHON`。缺少工具会明确停止作图；修复后只运行04，不需重跑分析。可单独执行`python scripts/bundle_MC19_figures.py --figure-dir "图运行目录"`，兼容旧版平铺和新supporting_files来源。

当前已完成结果与本轮验证详见[更新日志](../CHANGELOG.md)。

## MC19 GLMM完整数据推断（2026-10-27修改，待用户运行）

当前`GLMM.R`默认`analysis_mode="full_data_inference"`，续接已完成120255的原8模型。主分析始终使用98个coral样本/7TS，含TS4、TS6；无水温、Host、Repeat项，不平均。当前仍为MC10、MC6、MC95、MC376与硝酸盐，各有/无TS随机截距，并列报告，不根据LOTO翻转、诊断P或最小P选择模型。

在RStudio运行：

```r
source("C:/Users/xzy/Desktop/CUHK_Study/GNBF6010_Research_Projects/Research/Msc main project-time series/scripts/GLMM.R", encoding = "UTF-8")
```

本模式只需glmmTMB，复用原8个完整估计和原诊断，不重新拟合原8模型，不重跑56次LOTO。每MC/结构的检验目标均为硝酸盐系数β=0：完整公式分别为`cbind(MC_reads, Other_reads) ~ E_z`和`~ E_z + (1 | TS)`；相应零模型只去掉`E_z`，其它结构不变。用原98行拟合零模型，再在零模型下模拟响应，每次按同一98行设计重拟合零/完整模型；TS、硝酸盐和每样本的总Ruegeria reads分母固定，有TS结构每次重新抽7个随机截距。所有似然均由glmmTMB计算，LRT为`2*(logLik(full)-logLik(null))`。

默认`bootstrap_B=999L`、固定seed；最多新增8个零模型和15984次bootstrap重拟合，明显比LOTO耗时长。每25次保存检查点并显示本模型预计剩余时间；不按显著性提前停止。中断后在GLMM.R中把`resume_results_dir`填为控制台打印的本次`GLMM_MC19_nitrate_full_data_...`完整路径，保持B/seed、原始结果、helper及R/glmmTMB版本不变，再Source。已完成模型只读检查点；每模型最多重算最近未保存的24次draw。不要通过重跑不同seed来选择较小P。

- `beta_P_value`/明细`P_value`：检验β=0的双侧关联参数bootstrap LRT P，`(1+不小于实际LRT的次数)/(B+1)`。这不是诊断P，也不是有/无TS结构之间的检验。
- `beta_q_value`/明细`q_value`：本轮预定4MC×硝酸盐×2结构共8项P统一BH，失败模型仍计入n=8；这是本阶段q，未来扩展19MC×4指标后需按届时预定的完整检验族重新校正。
- bootstrap重拟合失败不替换、不从分母删除：明细给`P_numerical_lower/upper`，主表使用将未知失败结果全部计为超越的保守上界并明确标记。零模型本身不可用时保留该行，P/q留空。
- `MC95_lower/upper`是有限模拟次数引入的二项区间；`CI95_low/high_approx`是原完整模型的Wald近似效应区间。零模型bootstrap只校准P，没有反演为效应CI。模型化校准仍依赖Beta-binomial及随机项假设，不能分离遗漏的季节混杂。
- 高硝酸盐信息集中于TS4=.098、TS6=.110 mg N/L，写入主表、说明与图。全部原始点保留；LOTO汇总在配套文件中提供参考，不参与P/q、样本或模型筛选。

```text
results/GLMM_MC19_nitrate_full_data_时间戳/
├─ 00_result.csv                    # 8行：完整beta、SE、近似CI、幅度、beta_P/q与状态
├─ 00_readme.txt
└─ supporting_files/
   ├─ inference_details.csv        # 完整数值、P上下界/模拟精度、失败与诊断状态
   ├─ LOTO_supplement_only.csv     # 仅作补充，不改原LOTO结果
   ├─ bootstrap/                  # 每模型零拟合、模拟计数、逐draw记录与恢复检查点
   └─ 冻结输入、来源指纹、源码、manifest与日志

figures/GLMM_MC19_nitrate_full_data_同一时间戳/
├─ 00_full_data_inference.pdf      # 7页：说明、完整效应、4MC曲线、beta P/q汇总
└─ supporting_files/              # 曲线坐标、主表快照、源码及原诊断/LOTO图路径
```

full_data的result/figure latest独立于pilot与LOTO；只有用户运行结束才创建，不把本次源码交付标成已完成结果。仅重画新图运行`plot_GLMM_MC19_full_data.R`，或默认模式的`plot_GLMM_MC19.R`，均不拟合、不模拟。原诊断和LOTO图继续读取既有对应目录。

本轮仅parse语法检查和静态核对，未拟合、模拟、校准或绘图。旧源码已哈希备份至`archive/script_backups/GLMM_before_full_data_inference_20261027/`。方法依据：[零模型参数bootstrap LRT](https://hojsgaard.github.io/pbkrtest/reference/pb-modcomp.html)、[glmmTMB模拟](https://glmmtmb.github.io/glmmTMB/reference/simulate.glmmTMB.html)、[R BH校正](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/p.adjust.html)。代码使用glmmTMB自定义同类流程，不调用不支持该模型类的PBmodcomp。

## MC19 GLMM硝酸盐先导（2026-10-27，用户已完成120255）

当前批次`GLMM_MC19_nitrate_pilot_20261027_120255`：8模型均收敛/Hessian通过，但均有诊断标记，P/q未校准；效应/曲线PDF为7页，诊断为17页。QQ诊断页重复标题的排版源码已修正；刷新原先导图时，将plot_GLMM_MC19.R的plot_mode设为"pilot"，不拟合、不模拟。新的完整数据推断不改写本批原输出。

本轮仅4MC（MC10、MC6、MC95、MC376）×硝酸盐×BB_GLM/BB_TS，共8个先导模型。两种结构分别为`cbind(MC_reads, Other_reads) ~ E_z`及`~ E_z + (1|TS)`；均使用glmmTMB Beta-binomial/logit、`dispformula=~1`、`ziformula=~0`。无Host/Repeat、reads样条或偏态扰动。全部生物样本逐行保留，不平均。

从本项目RStudio的Console执行（也可打开对应文件后Source整个文件）：

```r
source("C:/Users/xzy/Desktop/CUHK_Study/GNBF6010_Research_Projects/Research/Msc main project-time series/scripts/GLMM.R", encoding = "UTF-8")
```

以下为从头重做先导的说明：须先将GLMM.R的analysis_mode改为"fit_pilot"，才先8个拟合，再对数值可用模型各模拟2000次诊断，最后生成两份多页PDF。该模式依赖readxl、glmmTMB、DHARMa；不自动安装。本模式不做LOTO、bootstrap或76组合全量分析。当前续接完整结果应保留默认full_data_inference，避免重做先导。

输入由`PERMANOVA_input_MC19_env7_20261026.xlsx`的Metadata/Community/Environment_TS与`mc/MC_abundance_absolute.tsv`按New_ID对应。使用全部MC计数之和作Ruegeria分母，运行时核对19MC比例；仅coral、硝酸盐完整样本，当前预期98样本/7TS，变化时停止提示。TS8缺硝酸盐而排除。其余营养盐原值和检出限信息保存在输入快照中；本轮只拟合硝酸盐。E_z由7个有效TS等权计算中心与SD，样本不求均值。

```text
results/GLMM_MC19_nitrate_pilot_时间戳/
├─ 00_result.csv
├─ 00_readme.txt
└─ supporting_files/       # 模型、原始样本输入、TS环境值、诊断、比较表和快照

figures/GLMM_MC19_nitrate_pilot_同一时间戳/
├─ 00_1_effects_and_curves.pdf  # 7页：说明、8模型效应、4MC曲线、TS偏移
├─ 00_2_diagnostics.pdf        # 全部模型可诊断时17页；失败模型保留说明页
└─ supporting_files/          # 曲线坐标、结果快照、日志
```

每个MC并列展示两种结构；保留失败模型。数值检查包括优化器、Hessian、梯度、系数尺度及TS方差边界；标记需人工核查。`P_value`与`q_value`留空，`inference_status=pilot_not_calibrated`；Wald95区间是未作TS校准的先导近似区间。诊断P只用于模型失配筛查，不是硝酸盐关联P。AIC仅作结构描述，不自动择模型。

图中实线为b_TS=0固定部分、虚线为对估计TS随机截距分布积分的边际占比，叉号为含TS偏移的拟合值；这三者分开。阴影仅为固定部分近似Wald区间。GLMM系数为条件logit斜率；比较表另给共同硝酸盐范围的边际预测占比变化，不能把斜率变动直接等同于环境作用被“扣掉”。

仅重做诊断（读取已有8模型，不重新拟合，默认同步刷新图）：

```r
source("C:/Users/xzy/Desktop/CUHK_Study/GNBF6010_Research_Projects/Research/Msc main project-time series/scripts/GLMM_reads_diagnostics_full.R", encoding = "UTF-8")
```

仅重画图（不拟合、不模拟）：

```r
source("C:/Users/xzy/Desktop/CUHK_Study/GNBF6010_Research_Projects/Research/Msc main project-time series/scripts/plot_GLMM_MC19.R", encoding = "UTF-8")
```

两个配套入口顶部的`results_dir=""`只读取专用`results/00_latest_GLMM_MC19_pilot.txt`，也可填本轮某次运行的完整路径；不会退回旧MC8指针。主入口拟合后才创建新指针，本轮修改不会伪造已运行结果。旧GLMM与诊断源码已核对备份于`archive/script_backups/GLMM_before_MC19_pilot_20261027/`。08_plot_GLMM_host_beta_heatmap.R仅为历史MC8宿主图。

## 逐TS留一敏感性（用户已完成140106）

当前结果`GLMM_MC19_nitrate_LOTO_20261027_140106`：56次拟合全部收敛/Hessian通过，MC10/BB_TS/留TS5有一个TS方差接近0的标记。主表8行、明细64行、对照PDF5页。四MC的关系对TS4/TS6敏感，P/q仍未校准。面板脚注与横轴重叠的排版源码已修正，原PDF未重画；仅刷新图时运行plot_GLMM_MC19_LOTO.R，不要重新执行56次拟合。

已有120255八模型后，无需重跑GLMM.R。单独执行：

```r
source("C:/Users/xzy/Desktop/CUHK_Study/GNBF6010_Research_Projects/Research/Msc main project-time series/scripts/GLMM_leave_one_TS_out.R", encoding = "UTF-8")
```

该入口默认固定读取`results/GLMM_MC19_nitrate_pilot_20261027_120255/`，可在顶部修改`source_results_dir`。它只读复用原8个完整模型及其冻结样本，每次从同一份98样本开始，分别留出TS1–TS7，4MC×2结构×7轮新增56次拟合。不重拟合完整模型、不读取或修改当前原始Excel、不模拟、不置换、不bootstrap。只需要glmmTMB；与现有模型共用Beta-binomial及优化设置，无新分布、Host或Repeat项。

| 留出TS | 训练样本数 | 训练TS数 |
|---|---:|---:|
| TS1 | 81 | 6 |
| TS2 | 87 | 6 |
| TS3 | 90 | 6 |
| TS4 | 89 | 6 |
| TS5 | 75 | 6 |
| TS6 | 86 | 6 |
| TS7 | 80 | 6 |

每轮`droplevels`后以剩余6TS拟合；E_z沿用完整7TS中心和SD，不重新标准化。β变化在同一结构内部与完整模型比较。预测占比变化使用7个训练范围的交集（本批0.007–0.098 mg N/L），完整模型也在该范围重新计算点预测；不把训练范围不同造成的差异误认为关系改变。边际预测对各轮TS随机截距分布积分，另记录TS_SD变化。

```text
results/GLMM_MC19_nitrate_LOTO_时间戳/
├─ 00_result.csv              # 8行：每MC×结构的完整估计和7轮变化汇总
├─ 00_readme.txt
└─ supporting_files/
   ├─ all_estimates.csv       # 64行：原8估计 + 56留一估计，含近似区间
   ├─ fold_plan.csv           # 留出TS、样本数、环境范围
   ├─ frozen_scaling_and_range.csv
   ├─ frozen_full_input.csv   # 原完整样本快照，不是平均数据
   ├─ models/                # 56个留一模型，失败情况同样保留
   └─ 其余原模型引用/指纹、源码快照、清单和状态

figures/GLMM_MC19_nitrate_LOTO_同一时间戳/
├─ 00_LOTO_stability.pdf      # 5页：说明 + 每MC一页的两结构β/占比变化对照
└─ supporting_files/         # 图数据、结果快照及作图源码
```

主表重点字段：`n_usable_folds`、`unusable_TS`、`numerical_review_TS`、`sign_reversal_count/TS`、`beta_min/max`、`largest_abs_beta_change`、`largest_abs_change_in_full_SE`、`largest_shift_TS`、`SE_ratio_min/max`、共同范围的边际占比变化。`prediction_unavailable_TS`单独标记预测计算问题。全部失败/边界信息可追溯；不按“全部同号”或固定百分比门槛自动判断通过。原β接近0时结合绝对变化与区间解释，不采用不稳定的百分比变化。

图中菱形为完整模型，圆点为留一模型，括号为训练样本数；β横线为未作TS校准的近似Wald95区间，下面两个面板为共同范围边际占比变化的点估计。方向翻转不等于显著性改变；P/q仍留空。主分析仍使用完整数据，不选最好的一轮替换原结果。此检查不承担正式P值校准或独立验证。

新结果和图使用LOTO专用latest，原8模型及其latest不变。仅重画LOTO图时运行：

```r
source("C:/Users/xzy/Desktop/CUHK_Study/GNBF6010_Research_Projects/Research/Msc main project-time series/scripts/plot_GLMM_MC19_LOTO.R", encoding = "UTF-8")
```

作图入口只读LOTO已保存CSV与清单，不拟合、不模拟；顶部`loto_results_dir=""`读取LOTO latest，也可指定LOTO目录。两个入口共用`glmm_MC19_loto_helpers.R`及先导helper。修改前源码已核对备份至`archive/script_backups/GLMM_before_LOTO_20261027/`。源码交付时仅parse/静态检查；56次拟合与图已由用户完成，助手本轮仅读结果并修排版代码。

以下为脚本导航，包含MC8历史模块：

# R script entry points

MC19先导入口已改写；历史MC8结果和旧源码备份保留。

| 任务 | 入口 |
|---|---|
| 完整PERMANOVA分析 | [run_permanova.R](run_permanova.R) |
| MC19硝酸盐完整数据beta推断；可切换重做先导 | [GLMM.R](GLMM.R) |
| 已保存完整数据推断的7页合并图 | [plot_GLMM_MC19_full_data.R](plot_GLMM_MC19_full_data.R) |
| 已保存MC19模型的模拟诊断 | [GLMM_reads_diagnostics_full.R](GLMM_reads_diagnostics_full.R) |
| MC19合并图；默认full_data，可切换pilot | [plot_GLMM_MC19.R](plot_GLMM_MC19.R) |
| 原8模型的逐TS留一：新增56拟合 | [GLMM_leave_one_TS_out.R](GLMM_leave_one_TS_out.R) |
| 已保存LOTO结果的对照合并图 | [plot_GLMM_MC19_LOTO.R](plot_GLMM_MC19_LOTO.R) |
| PERMANOVA总体图与交互3D HTML | [04_plot_permanova.R](04_plot_permanova.R) |
| 两两比较热图 | [06_plot_pairwise_permanova.R](06_plot_pairwise_permanova.R) |
| 宿主高亮分面 | [07_plot_host_highlight_panels.R](07_plot_host_highlight_panels.R) |
| 历史MC8 GLMM宿主系数热图 | [08_plot_GLMM_host_beta_heatmap.R](08_plot_GLMM_host_beta_heatmap.R) |
| HOST/Seawater/Sediment图 | [09_plot_sample_type.R](09_plot_sample_type.R) |
| DO组成预览 | [plot_DO_composition_preview.R](plot_DO_composition_preview.R) |
| PM11环境时间序列 | [plot_PM11_timeseries_5years.R](plot_PM11_timeseries_5years.R) |

其余编号脚本和 `permanova_*`、`pairwise_letters.R` 为上述流程使用的单步分析或辅助函数。

RStudio项目 `scripts.Rproj`、`.RData`、`.Rhistory` 和 `.Rproj.user/` 保留原位置。
历史 `.bak`、before版本及备份目录集中在 [archive/script_backups/](../archive/script_backups/)。
