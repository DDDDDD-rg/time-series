# MC19 × 4营养盐：TS随机截距Wald流程与代码解读

更新：2026-10-27。用户查看172226两结构结果后，确认后续所有当前模型统一采用TS随机截距。用户随后已完成233529的76模型批次；2026-10-28发布时核对Wald/BH、全部90个文件指纹及9份PDF页数，未重新拟合或模拟。172226的152项P/q和原图保持原样。

## 1. 从哪里运行

在RStudio打开scripts/GLMM.R并点击Source，或运行：

```r
source("C:/Users/xzy/Desktop/CUHK_Study/GNBF6010_Research_Projects/Research/Msc main project-time series/scripts/GLMM.R", encoding = "UTF-8")
```

主入口依次分析、诊断、作图。只重做诊断用GLMM_reads_diagnostics_full.R，只重画图用plot_GLMM_MC19.R。辅助文件只定义函数，不需要逐个运行。

| 文件 | 当前用途 |
|---|---|
| GLMM.R | 路径、复用/诊断/作图开关与总入口 |
| glmm_MC19_wald_helpers.R | 76组合、TS专用拟合、Wald/BH、结果保存与读取、诊断 |
| glmm_MC19_wald_plot_helpers.R | TS-only总览、森林图、单MC曲线及诊断PDF |
| glmm_MC19_pilot_helpers.R | 复用输入验证、提取和诊断基础函数；保留历史两结构定义 |

## 2. 研究范围与设置

```r
cfg$model_ids <- "BB_TS"
cfg$n_models <- cfg$q_family_size <- 76L
```

全部19MC分别与硝酸盐、亚硝酸盐、氨氮、DIP拟合，每个组合一个模型。98个coral原始样本，TS1–TS7含TS4/TS6；无Host、Repeat、水温或平均。不因q值、极端点或LOTO变号筛除样本/MC。

`reuse_pilot_fits=TRUE`时，核对120255先导中4个有TS硝酸盐模型，再新增72次拟合；FALSE时重新拟合76个。复用检查样本ID、TS、MC/其它/总reads、E_raw、E_z、分布和实际均值公式。旧先导无TS模型不会进入新任务。

`run_diagnostics`与`run_plots`控制后续步骤。`n_sim=2000L`只用于分布诊断，beta的P不使用这些模拟。主流程没有bootstrap、置换或LOTO。

## 3. 输入与标准化

读取data/PERMANOVA_input_MC19_env7_20261026.xlsx和data/mc/MC_abundance_absolute.tsv，按New_ID匹配。输入验证覆盖全部19MC比例、原始营养盐/检出限、98样本和七个TS。

```r
Other_reads <- Ruegeria_reads_total - MC_reads
MC_proportion <- MC_reads / Ruegeria_reads_total
center <- mean(values)
scale <- sd(values)
E_z <- (E_raw - center) / scale
```

分母加总原始表全部MC行（175行），不按19MC重新归一化。`values`仅含当前营养盐在七个TS的七个值，各TS等权计算中心/SD。98个生物样本仍逐行进入模型。7448行长表是19×4×98的分组展开，不增加独立样本数量。

输入中的MC10模板只用于取得共有样本元数据，不按MC10丰度筛选。低于检出限仍按L/2处理并保留标记：NO2四个TS，DIP三个TS。

## 4. 唯一的拟合方程

当前`mc19_wald_fit_one()`直接固定以下公式，不再调用历史两结构拟合器：

```r
formula <- cbind(MC_reads, Other_reads) ~ E_z + (1 | TS)
fit <- glmmTMB::glmmTMB(
  formula = formula,
  data = d,
  family = glmmTMB::betabinomial(link = "logit"),
  ziformula = ~0,
  dispformula = ~1,
  REML = FALSE
)
```

`cbind`两列为该MC与其它Ruegeria的reads。Beta-binomial描述计数比例及额外离散；logit为占比的链接函数；E_z系数为环境增加一个TS间标准差的logit变化；`(1|TS)`为同次采样共享的背景偏移。每模型一个恒定离散参数，没有额外零膨胀项。优化控制沿用单线程、iter.max/eval.max均10000。

TS随机截距方差接近零时显示边界审查标记，不按此自动删除模型。其假设及只有七次环境采样的信息限制仍需结合结果解释。

## 5. Wald beta检验与76项BH

```r
tab <- summary(entry$fit)$coefficients$cond
values <- unname(tab["E_z", c("Estimate", "Std. Error", "z value", "Pr(>|z|)")])
```

四个值依次为beta、SE、z和原始双侧Wald P，直接来自glmmTMB。CI95为beta ± qnorm(.975)×SE，OR为exp(beta)。数值不可用时保留该组合，P/q留NA。

```r
valid <- is.finite(results$beta_P_value) & results$numerical_ok %in% TRUE
results$beta_q_value[valid] <- stats::p.adjust(
  results$beta_P_value[valid], method = "BH", n = cfg$q_family_size
)
```

`cfg$q_family_size=76L`，全部76项只校正一次。`valid`不按P大小筛选；失败项仍计入计划族76。TS-only是用户查看152项结果后确定的后续分析范围；本轮未计算新q，旧批次152项q仍保留，不能直接拼接成76项q。

## 6. CSV重点参数前置

两个CSV共享以下22列前缀，完整明细再接其它技术字段：

```text
MC_ID, environment,
beta_per_TS_SD, CI95_low, CI95_high,
beta_P_value, beta_q_value,
marginal_change_pp, odds_ratio_per_TS_SD,
SE, TS_SD_logit,
numerical_ok, diagnostic_status,
n_samples, n_TS, model, z_value, censored_TS,
fit_status, inference_status, q_family_size, fit_origin
```

`marginal_change_pp`是营养盐从最低到最高观测值时，模型边际预测占比的百分点变化；不是相对百分比变化。含TS的边际预测对正态随机截距分布进行积分。模型/样本标识、方法状态和来源等信息完整保留，重点效应与证据先显示。

读取时先按MC_ID/environment/model键对应，主表允许1e-7绝对误差的显示舍入；诊断和作图始终使用全精度all_results.csv。删除/重复模型键或超出舍入范围的改值会提示检查。

## 7. 诊断与图

每个保存的有TS模型模拟2000次响应，重抽TS随机截距而不重新拟合。DHARMa与模拟摘要评价总体/分TS零值、中心、q95和最大值等；诊断P与环境beta P用途分开。review_flags用于人工审查，不自动否定模型或要求全部诊断P>.05。

曲线每MC一整页：原始98个样本点，实线为b_TS=0的固定部分，虚线为对TS分布积分的边际均值，叉号为条件TS拟合值，灰带为固定部分Wald均值区间。总览q采用三位有效数字，避免极小q显示为0.000。

## 8. 输出位置和独立指针

```text
results/GLMM_MC19_env4_TS_Wald_coral_时间戳/
├─ 00_result.csv                 # 76行，重点列前置
├─ 00_readme.txt
└─ supporting_files/
   ├─ all_results.csv            # 同样前置重点列，后接完整元数据
   ├─ analysis_input.csv         # 7448行长表
   ├─ input_scaling.csv
   ├─ models/                    # 76个记录，默认4复用+72新增
   ├─ diagnostics/
   ├─ scripts_used/
   └─ wald_manifest.rds

figures/GLMM_MC19_env4_TS_Wald_coral_同一时间戳/
├─ 00_summary.pdf                # 2页
├─ 01_NO3_effects_and_curves.pdf  # 21页
├─ 02_NO2_effects_and_curves.pdf  # 21页
├─ 03_NH4_effects_and_curves.pdf  # 21页
├─ 04_DIP_effects_and_curves.pdf  # 21页
└─ supporting_files/diagnostics/ # 诊断完整时四份PDF各39页
```

新机器指针为results和figures下的`00_latest_GLMM_MC19_env4_TS_Wald.txt`，拟合完成且指纹校验通过后才创建。当前配套入口严格读取新版本76个TS模型，不接旧152项批次。历史原图和目录不修改；旧源码在各批scripts_used以及本轮备份中保留。

本轮11份修改前源码/文档备份：archive/script_backups/GLMM_before_TS_only_Wald_20261027/。仅parse与静态审计，不source/eval生产分析，不执行拟合、P提取、BH、响应模拟或正式图生成。
