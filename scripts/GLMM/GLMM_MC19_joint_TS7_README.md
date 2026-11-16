# 七TS普通联合GLMM：运行与重要代码

发布核对2026-11-16：本五指标版本两范围38模型已完成，诊断/敏感性/先验预测亦有保存结果，具体状态见[当前审查](../../../../Research_Knowledge/questions/MC19_joint_TS7_results_review_v1.md)。下文“拟合数为0”是当时准备记录，非最新运行状态。

2026-11-16更新：当前版本同时放入NO3-N、NO2-N、总氨氮TAN、DIP和温度，不放DIN。已完成只读共线性和隔离验证，真实模型拟合数为0，正式fit由用户运行。旧四营养盐run001及其输入/源码快照保留。

主脚本：`GLMM_MC19_joint_TS7.R`。本轮仅准备普通最大似然版本；贝叶斯收缩版的先验尚未确定。本脚本是用户明确要求的七TS独立研究入口，不修改八TS主流程，也不把缺失TS8自动回退为七TS。

## 模型范围

- TS1–TS7，coral 98样本或ecosystem 310样本；每次运行一个范围。
- 19个MC全部纳入，各拟合一个模型。每模型同时估计NO3-N、NO2-N、总氨氮TAN、DIP、温度五个环境系数：共19模型、95个环境系数。
- Beta-binomial/logit，TS随机截距，常数离散参数；无额外零膨胀。
- 样本不平均，不加入Host/Repeat或DIN，不逐步筛选环境；仍为普通最大似然版，未设置斜率收缩。
- 全部零计数保留。分母使用原始计数表全部175行，包括Unassigned；19MC不重新归一化。

## 如何运行

在RStudio打开`../scripts.Rproj`，当前工作目录为scripts。

先做只读检查：

```r
options(mc19.analysis_scope = "coral", mc19.joint.stage = "check")
source("GLMM/GLMM_MC19_joint_TS7.R", encoding = "UTF-8")
```

检查会读取文件并计算标准化、相关性及VIF，不拟合、不创建结果目录。

准备好正式运行后：

```r
options(mc19.analysis_scope = "coral", mc19.joint.stage = "fit")
source("GLMM/GLMM_MC19_joint_TS7.R", encoding = "UTF-8")
```

另一个范围将`"coral"`改为`"ecosystem"`。两个范围分别运行、分别校正。脚本需要readxl；fit阶段另需glmmTMB，不自动安装软件包。

若要依次运行两范围，在同一scripts工作目录执行：

```r
for (scope in c("coral", "ecosystem")) {
  options(mc19.analysis_scope = scope, mc19.joint.stage = "fit")
  source("GLMM/GLMM_MC19_joint_TS7.R", encoding = "UTF-8")
}
```

## 重要代码的含义

### 1. 先筛七TS，再核对环境完整性

`env <- env[match(cfg$ts, env$TS), , drop=FALSE]`

输入仍读取当前`data/MC19_TS8_env5.xlsx`，但本轮明确只取TS1–7。随后只核对这七行原始报告值和实际水样月份，因此TS8缺测不阻止本轮，也不会被补零。Metadata与Community按New_ID核对，保留原样本。

### 2. 原始报告值和低于报告限的记录

`value[censored] <- value[censored]/2`

例如`<0.002`保留原始文本和删失标记，计算值沿用0.001；缺失、负值和无效格式会停止，不自动删样本。这是沿用原分析的L/2约定，并未估计检出限以下的真实浓度。

温度从Environment_TS的`SST_C_raw`读取，使用独立数值解析，不采用L/2。缺失、非有限数或<L温度直接停止。脚本逐TS匹配`data/water_quality/Bluff Island_PM11_marine_water_quality.csv`同一水样日期、PM11、Surface Water的唯一温度记录，要求与工作簿一致。七个值依次为25.8、31.4、28.6、16.7、24.6、26.6、25.2°C；它们是PM11水样温度。CSV、工作簿和计数表均记录指纹并保存输入快照。

### 3. 计数分母

`totals <- colSums(counts)`

counts包含原始175个MC行，故totals是每个样本的全部Ruegeria reads。随后用原始MC计数除以totals，与Community中的19MC比例逐项核对。

### 4. 七个TS等权标准化

`env[[cfg$z[j]]] <- (env[[field]]-mu)/s`

mu和s只从七个环境值求得，再映射到每个样本。TS5样本较多，不会因此在定义环境均值和SD时权重更高；生物计数仍是全部原始样本。

### 5. 共线性诊断

`X <- model.matrix(reformulate(cfg$z), data=env)`

在7行TS环境设计上检查截距加五个环境项是否满秩。VIF来自标准化环境相关矩阵的逆矩阵对角线。秩不足直接停止；VIF较高只提示，不自动删变量或按P选模型。

当前五指标输入只读核对：秩6/6；VIF为NO3=68.4221、NO2=8.5470、TAN=23.2286、DIP=99.9126、温度=1.6187，条件数24.2796。NO3–DIP相关0.975620，NO2–TAN相关0.890665；温度与四营养盐的相关依次为-0.554534、-0.149215、-0.219944、-0.492718。满秩表示形式上可估，营养盐各项信息仍高度重叠；只有7套环境组合估计截距和5个斜率，后续应重点检查逐TS稳定性。此处是输入检查，尚无加入温度后的联合拟合结果。

### 6. 模型公式

```r
cbind(MC_reads, Other_reads) ~
  Nitrate_z + Nitrite_z + Ammonia_z + DIP_z + Temperature_z + (1 | TS)
```

左边保留MC计数与其余Ruegeria计数。右边五个指标在同一模型，同时估计；各环境beta表示其他四项给定时，增加一个TS层SD对应的条件logit变化。`Temperature_z`由七个TS的`SST_C`等权标准化而来。`(1 | TS)`是TS共同背景的随机截距。

`family = glmmTMB::betabinomial(link="logit")`使用计数与分母建模；包括零计数。

`ziformula=~0`表示不另设结构零机制；`dispformula=~1`表示每个MC模型使用一个共同离散参数；`REML=FALSE`使用最大似然；`na.action=na.fail`防止拟合时悄悄删除缺失行。

### 7. 输出系数和数值检查

先检查优化器收敛、正定Hessian、有限的系数/标准误/似然以及合法的方差和离散参数，再给出Wald P及区间。数值失败项保留行和原因，不填造P值。

beta和SE在模型可提取时保留；95%区间为`beta ± qnorm(.975)*SE`。`exp(beta)`为MC相对比例的odds ratio，不是检出OR，也不是MC绝对数量倍数。

梯度>0.01、极大系数及接近零的TS方差另行标记，不把所有警告都自动视为失败。数值通过与残差分布适配是不同检查。

### 8. BH检验族

`p.adjust(results$P_value, method="BH", n=n_coefficients)`，其中`n_coefficients=19*5=95`。

同一范围19MC×5个联合条件系数统一校正，不包括截距。失败项保留NA，计划检验族保持95。不会复用旧四营养盐联合模型的系数、P或q。普通版P和CI是七TS下的探索性Wald近似，未做bootstrap、置换或跨TS验证；本脚本不按P或q给环境做因果重要性排序。

## 正式运行后会生成什么

自动创建`results/joint_glmm_TS7_env5_coral_run001/`或`results/joint_glmm_TS7_env5_ecosystem_run001/`，已存在时递增编号；与旧四营养盐`joint_glmm_TS7_{scope}_run001`区分，不覆盖、不发布到旧八TS结果指针。

- `00_coefficients.csv`：95行联合条件系数、SE、95%CI、OR、P/q及数值状态。
- `01_model_checks.csv`：19个模型的收敛、Hessian、TS方差、离散参数、警告和失败原因。
- `supporting_files/models/`：每MC的模型RDS。
- `supporting_files/`：样本和环境输入、标准化、VIF、相关矩阵、逐TS温度来源、输入与脚本快照、运行环境及来源指纹。
- `run_status.txt`：未完成、部分模型待复核，或拟合结束但残差诊断待做。

本脚本不自动开展DHARMa模拟、逐TS留一、重要性删项分析或作图；保存模型后可针对结果继续设计诊断。

## 准备阶段验证

```r
source("GLMM/test_GLMM_MC19_joint_TS7.R", encoding="UTF-8")
```

测试在只读check模式核对两个范围，逐样本比对旧七TS模型的原始计数、分母和四项营养盐，温度与已完成的六指标单环境流程及原CSV一致；检查温度无删失/数值解析、五项等权标准化、满秩和95项输出结构。另用临时CSV故意改变温度，确认check拒绝来源不一致。拟合函数被临时替换为哨兵并在退出时恢复；真实模型拟合数为0，未创建正式结果目录。

若在Windows命令行遇到已知的R启动`compiler namespace`报错，当前本机验证使用的进程级设置是`LC_ALL=English_United States.utf8`、`LANG=en_US.UTF-8`，并用`Rscript --vanilla --encoding=UTF-8`。未修改系统或RStudio配置。

方法依据：[glmmTMB](https://glmmtmb.github.io/glmmTMB/reference/glmmTMB.html)、[数值问题说明](https://glmmtmb.github.io/glmmTMB/articles/troubleshooting.html)、[BH实现](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/p.adjust.html)。
