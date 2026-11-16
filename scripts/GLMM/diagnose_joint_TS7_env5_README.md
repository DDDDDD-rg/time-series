# 七TS五指标联合GLMM：保存模型的残差诊断

入口为diagnose_joint_TS7_env5.R；默认check，仅核对输入。正式诊断读取五指标联合coral/ecosystem的run001，各19模型，不读取四营养盐旧模型、不重新拟合、不重算环境关联P/q。两范围一起处理，原结果只读。

## 运行

在scripts.Rproj对应的scripts目录运行：

```r
options(mc19.joint.diagnostics.stage="diagnostics",
        mc19.joint.diagnostics.n_sim=1000L,
        mc19.joint.diagnostics.seed=81270L)
source("GLMM/diagnose_joint_TS7_env5.R",encoding="UTF-8")
```

需已安装glmmTMB和DHARMa。n_sim默认1000、至少250，可改2000；两范围/各MC使用固定且不同种子。每个MC沿用真实环境、TS和原始Ruegeria reads分母，用glmmTMB::simulate生成计数；TS随机截距重新抽取，为unconditional模拟，无refit。该设置显式指定，不依赖DHARMa版本的默认模拟层级。

输出独立results/joint_glmm_TS7_env5_diagnostics_runNNN，每次递增且不覆盖。00_all_38_models.csv汇总38模型的uniformity/dispersion/zero-inflation诊断原始P、状态与错误；两范围子目录保存每MC测试、逐TS/总体零值与正值中位数/尾部的模拟95%区间、每样本残差和完整RDS。模拟摘要区间为描述性逐项包络，不是多重校正后的联合区间。

review_flags表示任一诊断P<.05、任一摘要超出模拟区间，或单项检验失败，供人工检查；不自动否定环境关联。诊断P不与95项环境关联q混合校正。单模型诊断失败会记录原因并继续，其余模型保留。开始时验证38份模型，正式执行时另设防重拟合哨兵；结束复核模型文件MD5。来源、种子、版本、源码快照随输出保存。

## 图与结果分开

先完成diagnostics；再将实际返回的新目录赋给run_dir：

```r
options(mc19.joint.diagnostics.stage="plots",
        mc19.joint.diagnostics.run_dir="填入刚生成的诊断目录绝对路径")
source("GLMM/diagnose_joint_TS7_env5.R",encoding="UTF-8")
```

plots只读取保存诊断，每范围一个diagnostics.pdf，包括总体残差图和残差对五项环境的图，不再模拟。PDF存在时停止以避免覆盖。逐TS比较表在结果阶段已生成；诊断不包含逐TS留一重拟合，也不处理共线性本身。

本次仅准备源码和语法检查，未运行38模型诊断、模拟或作图。
