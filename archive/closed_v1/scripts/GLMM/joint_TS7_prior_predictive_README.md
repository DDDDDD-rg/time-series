# 七TS联合GLMM：先验预测检查

运行入口（PowerShell，从研究目录）：

```powershell
& './scripts/GLMM/start_joint_TS7_prior_predictive.ps1'
```

先检查Beta-binomial生成器的理论矩和当前原始输入，再隐藏窗口启动独立R进程。
默认每项4000次独立先验抽样；19MC×2范围×3档斜率先验，共114项。无需Stan或brms；使用base R生成器，输入核查复用原脚本及readxl。不调用原脚本拟合入口。

## 候选模型

五个环境斜率独立Normal(0,tau²)，tau=0.5/1/2；alpha Normal(-4,2.5²)；sigma_TS Half-Normal(0,1²)；log(phi) Normal(log(50),2²)。
标准化沿用七TS等权的原输入，截距在环境z=0处定义；每次抽样为每个TS生成一个共享随机截距。按实际样本N生成 q~Beta(mu*phi,(1-mu)*phi)，再Y~Binomial(N,q)。对极端浮点概率做machine-epsilon边界保护并记录次数。

这是待校准的候选先验，不利用MC响应更新参数。观察值只用于叠加比较，不决定抽样。三档tau共用基础参数随机数，MC与范围间种子不同；各MC是独立one-versus-rest预测，不能叠成总和为1的组成图。ecosystem包含coral，范围比较不构成独立重复。

## 输出

新建`results/joint_TS7_prior_predictive_runNNN`，旧结果保留。最新路径在`results/00_latest_joint_TS7_prior_predictive.txt`。

- `run_status.txt`：running/complete/failed；`01_progress.csv`：逐项完成记录。
- `02_prior_predictive_summary.csv`：样本平均比例、TS等权平均、零值率、最大比例、TS均值范围及1%/10%/50%超过率的先验分布与观测对照。
- `03_prior_predictive_by_TS.csv`：逐TS的比例、零值、高值及平均计数。
- `tasks/*.rds`：参数抽样、潜在均值与预测统计量检查点；不保存巨大的逐样本完整抽样阵列。
- `figures/*.pdf`：每MC/范围一页、三档tau对比，原始期望比例尺度及预测零值/最大值；红线是观测参考。不是PPT结果图。
- `supporting_files`：输入上下文、指纹、标准化、环境、源码及R会话。

所有95%区间为先验/先验预测区间，不是后验可信区间；观测越界及阈值只是描述，不作通过/拒绝或自动选先验的规则。脚本不生成P/q，不拟合收缩模型，不修改数据与PPT。中断后保留检查点；再次运行创建新目录，不隐式续跑旧任务。

方法依据见知识库Gelman 2008、Wesner/Pomeranz 2021、Gabry 2019和联合GLMM问题记录；本脚本为项目实现，未宣称复现原文实验。
