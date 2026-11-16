# 联合GLMM：逐TS与极端样本敏感性入口

一键后台启动：运行同目录start_joint_TS7_sensitivity.ps1。默认2个并行模型任务，每个模型内部1线程，可用-Workers调整为1–8。启动前只读核对来源、原始计数和样本身份；启动日志/PID在archive/joint_TS7_sensitivity_launch/launchNNN。R入口默认只读check，命令行加--run才拟合。

范围固定为五环境联合run001（coral与ecosystem），保留NO3、NO2、总氨氮、DIP和温度，原七TS标准化值、原175行Ruegeria分母、全部原始零值和原模型分布/随机截距不变。不会重新读取变更后的工作簿：使用原运行冻结的输入快照，原拟合与结果只读。

- MC15、MC6、MC98：两范围各逐一留出TS1–7，共42模型，剩余TS因子drop未使用层级；检查每次设计秩/条件数。
- MC95、MC43：两范围逐样本核对冻结原始reads、175行分母、工作簿比例、Metadata样本ID/TS/样本类型；输出所有样本的审计表。按预设规则取比例最高3个与计数残差平方贡献最高3个的并集，分别单样本留一，再补一个并集联合留出压力测试。真实高值保留为主分析，这些留出不等于数据清洗。
- 新拟合使用原普通联合估计函数；记录数值失败、随机方差边界和完整系数。输出原系数、原q、差值/原SE、方向变化、系数范围等。敏感性Wald P仅作描述，不对选择性子集重新计算BH，不用各次是否显著作稳定标准。六TS六设计列的饱和状态单独标记。

结果独立存results/joint_TS7_sensitivity_runNNN，每任务一个RDS检查点，控制器每批更新00_task_progress.csv。01_all_coefficient_changes.csv含全部变化；02_original_significant_changes.csv只筛原模型q<.05的关联供优先阅读；03_LOTO_summary.csv汇总符号/幅度；04_extreme_sample_source_audit.csv是样本审计。数值失败/任务异常保持可见，不自动删点或更换模型。此入口不新增残差模拟、不作图、不修改PPT。

恢复中断任务：启动器加-Resume '该结果目录绝对路径'，会核对输入/源码指纹并跳过已保存任务（包含已记录失败项）。运行锁防止同目录并发；如果进程非正常结束，核对controller_pid.txt确认进程退出后，才手动移除该目录下.active_lock并恢复。完成结果不受恢复操作覆盖；如希望改变模型/任务规则，另起新运行。

RStudio只读检查可设options(mc19.joint.sensitivity.define_only=TRUE)，source入口后调用js_main(project,stage='check')；正式建议使用PS1启动器，关闭聊天不影响后台工作。
