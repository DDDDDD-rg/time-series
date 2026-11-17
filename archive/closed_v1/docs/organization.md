# Folder organization

更新：2026-11-17。

- data：原始资料及输入，原始来源文件名保留；derived放辅助派生工作簿。
- scripts：可运行入口、共享函数及必要的验证程序；GLMM/PERMANOVA目录保持RStudio调用结构。
- results：当前六指标单项、五指标联合、配套检查及仍被读取的历史基线。模型run编号与latest指针保持机器接口稳定。
- figures：六指标图包、可编辑图源及仍被调用的绘图缓存。
- presentations：只保留当前time_series.pptx。
- archive/superseded_results：旧温度增量批次及旧四营养盐联合模型。
- archive/superseded_figures：旧PPT包与温度五指标检查包。
- archive/execution：已结束的运行日志。
- 其余archive子目录：既有脚本备份、布局修订、审查与历史来源；部分脚本仍引用这些材料，保留内部结构。

命名：人工交付及独立检查采用简短英文、下划线和v1/v2编号，不再给新目录追加日期。既有模型ID、原始资料和程序指针不强行改名。ecosystem为当前范围名，历史all_hosts基线名称保留，避免混淆来源。根目录名称保留以兼容当前工程和外部引用。

路径对照和校验在[folder_cleanup_v1](../archive/folder_cleanup_v1)。未删除研究数据/模型/图像，未新建兼容快捷方式；移除了两个指向已归档旧图包的历史隐藏联接。其他119个联接中的117个保留并核对目标。归档RDS保留历史内部路径字段，用move_map.json追溯；恢复历史运行时应按路径对照读取，而不要直接运行内部保存的旧绝对路径。
