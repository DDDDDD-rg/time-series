# MC19输出命名与读取

2026-11-17：当前已完成七TS六指标单项及五指标联合模型；八TS为待补输入方案。既有runNNN、plotNNN及注册表不重命名。独立检查与人工成果使用简短英文_v1/_v2；下方为原8TS命名规范，按相应入口区分。

更新：2026-11-02。结果仍用递增run编号，只重画用递增plot编号；目录名不追加时间。当前分析仅8TS、五营养盐，方法版本另存manifest。

| 内容 | 示例 |
|---|---|
| 原GLMM | results/MC19_GLMM_coral_run001/ |
| PERMANOVA | results/MC19_PERMANOVA_ecosystem_run001/ |
| 两部分 | results/MC19_two_part_coral_run001/ |
| talk图 | figures/MC19_talk_coral_run001/ |
| 配套图修订 | figures/MC19_GLMM_coral_run001_plot001/ |

各类型、范围独立编号，已分配ID由.run_registry保留，归档后也不复用。run_record.csv保存来源、方法版本、完成状态和时间。失败或不完整的批次不会发布为最新完成结果。

未指定目录时，读取latest指针并检查编号目录中的完成批次，按run数字选择最新。统计版本、范围、8TS和模型/环境数量必须匹配；旧7TS批次不会自动选中，显式指定也会提示不符合当前版本。

通用指针为`00_latest_MC19_<kind>_<scope>.txt`；辅助指针带`env5_TS8`（两部分带`TS8_env5`）。本次准备不改动任何已有结果指针，首次8TS新批次完成后自动发布。

旧7TS结果和图保留原目录供查阅，不增加旧数据兼容分支。检查程序`test_mc19_run_io.R`只操作临时测试文件，不运行研究分析。
