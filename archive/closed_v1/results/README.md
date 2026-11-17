# Analysis results

更新：2026-11-17。以下为当前七TS输出。

| 分析 | Coral | Ecosystem |
|---|---|---|
| 六指标单项GLMM | [结果](MC19_GLMM_TS7_env6_coral_run001/00_result.csv) | [结果](MC19_GLMM_TS7_env6_ecosystem_run001/00_result.csv) |
| 六指标PERMANOVA | [目录](MC19_PERMANOVA_TS7_env6_coral_run001) | [目录](MC19_PERMANOVA_TS7_env6_ecosystem_run001) |
| two-part | [目录](MC19_two_part_TS7_env6_coral_run001) | [目录](MC19_two_part_TS7_env6_ecosystem_run001) |
| 五指标联合GLMM | [目录](joint_glmm_TS7_env5_coral_run001) | [目录](joint_glmm_TS7_env5_ecosystem_run001) |

- [联合残差诊断](joint_glmm_TS7_env5_diagnostics_run001)
- [联合敏感性检查](joint_TS7_sensitivity_run001)
- [先验预测检查](joint_TS7_prior_predictive_run001)：待进一步校准，尚未拟合收缩后验。
- [二维NMDS稳定性](nmds_2d_stability_v1)、[ecosystem三维稳定性](nmds_3d_ecosystem_stability_v1)
- [reads核对](reads_audit_v1)

旧温度增量批次、旧四营养盐联合模型已归入[历史结果](../archive/superseded_results)。运行日志在[archive/execution](../archive/execution)。

保留的glmm_coral_v1、glmm_all_hosts_v1、two_part_coral_v1、two_part_all_hosts_v1及permanova两个v1目录属于历史基线依赖，增量检查及已有图源读取仍使用它们，不能当作当前六指标结果。00_latest指针供程序定位批次，其既有文件名保持稳定。
