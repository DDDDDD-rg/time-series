# Time-series research project

更新：2026-11-17。当前已完成TS1–TS7；TS8环境资料待补。当前单指标分析使用六指标，联合GLMM使用五指标，两条流程分别保存。

| 内容 | 入口 |
|---|---|
| 当前PPT | [time_series.pptx](presentations/time_series.pptx) |
| 当前模型及检验 | [results/README.md](results/README.md) |
| 当前图、交互图与绘图来源 | [figures/README.md](figures/README.md) |
| 脚本与运行次序 | [scripts/README.md](scripts/README.md) |
| 数据与TS8状态 | [data/00_CURRENT_DATA.md](data/00_CURRENT_DATA.md) |
| 归档、命名与路径对应 | [docs/organization.md](docs/organization.md) |

七TS：coral 98、ecosystem 310个样本。单指标GLMM包括NO3、NO2、总氨氮、DIN、DIP、温度，每范围114项，BH显著10/18项；PERMANOVA六指标均未过BH；two-part跟进显著组合。五指标联合GLMM排除DIN，包含温度；诊断、逐TS敏感性与先验预测已完成，收缩后验拟合尚未开展。

NMDS：coral二维稳定；ecosystem二维整体相似但少数点不稳定，三维已通过重复布局检查（stress约0.156899）。本轮只整理文件和导航，没有重跑模型或修改PPT内容。
