# MC19运行编号与最新结果

当前GLMM、PERMANOVA、两部分分析和talk作图共用`mc19_run_io.R`。现有时间戳目录及其中的模型、P/q、图和快照保留原样。历史MC8、先导、LOTO和bootstrap独立入口保留原有规则。

## 输出目录

| 内容 | 新目录示例 |
|---|---|
| 原GLMM结果 | `results/MC19_GLMM_coral_run001/` |
| PERMANOVA结果 | `results/MC19_PERMANOVA_ecosystem_run001/` |
| 两部分结果 | `results/MC19_two_part_coral_run001/` |
| talk图 | `figures/MC19_talk_coral_run001/` |
| 同批结果的第1版配套图 | `figures/MC19_GLMM_coral_run001_plot001/` |
| 同批结果的第2版配套图 | `figures/MC19_GLMM_coral_run001_plot002/` |
| 旧时间戳结果的重新作图 | `figures/MC19_GLMM_coral_legacy_plot001/` |

coral/ecosystem及每种分析各自编号。结果重跑增加`run`，只重画配套图增加`plot`；目录内PDF/PNG/CSV按内容命名。旧结果重新作图的来源完整路径记录在`run_record.csv`，`legacy_plot`不使用旧时间尾缀。

每批`run_record.csv`记录类型、范围、方法版本、来源、UTC时间及完成状态；原有统计manifest和来源指纹继续保留。`v1/v2`表示分析方案版本，重新运行同一方案不会自动提高方法版本。

项目根`.run_registry/`保存已分配编号，随输出归档仍保留，因此失败批次和已归档批次的编号不会重用；编号允许有空缺。该目录保持本地忽略。分配使用独立目录的排他创建，支持同时启动两个会话。

默认自动递增，通常不用设置编号。确需指定时，可在运行前设置一个高于已分配编号的值，使用后清除：

```r
options(mc19.run_id = "run010")
# source相应分析或talk入口
options(mc19.run_id = NULL)
```

## 最新结果如何读取

1. 显式填写的结果目录优先；检查分析类型、范围和完成状态，填写错误时直接报错。
2. 未指定目录时，同时检查新旧latest指针和当前results目录；只选择已完成的对应分析/范围。
3. 新规则按运行编号的数值排序，`run1000`晚于`run999`；没有新编号批次时，按旧目录的时间戳选择。
4. 指针缺失、指向不存在的目录、失败或未完成批次时，自动从其它已完成候选中恢复读取。读取不会重写指针或改动结果。
5. 自动扫描限于当前results目录；已归档的结果可通过仍有效的指针或显式完整路径读取。不会把整个archive当作新结果扫描。

新机器指针采用`00_latest_MC19_<类型>_<范围>.txt`；只有结果保存完整后才发布，同时更新原来的类型/范围专属指针。GLMM结果完成与模拟诊断完成分开记录，诊断状态继续保留在原配套文件中。两部分已完成但有需复核分模型的批次仍可读取。

## 继续使用原入口

原GLMM、PERMANOVA、诊断和talk入口不变；默认仍由`mc19.analysis_scope`选择范围。保留`mc19.glmm.result_dir`和`mc19.result_dir`指定历史来源的用法。

两部分只重画最新完成结果，不用再手填目录：

```r
options(mc19.analysis_scope = "coral",
        mc19.two_part.stage = "plot",
        mc19.two_part.run_dir = "")
source(file.path(project_dir, "scripts/GLMM/GLMM_MC19_two_part.R"), encoding="UTF-8")
```

两部分恢复拟合仍须显式指定原运行目录，沿用原输入、代码和软件版本校验。旧批次作图可直接读取；原代码快照保留在原批次配套文件中。

talk默认每次新建一个图批次。若先画森林/NP，再向同批图补充NMDS，可设：

```r
options(mc19.talk.output_dir = "latest", mc19.talk.parts = c("NMDS2D", "NMDS3D"))
# source plot_MC19_talk_figures.R
options(mc19.talk.output_dir = NULL)
```

这会复用当前范围最新完成的talk目录及兼容坐标；GLMM来源不同时要求显式指定其原始来源。也可直接填写已有图目录。常规配套图默认新建`plot`修订；高级选项`mc19.glmm.figure_dir`、`mc19.figure_dir`、`mc19.two_part.figure_dir`可显式指定图目录。

## 本次验证

仅执行R语法检查、临时文件的编号/完成状态/指针/范围/归档保留编号测试，以及现有六批结果的只读加载。两个原GLMM明细读取、两部分筛选输入和两个PERMANOVA绘图上下文均只读通过。未运行拟合、PERMANOVA、NMDS、模拟诊断或正式绘图。

回归检查入口：`test_mc19_run_io.R`。代码修改前备份：`archive/code_before_run_names_v1/`。
