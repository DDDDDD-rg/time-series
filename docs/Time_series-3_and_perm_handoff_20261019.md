# Time series-3 与 Time series-3 perm 综合交接

整理日期：2026-10-19（Asia/Shanghai）。本文件用于 Time series-4 及后续对话接续。

后续研究重点已更新，先读 [2026-10-19 研究方向调整记录](Research_direction_update_20261019.md)：扩展 MC，先画时间趋势；主范围改为仅珊瑚和全部生物 hosts、排除海水沉积物；营养盐拆分，评估历史同月均值及不调整宿主的路线。下文保留此前已完成分析的实际状态。

## 1. 接续时先读

- 研究目标：解释 Ruegeria 内候选 MC 的组成与宿主、采样批次、季节及环境梯度的关联。群落层面使用 PERMANOVA，单个 MC 使用共同结构的 GLMM，各 MC 分别估参。
- 当前 8 个 MC：MC0、MC6、MC10、MC15、MC21、MC25、MC507、MC520。
- 当前完整 GLMM 结果为 `results/GLMM_Z_TS_SD_Q_MC8_20261013_185720/`；当前完整拟合脚本为 `scripts/GLMM/GLMM.R`。
- 当前完整 PERMANOVA 结果为 `results/PERMANOVA_MC8_20261014_183535/`，已包括 Seawater、Sediment；入口为 `scripts/PERMANOVA/run_permanova.R`。
- 当前 GLMM：7 个 MC 数值有效；MC15 的 Hessian 检查失败，不能报告有效 SE、CI、p、q。DO 对 MC6、MC507 为正关联，对 MC21 为负关联，均通过计划 8 项的 BH 校正。
- 当前 PERMANOVA：宿主类别与 TS 均显著；干湿季未检出显著差异；8 项环境筛查中仅 DO 在全样本与仅珊瑚范围均通过 BH 校正。
- 最近讨论转到环境 Excel、2015–2025 年 PM11 表层历史数据，以及如何按独立、配对、年份分块的设计选择检验。尚未对这份 11 年环境表开展新的显著性检验。
- 用户最近授权：将 Time series-3 与 Time series-3 perm 一并整理、保存、读取。当前不需要继续改分析脚本、重新拟合模型或运行新检验。
- 2026-10-12 的旧交接是历史起点；其中“惩罚尚未实施”“GLMM.R 还是旧模型”等状态，已被本文件第 4、5 节记录的后续工作取代。

## 2. 来源与核对范围

| 对话 | 对话 ID | 本次保存的文字 |
|---|---|---|
| Time series-3 | `01a12e8a-d27b-7010-bee3-dd29e9d0f5d1` | 110 条用户消息、84 条最终答复 |
| Time series-3 perm | `01a130ff-ec29-77d3-9f34-a81a3e7604cf` | 58 条用户消息、47 条最终答复 |

本次先通过对话工具确认标题与近期内容，再从两份本地会话记录抽取用户消息和最终答复。文字存档不包含内部推理、工具调用过程或图片二进制；原消息中的图片、附件路径引用保留。存档保留当时的表述，后续修订后的状态在本综合交接中说明。

本次重新读取了项目首页、结果导航、当前 GLMM 与 PERMANOVA 关键输出、分范围原始 CSV、相关脚本、PERMANOVA 方法说明，以及历史环境数据 audit JSON。核实了下文关键文件存在。Excel/PPT/Word 的内容细节来自原对话及其已记录的核对，本次未逐页重审这些文件；后续编辑前再读取实际版本。

完整文字存档与 manifest 位于研究根目录之外、项目工作区内：

```text
C:/Users/xzy/Desktop/CUHK_Study/GNBF6010_Research_Projects/Codex_Workspace/outputs/01a15454-d30e-77a3-a7aa-702e4b1dcd5f/context_handoff/
```

- [Time series-3 文字存档](../../../Codex_Workspace/outputs/01a15454-d30e-77a3-a7aa-702e4b1dcd5f/context_handoff/Time_series-3_transcript.md)
- [Time series-3 perm 文字存档](../../../Codex_Workspace/outputs/01a15454-d30e-77a3-a7aa-702e4b1dcd5f/context_handoff/Time_series-3_perm_transcript.md)
- [来源、消息数量与 SHA-256 清单](../../../Codex_Workspace/outputs/01a15454-d30e-77a3-a7aa-702e4b1dcd5f/context_handoff/manifest.json)

## 3. 项目位置与工作约定

研究根目录：

```text
C:/Users/xzy/Desktop/CUHK_Study/GNBF6010_Research_Projects/Research/Msc main project-time series
```

下文相对路径均相对于此研究根目录。项目首页为 [README.md](../README.md)。旧 `C:/Users/xzy/Documents/CUHK Project` 是历史兼容入口，正式研究数据以当前根目录为准。

- 中文交流；讲统计和代码时先给核心结论，再用直观例子、公式和真实结果解释。合适时配示意图，标明是假设参数还是实际拟合，并核对坐标轴与分布含义。
- 表格默认白底黑字，必要加粗和细线，不加花哨底色。模型结果表先放方程及各 MC 系数，检测放第二个工作表，解释与来源放后面。
- 图用英文、无图题，白底、高对比配色；字体、数字、显著星号要大且不拥挤。放大字体时固定画布，避免同步放大画布抵消视觉效果。
- 用户指定“只改脚本，不运行／我来运行”时，只编辑和静态核查。已授权的常规工作无需重复确认。
- 修改原始资料前保留备份并验证。不要为整理随意改变样本 ID、输入定义或分析范围。
- 不自动打开文件、图窗或成果侧栏；聊天内交付文件链接。不要自行创建、分叉或切换其他对话。
- 区分已实现结构、已保存结果、历史比较和未实施建议。当前结果以真实脚本、输出及运行日期核对。

## 4. GLMM：数据和实际模型

### 4.1 数据与分析单位

- 输入：[data/GLMM_input_MC8_20261011.xlsx](../data/GLMM_input_MC8_20261011.xlsx)，`Input_long`。
- 每个 MC 源记录 400 条，实际拟合 348 条、TS1–TS7；TS8 的 52 条记录缺 DO，按统一规则排除。
- 拟合单位是原始样本记录，不是 Repeat_group 均值；Repeat_group 随机截距描述组内共享偏移。历史交接记录当前拟合有 103 个 Repeat_group。
- `Y = MC_reads`，`N = Ruegeria_reads_total`，`Other_reads = N − Y`；二项响应为 `cbind(MC_reads, Other_reads)`。
- 分母为样本全部 Ruegeria reads。全部观测零值保留，不加伪计数，不按所选 8 个 MC 的总和重新归一化。
- 宿主参考组为 `Acropora sp.`。GLMM 保留的宿主范围不同于 PERMANOVA 的至少两个重复组筛选。
- `DO_c = DO_mg_L − 6.3`，单位仍为 mg/L；`x = log10(N) − c`，仅中心化，reads 增加 10 倍时 x 增加 1。各 MC 共用样本、中心和样条设计。

### 4.2 Z_TS_SD_Q 的共同结构

```text
Y_i ~ pi_i * point_mass_at_0 + (1 - pi_i) * Binomial(N_i, p_i)

logit(p_i) = beta0 + beta_DO * DO_c_i + beta_Host[i]
             + ns(x_i, df=3) + b_Repeat[i] + b_TS[i] + u_i

log sigma_i = a0 + a1*x_i + a2*x_i^2

u_i：均值为0的中心化偏正态扰动；标准差为sigma_i；偏斜参数alpha每MC为常数。

logit(pi_i) = gamma0 + gamma_reads*x_i + z_TS[i]
z_t ~ Normal(0, 1.5^2)
零值reads斜率惩罚：lambda * gamma_reads^2；lambda = 0.125
```

Repeat 和计数部分的 TS 随机截距各自的 SD 从数据估计。零分量 TS 偏移采用固定 1.5 尺度的正则；这与计数部分 TS 随机截距的作用位置不同。

零膨胀成分以概率 pi 直接生成零，普通二项成分也可以生成零。本次 MC6、MC10 实际采用 `exact_pi0` 边界，其余 MC 的保存模式为 `reads_dependent_pi`。因此描述共同框架时，应保留实际的边界处理。

偏正态描述固定 reads 位置下 u 的正负不对称；均值为 0 通过中心化实现。`ns(x, df=3)` 描述中心趋势，sigma 的二次项描述波动宽度，alpha 描述偏斜。实际 7 个有效模型的 sigma 点估计均呈 U 形，最小宽度位于约 6,300–15,100 reads；不能用早期任意参数示意图的“中间最宽”代替真实结果。

### 4.3 当前运行和推断

- 完整拟合入口：[scripts/GLMM/GLMM.R](../scripts/GLMM/GLMM.R)；已有结果复核入口：[scripts/GLMM/GLMM_reads_diagnostics_full.R](../scripts/GLMM/GLMM_reads_diagnostics_full.R)。两者职责已不同于 10 月 12 日旧交接。
- `GLMM.R` 默认最多 8 个独立进程，按 CPU 和可用内存限制实际数量，内部单线程；历史参数只作起点，当前完整运行重新拟合。
- 结果目录：[results/GLMM_Z_TS_SD_Q_MC8_20261013_185720](../results/GLMM_Z_TS_SD_Q_MC8_20261013_185720/)。本次从 latest 导航确认其仍为当前保存的完整运行；没有在 10 月 19 日重跑。
- 保留多起点、Hessian/梯度、积分精度、零值、正值和尾部核查。2,000 次模拟条件于拟合的 Repeat/TS/零值 TS 偏移，不是 2,000 次重新拟合。
- beta 的 SE 来自惩罚目标 Hessian 的逆矩阵近似；`z = beta/SE`，双侧正态近似 p，逐项 95% CI 为 `beta ± 1.96 SE`。
- DO 的 BH 族为计划 8 个 MC，MC15 虽无有效 p 仍计入计划数量。宿主全部对照系数跨 MC 为另一族，reads 样条系数和截距另分族。样条单个基系数的 p 不表示整体 reads 效应。
- 残差模拟用于模型检查，beta 显著性来自近似 Wald。未开展 CV、自动最优 lambda 选择或重拟合 bootstrap；TS 留组验证曾被用户明确暂缓。

## 5. GLMM：已核实结果与历史变化

### 5.1 DO 结果

来源：[00_DO_significance.csv](../results/GLMM_Z_TS_SD_Q_MC8_20261013_185720/00_DO_significance.csv)。本表 q 为 BH 校正值，CI 是逐项区间。

| MC | DO beta | 95% CI | BH q | 当前判断 |
|---|---:|---|---:|---|
| MC0 | 0.04596 | −0.02850～0.12042 | 0.301775 | 未通过 BH |
| MC6 | 0.13611 | 0.05773～0.21449 | 0.002661 | 正关联 |
| MC10 | 0.13871 | −0.05402～0.33144 | 0.275977 | 未通过 BH |
| MC15 | 暂不报告有效估计 | — | — | 数值检查失败 |
| MC21 | −0.22841 | −0.31804～−0.13878 | 0.000004716 | 负关联 |
| MC25 | −0.10897 | −0.26552～0.04758 | 0.275977 | 未通过 BH |
| MC507 | 0.22185 | 0.03902～0.40467 | 0.046385 | 正关联 |
| MC520 | −0.11093 | −0.46245～0.24059 | 0.612856 | 未通过 BH |

MC15 在原始系数表中有点估计（DO 约 0.71344），但 Hessian 未正定，且未由至少两个有效起点复现。旧对话进一步核对过：六个起点的目标值接近，主要问题是这些起点均未通过 Hessian；不能简单描述成各起点跑到完全不同的答案。本轮跳过其残差诊断；展示保留 MC15，写 `Not estimable`。

### 5.2 解释与剩余问题

- 宿主 beta 是相对 Acropora、调整 DO、reads、Repeat/TS 后的计数比例倾向；exp(beta) 为条件优势比。检测到某个 MC 的概率还同时取决于零分量与抽样读深。
- 旧对话核对过 MC21–Montipora 正对照，MC507 有 9 个正向宿主/环境类别对照通过 BH；宿主完整数值从 [00_beta_significance.csv](../results/GLMM_Z_TS_SD_Q_MC8_20261013_185720/00_beta_significance.csv) 读取。
- 环境梯度实际来自 7 个 TS。DO 和其他批次环境因素可能同步变化，解释需结合批次设计和结构敏感性。
- [00_model_checks.csv](../results/GLMM_Z_TS_SD_Q_MC8_20261013_185720/00_model_checks.csv) 中 7 个有效 MC 未触发当前精简诊断的主要警报；高比例校准仍有具体偏差：MC6、MC507 高尾低估，MC520 高预测组均值低估，MC0 极端最大值略超模拟范围。
- 最新精简诊断和早期完整诊断的项目不同，不能只按“全部诊断 p 大于 0.05”判断模型改善。
- 推荐已有展示方向：DO beta 森林图、宿主 beta 热图、条件预测曲线；补充展示 sigma(reads)、按 TS 零值校准、正值与尾部校准。预测曲线目前的参照条件和是否有置信带须从保存数据核对。
- 当前简洁汇总工作簿：[results/Z_TS_SD_Q_模型结果与主要结论_20261013_185720.xlsx](../archive/reports/Z_TS_SD_Q_模型结果与主要结论_20261013_185720.xlsx)。

### 5.3 模型推进历史

1. 10 月 12 日起点：未正则化 A/B/C/D 四候选。随后用户授权改脚本、自行运行，完成 D 模型的三档惩罚 `lambda=0.03125/0.125/0.5`，稳定了分离问题，但仍有残差缺口。
2. 第一轮结构筛查：D0、Z_DO、Z_TS、SD_Q、Z_DO_SD_Q、Z_TS_SD_Q、TAIL。支持批次零值结构与二次宽度，优先保留 Z_TS_SD_Q。
3. 后续尝试 SD_M、ALPHA_M、SD_ALPHA_M；部分 MC 改善，另一些数值不稳或结论敏感。用户认为模型过于复杂，认可先做减法、整理科学结论。
4. 精简复核比较 D0、Z_TS_SD_Q、SD_M，复用保存结果：MC21、MC507 的 DO 方向和预测差异较稳定；MC0、MC10、MC15 的效应幅度更依赖结构。历史 MC15 有可复用的有效拟合，不代表 10 月 13 日完整重拟合也有效。
5. 随后改写 GLMM.R 完整重拟合 Z_TS_SD_Q，并输出 beta 近似显著性；当前以第 5.1 节所列完整运行为准。

历史运行主要在 `archive/mc8_models/`，旧位置可能保留兼容联接。需要具体历史时按运行名定位：`GLMM_D_zero_penalty_MC8_20261012_160110`、`GLMM_structure_suite_MC8_20261012_185348`、`GLMM_structure_suite_MC8_20261012_235153`、`GLMM_saved_review_MC8_20261013_112414`。

## 6. PERMANOVA：数据范围、模型和置换

### 6.1 当前分析范围

输入：[data/PERMANOVA_input_MC8_20261011.xlsx](../data/PERMANOVA_input_MC8_20261011.xlsx)。先在 Repeat_group 内对原始 MC 比例取均值，再平方根转换计算 Bray–Curtis。所选 MC 不重新归一化为 100%。

| 范围 | 总体/季节重复组 | 类别 | 环境检验重复组 |
|---|---:|---|---:|
| all_hosts（当前为全样本） | 115 | 14 类，包括 Seawater、Sediment | 101，TS1–TS7 |
| coral | 39 | 3 类珊瑚 | 33，TS1–TS7 |
| sample_type | 117 | HOST=103、Seawater=7、Sediment=7 | 本层回答三类样本比较 |

逐宿主分析排除只有一个 Repeat_group 的 Cuttlefish、Sargassum；三类分析将其并入 HOST，故总数为 117。HOST 表示合并的生物宿主样本类别。

Seawater、Sediment 仅采于 TS2–TS8，TS1 缺采格子不补零、不造行。Host/TS/Season 使用具备所需信息的实际观测；各环境变量独立筛选有效批次，当前均为 TS1–TS7。完整 MC 响应矩阵须无缺失；不完整的“宿主×批次”布局与某条观测的环境协变量缺失分开处理。

### 6.2 公式、置换与多重校正

| 问题 | 模型/调整项 | 打乱什么 | BH 范围 |
|---|---|---|---|
| Host_type | 距离 ~ Host_type + TS，取宿主边际项 | TS 内置换 Repeat_group 的完整 MC 组成行 | 与 TS 总体检验合为 2 项，各分析范围分别校正 |
| TS | 同上，取 TS 边际项 | Host_type 内置换组成行 | 同上 |
| Season | 调整 Host_type | TS 之间整组重新分配季节标签，MC 行和宿主原位 | 单项，报告原始 p |
| 每项环境指标 | 调整 Host_type，逐变量检验 | TS 之间整组重新分配该环境值 | 每范围 8 项 |
| Sample_type 总体 | 距离 ~ Sample_type + TS，取样本类型项 | TS 内置换组成行 | 单项原始 p |
| Host/Sample_type 两两 | 只保留共同 TS，调整 TS | 共同 TS 内置换组成行 | 各比较族全部计划对数 |
| TS 两两 | 只保留共同宿主，调整 Host_type | 共同宿主内置换组成行 | 各范围计划 28 对 |

宿主/TS 用 `adonis2(..., by="margin")`。原始 Bray 是主结果，Lingoes 为敏感性检查。分类因素另有 PERMDISP，连续环境指标不人为分组做 PERMDISP。

环境批次整组置换保持同一 TS 的所有行共用一个值，不把同批次的多个宿主当成新的独立环境水平。7 个互异值有 5,040 种分配，有并列值时去除重复分配；当前 pH 为 1,260 种，DIN、浊度各 2,520 种。季节 TS1–TS8 的不同标签分配为 70 种。常规行置换设置上限 9,999，具体比较的实际置换空间需读输出。

DIP 的原始低于检出限记录保留，但未纳入当前 8 项环境主检验。后来 SPSS 表按用户要求将 DIP 填成检出限一半，未回写或重跑这轮 PERMANOVA。

### 6.3 R²、F、p、BH 和热图解释

- 总 SS 来自全部参与观测的距离关系；在通常的距离平方和表达下为 `sum(d_ij^2, i<j) / n`。
- 边际项 SS 为加入该项后、相对仅含其他调整项模型减少的残差 SS；R² 为该项 SS 除以本次分析的总 SS。分母保留总 SS，并非扣除 TS 后的残差 SS。
- pseudo-F 比较每自由度的项 SS 与残差 SS；p 根据所规定置换中 F 不小于观测 F 的比例得出。BH 用于控制同一比较族的多重检验，q 是调整后的 p。
- 两两热图每格来自该对实际样本和共同覆盖范围，R² 的分母因比较而异。较大的 R² 表示该对分析中更高的解释比例，可能与组间位置、组内波动及样本组成共同相关。
- 热图星号依据 BH q；`†` 表示两两分散度检验提示，NA 保留无法估计的比较。读取当轮图例确认具体标注。
- 同一组成距离可有旋转/翻转后的 NMDS 坐标。NMDS 展示距离排序，PERMANOVA 检验使用完整距离矩阵。图上点和椭圆未扣除 TS；调整 TS 的是统计检验。
- DO 组成柱图中 R1、R2 是“宿主×TS”内重新编号的 Repeat_group，每柱为组内平均组成；不同 TS 的 R1 不代表同一个持续追踪个体。

## 7. PERMANOVA：本次核实的结果

原始结果目录：[results/PERMANOVA_MC8_20261014_183535](../results/PERMANOVA_MC8_20261014_183535/)。数字核对分范围原始 CSV，汇总文件若有 NA 优先回到这些表。

| 检验 | 全样本 R² | 全样本 p/q | 珊瑚 R² | 珊瑚 p/q |
|---|---:|---:|---:|---:|
| Host_type，调整 TS | 0.192765 | q=0.0001 | 0.123567 | q=0.0005 |
| TS，调整宿主 | 0.298723 | q=0.0001 | 0.373633 | q=0.0002 |
| Season，调整宿主 | 0.043179 | p=0.428571 | 0.045732 | p=0.485714 |
| DO，调整宿主 | 0.157613 | q=0.020635 | 0.217644 | q=0.015873 |

其余 7 项环境指标均未通过各范围的 8 项 BH 校正。宿主 PERMDISP 在两个范围均 p=0.0001；TS PERMDISP 全样本 p=0.0041、珊瑚 p=0.0564。解释组成效应时结合这些分散度结果。

| 两两比较族 | 计划对数 | 可估计 | BH 显著 |
|---|---:|---:|---:|
| 全样本 Host_type | 91 | 86 | 15 |
| 全样本 TS | 28 | 28 | 24 |
| 珊瑚 Host_type | 3 | 3 | 2 |
| 珊瑚 TS | 28 | 28 | 6 |
| HOST/Seawater/Sediment | 3 | 3 | 2 |

跳过比较保留 NA，BH 的计划族大小不随跳过行缩小。共同区组不足、某组重复不足、模型秩不足或无残差自由度均可能导致跳过；具体原因查 `Status`。

三类样本总体：N=117，R²=0.030686，原始 p=0.0157。三类两两均仅用共同 TS2–TS8：

| 比较 | N1+N2 | R² | BH q |
|---|---|---:|---:|
| HOST–Seawater | 88+7 | 0.025147 | 0.0468 |
| HOST–Sediment | 88+7 | 0.013549 | 0.1292 |
| Seawater–Sediment | 7+7 | 0.098072 | 0.046875 |

珊瑚宿主两两：Acropora–Montipora q=0.0594；Acropora–Platygyra q=0.0054；Montipora–Platygyra q=0.0003。珊瑚 TS4–TS7 的 R² 约 0.455，具体 q 与共同宿主范围从该对原始行读取。

## 8. PERMANOVA：脚本、图和交互要求

- 完整分析：`scripts/PERMANOVA/run_permanova.R`；数据规则：`01_data_settings.R`；Host/TS：`02_test_host_time.R`；Season/环境：`03_test_season_environment_batch_labels.R`（旧 03 入口调用同版）。
- 两两：`05_test_pairwise_host_time.R`、`permanova_pairwise_helpers.R`；三类样本：`08_test_sample_type.R`。
- 总体图与 3D：`04_plot_permanova.R`；两两热图：`06_plot_pairwise_permanova.R`；逐宿主高亮：`07_plot_host_highlight_panels.R`；三类作图：`09_plot_sample_type.R`；DO 柱图：`plot_DO_composition_preview.R`。
- 当前图：`figures/PERMANOVA_MC8_20261014_183535/`，总览 `00_all_final_results.pdf`，五张独立两两热图在 `pairwise_heatmaps/`。
- 当前宿主 GLMM 热图脚本：`08_plot_GLMM_host_beta_heatmap.R`，当前图入口见 `figures/README.md`。
- 3D 仅制作 Host 和 TS，all_hosts、coral 各两份：`05_host_NMDS_3D.html`、`06_time_NMDS_3D.html`。只换展示可以复用匹配缓存；数据/维数改变须重新拟合 NMDS。
- 3D 悬停图例时放大目标组，其余保留原色；点击时可多选，未选组灰色淡化；移出悬停后恢复点击选择，`Show all` 清除。固定轴范围、轴位置和 XYZ，保留旋转缩放与正交投影。
- 分类宿主用离散高对比色，DO 用连续渐变。总览统一圆点、细边椭圆；早期逐宿主分面曾尝试背景椭圆，用户随后撤回。分面每行 3 列，当前 14 类自动为 5 行，早期 12 类才是 3×4。
- HTML 若有同名 `_files` 目录，保存分享时整组保留。
- `docs/README_permanova_readable.md` 是包含环境样本后的说明；`docs/README_batch_labels.md` 后半部保留旧排除环境样本时的 101/89 组数量，阅读时用本交接和当前结果更新这些数字。

## 9. 环境 Excel 与 PM11 历史数据

交付文件：[results/TS1-TS7_环境指标_SPSS.xlsx](../results/TS1-TS7_环境指标_SPSS.xlsx)。继续修改前先重新读取，用户可能手动编辑。

### 9.1 Environment_TS1_TS7

每 TS 一行，共 7 行；9 个环境指标在前，TS 在最后：`Salinity_psu, SST_C, pH, DO_mg_L, Chla_ug_L, DIN_mgN_L, DIP_mgP_L, PSM_mg_L, Turbidity_NTU, TS`。

来自 PERMANOVA 输入的 Metadata，同 TS 在多个宿主行的环境记录已核对一致并去重。TS8 无对应环境数据。TS1–TS3 的 DIP 原文 `<0.002 mg P/L`，用户指定取一半，填 `0.001`；真正缺失值留空。

| TS | 月份 | 季节 | DO，mg/L |
|---|---|---|---:|
| TS1 | 2024-04 | Wet | 6.8 |
| TS2 | 2024-07 | Wet | 5.5 |
| TS3 | 2024-10 | Dry | 6.7 |
| TS4 | 2025-02 | Dry | 8.1 |
| TS5 | 2025-04 | Wet | 7.9 |
| TS6 | 2025-08 | Wet | 3.8 |
| TS7 | 2025-11 | Dry | 5.3 |

湿季为 4–9 月，干季为 10 月至次年 3 月；当前 7 TS 为干 3、湿 4。

### 9.2 PM11_Surface_2015_2025

- 用户初称“十年”，明确选择后最终为 2015–2025 共 11 年，与原季节图一致。按年月顺序保留各年观测，没有先平均成 12 个气候月均值。
- 132 行=11×12，其中 127 个有观测年月，5 个缺测年月的环境值留空。前 9 列环境指标，接 `Season`、`Month`；用户明确要求无年份列。Season：1=Dry、2=Wet；Month 为数字 1–12，最后一个有效列。多余空白列不无故改动。
- 原始 CSV：`data/water_quality/Bluff Island_PM11_marine_water_quality.csv`，选 `Station=PM11`、`Depth=Surface Water`。
- 共 128 条表层采样记录。2020 年 11 月有两次，逐指标取算术月均形成 127 个年月观测。低于检出限的值先替换成一半，再月内平均。
- 缺月：2020-02、2020-03、2020-08、2021-01、2022-03。
- DIN 对应原始 `Total Inorganic Nitrogen`；DIP 对应 `Orthophosphate Phosphorus`；PSM 对应 `Suspended Solids`。
- 可追溯 JSON、提取/写入脚本和加历史表前备份位于 `Codex_Workspace/outputs/01a12e8a-d27b-7010-bee3-dd29e9d0f5d1/glmm_result_summary/`：`historical_environment.json`、`extract_history.py`、`edit_history.mjs`、`TS1-TS7_before_history_backup.xlsx`。该 JSON 保留年月 audit，可用于恢复年份。

## 10. 最近环境统计讨论和待继续的问题

- 当前每 TS 一个环境值，若把 7 个 TS 当作 7 组，常规单因素 ANOVA 残差自由度为 7−7=0。复制到宿主样本行不会产生新的独立环境重复。
- 环境指标与 MC 的连续梯度关联，无需先通过 TS 或干湿季差异检验；季节分组可能掩盖季节内部连续变化。研究问题分别是“组间环境是否不同”和“环境变化是否对应 MC 变化”。
- 用户自述 7 TS 的环境季节对比不显著，尚未提供具体 SPSS 输出和方法；这与已保存 PERMANOVA 的 MC 组成季节检验是两个问题。
- 干 3、湿 4 时，无并列值的常规双侧精确 Mann–Whitney 最小 p 为 2/35≈0.057；这是该检验约定下的分辨率，不泛化为所有置换检验。
- 最新问题：“不正态用非参数检验对吗？”既有讲解重点：先看组内分布/模型残差和设计，再选独立或配对方法；不要仅用所有月份混合后的正态检验作决定。
- 已讨论的备选：独立两组 Mann–Whitney，独立多组 Kruskal–Wallis；按年份区组比较 12 月可考虑 Friedman（常规需要完整区组）；每年干湿均值配对可考虑 Wilcoxon 符号秩（还需考虑差值对称性）。多年记录的年份结构和时间相关性应进入具体方案。
- 尚未对 11 年数据运行上述显著性检验，也未替用户决定最终方法。若下一步授权分析，先明确比较目标、年份/季节的配对单位、缺月处理、时间相关性与多指标 BH 范围，再运行。
- 可继续讨论的其他方向：环境变化图、指标间相关性、DO 与其他批次环境因素混杂、相对历史同期的单个 TS 异常评估；这些均是待选方向。

## 11. PPT、筛选理由与讲解文档

- 当前主要 PPT：[timeseries_xzy_20261015.pptx](../archive/old_outputs/presentations/timeseries_xzy_20261015.pptx)，备份 `timeseries_修改前备份_20261015.pptx`。用户曾由 timeseries.pptx 更名；编辑前确认指定文件，不误改 time+lasso.pptx。
- 旧对话最后核对主 PPT 为 15 页；涵盖营养盐、宿主/时间/季节 PERMANOVA、DO 梯度、MC 丰度、GLMM、宿主 beta 热图及最后 Lasso。
- 第 12 页已按要求加入 Selection rationale、保留 MC15 并写 Not estimable；用户曾要求删除整份 PPT 备注。后续核对以当前实际 PPT 为准。
- 明确给定的标签：MC10=NGP，MC15/MC0=NGP control，MC6=Denitrifying specialists。功能标签原始证据与生态筛选依据分别记录。
- 生态筛选来源：`data/source_tables/ISA_asv.xlsx`（IndVal>0.5）、`data/Cooccurence.xlsx`（共现对 n_TS≥4）、`data/source_tables/Network.xlsx`（betweenness 全库 top 5%，部分 top 1%）。这些表的单位为 ASV，通过 ASV–MC 对应提供候选证据。
- 既有核对：MC6、MC21 三表都有；MC507 有宿主关联与网络；MC520 有宿主关联；MC25 有共现；MC10、MC15 有网络；MC0 未列三表。部分匹配为 `WEAK_ownMC`，原始质量标记应保留。
- IndVal 入选阈值与 p 显著性分开；表中曾有 p=0.073 的入选记录。共现 `max_abs_r` 不保留方向；网络 top 5% 是排名，后续功能解释应追溯原始证据。
- 已交付 [timeseries_逐页讲解与初步结论_20261015.docx](../archive/old_outputs/presentations/timeseries_逐页讲解与初步结论_20261015.docx)，原对话记录为 35 页，逐页解释 PPT 的读图方法、知识点、初步结论与问答，含 GLMM 方程和来源。
- PERMANOVA 的 PPT 静态图与四份嵌入 3D HTML 曾核对更新；用户要求气泡图保持完整、已附 PDF，尽量不改原有内容。
- 最后 Lasso 的 selection rate 表示入选稳定性，区别于 p 值；候选生物学解释需结合注释来源和具体效应。

## 12. 文件整理记录与下一次接续顺序

10 月 14 日曾整理研究目录：当前脚本留 `scripts/`，当前结果与图留 `results/`、`figures/`，方法说明到 `docs/`，历史运行与备份到 `archive/`。旧路径的隐藏兼容联接指向同一份归档文件。原对话报告当时 9,437 个非占用标记文件内容一致、19 个 R 脚本语法通过；本次未重做该全目录校验。

移动记录：[archive/folder_layout_v1/move_plan.csv](../archive/folder_layout_v1/move_plan.csv)。项目后来新增的 PPT/Word 仍可在根目录，不能仅凭整理约定猜其位置。

后续接续：

1. 先读本文件第 1 节，明确用户本次实际问题。
2. 涉及当前结果时读相关原始 CSV 和 `scripts_used/` 或 `01_analysis_script_used.R`；解释当轮结果以运行快照为依据。
3. 涉及环境检验时先读取当前 Excel，并据 audit 恢复必要的年份信息；按用户新指令决定是否运行。
4. 涉及 PPT 时先确认当前目标文件，必要时读取逐页讲解文档；延续 MC15 保留、简洁表格、英文图及不自动打开侧栏的约定。
5. 需要更早的细节再定向检索两个文字存档；不把历史候选、建议或旧诊断当作当前结果。

本次仅新建交接文件和文字存档，并只读核对来源；未改变研究数据、模型脚本、既有分析结果或演示文件。
