19MC x 4营养盐，仅TS随机截距模型，共76项Wald beta检验，单一BH族n=76。
98个coral原始样本，TS1-TS7含TS4/TS6，不平均，无Host/Repeat/SST项。
响应为MC reads和其它Ruegeria reads，分母加总原表全部MC，不对19MC重新归一化。
公式固定为cbind(MC_reads, Other_reads) ~ E_z + (1 | TS)；Beta-binomial/logit，dispformula=~1，ziformula=~0。
每次单一营养盐，按7TS等权中心/SD标准化。
beta_P_value直接取summary(fit)$coefficients$cond[E_z, Pr(>|z|)]。不运行bootstrap、标签置换或LOTO。
Wald P/CI是渐近近似；环境信息来自7TS。高硝酸盐信息集中于TS4、TS6。
beta_q_value对全部76项统一BH，无法可靠估计项留NA但保持n=76，不按显著性选择MC。
数值问题、分布诊断标记分开；少数极端值或LOTO翻转不自动淘汰完整模型。
保留原始报告值及L/2删失规则：NO2有4/7、DIP有3/7低检出限，详见input_environment_TS.csv。
默认核对后仅复用120255的4个有TS硝酸盐拟合，新增72个；关闭复用则拟合全部76个。
TS-only范围为查看172226两结构结果后由用户确定的后续方案；历史152项P/q及输出保持原样。
默认另做已保存模型的2000次模拟诊断，不重新拟合；诊断P不作为beta关联P。
00_result.csv为主要结果；完整字段、模型、输入、诊断、源码和状态在supporting_files。
