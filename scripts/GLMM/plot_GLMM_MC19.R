# 当前MC19 x 5营养盐/TS-only Wald作图：只读保存的95个有TS模型及P/q，不拟合、不模拟。
# 主图：19x4效应热图、全部BH显著MC×营养盐曲线（分页PDF及PNG）。无图题，英文标签。
# 诊断PDF仍位于supporting_files/diagnostics；旧5份主PDF运行时移至supporting_files/previous_main_figures。
# 历史bootstrap图用plot_GLMM_MC19_full_data.R；历史LOTO图用plot_GLMM_MC19_LOTO.R。

local({
  analysis_scope <- getOption("mc19.analysis_scope", "coral")  # coral / ecosystem
  project_dir <- "C:/Users/xzy/Desktop/CUHK_Study/GNBF6010_Research_Projects/Research/Msc main project-time series"
  results_dir <- ""              # 空字符串读env5/TS8/TS-only Wald独立latest。

  engine <- new.env(parent = globalenv())
  for (helper in c("glmm_MC19_pilot_helpers.R", "glmm_MC19_wald_helpers.R",
      "glmm_MC19_wald_plot_helpers.R")) sys.source(file.path(project_dir, "scripts", "GLMM", helper), envir = engine)
  cfg <- engine$mc19_wald_config(project_dir, analysis_scope)
  engine$mc19_wald_plots(cfg, results_dir)
})

# 2026-10-29版式修订：同一MC的所有显著指标放在同一行、同一页（最多4列，列宽固定，页面宽度随列数扩展）；每页最多3个MC，单指标MC独占一行。底部预留1.4英寸、行间1英寸，轴标签和三行简短图注分开。相同MC主图纵轴一致；极端值压缩曲线时主图改为放大范围，标记Full range的小图显示全部原始点与完整纵轴，自动选择较少遮挡曲线的上方角落。移除原1.8倍纵轴留白。只修改源码、检查分组及版式坐标，不运行正式绘图。

# 2026-10-29当前版：用户确认大图必须full range，小图显示局部区域。大图覆盖全部原始点及95%CI，同MC统一纵轴，不加1.8倍留白。需要放大时，在大图右侧独立留出Zoom小图区（包含边框、刻度、文字的区域也不与大图相交），相同x范围、缩小y范围；不再向大图绘图区内叠放小图。列距自动增至7.8英寸，主绘图区宽度和16pt基准字号保持。显著组合/同行同页/模型/Pq/热图均不变。supporting_files/curve_axis_ranges.csv记录主图和放大图范围。此次只修改与只读/坐标检查，现有图由用户重画。
