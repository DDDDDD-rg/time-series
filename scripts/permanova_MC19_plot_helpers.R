# MC19 coral：只读取已完成的分析结果；不重跑PERMANOVA或BH。
# NMDS仅在没有兼容坐标缓存时拟合；每个点/热图列/组成柱均为原始样本。

mc19_plot_context <- function() {
  packages <- c("ggplot2", "vegan")
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly=TRUE)]
  if (length(missing)) stop("缺少作图R包：", paste(missing, collapse=", "))
  suppressPackageStartupMessages(library(ggplot2))
  pointer <- file.path(project_dir, "results", "00_latest_PERMANOVA_MC19_env7_coral.txt")
  result_dir <- getOption("mc19.result_dir", NULL)
  if (is.null(result_dir)) {
    if (!file.exists(pointer)) stop("尚无MC19结果。请先运行run_permanova.R，或指定mc19.result_dir。")
    result_dir <- readLines(pointer, n=1L, warn=FALSE, encoding="UTF-8")
  }
  status <- mc19_support_file(result_dir, "00_run_status.txt")
  if (!file.exists(status) || !identical(readLines(status, n=1L, warn=FALSE), "complete")) {
    stop("该分析尚未完整完成，不能据此出图：", result_dir)
  }
  source_rds <- mc19_support_file(result_dir, "coral_environment_analysis.rds")
  saved <- readRDS(source_rds)
  if (!is.character(saved$version) || length(saved$version) != 1L ||
      !saved$version %in% compatible_analysis_versions) stop("结果版本不符，不能使用旧MC8结果。")
  input <- saved$input
  stopifnot(identical(colnames(input$community), mc_cols),
    identical(rownames(input$community), as.character(input$metadata$New_ID)),
    nrow(input$metadata)==nrow(input$community), all(input$metadata$Host_niche=="coral"),
    !anyDuplicated(input$metadata$New_ID), identical(saved$results$Factor, env_cols),
    all(saved$results$N_samples==nrow(input$metadata)), all(saved$results$N_MC==19L),
    all(saved$results$Status=="ok"), all(is.finite(saved$results$P_BH_raw)))
  snapshot <- mc19_support_file(result_dir, "00_input_used_MC19_env7.xlsx")
  if (!file.exists(snapshot) || !identical(unname(tools::md5sum(snapshot)), input$input_md5)) {
    stop("保存的输入快照与分析指纹不符。")
  }
  figure_root <- getOption("mc19.figure_dir", file.path(project_dir, "figures", basename(result_dir)))
  figure_dir <- file.path(figure_root, "supporting_files")
  dir.create(figure_dir, recursive=TRUE, showWarnings=FALSE)
  labels <- c(Salinity_psu="Salinity (psu)", SST_C="Temperature (deg C)", pH="pH",
    Nitrate_mgN_L="Nitrate (mg N/L)", Nitrite_mgN_L="Nitrite (mg N/L)",
    Ammonia_mgN_L="Ammonia nitrogen (mg N/L)", DIP_mgP_L="DIP (mg P/L)")
  ts_colours <- setNames(c("#0072B2","#D55E00","#009E73","#CC79A7",
    "#E69F00","#332288","#56B4E9","#333333"), ts_levels)
  host_colours <- c("Acropora sp."="#0072B2", "Montipora sp."="#D55E00", "Platygyra sp."="#009E73")
  mc_colours <- setNames(grDevices::hcl.colors(length(mc_cols), "Dark 3"), mc_cols)
  mc_colours <- c(mc_colours, Other_Ruegeria="#DDDDDD")
  list(saved=saved, input=input, result_dir=result_dir, source_rds=source_rds,
       figure_root=figure_root, figure_dir=figure_dir, labels=labels,
       ts_colours=ts_colours, host_colours=host_colours, mc_colours=mc_colours)
}

mc19_theme <- function() {
  ggplot2::theme_classic(base_size=12, base_family="Arial") +
    ggplot2::theme(plot.title=ggplot2::element_blank(), plot.subtitle=ggplot2::element_blank(),
      axis.text=ggplot2::element_text(colour="black"),
      legend.key=ggplot2::element_rect(fill="white", colour=NA),
      plot.caption=ggplot2::element_text(hjust=0, size=9),
      plot.margin=ggplot2::margin(10,25,10,10))
}

mc19_percent <- function(x) paste0(format(round(100*x, 2), trim=TRUE), "%")

mc19_nmds <- function(context, k=2L, fit_if_missing=TRUE) {
  cache_file <- mc19_support_file(context$result_dir, paste0("plot_NMDS", k, "D.rds"))
  ids <- as.character(context$input$metadata$New_ID)
  if (file.exists(cache_file)) {
    cache <- readRDS(cache_file)
    compatible <- identical(cache$version, context$saved$version) &&
      identical(cache$input_md5, context$input$input_md5) && identical(cache$New_ID, ids) &&
      identical(cache$MC, mc_cols) && identical(cache$k, k)
    if (!compatible) stop("NMDS缓存版本不符，请移走旧缓存后重新运行04：", cache_file)
    return(cache)
  }
  if (!fit_if_missing) stop("缺少二维NMDS缓存，请先运行04_plot_permanova.R。")
  set.seed(20261026L + k)
  fit <- vegan::metaMDS(context$saved$distance, k=k, trymax=100,
    autotransform=FALSE, wascores=FALSE, trace=FALSE, plot=FALSE)
  coordinates <- as.matrix(vegan::scores(fit, display="sites", choices=seq_len(k)))
  if (is.null(rownames(coordinates)) || !setequal(rownames(coordinates), ids)) {
    stop("NMDS坐标的样本标签与分析样本不一致。")
  }
  coordinates <- coordinates[match(ids, rownames(coordinates)), , drop=FALSE]
  colnames(coordinates) <- paste0("NMDS", seq_len(k))
  if (any(!is.finite(coordinates)) || !is.finite(fit$stress)) stop("NMDS未产生有效坐标。")
  cache <- list(version=context$saved$version, input_md5=context$input$input_md5, New_ID=ids,
    MC=mc_cols, k=k, coordinates=coordinates, stress=fit$stress, fit=fit)
  dir.create(dirname(cache_file), recursive=TRUE, showWarnings=FALSE)
  saveRDS(cache, cache_file)
  cache
}

mc19_main_plots <- function(context) {
  d <- context$input$metadata
  cache <- mc19_nmds(context, 2L)
  scores <- cbind(d, as.data.frame(cache$coordinates))
  scores$TS <- factor(scores$TS, levels=ts_levels)
  write_csv_excel(scores, file.path(context$figure_dir, "NMDS2D_sample_coordinates.csv"))
  caption <- sprintf("Individual coral samples: n = %d; 19 MCs; NMDS stress = %.3f.", nrow(d), cache$stress)
  plots <- list()
  plots$NMDS_TS <- list(plot=ggplot(scores, aes(NMDS1, NMDS2, colour=TS)) +
    geom_point(size=2.5, alpha=0.85) + coord_equal() +
    scale_colour_manual(values=context$ts_colours, drop=TRUE) +
    labs(x="NMDS1", y="NMDS2", caption=caption) + mc19_theme(), width=8, height=6)
  for (variable in env_cols) {
    display <- scores; display$Env_value <- display[[variable]]
    result <- context$saved$results[context$saved$results$Factor==variable, , drop=FALSE]
    text <- paste(caption, sprintf("R2 = %.3f; BH q = %.4g (7 environment tests).", result$R2_raw, result$P_BH_raw), sep="\n")
    plots[[paste0("NMDS_environment_", variable)]] <- list(plot=
      ggplot(display, aes(NMDS1, NMDS2, colour=Env_value)) + geom_point(size=2.5, alpha=0.9) +
      coord_equal() + scale_colour_gradientn(colours=grDevices::hcl.colors(128, "Viridis")) +
      labs(x="NMDS1", y="NMDS2", colour=unname(context$labels[variable]), caption=text) + mc19_theme(), width=8, height=6)
  }
  results <- context$saved$results
  results$Indicator <- factor(unname(context$labels[results$Factor]), levels=rev(unname(context$labels[env_cols])))
  results$Label <- sprintf("R2 = %.3f; q = %.3g", results$R2_raw, results$P_BH_raw)
  plots$environment_R2_BH <- list(plot=ggplot(results, aes(R2_raw, Indicator)) +
    geom_col(fill="grey60", width=0.65) + geom_text(aes(label=Label), hjust=-0.05, size=3.2) +
    scale_x_continuous(expand=expansion(mult=c(0.02,0.65))) +
    labs(x="PERMANOVA R2", y=NULL, caption="Seven single-environment tests; BH q values read from saved results.") +
    mc19_theme(), width=10, height=6)
  batches <- context$input$batches
  environment_long <- do.call(rbind, lapply(env_cols, function(variable) {
    data.frame(TS=batches$TS, Date=as.Date(batches$Water_quality_date),
      Indicator=unname(context$labels[variable]), Value=batches[[variable]])
  }))
  plots$environment_time_series <- list(plot=ggplot(environment_long, aes(Date, Value)) +
    geom_line(linewidth=0.45) + geom_point(size=2) + facet_wrap(~Indicator, scales="free_y", ncol=3) +
    scale_x_date(date_labels="%Y-%m") + labs(x="Water-quality sampling date", y=NULL,
      caption="One observed environment record per TS; values below detection use L/2.") +
    mc19_theme() + theme(axis.text.x=element_text(angle=45,hjust=1)), width=12, height=9)
  raw <- do.call(rbind, lapply(mc_cols, function(mc) data.frame(New_ID=d$New_ID,
    TS=factor(d$TS, levels=ts_levels), MC=mc, Abundance=context$input$community[,mc])))
  raw$MC <- factor(raw$MC, levels=mc_cols)
  plots$MC_individual_samples_by_TS <- list(plot=ggplot(raw, aes(TS, Abundance, colour=TS)) +
    geom_point(position=position_jitter(width=0.15,height=0,seed=20261026), size=1.2, alpha=0.7) +
    facet_wrap(~MC, scales="free_y", ncol=4) + scale_colour_manual(values=context$ts_colours, guide="none") +
    scale_y_continuous(labels=mc19_percent) + labs(x="TS", y="Relative abundance within Ruegeria",
      caption="Each point is one original sample. Samples are not averaged or linked as repeated individuals.") +
    mc19_theme(), width=13, height=14)
  order <- order(match(d$TS, ts_levels), d$New_ID)
  sample_ids <- as.character(d$New_ID[order])
  raw$Sample <- factor(as.character(raw$New_ID), levels=sample_ids)
  plots$MC_sample_abundance_heatmap <- list(plot=ggplot(raw, aes(Sample, MC, fill=Abundance)) +
    geom_tile() + scale_fill_gradientn(colours=grDevices::hcl.colors(128,"Viridis"), labels=mc19_percent) +
    labs(x="Individual sample New_ID (ordered by TS)", y=NULL, fill="Relative abundance",
      caption="One column per sample; raw proportions retained, without renormalising the 19 MCs.") +
    mc19_theme() + theme(axis.text.x=element_text(angle=90,hjust=1,size=5),axis.ticks=element_blank()), width=16, height=7)
  write_csv_excel(d[order,c("New_ID","Sample_ID","TS","Host_type")], file.path(context$figure_dir,"heatmap_sample_order.csv"))
  plots
}

mc19_host_plot <- function(context, allow_nmds_fit=FALSE) {
  cache <- mc19_nmds(context, 2L, fit_if_missing=allow_nmds_fit)
  d <- cbind(context$input$metadata, as.data.frame(cache$coordinates))
  hosts <- unique(as.character(d$Host_type))
  foreground <- d; foreground$Panel <- as.character(foreground$Host_type)
  background <- do.call(rbind,lapply(hosts,function(host) {x<-d; x$Panel<-host; x}))
  p <- ggplot() + geom_point(data=background,aes(NMDS1,NMDS2),colour="grey85",size=1.8) +
    geom_point(data=foreground,aes(NMDS1,NMDS2,colour=Host_type),size=2.5) +
    facet_wrap(~Panel,ncol=3) + coord_equal() + scale_colour_manual(values=context$host_colours,guide="none") +
    labs(x="NMDS1",y="NMDS2",caption=sprintf("Descriptive host labels on the same %d sample coordinates; stress = %.3f. No host model or pairwise test.",nrow(d),cache$stress)) + mc19_theme()
  list(host_highlight_descriptive=list(plot=p,width=13,height=5))
}

mc19_composition_plots <- function(context) {
  d <- context$input$metadata; y <- context$input$community
  other <- 1-rowSums(y)
  if (any(other < -1e-10)) stop("19个MC的总占比超过1，无法计算其余Ruegeria。")
  y <- cbind(y,Other_Ruegeria=pmax(other,0))
  plots <- list(); orders <- list()
  for (variable in env_cols) {
    ix <- order(d[[variable]], match(d$TS,ts_levels), d$New_ID)
    ordered <- data.frame(Rank=seq_along(ix),New_ID=d$New_ID[ix],TS=d$TS[ix],Env_value=d[[variable]][ix],Factor=variable)
    orders[[variable]] <- ordered
    long <- do.call(rbind,lapply(colnames(y),function(mc) data.frame(Rank=ordered$Rank,MC=mc,Abundance=y[ix,mc])))
    long$MC <- factor(long$MC,levels=colnames(y))
    result <- context$saved$results[context$saved$results$Factor==variable,,drop=FALSE]
    caption <- sprintf("Each bar is one sample, ordered by environment value (ties: TS, New_ID). R2 = %.3f; BH q = %.4g. Other includes unselected Ruegeria reads.",result$R2_raw,result$P_BH_raw)
    plots[[paste0("MC_composition_",variable)]] <- list(plot=ggplot(long,aes(Rank,Abundance,fill=MC)) +
      geom_col(width=1) + scale_fill_manual(values=context$mc_colours,drop=FALSE) +
      scale_y_continuous(labels=mc19_percent,expand=expansion(mult=c(0,0.02))) +
      labs(x=paste("Sample rank by",unname(context$labels[variable])),y="Relative abundance within Ruegeria",fill="MC",caption=caption) +
      mc19_theme() + guides(fill=guide_legend(ncol=2)),width=15,height=6)
  }
  write_csv_excel(do.call(rbind,orders),file.path(context$figure_dir,"composition_sample_orders.csv"))
  plots
}

mc19_export_3d <- function(context) {
  required <- c("plotly","htmlwidgets")
  missing <- required[!vapply(required,requireNamespace,logical(1),quietly=TRUE)]
  if (length(missing)) stop("可选3D图缺少R包：",paste(missing,collapse=", "))
  cache <- mc19_nmds(context,3L)
  d <- cbind(context$input$metadata,as.data.frame(cache$coordinates))
  d$TS <- factor(d$TS,levels=ts_levels[ts_levels %in% d$TS])
  d$Hover <- paste0("New_ID: ",d$New_ID,"<br>Sample: ",d$Sample_ID,"<br>TS: ",d$TS,"<br>Host: ",d$Host_type)
  for (variable in env_cols) d$Hover <- paste0(d$Hover,"<br>",context$labels[variable],": ",d[[variable]])
  widget <- plotly::plot_ly(d,x=~NMDS1,y=~NMDS2,z=~NMDS3,type="scatter3d",mode="markers",
    color=~TS,colors=unname(context$ts_colours[levels(d$TS)]),text=~Hover,hoverinfo="text",marker=list(size=4,opacity=0.85))
  widget <- plotly::layout(widget,scene=list(aspectmode="data"),
    annotations=list(list(text=sprintf("coral samples; n = %d; 3D stress = %.3f",nrow(d),cache$stress),xref="paper",yref="paper",x=0,y=1.05,showarrow=FALSE)))
  htmlwidgets::saveWidget(widget,file.path(context$figure_dir,"NMDS3D_TS.html"),selfcontained=FALSE)
  write_csv_excel(d,file.path(context$figure_dir,"NMDS3D_sample_coordinates.csv"))
}

mc19_pdf_tools <- function() {
  profile <- Sys.getenv("USERPROFILE", unset="")
  bundled <- if (nzchar(profile)) file.path(profile, ".cache", "codex-runtimes",
    "codex-primary-runtime", "dependencies", "python", "python.exe") else character()
  candidates <- unique(c(getOption("mc19.python", character()), Sys.getenv("MC19_PYTHON", unset=""),
    unname(Sys.which(c("python", "python3"))), bundled))
  candidates <- candidates[nzchar(candidates) & file.exists(candidates)]
  script <- file.path(project_dir, "scripts", "bundle_MC19_figures.py")
  if (!file.exists(script)) stop("缺少PDF合并脚本：", script)
  for (python in candidates) {
    probe <- suppressWarnings(tryCatch(system2(python, c("-X", "utf8", "-c", shQuote("import pypdf")),
      stdout=FALSE, stderr=FALSE), error=function(e) 1L))
    if (identical(as.integer(probe), 0L)) return(list(python=python, script=script))
  }
  stop("合并PDF需要Python和pypdf。可先设置 options(mc19.python='解释器完整路径')，或环境变量MC19_PYTHON；修复后只运行04即可。")
}

mc19_bundle_figures <- function(context, tools, allow_partial=FALSE) {
  arguments <- c("-X", "utf8", shQuote(tools$script), "--figure-dir", shQuote(context$figure_root),
    if (allow_partial) "--allow-partial")
  output <- suppressWarnings(system2(tools$python, arguments, stdout=TRUE, stderr=TRUE))
  status <- attr(output, "status")
  if (!is.null(status) && status != 0L) stop("合并PDF失败：", paste(output, collapse="\n"))
  message(paste(output, collapse="\n"))
  invisible(output)
}

plot_mc19_suite <- function(parts=c("main","host","composition")) {
  if (!all(parts %in% c("main","host","composition"))) stop("未知作图分支。")
  pdf_tools <- mc19_pdf_tools()
  context <- mc19_plot_context()
  status <- file.path(context$figure_dir,"00_plot_status.txt")
  writeLines(c("running",context$result_dir),status)
  tryCatch({
    plots <- list()
    if ("main" %in% parts) plots <- c(plots,mc19_main_plots(context))
    if ("host" %in% parts) plots <- c(plots,mc19_host_plot(context,allow_nmds_fit="main" %in% parts))
    if ("composition" %in% parts) plots <- c(plots,mc19_composition_plots(context))
    source_names <- c("permanova_MC19_plot_helpers.R","permanova_MC19_data_settings.R","04_plot_permanova.R",
      "07_plot_host_highlight_panels.R","plot_MC_environment_composition.R","bundle_MC19_figures.py")
    snapshot <- file.path(context$figure_dir,"scripts_used")
    dir.create(snapshot,showWarnings=FALSE)
    if (!all(file.copy(file.path(project_dir,"scripts",source_names),snapshot,overwrite=TRUE))) stop("作图脚本快照失败。")
    pdf_device <- if (isTRUE(capabilities("cairo"))) grDevices::cairo_pdf else grDevices::pdf
    manifest <- list()
    for (name in names(plots)) {
      spec <- plots[[name]]
      png <- file.path(context$figure_dir,paste0(name,".png"))
      pdf <- file.path(context$figure_dir,paste0(name,".pdf"))
      ggplot2::ggsave(png,plot=spec$plot,width=spec$width,height=spec$height,dpi=300,bg="white")
      pdf_plot <- if (isTRUE(capabilities("cairo"))) spec$plot else spec$plot+ggplot2::theme(text=ggplot2::element_text(family="sans"))
      ggplot2::ggsave(pdf,plot=pdf_plot,device=pdf_device,width=spec$width,height=spec$height,bg="white")
      manifest[[name]] <- data.frame(Figure=name,Source_RDS=context$source_rds,
        Input_MD5=context$input$input_md5,Version=context$saved$version,N_samples=nrow(context$input$metadata),PNG=png,PDF=pdf)
    }
    manifest_path <- file.path(context$figure_dir,"figure_sources.csv")
    manifest <- do.call(rbind,manifest)
    if (file.exists(manifest_path)) {
      previous <- read.csv(manifest_path,check.names=FALSE,stringsAsFactors=FALSE,fileEncoding="UTF-8-BOM")
      previous <- previous[!previous$Figure %in% manifest$Figure,,drop=FALSE]
      manifest <- rbind(previous,manifest)
    }
    write_csv_excel(manifest,manifest_path)
    if ("main" %in% parts && isTRUE(getOption("mc19.interactive_3d",FALSE))) mc19_export_3d(context)
    writeLines(enc2utf8(c("MC19 coral figures",paste("Source:",context$result_dir),
      "All points/heatmap columns/composition bars are individual samples; no sample averaging.",
      "All seven environment tests are displayed, including nonsignificant tests; q comes from saved BH results.",
      "NMDS is a display of the saved Bray-Curtis distances; cached coordinates are reused.",
      "Host panels are descriptive only, with no Host_type or pairwise test.",
      "Environment lines connect observed TS records. MC-by-TS panels show raw sample points.",
      "Optional 3D: retain the HTML file together with its _files directory.")),
      file.path(context$figure_dir,"README.txt"),useBytes=TRUE)
    mc19_bundle_figures(context, pdf_tools, allow_partial=!setequal(parts,c("main","host","composition")))
    writeLines(c("complete",context$result_dir),status)
    writeLines(enc2utf8(context$figure_root),file.path(project_dir,"figures","00_latest_PERMANOVA_MC19_env7_coral.txt"),useBytes=TRUE)
    message("合并PDF已保存：",context$figure_root,"；单图和配套文件：",context$figure_dir)
    invisible(context$figure_root)
  },error=function(e) {
    writeLines(c("failed",conditionMessage(e)),status)
    stop(conditionMessage(e),call.=FALSE)
  })
}
