# 当前MC19 coral方案停用本历史入口；原代码留作查阅，原文件另有备份。
stop("本轮仅coral，不绘制旧HOST/Seawater/Sediment图。请使用run_permanova.R或04_plot_permanova.R。",call.=FALSE)

# HOST / Seawater / Sediment: 2D NMDS overview and pairwise heatmap.
# Read saved PERMANOVA results; never refit tests or recalculate BH.
# If no compatible coordinates exist, fit ONLY a 2D NMDS (4 workers).
# English labels, no titles/subtitles. Run this file when ready to export.
(function() {
  project_dir <- "C:/Users/xzy/Desktop/CUHK_Study/GNBF6010_Research_Projects/Research/Msc main project-time series"
  latest_pointer <- file.path(project_dir,"results","00_latest_PERMANOVA_MC8.txt")
  # Optional override: options(mc.results_root="C:/.../results/PERMANOVA_MC8_...")
  results_root <- getOption("mc.results_root",if(file.exists(latest_pointer))
    readLines(latest_pointer,n=1L,warn=FALSE,encoding="UTF-8") else "")
  workers <- 4L
  random_seed <- 20261014L
  trymax <- 200L
  show_ellipses <- TRUE
  ellipse_level <- 0.95
  reuse_saved_coordinates <- TRUE
  save_pdf <- TRUE
  font_family <- if(isTRUE(capabilities("cairo"))) "Arial" else "sans"
  dpi <- 300L
  output_dir <- file.path(project_dir,"figures",paste0("sample_type_",
    format(Sys.time(),"%Y%m%d_%H%M%S")))
  required <- c("ggplot2","vegan","viridisLite","scales")
  missing <- required[!vapply(required,requireNamespace,logical(1),quietly=TRUE)]
  if(length(missing)) stop("Install required packages: ",paste(missing,collapse=", "))
  suppressPackageStartupMessages(library(ggplot2))
  source(file.path(project_dir,"scripts", "PERMANOVA","permanova_plot_style.R"),encoding="UTF-8",local=TRUE)
  result_file <- file.path(results_root,"sample_type","sample_type_analysis.rds")
  if(!file.exists(result_file)) stop("Missing saved three-category result: ",result_file,
    ". Set mc.results_root to the completed run, or run 08_test_sample_type.R first.")
  saved <- readRDS(result_file)
  data <- as.data.frame(saved$data)
  result <- saved$result
  pairs <- saved$family$results
  mc <- c("MC6","MC507","MC21","MC520","MC10","MC15","MC0","MC25")
  groups <- c("HOST","Seawater","Sediment")
  needed <- c("Repeat_group_ID","Host_type","TS","Sample_type",mc)
  stopifnot(all(needed %in% names(data)),!anyNA(data[,needed]),
    !anyDuplicated(data$Repeat_group_ID),setequal(unique(as.character(data$Sample_type)),groups),
    all(is.finite(as.matrix(data[,mc]))),all(as.matrix(data[,mc])>=0),
    all(rowSums(data[,mc])>0),nrow(result)==1L,result$Status=="ok",
    result$Factor=="Sample_type",result$Scope=="sample_type",result$N_units==nrow(data),
    is.finite(result$P_raw),is.finite(result$R2_raw),
    nrow(pairs)==3L,all(pairs$Factor=="Sample_type"),all(pairs$Scope=="sample_type"),
    all(pairs$BH_family_size==3L),all(pairs$Input_MD5==result$Input_MD5))
  keys <- vapply(seq_len(nrow(pairs)),function(i)
    paste(sort(c(pairs$Group1[i],pairs$Group2[i])),collapse="|"),character(1))
  expected_keys <- apply(combn(groups,2),2,function(x) paste(sort(x),collapse="|"))
  stopifnot(setequal(keys,expected_keys),!anyDuplicated(keys))
  data <- data[order(data$Repeat_group_ID),,drop=FALSE]
  data$Sample_type <- factor(data$Sample_type,levels=groups)
  palette <- pv_subset(pv_sample_types,groups)
  fingerprint <- as.character(result$Input_MD5)
  distance_spec <- "Bray-Curtis; sqrt of original-denominator MC8 Repeat_group means; k=2"
  p_text <- function(p) ifelse(is.finite(p),formatC(p,format="g",digits=3),"NA")

  # Full 117-row coordinates are needed; the all_hosts ordination has 115 rows.
  # Only reuse this script's caches after matching IDs, proportions and metadata.
  find_coordinates <- function(data) {
    folders <- list.dirs(file.path(project_dir,"figures"),recursive=FALSE,full.names=TRUE)
    folders <- folders[grepl("^sample_type_",basename(folders))]
    files <- file.path(folders,"sample_type_plot_data.rds")
    files <- files[file.exists(files)]
    files <- files[order(file.info(files)$mtime,decreasing=TRUE)]
    for(file in files) {
      cached <- tryCatch(readRDS(file),error=function(e) NULL)
      if(is.null(cached) || !identical(cached$MC,mc) ||
          !identical(cached$input_md5,fingerprint) ||
          !identical(cached$distance_spec,distance_spec)) next
      ord <- cached$ordination
      d <- ord$data
      if(is.null(d) || !all(c(needed,"NMDS1","NMDS2") %in% names(d)) ||
          !identical(as.character(data$Repeat_group_ID),as.character(d$Repeat_group_ID))) next
      same <- vapply(c("Sample_type","Host_type","TS",mc),function(column)
        if(column %in% mc) isTRUE(all.equal(data[[column]],d[[column]],tolerance=1e-12)) else
          identical(as.character(data[[column]]),as.character(d[[column]])),logical(1))
      if(!all(same) || !identical(as.integer(ord$fit$ndim),2L) ||
          !is.finite(ord$fit$stress) || !all(is.finite(as.matrix(d[,c("NMDS1","NMDS2")])))) next
      ord$data <- cbind(data,d[,c("NMDS1","NMDS2"),drop=FALSE])
      ord$cache_source <- file
      return(ord)
    }
    NULL
  }
  fit_coordinates <- function(data) {
    matrix <- as.matrix(data[,mc]);rownames(matrix)<-data$Repeat_group_ID
    distance <- vegan::vegdist(sqrt(matrix),method="bray")
    cores <- parallel::detectCores(logical=FALSE)
    if(!is.finite(cores)) cores <- 1L
    n_workers <- max(1L,min(workers,as.integer(cores)))
    cluster <- NULL
    if(n_workers>1L) {
      cluster <- parallel::makePSOCKcluster(n_workers)
      on.exit(parallel::stopCluster(cluster),add=TRUE)
      parallel::clusterEvalQ(cluster,library(vegan))
      parallel::clusterSetRNGStream(cluster,iseed=random_seed)
    }
    set.seed(random_seed)
    message("Fitting 2D NMDS only: ",nrow(data)," groups; ",n_workers," workers.")
    fit <- vegan::metaMDS(distance,k=2,try=30,trymax=trymax,trace=FALSE,
      parallel=if(is.null(cluster)) 1L else cluster)
    xy <- as.data.frame(vegan::scores(fit,display="sites",choices=1:2))
    stopifnot(identical(rownames(xy),as.character(data$Repeat_group_ID)),
      all(is.finite(as.matrix(xy))),is.finite(fit$stress))
    list(fit=fit,data=cbind(data,xy),workers=n_workers)
  }
  ordination <- if(reuse_saved_coordinates) find_coordinates(data) else NULL
  if(is.null(ordination)) ordination <- fit_coordinates(data) else
    message("Reusing 2D NMDS coordinates: ",ordination$cache_source)
  if(is.null(ordination$fit$converged) || ordination$fit$converged<1L)
    message("The best NMDS solution was not repeated; inspect saved convergence diagnostics.")

  points <- ordination$data
  split_rows <- split(seq_len(nrow(points)),points$Sample_type)
  eligible <- vapply(split_rows,function(index) {
    xy <- as.matrix(points[index,c("NMDS1","NMDS2"),drop=FALSE])
    nrow(xy)>=4L && qr(scale(xy,center=TRUE,scale=FALSE))$rank==2L
  },logical(1))
  ellipse_points <- points[unlist(split_rows[eligible],use.names=FALSE),,drop=FALSE]
  legend_labels <- setNames(paste0(groups," (n = ",as.integer(table(points$Sample_type)),")"),groups)
  nmds_plot <- ggplot(points,aes(NMDS1,NMDS2))
  if(show_ellipses && nrow(ellipse_points)) nmds_plot <- nmds_plot +
    stat_ellipse(data=ellipse_points,aes(colour=Sample_type,group=Sample_type),
      type="norm",level=ellipse_level,linewidth=0.4,alpha=0.65,show.legend=FALSE)
  nmds_plot <- nmds_plot +
    geom_point(aes(fill=Sample_type),shape=21,colour="white",size=2.6,stroke=0.25) +
    scale_fill_manual(values=palette,limits=groups,labels=legend_labels,name="Sample type") +
    scale_colour_manual(values=palette,guide="none") +
    guides(fill=guide_legend(override.aes=list(shape=21,colour="white",size=3.4))) +
    coord_equal() + labs(title=NULL,subtitle=NULL,x="NMDS1",y="NMDS2",
      caption=sprintf("Stress = %.3f; n = %d; R\u00b2 = %.3f; p = %s",
        ordination$fit$stress,nrow(points),result$R2_raw,p_text(result$P_raw))) + pv_theme() +
    theme(text=element_text(family=font_family))

  # Pairwise statistics are copied directly from the saved three-pair BH family.
  pairs$Column <- groups[pmin(match(pairs$Group1,groups),match(pairs$Group2,groups))]
  pairs$Row <- groups[pmax(match(pairs$Group1,groups),match(pairs$Group2,groups))]
  valid <- pairs[pairs$Status=="ok" & is.finite(pairs$R2_raw) & is.finite(pairs$P_BH_raw),,drop=FALSE]
  unavailable <- pairs[!pairs$Pair_ID %in% valid$Pair_ID,,drop=FALSE]
  upper <- if(nrow(valid)) max(0.1,ceiling(max(valid$R2_raw)*20)/20) else 0.1
  lower <- if(nrow(valid)) min(0,floor(min(valid$R2_raw)*20)/20) else 0
  limits <- c(lower,upper)
  valid$Text_colour <- if(nrow(valid)) pv_text_colour(valid$R2_raw,limits) else character()
  valid$R2_label <- sprintf("%.3f",valid$R2_raw)
  valid$Q_label <- paste0("q = ",p_text(valid$P_BH_raw))
  valid$Stars <- ifelse(valid$P_BH_raw<=0.001,"***",ifelse(valid$P_BH_raw<=0.01,"**",
    ifelse(valid$P_BH_raw<=0.05,"*","")))
  valid$Dispersion <- ifelse((is.finite(valid$P_disp_BH_raw) & valid$P_disp_BH_raw<=0.05) |
    (is.finite(valid$P_disp_BH_lingoes) & valid$P_disp_BH_lingoes<=0.05),"\u2020","")
  heatmap <- ggplot() +
    geom_tile(data=data.frame(Column=groups,Row=groups),aes(Column,Row),fill="#F1F1F1",colour="white") +
    geom_tile(data=valid,aes(Column,Row,fill=R2_raw),colour="white") +
    geom_tile(data=unavailable,aes(Column,Row),fill="#DADADA",colour="white") +
    geom_text(data=valid,aes(Column,Row,label=R2_label,colour=Text_colour),size=4.5,
      position=position_nudge(y=0.02),family=font_family,show.legend=FALSE) +
    geom_text(data=valid,aes(Column,Row,label=Q_label,colour=Text_colour),size=3.5,
      position=position_nudge(y=-0.25),family=font_family,show.legend=FALSE) +
    geom_text(data=valid,aes(Column,Row,label=Stars,colour=Text_colour),size=3.5,
      position=position_nudge(y=0.30),family=font_family,show.legend=FALSE) +
    geom_text(data=valid,aes(Column,Row,label=Dispersion,colour=Text_colour),size=3.5,
      position=position_nudge(x=0.34,y=0.30),family=font_family,show.legend=FALSE) +
    geom_text(data=unavailable,aes(Column,Row),label="NA",family=font_family) +
    scale_colour_identity() + scale_fill_gradientn(colours=pv_continuous,limits=limits,
      name=expression(PERMANOVA~R^2),guide=guide_colourbar(title.position="top",
        barwidth=grid::unit(65,"mm"),barheight=grid::unit(6,"mm"))) +
    scale_x_discrete(limits=groups,drop=FALSE) +
    scale_y_discrete(limits=rev(groups),drop=FALSE) + coord_fixed() +
    labs(title=NULL,subtitle=NULL,x="Sample type",y="Sample type",
      caption=paste0("R\u00b2: marginal PERMANOVA effect; q: BH-adjusted P across 3 pairs.\n",
        "* q \u2264 0.05; ** q \u2264 0.01; *** q \u2264 0.001.\n",
        "\u2020: BH-adjusted PERMDISP P \u2264 0.05 (raw or Lingoes).\n",
        "Each pair adjusts for TS and uses shared batches; sample sizes may differ.")) +
    pv_theme() + theme(text=element_text(family=font_family),
      legend.position="top",legend.justification="right")

  # File creation starts here, after validation, coordinates and plots are ready.
  if(dir.exists(output_dir)) stop("Output folder already exists: ",output_dir)
  if(!dir.create(output_dir,recursive=TRUE)) stop("Cannot create output directory.")
  plots <- list("01_sample_type_NMDS"=nmds_plot,"02_sample_type_pairwise"=heatmap)
  png_device <- if(requireNamespace("ragg",quietly=TRUE)) ragg::agg_png else "png"
  for(name in names(plots)) {
    width <- if(grepl("NMDS",name)) 10 else 8
    height <- if(grepl("NMDS",name)) 8 else 8
    ggsave(file.path(output_dir,paste0(name,".png")),plots[[name]],device=png_device,
      width=width,height=height,dpi=dpi,bg="white")
    if(save_pdf) {
      cairo <- isTRUE(capabilities("cairo"))
      pdf_plot <- if(cairo) plots[[name]] else plots[[name]] + theme(text=element_text(family="sans"))
      ggsave(file.path(output_dir,paste0(name,".pdf")),pdf_plot,
        device=if(cairo) grDevices::cairo_pdf else grDevices::pdf,
        width=width,height=height,bg="white")
    }
  }
  saveRDS(list(ordination=ordination,MC=mc,input_md5=fingerprint,distance_spec=distance_spec,
    result=result,pairs=pairs,source=result_file),file.path(output_dir,"sample_type_plot_data.rds"))
  write.csv(points,file.path(output_dir,"sample_type_NMDS_coordinates.csv"),row.names=FALSE,fileEncoding="UTF-8")
  write.csv(pairs,file.path(output_dir,"sample_type_pairwise_used.csv"),row.names=FALSE,fileEncoding="UTF-8")
  write.csv(data.frame(Dimensions=2,N_units=nrow(points),Stress=ordination$fit$stress,
    Random_starts=ordination$fit$tries,Repeated_best=ordination$fit$converged,
    Stop_code=ordination$fit$icause),file.path(output_dir,"NMDS_diagnostics.csv"),row.names=FALSE)
  captions <- c("All labels are English; no figure titles or subtitles.",
    "One point represents one Repeat_group mean. HOST pools all sampled biological host types.",
    "Eight original-denominator MC proportions are square-root transformed before Bray-Curtis NMDS.",
    "All observed groups are retained; absent sample-type x TS combinations are not filled with zeros.",
    "NMDS shows unadjusted distances. The saved omnibus PERMANOVA adjusts for TS and permutes within TS.",
    "The omnibus caption reports the raw P of a single overall test, not BH-adjusted P.",
    "Pairwise q values are the saved BH adjustment over three planned comparisons; no recalculation is performed.",
    "Pairwise comparisons use shared TS. HOST comparisons with environmental samples exclude TS1 only for those pairs.",
    sprintf("Omnibus PERMDISP P (raw / Lingoes, unadjusted): %.4g / %.4g.",result$P_disp_raw,result$P_disp_lingoes),
    if(show_ellipses) "Thin outlines are 95% data ellipses (normal approximation), not confidence regions for group centres." else "Ellipses are omitted.",
    if(any(!eligible)) paste("No ellipse (<4 points or collinear):",paste(names(eligible)[!eligible],collapse=", ")),
    if(!is.null(ordination$cache_source)) paste("Reused coordinates:",ordination$cache_source) else "New 2D NMDS coordinates fitted; PERMANOVA was not rerun.",
    paste("Saved result:",result_file),paste("Source RDS MD5:",unname(tools::md5sum(result_file))))
  writeLines(enc2utf8(captions),file.path(output_dir,"figure_captions.txt"),useBytes=TRUE)
  message("Saved sample-type PNG/PDF figures: ",output_dir)
})()
