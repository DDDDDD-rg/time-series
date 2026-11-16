# Export saved-model figures only. No fitting, simulation, PERMANOVA, NMDS or PPT editing.
# The old nutrient pages are copied unchanged; temperature gets explicit review annotations.
(function() {
  project <- getOption("mc19.project_dir","C:/Users/xzy/Desktop/CUHK_Study/GNBF6010_Research_Projects/Research/Msc main project-time series")
  output <- getOption("mc19.review.output",file.path(project,"figures","temperature_review_v1"))
  if(dir.exists(output) && !isTRUE(getOption("mc19.review.redraw",FALSE)))
    stop("Review output already exists; choose a new directory or explicitly set mc19.review.redraw=TRUE.")
  if(dir.exists(output) && !file.exists(file.path(output,"supporting_files","export_status.txt")))
    stop("Only a completed review export can be redrawn.")
  prior <- options(mc19.pipeline="TS7_temperature",mc19.project_dir=project)
  on.exit(options(prior),add=TRUE)
  # Catch an accidental call to any statistical analysis during export.
  hooks <- list(c("glmmTMB","glmmTMB"),c("vegan","adonis2"),c("vegan","metaMDS"),c("stats","simulate"))
  for(h in hooks) suppressMessages(trace(h[[2]],where=asNamespace(h[[1]]),
    tracer=quote(stop("Statistical analysis is forbidden in this figure export")),print=FALSE))
  on.exit(for(h in hooks) untrace(h[[2]],where=asNamespace(h[[1]])),add=TRUE)
  dir.create(file.path(output,"supporting_files"),recursive=TRUE,showWarnings=FALSE)
  writeLines("running; do not package until export completes",file.path(output,"supporting_files","export_status.txt"))
  e <- new.env(parent=globalenv())
  for(f in c("glmm_MC19_pilot_helpers.R","glmm_MC19_wald_helpers.R","glmm_MC19_wald_plot_helpers.R",
    "mc19_talk_plot_helpers.R","mc19_two_part_helpers.R")) sys.source(file.path(project,"scripts","GLMM",f),e)
  # Track every existing result file and PPT before drawing. Export never publishes pointers.
  sources <- unlist(lapply(c("coral","ecosystem"),function(scope) unlist(lapply(c("GLMM","PERMANOVA","two_part"),function(kind)
    list.files(file.path(project,"results",paste0("MC19_",kind,"_TS7_temp_",scope,"_run001")),recursive=TRUE,full.names=TRUE)))))
  sources <- c(sources,list.files(file.path(project,"presentations"),pattern="[.]pptx$",recursive=TRUE,full.names=TRUE))
  before <- tools::md5sum(sources)
  manifest <- list();originals <- list()
  revise <- function(fun,from,to) {
    code <- paste(deparse(body(fun),width.cutoff=500L),collapse="\n")
    if(!grepl(from,code,fixed=TRUE)) stop("Expected plot annotation point is missing: ",from)
    body(fun) <- parse(text=gsub(from,to,code,fixed=TRUE))[[1]];fun
  }
  # These changes live only in this export environment; generic plot scripts stay unchanged.
  e$mc19_main_curves <- revise(e$mc19_main_curves,
    'sprintf("Beta = %.2f; q = %s",','sprintf("Beta = %.2f; q = %s [R]",')
  e$mc19_main_curves <- revise(e$mc19_main_curves,
    'nutrient or a DIN component below RL (L/2).','nutrient below RL (L/2).')
  e$mc19_main_curves <- revise(e$mc19_main_curves,
    'Main panels: full range; right-hand panels: y-axis zoom (same x-axis).',
    'Full ranges retained; right panels: zoom. [R]: temperature model requires diagnostic review.')
  e$mc19_heatmap_forest <- revise(e$mc19_heatmap_forest,
    'Filled points / *: saved BH q < 0.05 (%d tests).',
    'Filled points / *: saved BH q < 0.05 (%d tests).\n[R]: diagnostic review shown for Temperature only; see review PDF.')
  # Reuse only drawing definitions from the actual compact-PPT renderer, never its analysis/export body.
  style <- new.env(parent=e);style$e <- e;style$W <- 12.7;style$H <- 6;style$evidence <- list()
  template <- file.path(project,"archive","ppt_readability_v2","render_compact.R")
  for(x in parse(template)) if(is.call(x) && identical(x[[1]],as.name("<-")) && is.symbol(x[[2]]) &&
    as.character(x[[2]]) %in% c("canvas","compact_rationale","draw_heat","draw_curves")) eval(x,style)
  style$draw_heat <- revise(style$draw_heat,"heatwidth/4","heatwidth/length(env)")
  style$draw_heat <- revise(style$draw_heat,"heatwidth/8","heatwidth/(2*length(env))")
  style$draw_heat <- revise(style$draw_heat,'c("NO3", "NO2", "NH4", "DIP")[j]',
    'if(env[j]=="SST_C") "Temp\\n[R]" else unname(cfg$short_names[env[j]])')
  style$draw_heat <- revise(style$draw_heat,'DIP_mgP_L = "DIP"','DIP_mgP_L = "DIP", SST_C = "Temp [R]"')
  style$draw_heat <- revise(style$draw_heat,"5.81","5.69")
  style$draw_heat <- revise(style$draw_heat,
    ' | Beta per nutrient SD; bars: full 95% Wald CI; * / filled points: BH q < .05 (76 tests).',
    ' | Beta per TS SD; full 95% CI; BH q < .05 (95 tests). [R]: Temperature diagnostic review.')
  style$draw_heat <- revise(style$draw_heat,'text(sx, yy, rationale, adj = 0, cex = 0.65)',
    'text(sx, yy, rationale, adj = 0, cex = 0.65, font = 1)')
  style$draw_curves <- revise(style$draw_curves,
    'mainwidth, 3.75 - offset, 5.08 - offset','mainwidth, 1.35 - offset, 5.08 - offset')
  style$draw_curves <- revise(style$draw_curves,
    '3.96 - offset, 4.96 - offset','2.45 - offset, 4.96 - offset')
  style$draw_curves <- revise(style$draw_curves,
    'Open points: nutrient below LOD (L/2). Dashed lines: NO2-N LOD = 0.002 mg N/L; DIP LOD = 0.002 mg P/L.',
    '[R]: saved temperature GLMM has diagnostic review flags; see temperature review PDF.')
  ppt_save <- function(stem,draw,width=12.7,height=6) {
    grDevices::pdf(paste0(stem,".pdf"),width=width,height=height,pointsize=20,useDingbats=FALSE)
    tryCatch(draw(),finally=grDevices::dev.off())
    grDevices::png(paste0(stem,".png"),width=width,height=height,units="in",res=240,pointsize=20,type="cairo",bg="white")
    tryCatch(draw(),finally=grDevices::dev.off())
  }
  original_labels <- e$mc19_tp_forest_labels
  e$mc19_tp_forest_labels <- function(tab) {
    rows <- original_labels(tab)
    rows$label <- paste0(sub("MC2s1","MC2S1",rows$label,fixed=TRUE)," [R]")
    rows
  }
  e$mc19_tp_curve_page <- revise(e$mc19_tp_curve_page,
    'paste(mc, "/", axis_label)','paste(mc, "/", axis_label, "[R]")')
  e$mc19_tp_curve_page <- revise(e$mc19_tp_curve_page,
    'Positive abundance = MC reads / total Ruegeria reads; shown as %. Detection observations overlap at 0 and 1.',
    '[R]: parent temperature GLMM flagged; see review PDF. Positive abundance = MC / total Ruegeria reads (%).')
  copy_old <- function(src,dest) {
    stopifnot(file.exists(src),file.copy(src,dest,overwrite=TRUE),identical(unname(tools::md5sum(src)),unname(tools::md5sum(dest))))
    originals[[length(originals)+1L]] <<- data.frame(source=src,copy=dest,md5=unname(tools::md5sum(src)))
  }
  for(scope in c("coral","ecosystem")) {
    message("Exporting ",scope)
    folder <- file.path(output,scope);dir.create(folder,showWarnings=FALSE)
    cfg <- e$mc19_wald_config(project,scope)
    run <- e$mc19_wald_read_run(cfg,file.path(project,"results",paste0("MC19_GLMM_TS7_temp_",scope,"_run001")))
    r <- run$results
    diag <- read.csv(file.path(run$support,"diagnostics","diagnostic_overview.csv"))
    temp <- diag[diag$environment=="SST_C",]
    stopifnot(nrow(r)==95L,!any(r$environment=="DIN_mgN_L"),nrow(temp)==19L,all(temp$status=="review_flags"))
    cfg$short_names["SST_C"] <- "Temp [R]"
    cfg$axis_labels["SST_C"] <- "Temperature (deg C)"
    ppt_save(file.path(folder,"02_GLMM_heatmap_forest"),function() style$draw_heat(r,cfg))
    old_glmm <- file.path(project,"figures",if(scope=="coral")"glmm_coral_v1" else "glmm_all_hosts_v1")
    reference <- file.path(project,"data","archive","temperature_review_style_reference")
    page_images <- if(scope=="coral") 10:12 else 13:16
    for(i in seq_along(page_images)) copy_old(file.path(reference,paste0("image",page_images[i],".png")),
      file.path(folder,sprintf("03_GLMM_nutrients_PPT_page_%02d.png",i)))
    entries <- lapply(run$manifest$model_files,readRDS)
    temperature <- r[r$environment=="SST_C",]
    selected <- e$mc19_significant_curve_rows(temperature,cfg)
    selected$page <- seq_len(nrow(selected));selected$row <- selected$column <- selected$panel <- 1L
    grids <- list()
    for(key in selected$key) {
      g <- e$mc19_predict_grid(entries[[key]],run$manifest$scaling[run$manifest$scaling$environment=="SST_C",])
      g$environment <- "SST_C";grids[[key]] <- g
    }
    # One temperature MC per page, with a wide full-range panel and separate zoom.
    e$mc19_curve_panel_boxes <- function(row,column,width,height,rows=1L) {
      top <- 1-1/height
      list(main=c(1/width,12/width,top-6.4/height,top),
        zoom=c(13/width,18.7/width,top-3.1/height,top-.12/height))
    }
    env <- read.csv(file.path(run$support,"input_environment_TS.csv"))
    grDevices::pdf(file.path(folder,"04_GLMM_temperature_curves.pdf"),width=12.7,height=6,pointsize=20,useDingbats=FALSE)
    for(i in seq_len(nrow(selected))) style$draw_curves(temperature,entries,grids,cfg,env,selected[i,,drop=FALSE])
    grDevices::dev.off()
    for(i in seq_len(nrow(selected))) {
      grDevices::png(file.path(folder,sprintf("04_GLMM_temperature_curves_page_%02d.png",i)),
        width=12.7,height=6,units="in",res=240,pointsize=20,type="cairo",bg="white")
      style$draw_curves(temperature,entries,grids,cfg,env,selected[i,,drop=FALSE]);grDevices::dev.off()
    }
    e$mc19_talk_save <- function(p,stem,...) {
      p <- p+ggplot2::scale_x_date(breaks=as.Date(env$Water_quality_date),labels=paste(env$TS,format(as.Date(env$Water_quality_date),"%y-%m"),sep="\n"))+
        ggplot2::theme_bw(base_size=16)+ggplot2::theme(text=ggplot2::element_text(face="plain"),
          panel.grid=ggplot2::element_blank(),strip.background=ggplot2::element_rect(fill="white"),
          strip.text=ggplot2::element_text(face="bold",size=17),legend.position="top",legend.text=ggplot2::element_text(size=13),
          axis.text=ggplot2::element_text(size=14,colour="black"),axis.text.x=ggplot2::element_text(angle=0,hjust=.5),
          plot.caption=ggplot2::element_text(size=10,hjust=0,lineheight=1.05),plot.margin=ggplot2::margin(3,7,3,3))
      for(ext in c("pdf","png")) ggplot2::ggsave(paste0(stem,".",ext),p,width=12.7,height=7.5,dpi=240,bg="white")
    }
    e$mc19_np_plot(env,folder)
    for(ext in c("pdf","png")) file.rename(file.path(folder,paste0("NP_by_TS.",ext)),file.path(folder,paste0("05_environment_by_TS.",ext)))
    file.rename(file.path(folder,"NP_by_TS_data.csv"),file.path(output,"supporting_files",paste0(scope,"_environment_by_TS.csv")))
    old_tp <- file.path(project,"figures",if(scope=="coral")"two_part_coral_v1" else "two_part_all_hosts_v1")
    copy_old(file.path(reference,if(scope=="coral")"image8.png" else "image9.png"),
      file.path(folder,"06_two_part_nutrients_PPT_original.png"))
    copy_old(file.path(old_tp,"two_part_curves.pdf"),file.path(folder,"07_two_part_nutrients_curves_original.pdf"))
    tp_path <- file.path(project,"results",paste0("MC19_two_part_TS7_temp_",scope,"_run001"))
    tp <- e$mc19_tp_validate(tp_path)
    tab <- tp$tab[tp$tab$environment=="SST_C",]
    pred <- tp$pred[tp$pred$key %in% tab$key,];obs <- tp$obs[tp$obs$key %in% tab$key,]
    stopifnot(nrow(tab)==2*nrow(selected),setequal(unique(tab$MC_ID),selected$MC_ID))
    forest <- e$mc19_tp_forest(tab)
    forest <- forest+ggplot2::labs(caption=paste0("* det q < .05; # pos q < .05; *# both. Full 95% Wald CI; BH: ",unique(tab$q_family_size),
      " component tests.\n[R]: parent temperature GLMM flagged; see review PDF."))+
      ggplot2::theme(axis.text.y=ggplot2::element_text(size=21,face="bold"),axis.text.x=ggplot2::element_text(size=18,face="plain"),
        axis.title.x=ggplot2::element_text(size=19),strip.text=ggplot2::element_text(size=21),
        plot.caption=ggplot2::element_text(size=11,hjust=0,lineheight=1.1),plot.margin=ggplot2::margin(4,12,5,5))
    for(ext in c("pdf","png")) ggplot2::ggsave(file.path(folder,paste0("08_two_part_temperature_effects.",ext)),forest,
      width=12.7,height=if(scope=="coral")4 else 4.5,dpi=240,bg="white")
    grDevices::pdf(file.path(folder,"09_two_part_temperature_curves.pdf"),width=20,height=11.25,useDingbats=FALSE)
    ranges <- lapply(unique(tab$key),function(key)e$mc19_tp_curve_page(tab,pred,obs,key));grDevices::dev.off()
    for(i in seq_along(unique(tab$key))) {
      grDevices::png(file.path(folder,sprintf("09_two_part_temperature_curves_page_%02d.png",i)),
        width=20,height=11.25,units="in",res=180,type="cairo",bg="white")
      e$mc19_tp_curve_page(tab,pred,obs,unique(tab$key)[i]);grDevices::dev.off()
    }
    obs$MC_ID <- sub("__.*$","",obs$key)
    obs$component <- factor(obs$component,levels=c("detection","positive"),labels=c("Detection","Positive proportion"))
    residual <- ggplot2::ggplot(obs,ggplot2::aes(fitted,pearson))+ggplot2::geom_point(alpha=.4)+
      ggplot2::geom_hline(yintercept=0,linetype=3,colour="grey50")+
      ggplot2::facet_grid(rows=ggplot2::vars(MC_ID),cols=ggplot2::vars(component),scales="free")+
      ggplot2::labs(x="Fitted value",y="Pearson residual",caption="Temperature components only. Descriptive residual display; no new model or test.")+e$mc19_tp_theme()
    for(ext in c("pdf","png")) ggplot2::ggsave(file.path(folder,paste0("10_temperature_component_residuals.",ext)),residual,
      width=16,height=if(scope=="coral")8 else 11,dpi=180,bg="white")
    grDevices::pdf(file.path(folder,"11_temperature_model_review.pdf"),width=12,height=9,useDingbats=FALSE)
    for(key in selected$key) {
      entry <- entries[[key]];diagnosis <- readRDS(file.path(run$support,"diagnostics",paste0(key,".rds")))
      e$mc19_wald_diagnostic_page(entry,diagnosis,cfg)
      e$mc19_plot_simulation_summary(entry,diagnosis)
    }
    grDevices::dev.off()
    write.csv(temp,file.path(output,"supporting_files",paste0(scope,"_temperature_review_all_MC.csv")),row.names=FALSE)
    write.csv(selected,file.path(output,"supporting_files",paste0(scope,"_temperature_significant_MC.csv")),row.names=FALSE)
    write.csv(tab,file.path(output,"supporting_files",paste0(scope,"_temperature_two_part_current.csv")),row.names=FALSE)
    write.csv(do.call(rbind,ranges),file.path(output,"supporting_files",paste0(scope,"_temperature_curve_ranges.csv")),row.names=FALSE)
    write.csv(r,file.path(output,"supporting_files",paste0(scope,"_GLMM_current_95.csv")),row.names=FALSE)
    manifest[[scope]] <- data.frame(scope=scope,GLMM=run$run_dir,two_part=tp_path,
      temperature_MC=paste(selected$MC_ID,collapse=";"),GLMM_BH_family=95L,two_part_BH_family=unique(tab$q_family_size))
  }
  stopifnot(identical(before,tools::md5sum(sources)))
  write.csv(do.call(rbind,manifest),file.path(output,"supporting_files","analysis_sources.csv"),row.names=FALSE)
  write.csv(do.call(rbind,originals),file.path(output,"supporting_files","unchanged_original_figures.csv"),row.names=FALSE)
  write.csv(data.frame(path=names(before),md5=unname(before)),file.path(output,"supporting_files","protected_result_and_PPT_md5.csv"),row.names=FALSE)
  writeLines("complete; statistical result and PPT hashes unchanged; no fitting/simulation/PERMANOVA/NMDS",file.path(output,"supporting_files","export_status.txt"))
  message("Review export complete: ",output)
})()
