# Draw current six-indicator results in the actual PPT's graphic style.
# Saved fits, diagnostics and validated NMDS coordinates only; no statistical rerun.
(function() {
  project <- getOption("mc19.project_dir",normalizePath("Research/Msc main project-time series",winslash="/"))
  output <- getOption("mc19.env6.figure_dir",file.path(project,"figures","TS7_env6_PPT_v1"))
  if(dir.exists(output)) stop("Output exists; choose a new mc19.env6.figure_dir.")
  prior <- options(mc19.pipeline="TS7_env6",mc19.project_dir=project);on.exit(options(prior),add=TRUE)
  hooks <- list(c("glmmTMB","glmmTMB"),c("vegan","adonis2"),c("vegan","metaMDS"),c("stats","simulate"))
  for(h in hooks) suppressMessages(trace(h[2],where=asNamespace(h[1]),tracer=quote(stop("Statistical computation forbidden in saved-result export")),print=FALSE))
  on.exit(for(h in hooks) untrace(h[2],where=asNamespace(h[1])),add=TRUE)
  dir.create(file.path(output,"supporting_files"),recursive=TRUE)
  protected <- unlist(lapply(c("coral","ecosystem"),function(s) unlist(lapply(c("GLMM","PERMANOVA","two_part"),function(k)
    list.files(file.path(project,"results",paste0("MC19_",k,"_TS7_env6_",s,"_run001")),recursive=TRUE,full.names=TRUE)))))
  protected <- c(protected,list.files(file.path(project,"presentations"),pattern="[.]pptx$",full.names=TRUE))
  protected <- protected[!startsWith(basename(protected),"~$")]
  hashes <- tools::md5sum(protected)
  e <- new.env(parent=globalenv())
  for(f in c("glmm_MC19_pilot_helpers.R","glmm_MC19_wald_helpers.R","glmm_MC19_wald_plot_helpers.R","mc19_talk_plot_helpers.R","mc19_two_part_helpers.R"))
    sys.source(file.path(project,"scripts","GLMM",f),e)
  revise <- function(fun,a,b) {
    code <- paste(deparse(body(fun),width.cutoff=500L),collapse="\n")
    if(!grepl(a,code,fixed=TRUE)) stop("Template mismatch: ",a)
    body(fun) <- parse(text=gsub(a,b,code,fixed=TRUE))[[1]];fun
  }
  style <- new.env(parent=e);style$e <- e;style$W <- 12.7;style$H <- 6;style$evidence <- list()
  for(x in parse(file.path(project,"archive","ppt_readability_v2","render_compact.R")))
    if(is.call(x)&&identical(x[[1]],as.name("<-"))&&is.symbol(x[[2]])&&as.character(x[[2]]) %in% c("canvas","compact_rationale","draw_heat","draw_curves"))eval(x,style)
  style$draw_heat <- revise(style$draw_heat,"heatwidth/4","heatwidth/length(env)")
  style$draw_heat <- revise(style$draw_heat,"heatwidth/8","heatwidth/(2*length(env))")
  style$draw_heat <- revise(style$draw_heat,'c("NO3", "NO2", "NH4", "DIP")[j]',
    'if(env[j]=="SST_C") "Temp" else unname(cfg$short_names[env[j]])')
  style$draw_heat <- revise(style$draw_heat,'DIP_mgP_L = "DIP"','DIP_mgP_L = "DIP", DIN_mgN_L="DIN", SST_C="Temp"')
  style$draw_heat <- revise(style$draw_heat,"5.81","5.69")
  style$draw_heat <- revise(style$draw_heat,'text(sx, yy, rationale, adj = 0, cex = 0.65)',
    'text(sx, yy, rationale, adj = 0, cex = 0.65, font=1)')
  style$draw_heat <- revise(style$draw_heat,' | Beta per nutrient SD; bars: full 95% Wald CI; * / filled points: BH q < .05 (76 tests).',
    ' | Beta per TS SD; full 95% CI; BH q < .05 (114 tests). ')
  style$draw_curves <- revise(style$draw_curves,'cw <- (W - 0.2)/2','cw <- (W - 0.2)/max(2L,max(panels$column))')
  style$draw_curves <- revise(style$draw_curves,'Open points: nutrient below LOD (L/2). Dashed lines: NO2-N LOD = 0.002 mg N/L; DIP LOD = 0.002 mg P/L.',
    'Open: nutrient/component below RL (L/2); dashed: NO2/DIP RL = 0.002. ')
  labels <- c(Nitrate_mgN_L="NO3-N (mg N/L)",Nitrite_mgN_L="NO2-N (mg N/L)",Ammonia_mgN_L="Ammonia-N (mg N/L)",DIN_mgN_L="DIN (mg N/L)",DIP_mgP_L="DIP (mg P/L)",SST_C="Temperature (deg C)")
  savebase <- function(stem,draw,w=12.7,h=6) {
    grDevices::pdf(paste0(stem,".pdf"),width=w,height=h,pointsize=20,useDingbats=FALSE)
    tryCatch(draw(),finally=grDevices::dev.off())
    grDevices::png(paste0(stem,".png"),width=w,height=h,units="in",res=240,pointsize=20,type="cairo",bg="white")
    tryCatch(draw(),finally=grDevices::dev.off())
  }
  savegg <- function(p,stem,w=12.7,h=6) for(ext in c("pdf","png"))
    ggplot2::ggsave(paste0(stem,".",ext),p,width=w,height=h,dpi=240,bg="white")
  oldlabels <- e$mc19_tp_forest_labels
  e$mc19_tp_forest_labels <- function(tab) {
    z <- oldlabels(tab);i <- z$environment=="SST_C";z$label[i] <- paste0(z$label[i],"");z
  }
  manifest <- list();coordinate_audit <- list()
  for(scope in c("coral","ecosystem")) {
    message("Exporting current results: ",scope)
    folder <- file.path(output,scope);aux <- file.path(folder,"auxiliary")
    dir.create(aux,recursive=TRUE)
    cfg <- e$mc19_wald_config(project,scope)
    run <- e$mc19_wald_read_run(cfg,file.path(project,"results",paste0("MC19_GLMM_TS7_env6_",scope,"_run001")))
    r <- run$results;env <- read.csv(file.path(run$support,"input_environment_TS.csv"))
    cfg$short_names["SST_C"] <- "Temp"
    style$W <- 12.7
    savebase(file.path(folder,"02_GLMM_heatmap_forest"),function()style$draw_heat(r,cfg))
    entries <- lapply(run$manifest$model_files,readRDS);selected <- e$mc19_significant_curve_rows(r,cfg)
    grids <- setNames(lapply(selected$key,function(k) {
      a <- selected[selected$key==k,];g <- e$mc19_predict_grid(entries[[k]],run$manifest$scaling[run$manifest$scaling$environment==a$environment,]);g$environment <- a$environment;g
    }),selected$key)
    pages <- split(selected,selected$page)
    # A larger page only where one MC has three indicators; same font and column dimensions.
    for(pg in names(pages)) {
      panels <- pages[[pg]];style$W <- if(max(panels$column)>2)18.95 else 12.7
      savebase(file.path(folder,sprintf("03_GLMM_curves_page_%02d",as.integer(pg))),
        function()style$draw_curves(r,entries,grids,cfg,env,panels),style$W,6)
    }
    # The packaging script combines these pages with their native sizes.
    style$W <- 12.7
    e$mc19_talk_save <- function(p,stem,...) {
      p <- p+ggplot2::scale_x_date(breaks=as.Date(env$Water_quality_date),labels=paste(env$TS,format(as.Date(env$Water_quality_date),"%y-%m"),sep="\n"))+
        ggplot2::theme_bw(base_size=16)+ggplot2::theme(text=ggplot2::element_text(face="plain"),panel.grid=ggplot2::element_blank(),
          strip.background=ggplot2::element_rect(fill="white"),strip.text=ggplot2::element_text(face="bold",size=17),legend.position="top",
          axis.text=ggplot2::element_text(size=14,colour="black"),axis.text.x=ggplot2::element_text(angle=0,hjust=.5),
          plot.caption=ggplot2::element_text(size=10,hjust=0),plot.margin=ggplot2::margin(3,7,3,3))
      savegg(p,file.path(folder,"01_environment_by_TS"),12.7,7.5)
    }
    e$mc19_np_plot(env,aux)
    tp_path <- file.path(project,"results",paste0("MC19_two_part_TS7_env6_",scope,"_run001"))
    tp <- e$mc19_tp_validate(tp_path);tab <- tp$tab
    stopifnot(setequal(unique(tab$key),paste(selected$MC_ID,selected$environment,sep="__")),nrow(r)==114L)
    forest <- e$mc19_tp_forest(tab)+ggplot2::labs(caption=paste0(tools::toTitleCase(scope)," | * det q < .05; # pos q < .05; *# both. Full 95% CI. BH: ",unique(tab$q_family_size)," components.\n"))+
      ggplot2::theme(axis.text.y=ggplot2::element_text(size=if(scope=="coral")19 else 16),axis.text.x=ggplot2::element_text(size=18),
        axis.title.x=ggplot2::element_text(size=19),strip.text=ggplot2::element_text(size=21),plot.caption=ggplot2::element_text(size=11,hjust=0),plot.margin=ggplot2::margin(4,12,5,5))
    built <- ggplot2::ggplot_build(forest)
    stopifnot(nrow(built$data[[2]])==nrow(tab),isTRUE(all.equal(sort(built$data[[2]]$xmin),sort(tab$CI95_low))))
    savegg(forest,file.path(folder,"04_two_part_forest"),12.7,if(scope=="coral")6 else 7.5)
    grDevices::pdf(file.path(aux,"two_part_component_curves.pdf"),width=20,height=11.25,useDingbats=FALSE)
    ranges <- lapply(unique(tab$key),function(k)e$mc19_tp_curve_page(tab,tp$pred,tp$obs,k));grDevices::dev.off()
    write.csv(do.call(rbind,ranges),file.path(aux,"two_part_curve_ranges.csv"),row.names=FALSE)
    grDevices::pdf(file.path(aux,"GLMM_significant_model_diagnostics.pdf"),width=12,height=9,useDingbats=FALSE)
    for(k in selected$key) {
      diagnosis <- readRDS(file.path(run$support,"diagnostics",paste0(k,".rds")))
      e$mc19_wald_diagnostic_page(entries[[k]],diagnosis,cfg);e$mc19_plot_simulation_summary(entries[[k]],diagnosis)
    };grDevices::dev.off()
    # Validate complete compositions/IDs against the PPT ordination source before reusing coordinates.
    per_path <- file.path(project,"results",paste0("MC19_PERMANOVA_TS7_env6_",scope,"_run001"))
    current <- readRDS(file.path(per_path,"supporting_files",paste0(scope,"_environment_analysis.rds")))
    old_per <- file.path(project,"results",if(scope=="coral")"permanova_coral_v1" else "permanova_all_hosts_v1")
    old <- readRDS(file.path(old_per,"supporting_files",paste0(scope,"_environment_analysis.rds")))
    stopifnot(isTRUE(all.equal(old$input$community,current$input$community,tolerance=1e-12)),
      identical(as.character(old$input$metadata$New_ID),as.character(current$input$metadata$New_ID)))
    coord_dir <- file.path(project,"figures",if(scope=="coral")"overview_coral_v1" else "overview_all_hosts_v1")
    cache <- lapply(2:3,function(k)readRDS(file.path(coord_dir,paste0("NMDS",k,"D.rds"))))
    for(i in 1:2) stopifnot(identical(cache[[i]]$New_ID,as.character(current$input$metadata$New_ID)),identical(cache[[i]]$MC,colnames(current$input$community)))
    metadata <- current$input$metadata;xy <- cache[[1]]$coordinates
    stopifnot(identical(rownames(xy),as.character(metadata$New_ID)))
    d <- cbind(metadata,as.data.frame(xy));stats <- current$results
    nmplot <- function(field,a="NMDS1",b="NMDS2",dat=d,stress=cache[[1]]$fit$stress) {
      z <- stats[stats$Factor==field,];stopifnot(nrow(z)==1)
      ggplot2::ggplot(dat,ggplot2::aes(x=.data[[a]],y=.data[[b]],colour=.data[[field]]))+ggplot2::geom_point(size=2.5,alpha=.8)+
        ggplot2::scale_colour_viridis_c()+ggplot2::coord_equal()+ggplot2::labs(colour=unname(labels[field]),
          caption=sprintf("%s; n=%d; stress=%.4f\nPERMANOVA R2=%.3f; raw P=%.4f; BH q=%.4f",tools::toTitleCase(scope),nrow(dat),stress,z$R2_raw,z$P_raw,z$P_BH_raw))+
        ggplot2::theme_classic(base_size=18)+ggplot2::theme(text=ggplot2::element_text(face="bold"),axis.title=ggplot2::element_text(size=19),
          axis.text=ggplot2::element_text(size=16,colour="black"),legend.position="bottom",legend.title=ggplot2::element_text(size=15),legend.text=ggplot2::element_text(size=14),
          plot.caption=ggplot2::element_text(size=12,hjust=0,lineheight=1.1),plot.margin=ggplot2::margin(7,12,5,6))+
        ggplot2::guides(colour=ggplot2::guide_colourbar(direction="horizontal",title.position="top",barwidth=grid::unit(2.6,"in"),barheight=grid::unit(.13,"in")))
    }
    for(field in cfg$environments)savegg(nmplot(field),file.path(folder,paste0("05_NMDS2D_",field)),6.35,6)
    focus <- lapply(c("Nitrite_mgN_L","Ammonia_mgN_L"),nmplot)
    savebase(file.path(folder,"06_NMDS2D_nitrite_ammonia"),function() {
      grid::grid.newpage();for(j in 1:2)print(focus[[j]],vp=grid::viewport(x=(j-.5)/2,y=.5,width=.49,height=1))
    })
    summary <- ggplot2::ggplot(stats,ggplot2::aes(factor(Factor,levels=cfg$environments),R2_raw))+
      ggplot2::geom_col(fill="grey55",width=.65)+ggplot2::geom_text(ggplot2::aes(label=sprintf("P=%.3f\nq=%.3f",P_raw,P_BH_raw)),vjust=-.15,size=4)+
      ggplot2::scale_x_discrete(labels=unname(cfg$short_names[cfg$environments]))+ggplot2::scale_y_continuous(expand=ggplot2::expansion(mult=c(0,.28)))+
      ggplot2::labs(x=NULL,y="PERMANOVA R2",caption=paste(scope,"| Six separate environmental tests; BH across 6 tests."))+
      ggplot2::theme_classic(base_size=18)+ggplot2::theme(plot.caption=ggplot2::element_text(hjust=0,size=12))
    savegg(summary,file.path(folder,"07_PERMANOVA_summary"))
    grDevices::pdf(file.path(aux,"NMDS_Shepard.pdf"),width=12.7,height=6,useDingbats=FALSE)
    for(i in 1:2) {par(mar=c(5,5,3,2),font=2,font.lab=2,font.axis=2);vegan::stressplot(cache[[i]]$fit)};grDevices::dev.off()
    xyz <- cbind(metadata,as.data.frame(cache[[2]]$coordinates))
    for(pair in list(c(1,2),c(1,3),c(2,3))) {
      stem <- file.path(aux,paste0("NMDS3D_projection_axes",paste(pair,collapse="")))
      grDevices::pdf(paste0(stem,".pdf"),width=6.35,height=6,useDingbats=FALSE)
      for(field in cfg$environments)print(nmplot(field,paste0("NMDS",pair[1]),paste0("NMDS",pair[2]),xyz,cache[[2]]$fit$stress))
      grDevices::dev.off()
    }
    widget <- e$mc19_talk_3d_widget(xyz,stats,scope,cache[[2]]$fit$stress)
    e$mc19_save_widget(widget,file.path(aux,"NMDS3D_six_indicators.html"))
    write.csv(d,file.path(aux,"NMDS2D_coordinates_current_environment.csv"),row.names=FALSE)
    write.csv(xyz,file.path(aux,"NMDS3D_coordinates_current_environment.csv"),row.names=FALSE)
    coordinate_audit[[scope]] <- data.frame(scope=scope,k=2:3,stress=vapply(cache,function(x)x$fit$stress,numeric(1)),source=file.path(coord_dir,c("NMDS2D.rds","NMDS3D.rds")),composition_identical=TRUE)
    write.csv(selected,file.path(aux,"GLMM_significant_combinations.csv"),row.names=FALSE)
    write.csv(tab,file.path(aux,"two_part_current_results.csv"),row.names=FALSE)
    manifest[[scope]] <- data.frame(scope=scope,GLMM=run$run_dir,PERMANOVA=per_path,two_part=tp_path,GLMM_pairs=nrow(selected),two_part_components=nrow(tab))
  }
  stopifnot(identical(hashes,tools::md5sum(protected)))
  write.csv(data.frame(path=names(hashes),md5=unname(hashes)),file.path(output,"supporting_files","protected_result_PPT_hashes.csv"),row.names=FALSE)
  write.csv(do.call(rbind,manifest),file.path(output,"supporting_files","analysis_sources.csv"),row.names=FALSE)
  write.csv(do.call(rbind,coordinate_audit),file.path(output,"supporting_files","NMDS_coordinate_audit.csv"),row.names=FALSE)
  write.csv(unique(do.call(rbind,style$evidence)),file.path(output,"supporting_files","GLMM_curve_ranges.csv"),row.names=FALSE)
  writeLines("complete; saved results and PPT unchanged; no refit/permutation/simulation/NMDS computation",file.path(output,"supporting_files","export_status.txt"))
  message("PPT-style export completed: ",output)
})()
