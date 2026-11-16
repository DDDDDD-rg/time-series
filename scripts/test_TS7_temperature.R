# Input/routing/rendering regression checks. No real fits, simulations, tests or NMDS.
(function() {
  project <- getOption("mc19.project_dir",file.path(getwd(),"Research","Msc main project-time series"))
  project <- normalizePath(project,winslash="/",mustWork=TRUE)
  old <- options(mc19.pipeline="TS7_env6",mc19.project_dir=project,
    mc19.analysis_scope=getOption("mc19.analysis_scope"),mc19.temperature.stage=getOption("mc19.temperature.stage"))
  on.exit(options(old),add=TRUE)
  e <- new.env(parent=globalenv())
  for(f in c("glmm_MC19_pilot_helpers.R","glmm_MC19_wald_helpers.R",
    "glmm_MC19_wald_plot_helpers.R","mc19_talk_plot_helpers.R","mc19_two_part_helpers.R"))
    sys.source(file.path(project,"scripts","GLMM",f),e)
  # A tripwire proves that input/render checks cannot silently start analysis.
  hooked <- list()
  on.exit(for(h in hooked) untrace(h[[2]],where=asNamespace(h[[1]])),add=TRUE)
  for(h in list(c("glmmTMB","glmmTMB"),c("vegan","adonis2"),c("vegan","metaMDS"),c("DHARMa","simulateResiduals"))) {
    suppressMessages(trace(h[[2]],where=asNamespace(h[[1]]),tracer=quote(stop("Analysis tripwire triggered")),print=FALSE))
    hooked[[length(hooked)+1L]] <- h
  }
  fail <- function(x) stopifnot(inherits(tryCatch({force(x);NULL},error=identity),"error"))
  original <- tools::md5sum(file.path(project,"scripts","GLMM","GLMM_MC19_joint_TS7.R"))
  dirs_before <- lapply(c("results","figures",".run_registry"),function(a) list.dirs(file.path(project,a),recursive=FALSE))
  fixture <- tempfile(".mc19t-",tmpdir=normalizePath(file.path(project,"../..")));dir.create(fixture)
  on.exit(unlink(fixture,recursive=TRUE),add=TRUE) # Only this explicitly created temporary directory.
  for(scope in c("coral","ecosystem")) {
    options(mc19.analysis_scope=scope)
    cfg <- e$mc19_wald_config(project,scope)
    input <- e$mc19_wald_load_input(cfg)
    stopifnot(!e$mc19_profile()$incremental,cfg$n_models==114L,cfg$q_family_size==114L,cfg$expected_TS==7L,
      length(input$blocks)==114L,nrow(input$scaling)==6L,input$all_MC_count==175L,
      input$scaling$n_samples==if(scope=="coral")98L else 310L,
      identical(input$environment_TS$SST_C,c(25.8,31.4,28.6,16.7,24.6,26.6,25.2)))
    t <- input$scaling[input$scaling$environment=="SST_C",]
    stopifnot(abs(t$center-mean(input$environment_TS$SST_C))<1e-12,
      abs(t$scale-sd(input$environment_TS$SST_C))<1e-12,t$censored_TS==0)
    nf <- c("Nitrate_mgN_L","Nitrite_mgN_L","Ammonia_mgN_L")
    stopifnot(all(abs(input$environment_TS$DIN_mgN_L-rowSums(input$environment_TS[nf]))<1e-12),
      all(input$environment_TS$DIN_mgN_L_censored ==
        (rowSums(input$environment_TS[paste0(nf,"_censored")])>0)))
    legacy <- read.csv(file.path(project,"results",if(scope=="coral")"glmm_coral_v1" else "glmm_all_hosts_v1",
      "supporting_files","analysis_input.csv"),stringsAsFactors=FALSE)
    for(field in setdiff(cfg$environments,"DIN_mgN_L")) if(field!="SST_C") {
      a <- input$blocks[[paste("MC15",field,sep="__")]]
      b <- legacy[legacy$MC_ID=="MC15" & legacy$environment==field,]
      b <- b[match(a$New_ID,b$New_ID),]
      stopifnot(identical(as.character(a$New_ID),as.character(b$New_ID)),
        all(a$MC_reads==b$MC_reads),all(a$Ruegeria_reads_total==b$Ruegeria_reads_total),all(abs(a$E_raw-b$E_raw)<1e-12))
    }
    p <- new.env(parent=globalenv())
    for(f in c("permanova_MC19_data_settings.R","permanova_MC19_plot_helpers.R"))
      sys.source(file.path(project,"scripts","PERMANOVA",f),p)
    perma <- p$read_mc19_samples()
    stopifnot(identical(as.character(perma$metadata$New_ID),input$blocks[[1]]$New_ID),
      nrow(perma$batches)==7L,length(p$env_cols)==6L)
    orders <- p$enumerate_ts_orders(7L) # Enumerate labels only; do not fit their statistics.
    stopifnot(nrow(orders)==5040L,identical(as.integer(orders[1,]),1:7))
    bad <- input$environment_TS;bad$SST_C_raw[2] <- NA_character_;fail(e$mc19_prepare_environment(bad))
    bad <- input$environment_TS;bad$SST_C[2] <- bad$SST_C[2]+1;fail(e$mc19_verify_temperature(bad,project))
    bad <- input$environment_TS[-1,];fail(e$mc19_prepare_environment(bad))
    # Check independent output namespaces and incompatible manifests, using fixtures only.
    run <- e$mc19_new_run(fixture,"GLMM",scope,cfg$version)
    support <- file.path(run,"supporting_files");dir.create(support)
    saveRDS(list(analysis_scope=scope,version=cfg$version,n_TS=7L,fit_complete=TRUE,
      model_files=rep("fixture",114)),file.path(support,"wald_manifest.rds"))
    writeLines("fixture",file.path(run,"00_result.csv"));writeLines("fixture",file.path(support,"all_results.csv"))
    stopifnot(e$mc19_result_complete(run,"GLMM",scope),grepl("TS7_env6",basename(run)),
      length(e$mc19_io_aliases(fixture,"GLMM",scope))==0L)
    options(mc19.pipeline="TS8_env5");stopifnot(!e$mc19_result_complete(run,"GLMM",scope))
    options(mc19.pipeline="TS7_temperature")
    stopifnot(!e$mc19_result_complete(run,"GLMM",scope),e$mc19_profile()$n_models==95L,
      e$mc19_profile()$incremental,!"DIN_mgN_L" %in% e$mc19_profile()$fields)
    options(mc19.pipeline="TS7_env6")
    # Artificial coefficients exercise all six columns, forest panels and temperature curves.
    fake <- expand.grid(MC_ID=cfg$mc_ids,environment=cfg$environments,stringsAsFactors=FALSE)
    fake$model <- "BB_TS";fake$numerical_ok <- TRUE
    fake$beta_per_TS_SD <- seq(-1,1,length.out=nrow(fake));fake$SE <- .15
    fake$CI95_low <- fake$beta_per_TS_SD-.3;fake$CI95_high <- fake$beta_per_TS_SD+.3
    fake$beta_q_value <- .01
    grDevices::pdf(file.path(fixture,paste0(scope,"_SYNTHETIC_heatmaps.pdf")),width=44,height=18)
    e$mc19_main_heatmap(fake,cfg);e$mc19_heatmap_forest(fake,cfg);grDevices::dev.off()
    selected <- e$mc19_significant_curve_rows(fake,cfg)
    stopifnot(nrow(selected)==114L,any(selected$environment=="SST_C"),
      any(selected$environment=="DIN_mgN_L"),max(selected$column)==6L)
    fake$beta_q_value <- 1
    stopifnot(nrow(e$mc19_significant_curve_rows(fake,cfg))==0L)
    grDevices::pdf(file.path(fixture,paste0(scope,"_SYNTHETIC_nonsignificant.pdf")),width=32,height=18)
    e$mc19_heatmap_forest(fake,cfg);grDevices::dev.off()
    cat("PASS:",scope,"input, 114 fresh-fit tasks, 7! label enumeration, source guards and six-column render fixtures.\n")
  }
  # Inspect the timeline and 3D control structure without NMDS or statistical fitting.
  plots <- list();e$mc19_talk_save <- function(p,stem,...) plots[[basename(stem)]] <<- p
  e$mc19_np_plot(input$environment_TS,fixture)
  timeline <- read.csv(file.path(fixture,"NP_by_TS_data.csv"))
  stopifnot(nrow(timeline)==42L,sum(grepl("Temperature",timeline$Indicator))==7L,
    length(ggplot2::ggplot_build(plots$NP_by_TS)$layout$layout$PANEL)==6L)
  d <- perma$metadata
  d$NMDS1 <- seq_len(nrow(d))/nrow(d);d$NMDS2 <- sin(seq_len(nrow(d)));d$NMDS3 <- cos(seq_len(nrow(d)))
  mock <- data.frame(Factor=cfg$environments,R2_raw=.01,P_raw=.5,P_BH_raw=.7)
  w <- plotly::plotly_build(e$mc19_talk_3d_widget(d,mock,"SYNTHETIC",.2,"SST_C"))
  stopifnot(length(w$x$layout$updatemenus[[1]]$buttons)==6L,
    grepl("deg C",w$x$data[[6]]$marker$colorbar$title$text,fixed=TRUE),
    all(grepl("deg C",w$x$data[[6]]$text,fixed=TRUE)))
  key <- "MC15__SST_C"
  tab <- data.frame(key=key,MC_ID="MC15",environment="SST_C",component=c("detection","positive"),beta=.2,q=.1)
  pred <- data.frame(key=key,component=rep(c("detection","positive"),each=7),E_raw=rep(input$environment_TS$SST_C,2),mean=.1,low=.05,high=.2)
  obs <- data.frame(key=key,component=rep(c("detection","positive"),each=7),E_raw=rep(input$environment_TS$SST_C,2),observed=c(1,1,1,0,1,1,0,rep(.1,7)))
  grDevices::pdf(file.path(fixture,"SYNTHETIC_temperature_two_part.pdf"),width=22,height=11.25)
  check <- e$mc19_tp_curve_page(tab,pred,obs,key);grDevices::dev.off()
  stopifnot(check$detection_points==7,check$positive_prediction_points==7)
  # The merged PDFs must include the temperature pages, and fail when one is missing.
  dir.create(file.path(fixture,"supporting_files"))
  stems <- c("environment_R2_BH","NMDS_TS",paste0("NMDS_environment_",cfg$environments),"host_highlight_descriptive",
    "environment_time_series","MC_individual_samples_by_TS","MC_sample_abundance_heatmap",paste0("MC_composition_",cfg$environments))
  for(stem in stems) {
    grDevices::pdf(file.path(fixture,"supporting_files",paste0(stem,".pdf")))
    graphics::plot.new();graphics::text(.5,.5,paste("SYNTHETIC TEST PAGE",stem));grDevices::dev.off()
  }
  output <- p$mc19_bundle_figures(list(figure_root=fixture),p$mc19_pdf_tools())
  record <- jsonlite::fromJSON(paste(output,collapse="\n"))
  stopifnot(identical(record$bundles$pages,c(9L,3L,6L)))
  missing <- file.path(fixture,"supporting_files","MC_composition_SST_C.pdf")
  stopifnot(file.remove(missing))
  fail(suppressWarnings(p$mc19_bundle_figures(list(figure_root=fixture),p$mc19_pdf_tools())))

  # Exercise the all-stage controller with a recording source(), never the analysis entries.
  sandbox <- new.env(parent=globalenv());calls <- list()
  sandbox$source <- function(file,...) {
    calls[[length(calls)+1L]] <<- list(file=basename(file),profile=getOption("mc19.pipeline"),
      two_part_stage=getOption("mc19.two_part.stage"),plots=getOption("mc19.run_plots",TRUE))
    invisible(NULL)
  }
  options(mc19.temperature.stage="all")
  sys.source(file.path(project,"scripts","TS7_temperature.R"),sandbox)
  stopifnot(identical(vapply(calls,`[[`,character(1),"file"),
    c("GLMM.R","run_permanova.R","plot_MC19_talk_figures.R","GLMM_MC19_two_part.R")),
    all(vapply(calls,function(x) identical(x$profile,"TS7_env6"),logical(1))),
    identical(calls[[4]]$two_part_stage,"fit"))
  before_plot_option <- getOption("mc19.run_plots")
  calls <- list();options(mc19.temperature.stage="results")
  sys.source(file.path(project,"scripts","TS7_temperature.R"),sandbox)
  stopifnot(identical(vapply(calls,`[[`,character(1),"file"),
    c("GLMM.R","run_permanova.R","GLMM_MC19_two_part.R")),
    all(vapply(calls,function(x) identical(x$plots,FALSE),logical(1))),
    identical(getOption("mc19.run_plots"),before_plot_option))
  calls <- list();options(mc19.temperature.stage="plots")
  sys.source(file.path(project,"scripts","TS7_temperature.R"),sandbox)
  stopifnot(identical(vapply(calls,`[[`,character(1),"file"),
    c("plot_GLMM_MC19.R","04_plot_permanova.R","plot_MC19_talk_figures.R","GLMM_MC19_two_part.R")),
    identical(calls[[4]]$two_part_stage,"plot"))
  # Exercise the PERMANOVA no-plot gate with a fake analysis source.
  gate <- new.env(parent=globalenv());seen <- character()
  gate$source <- function(file,...,local) {
    seen <<- c(seen,basename(file))
    local$run_success <- TRUE;local$support_dir <- fixture;local$output_dir <- fixture
  }
  gate_old <- options(mc19.run_plots=FALSE,mc19.two_part.stage="fit",mc19.two_part.run_dir="")
  on.exit(options(gate_old),add=TRUE)
  sys.source(file.path(project,"scripts","PERMANOVA","run_permanova.R"),gate)
  stopifnot(identical(seen,"03_test_season_environment.R"),
    readLines(file.path(fixture,"00_pipeline_status.txt"))=="analysis_only_complete")
  # Exercise the two-part no-plot gate with fit/plot recorders; no real models.
  gate$sys.source <- function(file,envir,...) {
    envir$mc19_tp_input <- function(...) list(selected=data.frame(MC_ID="MC15",environment="SST_C",beta_q_value=.01),source=fixture)
    envir$mc19_new_run <- function(...) fixture
    envir$mc19_method_version <- function(...) "fixture"
    envir$mc19_tp_run <- function(...) seen <<- c(seen,"fit")
    envir$mc19_tp_plot <- function(...) stop("Plot gate failed")
  }
  seen <- character()
  sys.source(file.path(project,"scripts","GLMM","GLMM_MC19_two_part.R"),gate)
  stopifnot(identical(seen,"fit"))
  options(gate_old)
  cat("PASS: results/plots routing, restored options, PERMANOVA and two-part no-plot gates.\n")
  options(mc19.temperature.stage=NULL)
  stopifnot(identical(original,tools::md5sum(names(original))),
    identical(dirs_before,lapply(c("results","figures",".run_registry"),function(a) list.dirs(file.path(project,a),recursive=FALSE))))
  cat("PASS: timeline 42 rows/6 panels; 3D six buttons and degrees C; joint script and production output directories unchanged.\n")
})()
