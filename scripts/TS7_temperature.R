# Seven TS, DIN + NO3 + NO2 + TAN + DIP + temperature; fresh full fits/tests.
# Default check is read-only. The user selects each computational stage explicitly.
# options(mc19.analysis_scope="coral", mc19.temperature.stage="check")
(function() {
  starts <- getwd()
  for(f in sys.frames()) if(exists("ofile",f,inherits=FALSE)) {
    ofile <- get("ofile",f,inherits=FALSE)
    if(is.character(ofile) && length(ofile)==1L) starts <- c(starts,dirname(ofile))
  }
  candidates <- getOption("mc19.project_dir",character())
  for(start in starts) {
    current <- normalizePath(start,winslash="/",mustWork=FALSE)
    for(i in 1:6) {
      candidates <- c(candidates,current,file.path(current,"Research","Msc main project-time series"))
      current <- dirname(current)
    }
  }
  valid <- vapply(candidates,function(p) file.exists(file.path(p,"scripts","mc19_profile.R")) &&
    file.exists(file.path(p,"data","MC19_TS8_env5.xlsx")),logical(1))
  if(!any(valid)) stop("Open scripts.Rproj or set options(mc19.project_dir='research project path').")
  project <- normalizePath(candidates[which(valid)[1]],winslash="/",mustWork=TRUE)
  scope <- getOption("mc19.analysis_scope","coral")
  stage <- getOption("mc19.temperature.stage","check")
  stages <- c("check","results","plots","GLMM","PERMANOVA","talk","two_part","all",
    "plot_GLMM","plot_PERMANOVA","plot_two_part")
  stopifnot(length(scope)==1L,scope %in% c("coral","ecosystem"),length(stage)==1L,stage %in% stages)
  old <- options(mc19.pipeline="TS7_env6",mc19.project_dir=project)
  on.exit(options(old),add=TRUE)
  e <- new.env(parent=globalenv())
  for(f in c("glmm_MC19_pilot_helpers.R","glmm_MC19_wald_helpers.R"))
    sys.source(file.path(project,"scripts","GLMM",f),e)
  if(stage=="check") {
    cfg <- e$mc19_wald_config(project,scope)
    input <- e$mc19_wald_load_input(cfg)
    p <- new.env(parent=globalenv())
    sys.source(file.path(project,"scripts","PERMANOVA","permanova_MC19_data_settings.R"),p)
    perma <- p$read_mc19_samples()
    first <- input$blocks[[1]]
    stopifnot(identical(as.character(first$New_ID),as.character(perma$metadata$New_ID)),
      !e$mc19_profile()$incremental,length(input$blocks)==114L,nrow(input$scaling)==6L,
      nrow(perma$batches)==7L,length(p$env_cols)==6L)
    for(field in cfg$environments) stopifnot(all(abs(input$environment_TS[[field]]-perma$batches[[field]])<1e-10))
    n_fields <- c("Nitrate_mgN_L","Nitrite_mgN_L","Ammonia_mgN_L")
    stopifnot(all(abs(input$environment_TS$DIN_mgN_L-rowSums(input$environment_TS[n_fields]))<1e-12))
    cat("Fresh GLMM fits planned: 114; fresh PERMANOVA tests planned: 6 (7! whole-TS label orders each). No legacy model/test reuse.\n")
    print(input$scaling[c("environment","center","scale","n_TS","n_samples")],row.names=FALSE)
    print(e$mc19_verify_temperature(input$environment_TS,project),row.names=FALSE)
    cat("CHECK PASSED: 7 TS;",nrow(first),"samples; 19 MC; 6 indicators; 114 GLMM rows; BH n=114 / PERMANOVA BH n=6.\n",
      "Read-only: no models, residual simulations, permutations, NMDS, plots or output directories created.\n")
    return(invisible(list(config=cfg,GLMM=input,PERMANOVA=perma)))
  }
  run <- function(file) source(file.path(project,"scripts",file),encoding="UTF-8",local=new.env(parent=globalenv()))
  if(stage %in% c("results","GLMM","PERMANOVA","two_part")) {
    prior_plots <- options(mc19.run_plots=FALSE)
    on.exit(options(prior_plots),add=TRUE)
  }
  if(stage %in% c("GLMM","all","results")) run("GLMM/GLMM.R")
  if(stage %in% c("PERMANOVA","all","results")) run("PERMANOVA/run_permanova.R")
  if(stage=="plots") {
    run("GLMM/plot_GLMM_MC19.R")
    run("PERMANOVA/04_plot_permanova.R")
  }
  if(stage %in% c("talk","all","plots")) run("GLMM/plot_MC19_talk_figures.R")
  if(stage %in% c("two_part","all","results")) {
    prior <- options(mc19.two_part.stage="fit")
    on.exit(options(prior),add=TRUE)
    run("GLMM/GLMM_MC19_two_part.R")
  }
  if(stage=="plot_GLMM") run("GLMM/plot_GLMM_MC19.R")
  if(stage=="plot_PERMANOVA") run("PERMANOVA/04_plot_permanova.R")
  if(stage %in% c("plot_two_part","plots")) {
    prior <- options(mc19.two_part.stage="plot")
    on.exit(options(prior),add=TRUE)
    run("GLMM/GLMM_MC19_two_part.R")
  }
  invisible(stage)
})()
