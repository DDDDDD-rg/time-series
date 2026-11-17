# Reuse/output regression checks only. Real fit/permutation/NMDS calls are blocked.
(function() {
  project <- getOption("mc19.project_dir",file.path(getwd(),"Research","Msc main project-time series"))
  project <- normalizePath(project,winslash="/",mustWork=TRUE)
  old <- options(mc19.pipeline="TS7_temperature",mc19.project_dir=project,
    mc19.temperature.reuse_saved=getOption("mc19.temperature.reuse_saved"))
  on.exit(options(old),add=TRUE)
  e <- new.env(parent=globalenv())
  for(f in c("glmm_MC19_pilot_helpers.R","glmm_MC19_wald_helpers.R","mc19_two_part_helpers.R"))
    sys.source(file.path(project,"scripts","GLMM",f),e)
  hooks <- list(c("glmmTMB","glmmTMB"),c("vegan","adonis2"),c("vegan","metaMDS"))
  for(h in hooks) suppressMessages(trace(h[[2]],where=asNamespace(h[[1]]),
    tracer=quote(stop("Forbidden real analysis in reuse test")),print=FALSE))
  on.exit(for(h in hooks) untrace(h[[2]],where=asNamespace(h[[1]])),add=TRUE)
  fixture <- tempfile(".mc19inc-",tmpdir=normalizePath(file.path(project,"../..")))
  dir.create(fixture);on.exit(unlink(fixture,recursive=TRUE),add=TRUE)
  actual_project <- project
  before <- lapply(c("results","figures",".run_registry"),function(a) list.dirs(file.path(actual_project,a),recursive=FALSE))
  # Use an isolated project copy, so helper reloads cannot redirect writes to production.
  sources <- list.files(file.path(actual_project,"scripts"),recursive=TRUE,full.names=TRUE)
  sources <- sources[file.info(sources)$isdir %in% FALSE]
  data_files <- c("data/MC19_TS8_env5.xlsx","data/mc/MC_abundance_absolute.tsv",
    "data/water_quality/Bluff Island_PM11_marine_water_quality.csv")
  for(src in c(sources,file.path(actual_project,data_files))) {
    rel <- substring(src,nchar(actual_project)+2L);dest <- file.path(fixture,rel)
    dir.create(dirname(dest),recursive=TRUE,showWarnings=FALSE)
    stopifnot(file.copy(src,dest))
  }
  candidate <- file.path(actual_project,"results","MC19_GLMM_TS7_temperature_coral_run001","supporting_files")
  for(name in c("wald_manifest.rds","all_results.csv")) if(file.exists(file.path(candidate,name))) {
    dest <- file.path(fixture,"results","MC19_GLMM_TS7_temperature_coral_run001","supporting_files",name)
    dir.create(dirname(dest),recursive=TRUE,showWarnings=FALSE);stopifnot(file.copy(file.path(candidate,name),dest))
  }
  old_sources <- options(mc19.project_dir=fixture,
    mc19.temperature.legacy_glmm=file.path(actual_project,"results","glmm_coral_v1"),
    mc19.temperature.legacy_permanova=file.path(actual_project,"results","permanova_coral_v1"),
    mc19.temperature.legacy_two_part=file.path(actual_project,"results","two_part_coral_v1"))
  on.exit(options(old_sources),add=TRUE)
  project <- normalizePath(fixture,winslash="/",mustWork=TRUE)
  cfg <- e$mc19_wald_config(project,"coral")
  input <- e$mc19_wald_load_input(cfg)
  plan <- e$mc19_incremental_glmm_plan(cfg,input)
  stopifnot(length(plan$nutrients$entries)==76L,length(cfg$environments)==5L,!"DIN_mgN_L" %in% cfg$environments)
  first <- plan$nutrients$entries[[1]]
  changed <- first$data;changed$E_z[1] <- changed$E_z[1]+.1
  stopifnot(inherits(tryCatch(e$mc19_incremental_entry(first,changed),error=identity),"error"))
  # All production writers now receive the isolated project root.
  calls <- character()
  e$mc19_wald_fit_one <- function(mc,d,cfg) {
    stopifnot(identical(unique(d$environment),"SST_C"))
    calls <<- c(calls,mc)
    # Deliberately unavailable synthetic fit exercises NA/BH family handling.
    list(MC_ID=mc,model="BB_TS",data=d,fit=NULL,result=e$mc19_empty_result(mc,"BB_TS",d),
      formula="cbind(MC_reads, Other_reads) ~ E_z + (1 | TS)",version=cfg$version)
  }
  options(mc19.temperature.reuse_saved=TRUE)
  run <- e$mc19_incremental_glmm_run(cfg)
  m <- readRDS(file.path(run,"supporting_files","wald_manifest.rds"))
  tab <- read.csv(file.path(run,"supporting_files","all_results.csv"))
  stopifnot(nrow(tab)==95L,all(tab$q_family_size==95L),!any(tab$environment=="DIN_mgN_L"),
    m$new_fit_count==length(calls),m$reused_fit_count+m$new_fit_count==95L)
  expected <- p.adjust(tab$beta_P_value,"BH",n=95L)
  stopifnot(e$mc19_same(expected,tab$beta_q_value))
  for(key in names(plan$nutrients$entries)) {
    original <- plan$nutrients$entries[[key]]$result
    row <- tab[tab$MC_ID==original$MC_ID & tab$environment==original$environment,]
    stopifnot(e$mc19_same(row$beta_P_value,original$beta_P_value),
      e$mc19_same(row$beta_per_TS_SD,original$beta_per_TS_SD),
      e$mc19_same(row$CI95_low,original$CI95_low))
  }
  readback <- e$mc19_wald_read_run(cfg,run)
  stopifnot(nrow(readback$results)==95L)
  # Reuse one old two-part pair in a fixture; validate model response/predictions, recompute its q.
  legacy <- file.path(actual_project,"results","two_part_coral_v1","supporting_files","manifest.rds")
  key <- readRDS(legacy)$keys[1]
  split <- strsplit(key,"__",fixed=TRUE)[[1]]
  selected <- tab[tab$MC_ID==split[1] & tab$environment==split[2],];selected$key <- key
  d <- input$blocks[[key]];d$key <- key
  ti <- list(selected=selected,data=d,source=run,scope="coral",
    fingerprints=tools::md5sum(file.path(run,"supporting_files",c("analysis_input.csv","all_results.csv"))))
  for(component in c("detection","positive")) {
    cache <- e$mc19_incremental_tp_cache(ti,project,key,component)
    stopifnot(!is.null(cache),!is.null(cache$reused_from))
  }
  tr <- e$mc19_new_run(project,"two_part","coral",e$mc19_method_version("coral","two_part"),run)
  e$mc19_tp_run(ti,project,tr)
  audit <- read.csv(file.path(tr,"supporting_files","reuse_audit.csv"))
  stopifnot(nrow(audit)==2L,all(audit$origin=="reused_validated_component"))
  # Force missing temperature cache: exactly 19 temperature calls, never a nutrient call.
  options(mc19.temperature.reuse_saved=FALSE);calls <- character()
  options(mc19.temperature.legacy_glmm=file.path(actual_project,"results","glmm_all_hosts_v1"))
  cfg2 <- e$mc19_wald_config(project,"ecosystem")
  out2 <- e$mc19_incremental_glmm_run(cfg2)
  mm <- readRDS(file.path(out2,"supporting_files","wald_manifest.rds"))
  stopifnot(length(calls)==19L,setequal(calls,cfg2$mc_ids),mm$new_fit_count==19L,mm$reused_fit_count==76L)
  # PERMANOVA reuses all five saved tests and coordinates without calling adonis2/metaMDS.
  opts <- options(mc19.analysis_scope="coral");on.exit(options(opts),add=TRUE)
  p <- new.env(parent=globalenv())
  sys.source(file.path(project,"scripts","PERMANOVA","permanova_MC19_data_settings.R"),p)
  sys.source(file.path(project,"scripts","mc19_temperature_incremental.R"),p)
  inp <- p$read_mc19_samples();pp <- p$mc19_incremental_permanova_plan(project,"coral",inp)
  pr <- p$mc19_incremental_permanova_run(project,"coral",inp,p$input_file)
  pt <- read.csv(file.path(pr,"supporting_files","environment_results_full.csv"),fileEncoding="UTF-8-BOM")
  stopifnot(nrow(pt)==5L,identical(pt$Factor,cfg$environments),
    e$mc19_same(pt$P_raw,pp$old$results$P_raw[pp$indices]),
    e$mc19_same(pt$P_BH_raw,p.adjust(pt$P_raw,"BH",n=5L)),
    file.exists(file.path(pr,"supporting_files","plot_NMDS2D.rds")))
  stopifnot(identical(before,lapply(c("results","figures",".run_registry"),function(a) list.dirs(file.path(actual_project,a),recursive=FALSE))))
  cat("PASS: 76 old GLMMs retained, only 19 temperature calls when needed, 95-test BH, reused two-part models, five reused PERMANOVA tests, cached NMDS and isolated fixture outputs. Real analyses: zero.\n")
})()
