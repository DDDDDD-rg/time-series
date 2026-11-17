# User-run only. Select scope / parts before source(); no PERMANOVA or GLMM refitting.
# parts: forest, NP, NMDS2D, NMDS3D. Existing ecosystem 2D may need new NMDS starts.
local({
  project_dir <- getOption("mc19.project_dir",'C:/Users/xzy/Desktop/CUHK_Study/GNBF6010_Research_Projects/Research/Msc main project-time series')
  scope <- getOption('mc19.analysis_scope','coral')
  parts <- getOption('mc19.talk.parts',c('forest','NP','NMDS2D','NMDS3D'))
  stopifnot(scope %in% c('coral','ecosystem'),all(parts %in% c('forest','NP','NMDS2D','NMDS3D')))
  e <- new.env(parent=globalenv())
  sys.source(file.path(project_dir,'scripts','mc19_run_io.R'),e)
  out <- getOption('mc19.talk.output_dir','')
  if(out=='latest') out <- e$mc19_latest_talk(project_dir,scope)
  for(f in c('glmm_MC19_pilot_helpers.R','glmm_MC19_wald_helpers.R','glmm_MC19_wald_plot_helpers.R','mc19_talk_plot_helpers.R'))
    sys.source(file.path(project_dir,'scripts','GLMM',f),e)
  cfg <- e$mc19_wald_config(project_dir,scope)
  run <- e$mc19_wald_read_run(cfg,getOption('mc19.glmm.result_dir',''))
  if(!nzchar(out)) out <- e$mc19_new_run(project_dir,'talk',scope,e$mc19_method_version(scope,'talk'),run$run_dir,
    requested=getOption('mc19.run_id',NULL)) else {
    if(!dir.exists(out)) dir.create(out,recursive=TRUE)
    record <- e$mc19_io_record(out)
    if(!is.null(record)) {
      if(!identical(record$scope,scope) || !identical(record$source,run$run_dir)) stop('Talk directory source/scope mismatch.')
      e$mc19_io_update(out,state='running')
    } else {
      old_readme <- file.path(out,'README.txt')
      if(file.exists(old_readme) && !paste('GLMM source:',run$run_dir) %in% readLines(old_readme,warn=FALSE))
        stop('Legacy talk directory source mismatch; select its original GLMM source explicitly.')
      write.csv(data.frame(kind='talk',scope=scope,method_version=e$mc19_method_version(scope,'talk'),source=run$run_dir,
        created_at=e$mc19_io_time(),finished_at='',state='running'),file.path(out,'run_record.csv'),row.names=FALSE)
    }
  }
  if('forest' %in% parts) e$mc19_export_heatmap_forest(file.path(out,'MC19_heatmap_forest'),run$results,cfg)
  if('NP' %in% parts) e$mc19_np_plot(read.csv(file.path(run$support,'input_environment_TS.csv')),out)
  if(any(c('NMDS2D','NMDS3D') %in% parts)) {
    for(f in c('permanova_MC19_data_settings.R','permanova_MC19_plot_helpers.R')) sys.source(file.path(project_dir,'scripts','PERMANOVA',f),e)
    ctx <- e$mc19_plot_context()
    gi <- read.csv(file.path(run$support,'input_environment_TS.csv'),stringsAsFactors=FALSE)
    stopifnot(identical(as.character(gi$TS),as.character(ctx$input$batches$TS)))
    for(field in cfg$environments) if(any(abs(gi[[field]]-ctx$input$batches[[field]])>1e-10))
      stop('GLMM/PERMANOVA environment inputs differ: ',field)
    e$mc19_talk_nmds(ctx,out,c(if('NMDS2D' %in% parts) 2L,if('NMDS3D' %in% parts) 3L),as.integer(getOption('mc19.nmds.trymax',500L)))
  }
  writeLines(c(paste('GLMM source:',run$run_dir),paste('Parts:',paste(parts,collapse=',')),
    sprintf('%d environmental panels use saved %d-TS PERMANOVA P/BH results.',e$mc19_profile()$n_env,e$mc19_profile()$n_TS),
    '3D static figures are projections; retain HTML together with companion asset directories.'),file.path(out,'README.txt'))
  record <- e$mc19_io_record(out)
  previous_parts <- if('parts' %in% names(record)) strsplit(record$parts,',',fixed=TRUE)[[1]] else character()
  e$mc19_io_update(out,parts=paste(unique(c(previous_parts,parts)),collapse=','))
  e$mc19_publish_plot(project_dir,'talk',scope,out)
  message('Talk figures: ',out)
})
