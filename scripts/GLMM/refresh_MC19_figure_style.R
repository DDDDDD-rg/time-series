# Style-only refresh of current saved MC19 outputs. Never fits or simulates.
# Set mc19.style.output_root to a staging figures directory for visual review.
local({
  project_dir <- 'C:/Users/xzy/Desktop/CUHK_Study/GNBF6010_Research_Projects/Research/Msc main project-time series'
  destination <- getOption('mc19.style.output_root',file.path(project_dir,'figures'))
  e <- new.env(parent=globalenv())
  sys.source(file.path(project_dir,'scripts','mc19_run_io.R'),e)
  for(f in c('glmm_MC19_pilot_helpers.R','glmm_MC19_wald_helpers.R',
    'glmm_MC19_wald_plot_helpers.R','mc19_talk_plot_helpers.R','mc19_two_part_helpers.R'))
    sys.source(file.path(project_dir,'scripts','GLMM',f),e)
  # Refreshing existing paths must not redirect the normal latest pointers.
  e$mc19_publish_plot <- function(project_dir,kind,scope,path) {
    e$mc19_io_update(path,state='complete',finished_at=e$mc19_io_time())
    invisible(path)
  }
  original_options <- list(mc19.analysis_scope=getOption('mc19.analysis_scope'),
    mc19.two_part.figure_dir=getOption('mc19.two_part.figure_dir'),mc19.result_dir=getOption('mc19.result_dir'))
  on.exit(options(original_options),add=TRUE)
  options(mc19.result_dir='')
  for(scope in c('coral','ecosystem')) {
    message('Style refresh: ',scope)
    options(mc19.analysis_scope=scope)
    cfg <- e$mc19_wald_config(project_dir,scope)
    run <- e$mc19_wald_read_run(cfg)
    old_glmm <- file.path(project_dir,'figures',basename(run$run_dir))
    glmm_fig <- file.path(destination,basename(run$run_dir))
    dir.create(glmm_fig,recursive=TRUE,showWarnings=FALSE)
    entries <- lapply(run$manifest$model_files,readRDS)
    # Reuse saved prediction values exactly; do not call predict or refit.
    curves <- read.csv(file.path(old_glmm,'supporting_files','prediction_curves.csv'))
    grids <- split(curves,e$mc19_wald_key(curves$MC_ID,curves$environment,curves$model))
    environment <- read.csv(file.path(run$support,'input_environment_TS.csv'))
    selected <- e$mc19_significant_curve_rows(run$results,cfg)
    dir.create(file.path(glmm_fig,'supporting_files'),recursive=TRUE,showWarnings=FALSE)
    e$mc19_write_csv(selected,file.path(glmm_fig,'supporting_files','significant_curve_index.csv'))
    e$mc19_export_significant_curves(file.path(glmm_fig,'02_all_significant_nutrient_curves'),
      run$results,entries,grids,cfg,environment,selected)
    # Re-render existing supplemental diagnostics, using stored residuals only.
    if(isTRUE(getOption('mc19.style.diagnostics',TRUE))) for(env in cfg$environments) {
      short <- unname(cfg$short_names[env])
      target <- file.path(glmm_fig,'supporting_files','diagnostics',paste0(short,'_diagnostics.pdf'))
      dir.create(dirname(target),recursive=TRUE,showWarnings=FALSE)
      e$mc19_pdf(target,function() {
        e$mc19_wald_note(paste(short,'- diagnostic supplement'),
          c('Saved diagnostic residuals and simulation summaries; no new simulation.',
            'Diagnostic P values concern distributional fit, not nutrient associations.'))
        for(mc in cfg$mc_ids) {
          key <- e$mc19_wald_key(mc,env,'BB_TS')
          path <- file.path(run$support,'diagnostics',paste0(key,'.rds'))
          if(file.exists(path)) {
            diagnosis <- readRDS(path)
            e$mc19_wald_diagnostic_page(entries[[key]],diagnosis,cfg)
            e$mc19_plot_simulation_summary(entries[[key]],diagnosis)
          } else e$mc19_wald_note(key,'Saved diagnostics unavailable.')
        }
      },width=15,height=12)
    }
    talk_name <- paste0('MC19_talk_',scope,'_20261031_153730')
    talk <- file.path(destination,talk_name)
    dir.create(talk,recursive=TRUE,showWarnings=FALSE)
    e$mc19_export_heatmap_forest(file.path(talk,'MC19_heatmap_forest'),run$results,cfg)
    if(scope=='coral') e$mc19_np_plot(environment,talk)
    for(f in c('permanova_MC19_data_settings.R','permanova_MC19_plot_helpers.R'))
      sys.source(file.path(project_dir,'scripts','PERMANOVA',f),e)
    # data_settings loads shared I/O again, so restore the nonpublishing wrapper.
    e$mc19_publish_plot <- function(project_dir,kind,scope,path) {
      e$mc19_io_update(path,state='complete',finished_at=e$mc19_io_time());invisible(path)
    }
    ctx <- e$mc19_plot_context()
    # Copy only when a new staging folder has no cache. Cache-only mode fails
    # rather than computing any new NMDS solution, even if best repeats is zero.
    for(k in 2:3) {
      cache <- file.path(talk,paste0('NMDS',k,'D.rds'))
      if(!file.exists(cache)) stopifnot(file.copy(file.path(project_dir,'figures',talk_name,basename(cache)),cache))
    }
    e$mc19_talk_nmds(ctx,talk,c(2L,3L),cache_only=TRUE)
    two <- e$mc19_latest_result(project_dir,'two_part',scope)
    options(mc19.two_part.figure_dir=file.path(destination,basename(two)))
    e$mc19_tp_plot(two,project_dir,parts=getOption('mc19.style.two_part.parts',c('forest','curves','residual')))
    # Keep the most recent standalone forest revision in sync as well.
    latest_forest <- file.path(destination,paste0('MC19_two_part_',scope,'_legacy_plot001'))
    options(mc19.two_part.figure_dir=latest_forest)
    e$mc19_tp_plot(two,project_dir,parts='forest')
  }
  message('Style-only refresh complete: ',destination)
})
