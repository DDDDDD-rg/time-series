# User-run entry. Default check is read-only. Choose stage="fit" explicitly for 44 formal fits.
# options(mc19.analysis_scope="coral", mc19.two_part.stage="check")
# For resumption or replot, set mc19.two_part.run_dir to the existing run directory.
local({
  project_dir <- 'C:/Users/xzy/Desktop/CUHK_Study/GNBF6010_Research_Projects/Research/Msc main project-time series'
  scope <- getOption('mc19.analysis_scope','coral')
  stage <- getOption('mc19.two_part.stage','check') # check / fit / plot
  stopifnot(scope %in% c('coral','ecosystem'),stage %in% c('check','fit','plot'))
  e <- new.env(parent=globalenv())
  sys.source(file.path(project_dir,'scripts','mc19_run_io.R'),e)
  sys.source(file.path(project_dir,'scripts','GLMM','mc19_two_part_helpers.R'),e)
  out <- getOption('mc19.two_part.run_dir','')
  if(stage=='plot') {
    out <- e$mc19_latest_result(project_dir,'two_part',scope,out)
    e$mc19_tp_plot(out,project_dir,parts=getOption('mc19.two_part.plot_parts',c('forest','curves','residual')))
  } else {
    input <- e$mc19_tp_input(project_dir,scope,getOption('mc19.glmm.result_dir',''))
    print(input$selected[c('MC_ID','environment','beta_q_value')],row.names=FALSE)
    message('Scope: ',scope,'; selected pairs: ',nrow(input$selected),'; component tasks: ',2*nrow(input$selected))
    if(stage=='fit') {
      if(!nzchar(out)) out <- e$mc19_new_run(project_dir,'two_part',scope,'MC19_two_part_v1',input$source,
        requested=getOption('mc19.run_id',NULL))
      if(out=='latest') stop('For resumption supply an explicit run directory; latest is for plot only.')
      e$mc19_tp_run(input,project_dir,out)
      e$mc19_tp_plot(out,project_dir)
      message('Two-part output: ',out)
    }
  }
})
