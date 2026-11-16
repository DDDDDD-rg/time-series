# Saved five-environment joint GLMMs only. Default check; never refits.
joint_diag_config <- function() {
  entry <- getOption('mc19.joint.diagnostics.project_dir','')
  starts <- c(getwd(),entry)
  for(f in sys.frames()) if(exists('ofile',envir=f,inherits=FALSE))
    starts <- c(starts,dirname(get('ofile',envir=f,inherits=FALSE)))
  candidates <- entry[nzchar(entry)]
  for(s in starts[nzchar(starts)]) for(i in 1:7){
    candidates<-c(candidates,s,file.path(s,'Research','Msc main project-time series'));s<-dirname(s)
  }
  ok<-vapply(candidates,function(x)file.exists(file.path(x,'scripts/GLMM/GLMM_MC19_joint_TS7.R')),logical(1))
  if(!any(ok))stop('Cannot locate project; set mc19.joint.diagnostics.project_dir.')
  project<-normalizePath(candidates[which(ok)[1]],winslash='/',mustWork=TRUE)
  stage<-getOption('mc19.joint.diagnostics.stage','check')
  stopifnot(length(stage)==1L,stage %in% c('check','diagnostics','plots'))
  nsim<-getOption('mc19.joint.diagnostics.n_sim',1000L)
  seed<-getOption('mc19.joint.diagnostics.seed',81270L)
  stopifnot(length(nsim)==1L,is.finite(nsim),nsim>=250,nsim==as.integer(nsim),
    length(seed)==1L,is.finite(seed),seed>=0,seed<.Machine$integer.max-100)
  source_runs<-setNames(file.path(project,'results',paste0('joint_glmm_TS7_env5_',c('coral','ecosystem'),'_run001')),c('coral','ecosystem'))
  list(project=project,stage=stage,n_sim=as.integer(nsim),seed=as.integer(seed),sources=source_runs,
    plot_run=getOption('mc19.joint.diagnostics.run_dir',''))
}

joint_diag_helpers <- function(project) {
  e<-new.env(parent=globalenv())
  wanted<-c('mc19_capture','mc19_simulated_counts','mc19_stat_vector','mc19_simulation_summaries')
  for(x in parse(file.path(project,'scripts/GLMM/glmm_MC19_pilot_helpers.R')))
    if(is.call(x)&&identical(x[[1]],as.name('<-'))&&is.symbol(x[[2]])&&as.character(x[[2]]) %in% wanted)eval(x,e)
  stopifnot(all(vapply(wanted,exists,logical(1),envir=e,inherits=FALSE)));e
}

joint_diag_read <- function(root) {
  m<-readRDS(file.path(root,'supporting_files/manifest.rds'))
  stopifnot(identical(m$config$version,'MC19_joint_TS7_env5_v1'),m$n_models==19L,m$n_coefficients==95L,
    identical(m$config$z,c('Nitrate_z','Nitrite_z','Ammonia_z','DIP_z','Temperature_z')))
  tab<-read.csv(file.path(root,'00_coefficients.csv'))
  stopifnot(nrow(tab)==95L,all(tab$numerical_ok),!anyDuplicated(paste(tab$MC_ID,tab$term)),
    setequal(unique(tab$MC_ID),m$config$mc),all(tab$q_family_size==95L))
  files<-setNames(file.path(root,'supporting_files/models',paste0(m$config$mc,'.rds')),m$config$mc)
  stopifnot(all(file.exists(files)))
  # Read each fit before creating any diagnostic output; identity/count checks only.
  for(mc in names(files)){
    entry<-readRDS(files[[mc]]);fit<-entry$fit;d<-fit$frame;r<-tab[tab$MC_ID==mc,]
    stopifnot(!is.null(fit),isTRUE(entry$model_summary$numerical_ok),entry$MC_ID==mc,
      nrow(d)==m$config$expected_n,nlevels(d$TS)==7L,
      all(m$config$z %in% all.vars(formula(fit))),!'DIN_mgN_L' %in% all.vars(formula(fit)),
      family(fit)$family=='betabinomial',family(fit)$link=='logit')
    stopifnot(max(abs(glmmTMB::fixef(fit)$cond[r$term]-r$beta_per_TS_SD))<1e-10)
  }
  input<-read.csv(file.path(root,'supporting_files/analysis_input.csv'),stringsAsFactors=FALSE)
  list(manifest=m,files=files,input=input,coefficients=tab)
}

joint_diag_one <- function(entry,d,cfg,seed,h) {
  fit<-entry$fit
  # glmmTMB simulate(): draws new TS random intercepts, no parameter refits.
  sims<-h$mc19_simulated_counts(stats::simulate(fit,nsim=cfg$n_sim,seed=seed),d)
  dh<-DHARMa::createDHARMa(simulatedResponse=sims,observedResponse=d$MC_reads,
    fittedPredictedResponse=rowMeans(sims),integerResponse=TRUE,seed=seed)
  funs<-list(uniformity=function()DHARMa::testUniformity(dh,plot=FALSE),
    dispersion=function()DHARMa::testDispersion(dh,plot=FALSE),
    zero_inflation=function()DHARMa::testZeroInflation(dh,plot=FALSE))
  tests<-do.call(rbind,lapply(names(funs),function(n){
    a<-h$mc19_capture(funs[[n]])
    data.frame(test=n,p_diagnostic_only=if(is.null(a$value))NA_real_ else a$value$p.value,
      statistic=if(is.null(a$value))NA_real_ else as.numeric(a$value$statistic)[1],
      error=a$error,warnings=paste(a$warnings,collapse=' | '))
  }))
  summaries<-h$mc19_simulation_summaries(d,sims)
  residual_table<-d[,c('New_ID','TS','MC_reads','Ruegeria_reads_total','MC_proportion',
    'Nitrate_z','Nitrite_z','Ammonia_z','DIP_z','Temperature_z')]
  residual_table$scaled_residual<-dh$scaledResiduals
  residual_table$conditional_proportion<-as.numeric(predict(fit,type='response',re.form=NULL))
  residual_table$simulation_mean_proportion<-rowMeans(sims)/d$Ruegeria_reads_total
  flagged<-any(tests$p_diagnostic_only<.05,na.rm=TRUE)||any(summaries$outside_simulated95 %in% TRUE)||any(nzchar(tests$error))
  list(status=if(flagged)'review_flags' else 'completed_no_automatic_flags',tests=tests,
    summaries=summaries,residual_table=residual_table,dharma=dh,simulated_counts=sims,
    method='Unconditional; TS random intercepts redrawn; no refits',seed=seed,n_sim=cfg$n_sim)
}

joint_diag_main <- function(cfg=joint_diag_config()) {
  for(pkg in c('glmmTMB','DHARMa'))if(!requireNamespace(pkg,quietly=TRUE))stop('Missing R package: ',pkg)
  if(cfg$stage=='plots'){
    if(!nzchar(cfg$plot_run))stop('Set mc19.joint.diagnostics.run_dir to a completed diagnostic run.')
    root<-normalizePath(cfg$plot_run,winslash='/',mustWork=TRUE)
    for(scope in names(cfg$sources)){
      folder<-file.path(root,scope);overview<-read.csv(file.path(folder,'00_diagnostic_overview.csv'))
      pdf_path<-file.path(folder,'diagnostics.pdf')
      if(file.exists(pdf_path))stop('Plot file already exists: ',pdf_path)
      grDevices::pdf(pdf_path,width=12,height=8)
      tryCatch(for(i in seq_len(nrow(overview))){
        if(overview$status[i]=='failed')next
        a<-readRDS(file.path(folder,'models',paste0(overview$MC_ID[i],'.rds')))
        plot(a$dharma);graphics::mtext(paste(scope,overview$MC_ID[i]),outer=TRUE,line=-1)
        for(field in c('Nitrate_z','Nitrite_z','Ammonia_z','DIP_z','Temperature_z'))
          DHARMa::plotResiduals(a$dharma,form=a$residual_table[[field]],main=paste(scope,overview$MC_ID[i],field))
      },finally=grDevices::dev.off())
    }
    return(invisible(root))
  }
  runs<-lapply(cfg$sources,joint_diag_read)
  cat('Verified 38 saved fits. Stage:',cfg$stage,'; simulations per model:',cfg$n_sim,'\n')
  if(cfg$stage=='check')return(invisible(runs))
  # Guard against an accidental call to fitting anywhere below.
  trace('glmmTMB',where=asNamespace('glmmTMB'),tracer=quote(stop('Refitting forbidden in diagnostics')),print=FALSE)
  on.exit(untrace('glmmTMB',where=asNamespace('glmmTMB')),add=TRUE)
  root<-''
  for(i in 1:9999){p<-file.path(cfg$project,'results',sprintf('joint_glmm_TS7_env5_diagnostics_run%03d',i));if(!dir.exists(p)&&dir.create(p)){root<-p;break}}
  if(!nzchar(root))stop('Cannot create an independent diagnostic directory.')
  status<-'incomplete';on.exit(writeLines(status,file.path(root,'run_status.txt')),add=TRUE)
  h<-joint_diag_helpers(cfg$project);write_csv<-function(x,p)write.csv(x,p,row.names=FALSE,fileEncoding='UTF-8',na='')
  fingerprints<-unlist(lapply(runs,function(r)r$files),use.names=FALSE)
  hashes<-tools::md5sum(fingerprints)
  saveRDS(list(config=cfg,source_models=data.frame(path=fingerprints,md5=unname(hashes)),
    method='Unconditional simulation; no refits; diagnostic P values are not association P/q'),file.path(root,'manifest.rds'))
  file.copy(file.path(cfg$project,'scripts/GLMM/diagnose_joint_TS7_env5.R'),file.path(root,'diagnostic_script.R'))
  file.copy(file.path(cfg$project,'scripts/GLMM/glmm_MC19_pilot_helpers.R'),file.path(root,'summary_helpers.R'))
  writeLines(capture.output(sessionInfo()),file.path(root,'sessionInfo.txt'))
  all_rows<-list();has_failure<-FALSE
  for(scope in names(runs)){
    run<-runs[[scope]];folder<-file.path(root,scope);dir.create(folder);dir.create(file.path(folder,'models'))
    rows<-list()
    for(i in seq_along(run$files)){
      mc<-names(run$files)[i];seed<-cfg$seed+i+if(scope=='ecosystem')100L else 0L
      message(scope,' ',i,'/19: ',mc)
      a<-h$mc19_capture(function(){
        entry<-readRDS(run$files[[i]]);d<-run$input[run$input$MC_ID==mc,,drop=FALSE]
        stopifnot(nrow(d)==run$manifest$config$expected_n,
          identical(as.character(d$TS),as.character(entry$fit$frame$TS)),
          all(d$MC_reads==entry$fit$frame[[1]][,1]),
          all(d$Ruegeria_reads_total==rowSums(entry$fit$frame[[1]])))
        joint_diag_one(entry,d,cfg,seed,h)
      })
      result<-a$value
      row<-data.frame(scope=scope,MC_ID=mc,status=if(is.null(result))'failed' else result$status,
        seed=seed,n_sim=cfg$n_sim,error=a$error,warnings=paste(a$warnings,collapse=' | '),
        source_model_md5=unname(hashes[run$files[[i]]]))
      if(!is.null(result)){
        result$source_model<-run$files[[i]];result$source_model_md5<-row$source_model_md5
        saveRDS(result,file.path(folder,'models',paste0(mc,'.rds')))
        write_csv(result$tests,file.path(folder,paste0(mc,'_tests.csv')))
        write_csv(result$summaries,file.path(folder,paste0(mc,'_TS_summary.csv')))
        write_csv(result$residual_table,file.path(folder,paste0(mc,'_residuals.csv')))
        for(test in result$tests$test)row[[paste0(test,'_P')]]<-result$tests$p_diagnostic_only[result$tests$test==test]
      }else has_failure<-TRUE
      for(test in c('uniformity_P','dispersion_P','zero_inflation_P'))if(is.null(row[[test]]))row[[test]]<-NA_real_
      rows[[i]]<-row;write_csv(do.call(rbind,rows),file.path(folder,'00_diagnostic_overview.csv'))
    }
    all_rows[[scope]]<-do.call(rbind,rows)
  }
  stopifnot(identical(unname(tools::md5sum(fingerprints)),unname(hashes)))
  overview<-do.call(rbind,all_rows);write_csv(overview,file.path(root,'00_all_38_models.csv'))
  status<-if(has_failure)'complete_with_failed_diagnostics' else 'complete'
  cat('Saved:',root,'\nModel refits: 0. Original association P/q unchanged.\n')
  invisible(root)
}
if(!isTRUE(getOption('mc19.joint.diagnostics.define_only',FALSE)))joint_diag_main()
