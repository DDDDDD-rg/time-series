# One entry: source audit -> 42 leave-TS fits -> targeted sample influence.
# CLI: Rscript joint_TS7_sensitivity.R --run --project=<root> --workers=2
# Resume: add --resume=<existing output directory>. Default is read-only check.
js_save <- function(x,path) {
  tmp<-paste0(path,'.tmp-',Sys.getpid());saveRDS(x,tmp)
  if(file.exists(path)||!file.rename(tmp,path))stop('Checkpoint exists or cannot publish: ',path)
}
js_csv <- function(x,path)write.csv(x,path,row.names=FALSE,fileEncoding='UTF-8',na='')
js_engine <- function(path) {
  e<-new.env(parent=globalenv());e$options_unused<-NULL
  old<-options(mc19.joint.define_only=TRUE);on.exit(options(old))
  sys.source(path,e);e
}
js_prepare <- function(project) {
  if(!requireNamespace('readxl',quietly=TRUE))stop('Missing readxl')
  runs<-list();audits<-list();tasks<-list();fingerprints<-character()
  add<-function(scope,mc,kind,label,ids=character(),ts=''){
    tasks[[length(tasks)+1L]]<<-list(id=sprintf('task_%03d',length(tasks)+1L),scope=scope,mc=mc,
      kind=kind,label=label,drop_ids=as.character(ids),drop_TS=ts)
  }
  for(scope in c('coral','ecosystem')){
    root<-file.path(project,'results',paste0('joint_glmm_TS7_env5_',scope,'_run001'))
    manifest<-readRDS(file.path(root,'supporting_files/manifest.rds'))
    stopifnot(manifest$n_models==19L,manifest$n_coefficients==95L,manifest$config$version=='MC19_joint_TS7_env5_v1')
    d<-read.csv(file.path(root,'supporting_files/analysis_input.csv'),stringsAsFactors=FALSE)
    d$New_ID<-as.character(d$New_ID)
    baseline<-read.csv(file.path(root,'00_coefficients.csv'),stringsAsFactors=FALSE)
    frozen<-file.path(root,'supporting_files/source_inputs')
    book<-file.path(frozen,'MC19_TS8_env5.xlsx');raw_path<-file.path(frozen,'MC_abundance_absolute.tsv')
    meta<-as.data.frame(readxl::read_excel(book,sheet='Metadata'));meta$New_ID<-as.character(meta$New_ID)
    community<-as.data.frame(readxl::read_excel(book,sheet='Community'));community$New_ID<-as.character(community$New_ID)
    raw<-read.delim(raw_path,check.names=FALSE,quote='',comment.char='');stopifnot(nrow(raw)==175L,!anyDuplicated(raw$MC))
    models<-setNames(file.path(root,'supporting_files/models',paste0(c('MC15','MC6','MC98','MC95','MC43'),'.rds')),c('MC15','MC6','MC98','MC95','MC43'))
    files<-c(book,raw_path,file.path(root,'supporting_files/manifest.rds'),file.path(root,'supporting_files/analysis_input.csv'),file.path(root,'00_coefficients.csv'),models)
    fingerprints<-c(fingerprints,files)
    for(mc in names(models)){
      rows<-d[d$MC_ID==mc,,drop=FALSE];entry<-readRDS(models[[mc]])
      ids<-rows$New_ID;stopifnot(nrow(rows)==manifest$config$expected_n,!anyDuplicated(ids),all(ids %in% names(raw)))
      counts<-as.matrix(raw[,ids,drop=FALSE]);N<-colSums(counts);Y<-as.numeric(counts[match(mc,raw$MC),])
      im<-match(ids,meta$New_ID);ic<-match(ids,community$New_ID)
      stopifnot(!anyNA(im),!anyNA(ic),all(Y==rows$MC_reads),all(N==rows$Ruegeria_reads_total),
        all(abs(Y/N-as.numeric(community[[mc]][ic]))<1e-10),
        identical(as.character(meta$TS[im]),as.character(rows$TS)),all(meta$Sample_ID[im]==rows$Sample_ID),
        all(entry$fit$frame[[1]][,1]==Y),all(rowSums(entry$fit$frame[[1]])==N))
      if(mc %in% c('MC15','MC6','MC98'))for(ts in paste0('TS',1:7))add(scope,mc,'leave_TS',ts,ts=ts)
      if(mc %in% c('MC95','MC43')){
        diag_path<-file.path(project,'results/joint_glmm_TS7_env5_diagnostics_run001',scope,'models',paste0(mc,'.rds'))
        a<-readRDS(diag_path);fingerprints<-c(fingerprints,diag_path)
        stopifnot(a$source_model_md5==unname(tools::md5sum(models[[mc]])),identical(as.character(a$residual_table$New_ID),ids))
        squared<-(Y-rowMeans(a$simulated_counts))^2;share<-squared/sum(squared)
        prop_rank<-rank(-rows$MC_proportion,ties.method='min');res_rank<-rank(-share,ties.method='min')
        selected<-sort(unique(c(order(-rows$MC_proportion,ids)[1:3],order(-share,ids)[1:3])))
        audit<-rows[,c('New_ID','Sample_ID','TS','MC_reads','Ruegeria_reads_total','MC_proportion')]
        audit$scope<-scope;audit$MC_ID<-mc;audit$Host_type<-meta$Host_type[im];audit$Host_niche<-meta$Host_niche[im]
        audit$raw_MC_reads<-Y;audit$raw_total_reads<-N;audit$workbook_proportion<-as.numeric(community[[mc]][ic])
        audit$raw_identity_checks_passed<-TRUE;audit$count_residual_square_share<-share
        audit$proportion_rank<-prop_rank;audit$count_residual_rank<-res_rank;audit$selected_for_influence<-seq_len(nrow(rows)) %in% selected
        audits[[paste(scope,mc)]]<-audit
        for(j in selected)add(scope,mc,'leave_sample',ids[j],ids=ids[j])
        add(scope,mc,'leave_selected_bundle',paste(ids[selected],collapse=';'),ids=ids[selected])
      }
    }
    runs[[scope]]<-list(data=d,baseline=baseline,config=manifest$config)
  }
  fingerprints<-unique(c(fingerprints,file.path(project,'scripts/GLMM/GLMM_MC19_joint_TS7.R'),
    file.path(project,'scripts/GLMM/joint_TS7_sensitivity.R')))
  list(runs=runs,tasks=tasks,audit=do.call(rbind,audits),fingerprints=data.frame(path=fingerprints,md5=unname(tools::md5sum(fingerprints))))
}
js_task <- function(task,root) {
  path<-file.path(root,'tasks',paste0(task$id,'.rds'))
  if(file.exists(path))return(readRDS(path))
  answer<-tryCatch({
    if(!exists('.js_context_root',envir=.GlobalEnv,inherits=FALSE)||get('.js_context_root',envir=.GlobalEnv)!=root){
      assign('.js_context',readRDS(file.path(root,'context.rds')),envir=.GlobalEnv)
      assign('.js_context_root',root,envir=.GlobalEnv)
    }
    ctx<-get('.js_context',envir=.GlobalEnv);run<-ctx$runs[[task$scope]]
    e<-js_engine(file.path(root,'joint_fit_source.R'))
    full<-run$data[run$data$MC_ID==task$mc,,drop=FALSE]
    keep<-if(task$kind=='leave_TS')full$TS!=task$drop_TS else !full$New_ID %in% task$drop_ids
    d<-full[keep,,drop=FALSE];d$TS<-droplevels(factor(d$TS,levels=paste0('TS',1:7)))
    stopifnot(nrow(d)>0)
    env<-d[!duplicated(d$TS),run$config$z,drop=FALSE];X<-model.matrix(reformulate(run$config$z),env)
    rank<-qr(X)$rank
    if(rank<ncol(X))stop('Omitted-data environmental design is rank deficient')
    fit<-e$joint_TS7_fit_one(task$mc,d,run$config)
    out<-fit$coefficients;base<-run$baseline[run$baseline$MC_ID==task$mc,]
    ix<-match(out$term,base$term);stopifnot(!anyNA(ix))
    out$original_beta<-base$beta_per_TS_SD[ix];out$original_SE<-base$SE[ix];out$original_q<-base$q_value[ix]
    out$delta_beta<-out$beta_per_TS_SD-out$original_beta
    out$delta_in_original_SE<-out$delta_beta/out$original_SE
    out$sign_changed<-ifelse(out$numerical_ok,sign(out$beta_per_TS_SD)!=sign(out$original_beta),NA)
    out$absolute_beta_ratio<-abs(out$beta_per_TS_SD)/ifelse(abs(out$original_beta)>1e-8,abs(out$original_beta),NA_real_)
    out$q_value<-NA_real_;out$q_family_size<-NA_integer_
    out$inference_status<-'sensitivity_only_Wald_descriptive_no_new_BH'
    out$task_id<-task$id;out$kind<-task$kind;out$omitted<-task$label
    out$design_rank<-rank;out$design_columns<-ncol(X);out$condition_number<-kappa(X,exact=TRUE)
    out$saturated_TS_design<-nrow(env)==ncol(X)
    list(task=task,status=if(fit$model_summary$numerical_ok)'complete' else 'numerical_review',
      coefficients=out,model_check=fit$model_summary,fit=fit$fit,error='')
  },error=function(err)list(task=task,status='failed',error=conditionMessage(err)))
  js_save(answer,path);answer
}
js_collect <- function(root,tasks) {
  files<-file.path(root,'tasks',paste0(vapply(tasks,`[[`,character(1),'id'),'.rds'))
  entries<-lapply(files[file.exists(files)],readRDS)
  if(!length(entries))return(invisible(NULL))
  overview<-do.call(rbind,lapply(entries,function(x)data.frame(task_id=x$task$id,scope=x$task$scope,MC_ID=x$task$mc,
    kind=x$task$kind,omitted=x$task$label,status=x$status,error=x$error)))
  js_csv(overview,file.path(root,'00_task_progress.csv'))
  rows<-lapply(entries,function(x)x$coefficients);rows<-rows[!vapply(rows,is.null,logical(1))]
  if(length(rows)){
    all<-do.call(rbind,rows);js_csv(all,file.path(root,'01_all_coefficient_changes.csv'))
    js_csv(all[all$original_q<.05,,drop=FALSE],file.path(root,'02_original_significant_changes.csv'))
    loto<-all[all$kind=='leave_TS',,drop=FALSE]
    if(nrow(loto)){
      summary<-do.call(rbind,lapply(split(loto,paste(loto$analysis_scope,loto$MC_ID,loto$term)),function(d){
        ok<-d$numerical_ok & is.finite(d$beta_per_TS_SD);v<-d$beta_per_TS_SD[ok]
        data.frame(scope=d$analysis_scope[1],MC_ID=d$MC_ID[1],term=d$term[1],original_beta=d$original_beta[1],original_q=d$original_q[1],
          completed_omissions=nrow(d),numerically_usable=sum(ok),beta_min=if(length(v))min(v) else NA,
          beta_max=if(length(v))max(v) else NA,sign_changes=sum(d$sign_changed[ok]),
          max_abs_delta_original_SE=if(length(v))max(abs(d$delta_in_original_SE[ok])) else NA)
      }));js_csv(summary,file.path(root,'03_LOTO_summary.csv'))
    }
  }
  cat('Completed checkpoints:',nrow(overview),'/',length(tasks),'\n');flush.console()
  invisible(overview)
}
js_main <- function(project,stage='check',workers=2L,resume='') {
  project<-normalizePath(project,winslash='/',mustWork=TRUE)
  stopifnot(stage %in% c('check','run'),workers>=1,workers<=8)
  if(!requireNamespace('glmmTMB',quietly=TRUE))stop('Missing glmmTMB')
  if(nzchar(resume)){
    root<-normalizePath(resume,winslash='/',mustWork=TRUE)
    if(!startsWith(root,paste0(project,'/results/joint_TS7_sensitivity_run')))stop('Resume path outside this workflow')
    ctx<-readRDS(file.path(root,'context.rds'))
    stopifnot(identical(unname(tools::md5sum(ctx$fingerprints$path)),ctx$fingerprints$md5))
  }else{ctx<-js_prepare(project);root<-''}
  cat('Validated sources; task count:',length(ctx$tasks),'(42 leave-TS plus targeted sample influence).\n');flush.console()
  if(stage=='check')return(invisible(ctx))
  if(!nzchar(root)){
    for(i in 1:9999){p<-file.path(project,'results',sprintf('joint_TS7_sensitivity_run%03d',i));if(!dir.exists(p)&&dir.create(p)){root<-p;break}}
    if(!nzchar(root))stop('Cannot create output')
    dir.create(file.path(root,'tasks'));js_save(ctx,file.path(root,'context.rds'))
    file.copy(file.path(project,'scripts/GLMM/GLMM_MC19_joint_TS7.R'),file.path(root,'joint_fit_source.R'))
    file.copy(file.path(project,'scripts/GLMM/joint_TS7_sensitivity.R'),file.path(root,'sensitivity_source.R'))
    js_csv(ctx$audit,file.path(root,'04_extreme_sample_source_audit.csv'))
    js_csv(ctx$fingerprints,file.path(root,'source_fingerprints.csv'))
    js_csv(do.call(rbind,lapply(ctx$tasks,function(t)data.frame(task_id=t$id,scope=t$scope,MC_ID=t$mc,kind=t$kind,omitted=t$label))),file.path(root,'task_manifest.csv'))
  }
  lock<-file.path(root,'.active_lock')
  if(!dir.create(lock,showWarnings=FALSE))stop('Run appears active; verify PID before removing .active_lock manually')
  stopifnot(identical(dirname(normalizePath(lock,winslash='/')),root))
  on.exit(unlink(lock,recursive=TRUE),add=TRUE)
  writeLines(as.character(Sys.getpid()),file.path(root,'controller_pid.txt'))
  writeLines('running',file.path(root,'run_status.txt'))
  status<-'interrupted_or_failed';on.exit(writeLines(status,file.path(root,'run_status.txt')),add=TRUE)
  writeLines(root,file.path(project,'results','00_latest_joint_TS7_sensitivity.txt'))
  writeLines(capture.output(sessionInfo()),file.path(root,'sessionInfo.txt'))
  writeLines(c('# Joint sensitivity checks','',
    '42 leave-TS checks: MC15/MC6/MC98, coral and ecosystem, all seven TS.',
    'MC95/MC43: top-three proportions union top-three count-residual contributions; individual and bundled omissions.',
    'Full original data and fits remain the primary reference. Omissions are sensitivity analyses, not data cleaning.',
    'All five predictors and original seven-TS scaling retained; no DIN. Six-TS design may be saturated.',
    'No new BH for this selected sensitivity set. Original q is copied for reference only.',
    'Read 03_LOTO_summary.csv for sign/range changes and numerical usability; inspect 04_extreme_sample_source_audit.csv.',
    'Task RDS checkpoints enable resume; failed or non-convergent tasks remain explicitly visible.'),file.path(root,'00_README.md'))
  pending<-ctx$tasks[!vapply(ctx$tasks,function(t)file.exists(file.path(root,'tasks',paste0(t$id,'.rds'))),logical(1))]
  Sys.setenv(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
  if(length(pending)&&workers>1L){
    cl<-parallel::makePSOCKcluster(min(workers,length(pending)),outfile=file.path(root,'worker.log'))
    on.exit(parallel::stopCluster(cl),add=TRUE,after=FALSE)
    parallel::clusterExport(cl,c('js_task','js_engine','js_save'),envir=environment(js_main))
    parallel::clusterEvalQ(cl,{Sys.setenv(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1');NULL})
  }else cl<-NULL
  batches<-if(length(pending))split(pending,ceiling(seq_along(pending)/workers)) else list()
  for(batch in batches){
    if(is.null(cl))lapply(batch,js_task,root=root) else parallel::parLapplyLB(cl,batch,js_task,root=root)
    js_collect(root,ctx$tasks)
  }
  overview<-js_collect(root,ctx$tasks)
  stopifnot(identical(unname(tools::md5sum(ctx$fingerprints$path)),ctx$fingerprints$md5))
  status<-if(all(overview$status=='complete'))'complete' else 'complete_with_review_or_failures'
  cat('Run:',root,'; status:',status,'\n')
  invisible(root)
}
if(!isTRUE(getOption('mc19.joint.sensitivity.define_only',FALSE))){
  args<-commandArgs(trailingOnly=TRUE)
  get_arg<-function(key,default=''){v<-args[startsWith(args,paste0('--',key,'='))];if(length(v))sub(paste0('^--',key,'='),'',v[1]) else default}
  project<-get_arg('project',getOption('mc19.joint.sensitivity.project_dir',''))
  if(!nzchar(project)){
    starts<-c(getwd(),dirname(sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1])))
    for(s in starts[!is.na(starts)])for(i in 1:7){
      candidates<-c(s,file.path(s,'Research','Msc main project-time series'))
      hit<-candidates[file.exists(file.path(candidates,'scripts/GLMM/GLMM_MC19_joint_TS7.R'))]
      if(length(hit)&&!nzchar(project))project<-hit[1];s<-dirname(s)
    }
  }
  js_main(project,if('--run' %in% args)'run' else 'check',as.integer(get_arg('workers','2')),get_arg('resume'))
}
