# Prior predictive simulation only: no likelihood updates, optimization or MCMC.
# CLI: Rscript joint_TS7_prior_predictive.R --run --draws=4000
# Default: validate current raw inputs against the completed env5 analysis.
pp_csv <- function(x, p) write.csv(x, p, row.names=FALSE, fileEncoding='UTF-8', na='')
pp_q <- function(x) unname(quantile(x, c(.025,.25,.5,.75,.975)))
pp_bb <- function(N, mu, phi) {
  # Guard only floating-point rounding of logistic probabilities to exactly 0/1.
  mu <- pmax(.Machine$double.eps, pmin(1-.Machine$double.eps, mu))
  rbinom(length(mu), N, rbeta(length(mu), mu*phi, (1-mu)*phi))
}
pp_stats <- function(x, ts) {
  tm <- vapply(paste0('TS',1:7), function(t) mean(x[ts==t]), numeric(1))
  c(mean_proportion=mean(x), TS_equal_mean=mean(tm), zero_fraction=mean(x==0),
    fraction_above_1pct=mean(x>.01), fraction_above_10pct=mean(x>.1),
    fraction_above_50pct=mean(x>.5), maximum_proportion=max(x),
    TS_mean_range=diff(range(tm)))
}
pp_prepare <- function(project) {
  e <- new.env(parent=globalenv())
  old <- options(mc19.joint.define_only=TRUE,mc19.joint.project_dir=project)
  on.exit(options(old)); sys.source(file.path(project,'scripts/GLMM/GLMM_MC19_joint_TS7.R'),e)
  contexts <- list(); files <- character()
  for(scope in c('coral','ecosystem')) {
    options(mc19.analysis_scope=scope)
    cfg <- e$joint_TS7_config(); a <- e$joint_TS7_load(cfg)
    baseline_path <- file.path(project,'results',paste0('joint_glmm_TS7_env5_',scope,'_run001'),
                               'supporting_files/analysis_input.csv')
    baseline <- read.csv(baseline_path,stringsAsFactors=FALSE)
    for(mc in cfg$mc) {
      b <- baseline[baseline$MC_ID==mc,,drop=FALSE]
      ix <- match(as.character(a$samples$New_ID),as.character(b$New_ID))
      stopifnot(nrow(b)==cfg$expected_n,!anyNA(ix),!anyDuplicated(b$New_ID),
        all(b$MC_reads[ix]==a$counts[mc,]),all(b$Ruegeria_reads_total[ix]==a$samples$Ruegeria_reads_total),
        all(as.character(b$TS[ix])==as.character(a$samples$TS)))
      for(z in cfg$z) stopifnot(max(abs(b[[z]][ix]-a$samples[[z]]))<1e-12)
    }
    contexts[[scope]] <- list(config=cfg,input=a)
    files <- c(files,a$fingerprints$path,baseline_path)
  }
  files <- unique(c(files,file.path(project,'scripts/GLMM/GLMM_MC19_joint_TS7.R'),
                    file.path(project,'scripts/GLMM/joint_TS7_prior_predictive.R')))
  list(contexts=contexts,fingerprints=data.frame(path=files,md5=unname(tools::md5sum(files))))
}
pp_simulate <- function(a,mc,tau,draws,seed) {
  set.seed(seed)
  s <- a$input$samples; ts <- as.character(s$TS); ti <- match(ts,paste0('TS',1:7))
  env <- a$input$environment_TS
  X <- as.matrix(env[match(paste0('TS',1:7),env$TS),a$config$z])
  N <- s$Ruegeria_reads_total
  alpha <- rnorm(draws,-4,2.5)
  sigma <- abs(rnorm(draws,0,1))
  phi <- exp(rnorm(draws,log(50),2))
  beta <- matrix(rnorm(draws*5),draws,5)*tau
  u <- matrix(rnorm(draws*7),draws,7)*sigma
  eta <- sweep(beta %*% t(X)+u,1,alpha,'+')
  mu <- plogis(eta)
  obs <- as.numeric(a$input$counts[mc,])/N
  stat_names <- names(pp_stats(obs,ts))
  stats <- matrix(NA_real_,draws,length(stat_names),dimnames=list(NULL,stat_names))
  ts_stats <- array(NA_real_,c(draws,7,4),dimnames=list(NULL,paste0('TS',1:7),
    c('mean_proportion','zero_fraction','maximum_proportion','mean_count')))
  for(k in seq_len(draws)) {
    y <- pp_bb(N,mu[k,ti],phi[k]); prop <- y/N
    stopifnot(all(is.finite(prop)),all(y>=0 & y<=N),all(y==floor(y)))
    stats[k,] <- pp_stats(prop,ts)
    for(t in 1:7) {
      keep <- ti==t
      ts_stats[k,t,] <- c(mean(prop[keep]),mean(y[keep]==0),max(prop[keep]),mean(y[keep]))
    }
  }
  pars <- data.frame(draw=seq_len(draws),alpha=alpha,sigma_TS=sigma,phi=phi,beta)
  names(pars)[5:9] <- a$config$z
  # Each MC is a separate one-versus-rest model, not a joint compositional draw.
  list(scope=a$config$scope,MC_ID=mc,tau=tau,seed=seed,draws=draws,parameters=pars,
       expected_p_TS=mu,simulated_statistics=stats,simulated_TS_statistics=ts_stats,
       observed_statistics=pp_stats(obs,ts),observed_proportion=obs,TS=ts,N=N,
       probability_rounding_guard=sum(mu<.Machine$double.eps | mu>1-.Machine$double.eps))
}
pp_summary <- function(x) do.call(rbind,lapply(colnames(x$simulated_statistics),function(stat) {
  v <- x$simulated_statistics[,stat]; q <- pp_q(v); obs <- x$observed_statistics[[stat]]
  data.frame(scope=x$scope,MC_ID=x$MC_ID,tau=x$tau,statistic=stat,
    observed=obs,prior_q025=q[1],prior_q25=q[2],prior_median=q[3],prior_q75=q[4],prior_q975=q[5],
    observed_outside_prior_95=obs<q[1]||obs>q[5],draws=x$draws,
    interpretation='descriptive_prior_check_not_a_test_or_acceptance_rule')
}))
pp_TS_summary <- function(x) do.call(rbind,lapply(1:7,function(t) {
  keep <- x$TS==paste0('TS',t)
  obs <- c(mean(x$observed_proportion[keep]),mean(x$observed_proportion[keep]==0),
           max(x$observed_proportion[keep]),mean(x$observed_proportion[keep]*x$N[keep]))
  do.call(rbind,lapply(1:4,function(j) {
    q <- pp_q(x$simulated_TS_statistics[,t,j])
    data.frame(scope=x$scope,MC_ID=x$MC_ID,tau=x$tau,TS=paste0('TS',t),
      statistic=dimnames(x$simulated_TS_statistics)[[3]][j],observed=obs[j],
      prior_q025=q[1],prior_q25=q[2],prior_median=q[3],prior_q75=q[4],prior_q975=q[5])
  }))
}))
pp_plot <- function(objects,path) {
  pdf(path,width=12,height=10); on.exit(dev.off())
  par(mfrow=c(3,3),mar=c(3.8,4,2.6,1),oma=c(0,0,3,0))
  for(x in objects) {
    plot(1:7,rep(0,7),type='n',ylim=c(0,1),xaxt='n',xlab='TS',ylab='Expected proportion',
         main=paste('tau =',x$tau,'| latent mean p'))
    axis(1,1:7,paste0('TS',1:7))
    for(t in 1:7) {
      q <- pp_q(x$expected_p_TS[,t]);segments(t,q[1],t,q[5],col='grey55');
      segments(t,q[2],t,q[4],lwd=5,col='#72A3BB');points(t,q[3],pch=16)
    }
    for(stat in c('zero_fraction','maximum_proportion')) {
      hist(x$simulated_statistics[,stat],breaks=30,xlim=c(0,1),col='grey85',border='white',
           main=paste('tau =',x$tau,'|',stat),xlab='Prior predictive dataset statistic')
      abline(v=x$observed_statistics[[stat]],col='#B53A35',lwd=2)
    }
  }
  mtext(paste(objects[[1]]$scope,objects[[1]]$MC_ID,
    '| prior only; bars 50% / 95%; red = observed reference'),outer=TRUE,cex=1)
}
pp_main <- function() {
  args <- commandArgs(TRUE)
  value <- function(key,default) {v<-grep(paste0('^--',key,'='),args,value=TRUE);if(length(v))sub(paste0('^--',key,'='),'',v[1]) else default}
  if('--self-test' %in% args) {
    set.seed(1931); y <- pp_bb(100,rep(.2,100000),20)
    stopifnot(abs(mean(y)-20)<.2,abs(var(y)/(100*.2*.8*120/21)-1)<.04)
    stopifnot(all(is.finite(pp_bb(100,c(0,1),20))))
    cat('SELF_TEST_OK: beta-binomial mean/variance and boundary handling\n');return(invisible(NULL))
  }
  project <- normalizePath(value('project',getwd()),winslash='/',mustWork=TRUE)
  if(!file.exists(file.path(project,'scripts/GLMM/GLMM_MC19_joint_TS7.R'))) {
    f <- sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1])
    project <- normalizePath(file.path(dirname(f),'../..'),winslash='/',mustWork=TRUE)
  }
  draws <- as.integer(value('draws','4000'));stopifnot(is.finite(draws),draws>=100L)
  ctx <- pp_prepare(project)
  cat('INPUT_CHECK_OK: coral 98 / ecosystem 310; 19 MC; five predictors; seven TS; no DIN\n')
  if(!'--run' %in% args) return(invisible(ctx))
  for(i in 1:9999) {
    root <- file.path(project,'results',sprintf('joint_TS7_prior_predictive_run%03d',i))
    if(dir.create(root,showWarnings=FALSE)) break
    if(i==9999)stop('No available output directory')
  }
  writeLines('running',file.path(root,'run_status.txt'))
  final_status <- 'failed'
  on.exit(writeLines(final_status,file.path(root,'run_status.txt')),add=TRUE)
  dir.create(file.path(root,'tasks'));dir.create(file.path(root,'figures'));dir.create(file.path(root,'supporting_files'))
  saveRDS(ctx,file.path(root,'supporting_files/input_context.rds'))
  pp_csv(ctx$fingerprints,file.path(root,'supporting_files/input_fingerprints.csv'))
  pp_csv(data.frame(tau=c(.5,1,2),alpha_mean=-4,alpha_SD=2.5,sigma_TS_halfnormal_SD=1,
    log_phi_mean=log(50),log_phi_SD=2,draws=draws),file.path(root,'00_candidate_priors.csv'))
  for(scope in names(ctx$contexts)) {
    a<-ctx$contexts[[scope]]
    pp_csv(a$input$scaling,file.path(root,'supporting_files',paste0(scope,'_scaling.csv')))
    pp_csv(a$input$environment_TS,file.path(root,'supporting_files',paste0(scope,'_environment_TS.csv')))
  }
  file.copy(file.path(project,'scripts/GLMM/joint_TS7_prior_predictive.R'),file.path(root,'supporting_files'),overwrite=FALSE)
  capture.output(sessionInfo(),file=file.path(root,'supporting_files/sessionInfo.txt'))
  writeLines(root,file.path(project,'results/00_latest_joint_TS7_prior_predictive.txt'))
  writeLines(c('Prior predictive simulation; no posterior fitting or new P/q.',
    'Candidate priors are uncalibrated. Observations are reference overlays only.',
    'Beta-binomial: q ~ Beta(mu*phi,(1-mu)*phi); Y ~ Binomial(N,q).',
    'One shared TS intercept per draw and TS, then independent sample-level Beta draws.',
    'Alpha is defined at the equal-TS standardized environmental zero.',
    'Separate MC models are NOT jointly normalized composition simulations.',
    '1%, 10%, 50% thresholds and observed-outside-95 flags are descriptive, not rejection rules.'),
    file.path(root,'README.txt'))
  summaries<-list();ts_summaries<-list();progress<-list();idx<-0L
  for(scope in names(ctx$contexts)) for(mc in ctx$contexts[[scope]]$config$mc) {
    a<-ctx$contexts[[scope]];objects<-list()
    # Common parameter random numbers across tau; independent draws across MC/scope.
    seed<-271828L+match(mc,a$config$mc)*100L+match(scope,names(ctx$contexts))*10000L
    for(tau in c(.5,1,2)) {
      idx<-idx+1L;cat(sprintf('[%d/114] %s %s tau=%s\n',idx,scope,mc,tau));flush.console()
      x<-pp_simulate(a,mc,tau,draws,seed)
      task<-sprintf('%s_%s_tau%s',scope,mc,tau)
      saveRDS(x,file.path(root,'tasks',paste0(task,'.rds')))
      summaries[[idx]]<-pp_summary(x);ts_summaries[[idx]]<-pp_TS_summary(x)
      objects[[length(objects)+1L]]<-x
      progress[[idx]]<-data.frame(task=task,status='complete',draws=draws,
                                  probability_rounding_guard=x$probability_rounding_guard)
      pp_csv(do.call(rbind,progress),file.path(root,'01_progress.csv'))
    }
    pp_plot(objects,file.path(root,'figures',paste0(scope,'_',mc,'.pdf')))
  }
  stopifnot(idx==114L,identical(unname(tools::md5sum(ctx$fingerprints$path)),ctx$fingerprints$md5))
  pp_csv(do.call(rbind,summaries),file.path(root,'02_prior_predictive_summary.csv'))
  pp_csv(do.call(rbind,ts_summaries),file.path(root,'03_prior_predictive_by_TS.csv'))
  final_status<-'complete';cat('COMPLETE: ',root,'\n',sep='')
  invisible(root)
}
if(!isTRUE(getOption('mc19.prior.define_only',FALSE)))pp_main()
