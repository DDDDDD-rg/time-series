# Independent two-part exploratory follow-up. Sourcing this helper never runs models.
mc19_tp_input <- function(project_dir,scope,source_dir='') {
  e <- new.env(parent=globalenv())
  for(f in c('glmm_MC19_pilot_helpers.R','glmm_MC19_wald_helpers.R')) sys.source(file.path(project_dir,'scripts','GLMM',f),e)
  run <- e$mc19_wald_read_run(e$mc19_wald_config(project_dir,scope),source_dir)
  s <- run$results; s <- s[s$numerical_ok %in% TRUE & is.finite(s$beta_q_value) & s$beta_q_value<.05,]
  s$key <- paste(s$MC_ID,s$environment,sep='__');stopifnot(nrow(s)>0,!anyDuplicated(s$key))
  path <- file.path(run$support,'analysis_input.csv');d <- read.csv(path,stringsAsFactors=FALSE)
  d$key <- paste(d$MC_ID,d$environment,sep='__');d <- d[d$key %in% s$key,]
  stopifnot(!anyDuplicated(d[c('New_ID','MC_ID','environment')]),all(is.finite(d$E_z)),
    all(d$MC_reads>=0),all(d$MC_reads<=d$Ruegeria_reads_total),all(d$Ruegeria_reads_total>0),
    all(abs(d$MC_proportion-d$MC_reads/d$Ruegeria_reads_total)<1e-10))
  stopifnot(all(d$MC_proportion[d$MC_reads>0]>0 & d$MC_proportion[d$MC_reads>0]<1))
  stopifnot(all(table(d$key)==unique(s$n_samples)),length(unique(d$TS))==7)
  if(scope=='coral') stopifnot(all(tolower(d$Host_niche)=='coral')) else
    stopifnot(!any(tolower(d$Host_niche) %in% c('env','seawater','sediment')))
  for(k in s$key) {
    z <- d[d$key==k,];u <- unique(z[c('TS','E_raw','E_z')]);stopifnot(nrow(u)==7,!anyDuplicated(u$TS))
  }
  fingerprints <- tools::md5sum(c(path,file.path(run$support,'all_results.csv')))
  list(selected=s,data=d,source=run$run_dir,fingerprints=fingerprints,scope=scope)
}
mc19_tp_fit <- function(d,component) {
  d$TS <- factor(d$TS); d$Detected <- as.integer(d$MC_reads>0)
  if(component=='positive') d <- droplevels(d[d$Detected==1,])
  ans <- data.frame(component=component,n=nrow(d),n_TS=length(unique(d$TS)),n_zero=sum(d$MC_reads==0),
    beta=NA_real_,SE=NA_real_,CI95_low=NA_real_,CI95_high=NA_real_,P=NA_real_,q=NA_real_,
    numerical_ok=FALSE,TS_SD=NA_real_,gradient=NA_real_,boundary=FALSE,status='not_fitted',warnings='',error='')
  warnings <- character(); fit <- NULL; predictions <- residual_data <- NULL
  tryCatch(withCallingHandlers({
    if(nrow(d)<4 || length(unique(d$TS))<2 || length(unique(d$E_z))<2) stop('Insufficient positive/sample/TS/environment support')
    if(component=='detection' && length(unique(d$Detected))<2) stop('Detection outcome has no variation')
    formula <- if(component=='detection') Detected ~ E_z+(1|TS) else MC_proportion ~ E_z+(1|TS)
    family <- if(component=='detection') stats::binomial('logit') else glmmTMB::beta_family('logit')
    fit <- glmmTMB::glmmTMB(formula,data=d,family=family,ziformula=~0,dispformula=~1)
    co <- summary(fit)$coefficients$cond['E_z',]
    ans$beta <- co[1];ans$SE <- co[2];ans$CI95_low <- co[1]-qnorm(.975)*co[2];ans$CI95_high <- co[1]+qnorm(.975)*co[2]
    ans$gradient <- tryCatch(max(abs(fit$obj$gr(fit$fit$par))),error=function(e) NA_real_)
    ans$TS_SD <- unname(attr(glmmTMB::VarCorr(fit)$cond$TS,'stddev')[1])
    ans$boundary <- ans$TS_SD<1e-4
    ans$numerical_ok <- isTRUE(fit$fit$convergence==0) && isTRUE(fit$sdr$pdHess) && all(is.finite(co)) && co[2]>0
    if(ans$numerical_ok) ans$P <- co[4]
    ans$status <- if(ans$numerical_ok) 'fitted' else 'numerical_review'
    ez <- seq(min(d$E_z),max(d$E_z),length.out=101)
    scale_rows <- unique(d[c('E_raw','E_z')]);scale_rows <- scale_rows[order(scale_rows$E_z),]
    g <- data.frame(E_z=ez,TS=factor(rep(levels(d$TS)[1],length(ez)),levels=levels(d$TS)))
    pred <- predict(fit,newdata=g,re.form=NA,type='link',se.fit=TRUE)
    predictions <- data.frame(E_z=ez,E_raw=approx(scale_rows$E_z,scale_rows$E_raw,ez,rule=2)$y,
      mean=plogis(pred$fit),low=plogis(pred$fit-qnorm(.975)*pred$se.fit),high=plogis(pred$fit+qnorm(.975)*pred$se.fit))
    residual_data <- data.frame(New_ID=d$New_ID,TS=d$TS,E_raw=d$E_raw,Ruegeria_reads_total=d$Ruegeria_reads_total,
      observed=if(component=='detection') d$Detected else d$MC_proportion,
      fitted=as.numeric(predict(fit,type='response')),pearson=as.numeric(residuals(fit,type='pearson')))
    if(component=='detection' && any(residual_data$fitted<1e-6 | residual_data$fitted>1-1e-6)) warnings <- c(warnings,'Extreme fitted detection probabilities; review separation')
    if(abs(ans$beta)>10 || ans$SE>10) warnings <- c(warnings,'Large coefficient/SE; review separation or weak identification')
    if(is.finite(ans$gradient) && ans$gradient>.001) warnings <- c(warnings,'Gradient review')
  },warning=function(w){warnings <<- c(warnings,conditionMessage(w));invokeRestart('muffleWarning')}),
  error=function(e){ans$error <<- conditionMessage(e);ans$status <<- 'failed'})
  ans$warnings <- paste(unique(warnings),collapse=' | ')
  if(nzchar(ans$error)) {ans$numerical_ok <- FALSE;ans$P <- NA_real_}
  list(result=ans,fit=fit,predictions=predictions,residuals=residual_data)
}
mc19_tp_forest_labels <- function(tab) {
  stopifnot(!anyDuplicated(tab[c('key','component')]),
    all(tab$component %in% c('detection','positive')))
  rows <- unique(tab[c('key','MC_ID','environment')])
  for(component in c('detection','positive')) {
    z <- tab[tab$component==component,]
    j <- match(rows$key,z$key)
    rows[[paste0(component,'_q')]] <- z$q[j]
    rows[[paste0(component,'_significant')]] <- z$numerical_ok[j] %in% TRUE &
      is.finite(z$q[j]) & z$q[j]<.05
  }
  rows$symbol <- paste0(ifelse(rows$detection_significant,'*',''),
    ifelse(rows$positive_significant,'#',''))
  labels <- c(Nitrate_mgN_L='NO3-N',Nitrite_mgN_L='NO2-N',Ammonia_mgN_L='Ammonia-N',DIP_mgP_L='DIP')
  env <- unname(labels[rows$environment]);env[is.na(env)] <- rows$environment[is.na(env)]
  rows$label <- paste0(rows$MC_ID,ifelse(nzchar(rows$symbol),paste0(' ',rows$symbol),''),' / ',env)
  rows[order(rows$MC_ID,match(rows$environment,names(labels))),]
}
mc19_tp_forest <- function(tab) {
  rows <- mc19_tp_forest_labels(tab)
  tab$label <- factor(rows$label[match(tab$key,rows$key)],levels=rev(rows$label))
  dodge <- ggplot2::position_dodge(width=.55)
  family <- unique(tab$q_family_size);stopifnot(length(family)==1L)
  ggplot2::ggplot(tab,ggplot2::aes(beta,label,colour=component,group=component))+
    ggplot2::geom_vline(xintercept=0,linetype=3,linewidth=.65,colour='grey40')+
    ggplot2::geom_errorbar(ggplot2::aes(xmin=CI95_low,xmax=CI95_high),
      orientation='y',width=.18,linewidth=.9,position=dodge,na.rm=TRUE)+
    ggplot2::geom_point(size=3.3,position=dodge,na.rm=TRUE)+
    ggplot2::scale_colour_discrete(name=NULL,breaks=c('detection','positive'),
      labels=c('Detection (det)','Positive proportion (pos)'))+
    ggplot2::scale_y_discrete(expand=ggplot2::expansion(add=.7))+
    ggplot2::labs(x='Log-odds slope per TS environment SD',y=NULL,
      caption=paste0('MC symbols: * det q < 0.05; # pos q < 0.05; *# both. No symbol: neither.\n',
        'Bars: pointwise 95% Wald CI. BH q across ',family,' selected component tests.'))+
    ggplot2::theme_classic(base_size=24)+
    ggplot2::theme(text=ggplot2::element_text(face='bold'),
      axis.text=ggplot2::element_text(size=24,face='bold',colour='black'),
      axis.title.x=ggplot2::element_text(size=24,face='bold',margin=ggplot2::margin(t=12)),
      legend.position='top',legend.text=ggplot2::element_text(size=20,face='bold'),
      plot.caption=ggplot2::element_text(size=17,face='bold',hjust=0,lineheight=1.2,margin=ggplot2::margin(t=16)),
      plot.caption.position='plot',plot.margin=ggplot2::margin(16,24,18,18))
}
mc19_tp_theme <- function() {
  ggplot2::theme_classic(base_size=21)+ggplot2::theme(
    text=ggplot2::element_text(face='bold'),axis.text=ggplot2::element_text(size=18,face='bold',colour='black'),
    axis.title=ggplot2::element_text(size=21,face='bold'),strip.text=ggplot2::element_text(size=20,face='bold'),
    plot.subtitle=ggplot2::element_text(size=20,face='bold',margin=ggplot2::margin(b=14)),
    plot.caption=ggplot2::element_text(size=15,face='bold',hjust=0,lineheight=1.2,margin=ggplot2::margin(t=14)),
    plot.caption.position='plot',panel.spacing=grid::unit(1.2,'lines'),plot.margin=ggplot2::margin(18,24,18,18))
}
mc19_tp_plot <- function(out,project_dir=dirname(dirname(out)),parts=c('forest','curves','residual')) {
  stopifnot(requireNamespace('ggplot2',quietly=TRUE))
  stopifnot(length(parts)>0,all(parts %in% c('forest','curves','residual')))
  tab <- read.csv(file.path(out,'supporting_files','all_results.csv'))
  pred <- readRDS(file.path(out,'supporting_files','predictions.rds'))
  obs <- readRDS(file.path(out,'supporting_files','observations.rds'))
  scope <- readRDS(file.path(out,'supporting_files','manifest.rds'))$scope
  if(!mc19_result_complete(out,'two_part',scope)) stop('Two-part run is incomplete.')
  fig <- mc19_plot_directory(project_dir,out,'two_part',scope,getOption('mc19.two_part.figure_dir',NULL))
  if('forest' %in% parts) {
    p <- mc19_tp_forest(tab)
    height <- 11.25
    ggplot2::ggsave(file.path(fig,'two_part_effects.pdf'),p,width=20,height=height)
    ggplot2::ggsave(file.path(fig,'two_part_effects.png'),p,width=20,height=height,dpi=180,bg='white')
    utils::write.csv(mc19_tp_forest_labels(tab),file.path(fig,'two_part_effects_labels.csv'),row.names=FALSE)
  }
  if('curves' %in% parts) {
  grDevices::pdf(file.path(fig,'two_part_curves.pdf'),width=15,height=8)
  for(k in unique(tab$key)) {
    z <- pred[pred$key==k,];d <- obs[obs$key==k,]
    if(!nrow(z)) next
    p <- ggplot2::ggplot(z,ggplot2::aes(E_raw,mean))+ggplot2::geom_ribbon(ggplot2::aes(ymin=low,ymax=high),alpha=.18)+
      ggplot2::geom_line()+ggplot2::geom_point(data=d,ggplot2::aes(E_raw,observed),alpha=.25,inherit.aes=FALSE)+
      ggplot2::facet_wrap(~component,scales='free_y',labeller=ggplot2::as_labeller(c(detection='Detection (det)',positive='Positive proportion (pos)')))+ggplot2::labs(x=unique(tab$environment[tab$key==k]),y='Probability / positive proportion',
        subtitle=gsub('__',' / ',k,fixed=TRUE),caption='Curves at TS random intercept zero; ribbon: pointwise 95% Wald CI.\nAll relevant raw observations retained.')+mc19_tp_theme()
    print(p)
  }
  grDevices::dev.off()
  }
  if('residual' %in% parts && nrow(obs)) {
    p <- ggplot2::ggplot(obs,ggplot2::aes(fitted,pearson))+ggplot2::geom_point(alpha=.25)+
      ggplot2::facet_wrap(~component,scales='free',labeller=ggplot2::as_labeller(c(detection='Detection (det)',positive='Positive proportion (pos)')))+mc19_tp_theme()
    ggplot2::ggsave(file.path(fig,'residual_review.pdf'),p,width=14,height=8)
  }
  mc19_io_update(fig,parts=paste(parts,collapse=','))
  mc19_publish_plot(project_dir,'two_part',scope,fig)
  invisible(fig)
}
mc19_tp_run <- function(input,project_dir,out) {
  stopifnot(requireNamespace('glmmTMB',quietly=TRUE))
  dir.create(file.path(out,'supporting_files','models'),recursive=TRUE,showWarnings=FALSE)
  support <- file.path(out,'supporting_files')
  script <- file.path(project_dir,'scripts','GLMM','mc19_two_part_helpers.R')
  signature <- list(scope=input$scope,source=input$source,hashes=input$fingerprints,keys=input$selected$key,
    helper_hash=unname(tools::md5sum(script)),
    io_hash=unname(tools::md5sum(file.path(project_dir,'scripts','mc19_run_io.R'))),R=R.version.string,glmmTMB=as.character(utils::packageVersion('glmmTMB')))
  manifest <- file.path(support,'manifest.rds')
  if(file.exists(manifest)) stopifnot(identical(readRDS(manifest),signature)) else saveRDS(signature,manifest)
  writeLines('running',file.path(out,'run_status.txt'))
  write.csv(input$selected,file.path(support,'selected_original_results.csv'),row.names=FALSE)
  write.csv(input$data,file.path(support,'analysis_input.csv'),row.names=FALSE)
  file.copy(script,file.path(support,basename(script)),overwrite=TRUE)
  file.copy(file.path(project_dir,'scripts','mc19_run_io.R'),file.path(support,'mc19_run_io.R'),overwrite=TRUE)
  writeLines(capture.output(sessionInfo()),file.path(support,'sessionInfo.txt'))
  ans <- preds <- observations <- list();tasks <- 2L*nrow(input$selected)
  for(i in seq_len(nrow(input$selected))) for(component in c('detection','positive')) {
    key <- input$selected$key[i];cache <- file.path(support,'models',paste0(key,'__',component,'.rds'))
    if(file.exists(cache)) fit <- readRDS(cache) else {
      message(i,'/',nrow(input$selected),' ',key,' ',component)
      fit <- mc19_tp_fit(input$data[input$data$key==key,],component)
      saveRDS(fit,cache)
    }
    row <- fit$result;row$key <- key;row$MC_ID <- input$selected$MC_ID[i];row$environment <- input$selected$environment[i]
    ans[[length(ans)+1L]] <- row
    for(type in c('predictions','residuals')) if(!is.null(fit[[type]])) {
      z <- fit[[type]];z$key <- key;z$component <- component
      if(type=='predictions') preds[[length(preds)+1L]] <- z else observations[[length(observations)+1L]] <- z
    }
  }
  tab <- do.call(rbind,ans);tab$q <- p.adjust(tab$P,'BH',n=tasks);tab$q_family_size <- tasks
  tab$OR <- exp(tab$beta)
  write.csv(tab,file.path(support,'all_results.csv'),row.names=FALSE)
  # Retain every selected combination even when one component fails.
  wide <- input$selected[c('key','MC_ID','environment','beta_per_TS_SD','CI95_low','CI95_high','beta_P_value','beta_q_value')]
  for(component in c('detection','positive')) {
    z <- tab[tab$component==component,];cols <- c('beta','CI95_low','CI95_high','P','q','status','warnings','error')
    for(n in cols) wide[[paste(component,n,sep='_')]] <- z[[n]][match(wide$key,z$key)]
  }
  write.csv(wide,file.path(out,'00_result.csv'),row.names=FALSE)
  bind <- function(x) if(length(x)) do.call(rbind,x) else data.frame()
  saveRDS(bind(preds),file.path(support,'predictions.rds'));saveRDS(bind(observations),file.path(support,'observations.rds'))
  byts <- aggregate(Detected ~ MC_ID+TS,data=unique(transform(input$data,Detected=as.integer(MC_reads>0))[c('New_ID','MC_ID','TS','Detected')]),
    FUN=function(x) c(detected=sum(x),total=length(x)))
  write.csv(byts,file.path(support,'detection_by_TS.csv'),row.names=FALSE)
  writeLines(c(paste('Source:',input$source),paste('Fixed family size:',tasks),
    'Exploratory decomposition after selection in the same data; within-selected BH is not whole-pipeline FDR control.',
    'Detection: Bernoulli/logit; positive: Beta/logit; both E_z + (1|TS); original denominator and E_z retained.',
    'No depth/Host/Repeat/SST covariate added. Inspect depth and residuals. Detection is observed sequencing detection.',
    'Predictions use TS random intercept zero, not population-integrated means.'),file.path(out,'README.txt'))
  writeLines(if(all(tab$numerical_ok)) 'fits_complete_review_diagnostics' else 'fits_complete_some_require_review',file.path(out,'run_status.txt'))
  mc19_publish_result(project_dir,'two_part',input$scope,out)
  invisible(out)
}
