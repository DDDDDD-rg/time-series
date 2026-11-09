# Independent two-part exploratory follow-up. Sourcing this helper never runs models.
mc19_tp_input <- function(project_dir,scope,source_dir='') {
  e <- new.env(parent=globalenv())
  for(f in c('glmm_MC19_pilot_helpers.R','glmm_MC19_wald_helpers.R')) sys.source(file.path(project_dir,'scripts','GLMM',f),e)
  run <- e$mc19_wald_read_run(e$mc19_wald_config(project_dir,scope),source_dir)
  s <- run$results; s <- s[s$numerical_ok %in% TRUE & is.finite(s$beta_q_value) & s$beta_q_value<.05,]
  if(!nrow(s)) stop('新的8TS原GLMM没有BH q<0.05组合，无需运行two-part。',call.=FALSE)
  s$key <- paste(s$MC_ID,s$environment,sep='__');stopifnot(!anyDuplicated(s$key))
  path <- file.path(run$support,'analysis_input.csv');d <- read.csv(path,stringsAsFactors=FALSE)
  d$key <- paste(d$MC_ID,d$environment,sep='__');d <- d[d$key %in% s$key,]
  stopifnot(!anyDuplicated(d[c('New_ID','MC_ID','environment')]),all(is.finite(d$E_z)),
    all(d$MC_reads>=0),all(d$MC_reads<=d$Ruegeria_reads_total),all(d$Ruegeria_reads_total>0),
    all(abs(d$MC_proportion-d$MC_reads/d$Ruegeria_reads_total)<1e-10))
  stopifnot(all(d$MC_proportion[d$MC_reads>0]>0 & d$MC_proportion[d$MC_reads>0]<1))
  stopifnot(all(table(d$key)==unique(s$n_samples)),length(unique(d$TS))==8)
  if(scope=='coral') stopifnot(all(tolower(d$Host_niche)=='coral')) else
    stopifnot(!any(tolower(d$Host_niche) %in% c('env','seawater','sediment')))
  for(k in s$key) {
    z <- d[d$key==k,];u <- unique(z[c('TS','E_raw','E_z')]);stopifnot(nrow(u)==8,!anyDuplicated(u$TS))
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
  labels <- c(Nitrate_mgN_L='NO3-N',Nitrite_mgN_L='NO2-N',Ammonia_mgN_L='Ammonia-N',DIN_mgN_L='DIN',DIP_mgP_L='DIP')
  env <- unname(labels[rows$environment]);env[is.na(env)] <- rows$environment[is.na(env)]
  rows$label <- paste0(rows$MC_ID,ifelse(nzchar(rows$symbol),paste0(' ',rows$symbol),''),' / ',env)
  # Keep the focal MCs together at the top; include only rows actually analysed.
  priority <- match(rows$MC_ID,c('MC15','MC98','MC80','MC0','MC10'))
  priority[is.na(priority)] <- 6L
  rows[order(priority,rows$MC_ID,match(rows$environment,names(labels))),]
}
mc19_tp_forest <- function(tab) {
  rows <- mc19_tp_forest_labels(tab)
  tab$label <- factor(rows$label[match(tab$key,rows$key)],levels=rev(rows$label))
  family <- unique(tab$q_family_size);stopifnot(length(family)==1L)
  tab$component <- factor(tab$component,levels=c('detection','positive'))
  ggplot2::ggplot(tab,ggplot2::aes(beta,label,colour=component,group=component))+
    ggplot2::geom_vline(xintercept=0,linetype=3,linewidth=.65,colour='grey40')+
    ggplot2::geom_errorbar(ggplot2::aes(xmin=CI95_low,xmax=CI95_high),
      orientation='y',width=.18,linewidth=.9,na.rm=TRUE)+
    ggplot2::geom_point(size=3.3,na.rm=TRUE)+
    ggplot2::facet_grid(cols=ggplot2::vars(component),scales='free_x',
      labeller=ggplot2::as_labeller(c(detection='Detection (det)',positive='Positive proportion (pos)')))+
    ggplot2::scale_colour_manual(values=c(detection='#DF6C68',positive='#00A6AE'),guide='none')+
    ggplot2::scale_x_continuous(expand=ggplot2::expansion(mult=.07))+
    ggplot2::scale_y_discrete(expand=ggplot2::expansion(add=.7))+
    ggplot2::labs(x='Log-odds slope per TS environment SD',y=NULL,
      caption=paste0('MC symbols: * det q < 0.05; # pos q < 0.05; *# both. No symbol: neither.\n',
        'Bars: pointwise 95% Wald CI. BH q across ',family,' selected component tests.'))+
    ggplot2::theme_classic(base_size=24)+
    ggplot2::theme(text=ggplot2::element_text(face='bold'),
      axis.text=ggplot2::element_text(size=24,face='bold',colour='black'),
      axis.title.x=ggplot2::element_text(size=24,face='bold',margin=ggplot2::margin(t=12)),
      strip.text=ggplot2::element_text(size=22,face='bold'),strip.background=ggplot2::element_rect(fill='white'),
      panel.spacing=grid::unit(1.2,'lines'),
      plot.caption=ggplot2::element_text(size=17,face='bold',hjust=0,lineheight=1.2,margin=ggplot2::margin(t=16)),
      plot.caption.position='plot',plot.margin=ggplot2::margin(16,24,18,18))
}
mc19_tp_curve_page <- function(tab,pred,obs,key) {
  rows <- tab[tab$key==key,];z <- pred[pred$key==key,];d <- obs[obs$key==key,]
  field <- unique(rows$environment);mc <- unique(rows$MC_ID)
  stopifnot(length(field)==1L,length(mc)==1L)
  axis_label <- unname(c(Nitrate_mgN_L='Nitrate-N (mg N/L)',Nitrite_mgN_L='Nitrite-N (mg N/L)',
    Ammonia_mgN_L='Ammonia-N (mg N/L)',DIN_mgN_L='DIN (mg N/L)',DIP_mgP_L='DIP (mg P/L)')[field])
  # A common full/zoom range for the same MC across its nutrient pages.
  mc_keys <- tab$key[tab$MC_ID==mc]
  op <- obs[obs$key %in% mc_keys & obs$component=='positive',]
  pp <- pred[pred$key %in% mc_keys & pred$component=='positive',]
  full <- max(100*c(op$observed,pp$high))*1.08
  focus <- max(100*c(stats::quantile(op$observed,.85),pp$high))
  zoom <- min(full,focus*1.12)
  use_zoom <- full>2.5*focus && is.finite(zoom) && zoom>0
  fmt <- function(x) if(!is.finite(x)) 'NA' else formatC(x,format='g',digits=3)
  panel <- function(component,local_zoom=FALSE) {
    zz <- z[z$component==component,];dd <- d[d$component==component,]
    result <- rows[rows$component==component,];stopifnot(nrow(result)==1L)
    multiplier <- if(component=='positive') 100 else 1
    for(n in c('mean','low','high')) zz[[n]] <- multiplier*zz[[n]]
    dd$observed <- multiplier*dd$observed
    ymax <- if(component=='detection') 1 else if(local_zoom) zoom else full
    if(!local_zoom) stopifnot(all(dd$observed<=ymax),all(zz$high<=ymax))
    title <- if(component=='detection') 'Detection (det)' else if(local_zoom) 'Positive proportion: zoom' else 'Positive proportion: full range'
    subtitle <- if(local_zoom) 'Same data and curve; expanded y-axis' else
      paste0('Beta = ',fmt(result$beta),'; BH q = ',fmt(result$q),'; n = ',nrow(dd))
    p <- ggplot2::ggplot(zz,ggplot2::aes(E_raw,mean))
    if(field %in% c('Nitrite_mgN_L','DIP_mgP_L'))
      p <- p+ggplot2::geom_vline(xintercept=.002,linetype=2,colour='grey40',linewidth=.75)
    p <- p+ggplot2::geom_ribbon(ggplot2::aes(ymin=low,ymax=high),fill='#C9CED4',alpha=.70)+
      ggplot2::geom_point(data=dd,ggplot2::aes(E_raw,observed),inherit.aes=FALSE,alpha=.35,size=2.2)+
      ggplot2::geom_line(linewidth=.95,colour='#20272D')+
      ggplot2::coord_cartesian(ylim=if(component=='detection') c(-.035,1.035) else c(0,ymax))+
      ggplot2::labs(x=axis_label,y=if(component=='detection') 'Detection probability' else 'Positive relative abundance (%)',
        title=title,subtitle=subtitle)+mc19_tp_theme()+
      ggplot2::theme(plot.title=ggplot2::element_text(size=21,face='bold'),
        plot.subtitle=ggplot2::element_text(size=17,face='bold',margin=ggplot2::margin(b=12)),
        plot.margin=ggplot2::margin(12,22,18,18),axis.title=ggplot2::element_text(size=20),
        axis.title.y=ggplot2::element_text(margin=ggplot2::margin(r=12)))
    if(component=='detection') p <- p+ggplot2::scale_y_continuous(breaks=seq(0,1,.25))
    p
  }
  columns <- if(use_zoom) 3L else 2L
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout=grid::grid.layout(3,columns,
    heights=grid::unit(c(.65,8.90,1.70),'in'))))
  grid::grid.text(paste(mc,'/',axis_label),x=grid::unit(.35,'in'),y=grid::unit(11.0,'in'),
    just='left',gp=grid::gpar(fontsize=25,fontface='bold'))
  print(panel('detection'),newpage=FALSE,vp=grid::viewport(layout.pos.row=2,layout.pos.col=1))
  print(panel('positive'),newpage=FALSE,vp=grid::viewport(layout.pos.row=2,layout.pos.col=2))
  if(use_zoom) print(panel('positive',TRUE),newpage=FALSE,vp=grid::viewport(layout.pos.row=2,layout.pos.col=3))
  footer <- paste0('Curves: TS random intercept = 0; ribbons: pointwise 95% Wald CI. All observations and complete CI retained in full-range panels.\n',
    'Positive abundance = MC reads / total Ruegeria reads; shown as %. Detection observations overlap at 0 and 1.',
    if(field %in% c('Nitrite_mgN_L','DIP_mgP_L')) paste0('\nDashed line: reporting limit = 0.002 ',if(field=='DIP_mgP_L') 'mg P/L.' else 'mg N/L.') else '')
  grid::grid.text(footer,x=grid::unit(.35,'in'),y=grid::unit(.65,'in'),just=c('left','center'),
    gp=grid::gpar(fontsize=17,fontface='bold',lineheight=1.25))
  grid::popViewport()
  data.frame(key=key,positive_full_ymax_percent=full,positive_zoom_ymax_percent=zoom,zoom_enabled=use_zoom,
    detection_points=sum(d$component=='detection'),positive_points=sum(d$component=='positive'),
    detection_prediction_points=sum(z$component=='detection'),positive_prediction_points=sum(z$component=='positive'))
}
mc19_tp_theme <- function() {
  ggplot2::theme_classic(base_size=21)+ggplot2::theme(
    text=ggplot2::element_text(face='bold'),axis.text=ggplot2::element_text(size=18,face='bold',colour='black'),
    axis.title=ggplot2::element_text(size=21,face='bold'),strip.text=ggplot2::element_text(size=20,face='bold'),
    plot.subtitle=ggplot2::element_text(size=20,face='bold',margin=ggplot2::margin(b=14)),
    plot.caption=ggplot2::element_text(size=15,face='bold',hjust=0,lineheight=1.2,margin=ggplot2::margin(t=14)),
    plot.caption.position='plot',panel.spacing=grid::unit(1.2,'lines'),plot.margin=ggplot2::margin(18,24,18,18))
}
mc19_tp_analysis_hash <- function(path=NULL) {
  funcs <- c('mc19_tp_input','mc19_tp_fit')
  if(is.null(path)) bodies <- lapply(funcs,function(n) body(get(n,envir=environment(mc19_tp_analysis_hash)))) else {
    parsed <- parse(path); bodies <- setNames(vector('list',length(funcs)),funcs)
    for(expr in parsed) if(is.call(expr) && identical(expr[[1]],as.name('<-')) &&
        as.character(expr[[2]]) %in% funcs) bodies[[as.character(expr[[2]])]] <- expr[[3]][[3]]
    if(any(vapply(bodies,is.null,logical(1)))) stop('Missing analysis definitions in saved helper.')
  }
  canonical <- lapply(bodies,function(x) paste(deparse(x,width.cutoff=500L),collapse='\n'))
  tmp <- tempfile('mc19-signature-');on.exit(unlink(tmp));saveRDS(unname(canonical),tmp,version=2)
  unname(tools::md5sum(tmp))
}

# Deep, read-only consistency check, including legacy runs without seal hashes.
mc19_tp_validate <- function(out) {
  s <- file.path(out,'supporting_files');m <- readRDS(file.path(s,'manifest.rds'))
  tab <- read.csv(file.path(s,'all_results.csv'));wide <- read.csv(file.path(out,'00_result.csv'))
  input <- read.csv(file.path(s,'analysis_input.csv'));pred <- readRDS(file.path(s,'predictions.rds'));obs <- readRDS(file.path(s,'observations.rds'))
  same <- function(a,b,tol=1e-10) length(a)==length(b) && all((is.na(a)&is.na(b)) | (is.finite(a)&is.finite(b)&abs(a-b)<=tol*pmax(1,abs(b))))
  require <- function(ok,what) if(!isTRUE(ok)) stop('Two-part integrity: ',what,call.=FALSE)
  require(!anyDuplicated(m$keys) && setequal(unique(input$key),m$keys),'selected/input keys')
  require(all(tab$key==paste(tab$MC_ID,tab$environment,sep='__')) &&
    all(input$key==paste(input$MC_ID,input$environment,sep='__')) &&
    !anyDuplicated(input[c('key','New_ID')]),'MC/environment/sample identities')
  expected <- as.vector(outer(m$keys,c('detection','positive'),paste,sep='::'))
  require(!anyDuplicated(paste(tab$key,tab$component)) && setequal(paste(tab$key,tab$component,sep='::'),expected),'component task keys')
  require(!anyDuplicated(wide$key) && setequal(wide$key,m$keys),'primary keys')
  require(identical(unname(tools::md5sum(names(m$hashes))),unname(m$hashes)),'source input fingerprints')
  source_input <- read.csv(names(m$hashes)[basename(names(m$hashes))=='analysis_input.csv'])
  source_keys <- paste(source_input$New_ID,source_input$MC_ID,source_input$environment)
  ix <- match(paste(input$New_ID,input$MC_ID,input$environment),source_keys)
  require(!anyNA(ix),'original sample mapping')
  for(n in c('MC_reads','Other_reads','Ruegeria_reads_total','MC_proportion','E_raw','E_z'))
    require(same(input[[n]],source_input[[n]][ix]),paste('original input',n))
  require(identical(as.character(input$TS),as.character(source_input$TS[ix])),'original TS')
  seal <- file.path(s,'artifact_fingerprints.csv')
  if(file.exists(seal)) {
    hashes <- read.csv(seal);require(all(file.exists(file.path(out,hashes$file))) &&
      identical(unname(tools::md5sum(file.path(out,hashes$file))),hashes$md5),'saved artifact hashes')
  }
  require(same(tab$q,p.adjust(tab$P,'BH',n=length(expected))), 'fixed-family BH')
  require(all(tab$q_family_size==length(expected)), 'family size')
  require(setequal(paste(pred$key,pred$component,sep='::'),expected[expected %in% paste(pred$key,pred$component,sep='::')]) &&
    all(paste(obs$key,obs$component,sep='::') %in% expected), 'prediction/observation keys')
  for(i in seq_len(nrow(tab))) {
    r <- tab[i,];file <- file.path(s,'models',paste0(r$key,'__',r$component,'.rds'))
    require(file.exists(file),paste(r$key,r$component,'missing model'));cache <- readRDS(file)
    for(n in c('beta','SE','CI95_low','CI95_high','P')) require(same(r[[n]],cache$result[[n]]),paste(r$key,r$component,n))
    require(identical(r$numerical_ok,cache$result$numerical_ok), 'numerical state')
    w <- wide[match(r$key,wide$key),]
    require(identical(r$MC_ID,w$MC_ID) && identical(r$environment,w$environment),'MC/environment mapping')
    for(n in c('beta','CI95_low','CI95_high','P','q')) require(same(r[[n]],w[[paste(r$component,n,sep='_')]],1e-7),'primary/full precision mismatch')
    for(pair in list(c('predictions','pred'),c('residuals','obs'))) {
      stored <- if(pair[2]=='pred') pred else obs
      stored <- stored[stored$key==r$key & stored$component==r$component,,drop=FALSE]
      original <- cache[[pair[1]]]
      require(nrow(stored)==if(is.null(original))0L else nrow(original),'cached row count')
      if(!is.null(original)) for(n in names(original)) {
        a <- original[[n]];b <- stored[[n]]
        require(if(is.numeric(a)) same(a,b) else identical(as.character(a),as.character(b)),paste('cache field',n))
      }
    }
    if(isTRUE(r$numerical_ok)) {
      requireNamespace('glmmTMB',quietly=TRUE)
      d <- input[input$key==r$key,];if(r$component=='positive') d<-d[d$MC_reads>0,]
      co <- summary(cache$fit)$coefficients$cond['E_z',]
      require(same(unname(co[c(1,2,4)]),c(r$beta,r$SE,r$P)),'model coefficients')
      frame <- model.frame(cache$fit)
      require(same(model.response(frame),if(r$component=='positive')d$MC_proportion else as.integer(d$MC_reads>0)) &&
        same(frame$E_z,d$E_z) && identical(as.character(frame$TS),as.character(d$TS)),'model response/design')
      require(identical(as.character(cache$residuals$New_ID),as.character(d$New_ID)),'observation sample IDs')
      p <- cache$predictions;newdata<-d[rep(1,nrow(p)),];newdata$E_z<-p$E_z
      scale <- diff(range(d$E_raw))/diff(range(d$E_z));center<-d$E_raw[1]-d$E_z[1]*scale
      require(same(p$E_raw,p$E_z*scale+center),'prediction environment axis')
      q <- predict(cache$fit,newdata=newdata,re.form=NA,type='link',se.fit=TRUE)
      require(same(p$mean,plogis(q$fit)) && same(p$low,plogis(q$fit-qnorm(.975)*q$se.fit)) &&
        same(p$high,plogis(q$fit+qnorm(.975)*q$se.fit)),'model prediction values')
    }
  }
  list(tab=tab,pred=pred,obs=obs,scope=m$scope)
}

mc19_tp_seal <- function(out) {
  files <- c('00_result.csv',file.path('supporting_files',c('manifest.rds','analysis_input.csv','selected_original_results.csv','all_results.csv','predictions.rds','observations.rds')),
    list.files(file.path(out,'supporting_files','models'),pattern='[.]rds$',full.names=FALSE))
  model <- !grepl('[/\\]',files)&grepl('[.]rds$',files)
  files[model] <- file.path('supporting_files','models',files[model])
  write.csv(data.frame(file=files,md5=unname(tools::md5sum(file.path(out,files)))),file.path(out,'supporting_files','artifact_fingerprints.csv'),row.names=FALSE)
}

mc19_tp_plot <- function(out,project_dir=dirname(dirname(out)),parts=c('forest','curves','residual')) {
  stopifnot(requireNamespace('ggplot2',quietly=TRUE))
  stopifnot(length(parts)>0,all(parts %in% c('forest','curves','residual')))
  verified <- mc19_tp_validate(out)
  tab <- verified$tab;pred <- verified$pred;obs <- verified$obs;scope <- verified$scope
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
  grDevices::pdf(file.path(fig,'two_part_curves.pdf'),width=20,height=11.25,useDingbats=FALSE)
  ranges <- list()
  for(k in unique(tab$key)) {
    if(!any(pred$key==k)) next
    ranges[[k]] <- mc19_tp_curve_page(tab,pred,obs,k)
  }
  grDevices::dev.off()
  utils::write.csv(do.call(rbind,ranges),file.path(fig,'two_part_curve_display_ranges.csv'),row.names=FALSE)
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
    analysis_hash=mc19_tp_analysis_hash(),R=R.version.string,glmmTMB=as.character(utils::packageVersion('glmmTMB')))
  manifest <- file.path(support,'manifest.rds')
  if(file.exists(manifest)) {
    old <- readRDS(manifest)
    if(is.null(old$analysis_hash)) {
      saved <- file.path(support,'mc19_two_part_helpers.R')
      stopifnot(file.exists(saved),identical(unname(tools::md5sum(saved)),old$helper_hash))
      old$analysis_hash <- mc19_tp_analysis_hash(saved);old$helper_hash<-old$io_hash<-NULL
    }
    if(!setequal(names(old),names(signature)) || !identical(old[names(signature)],signature)) stop('Two-part analysis signature differs; use the matching analysis version or a new run.')
    saveRDS(signature,manifest)
  } else saveRDS(signature,manifest)
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
  mc19_tp_validate(out)
  mc19_tp_seal(out)
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
