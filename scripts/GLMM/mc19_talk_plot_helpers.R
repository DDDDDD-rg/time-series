# MC19 talk additions: saved-result plots and optional NMDS only.
mc19_talk_save <- function(p, stem, width=10, height=8) {
  ggplot2::ggsave(paste0(stem,'.pdf'),p,width=width,height=height,device=grDevices::cairo_pdf,bg='white')
  ggplot2::ggsave(paste0(stem,'.png'),p,width=width,height=height,dpi=180,bg='white')
}
mc19_talk_theme <- function() {
  ggplot2::theme_classic(base_size=21)+ggplot2::theme(
    text=ggplot2::element_text(face='bold'),
    axis.text=ggplot2::element_text(size=18,face='bold',colour='black'),
    axis.title=ggplot2::element_text(size=21,face='bold'),
    strip.text=ggplot2::element_text(size=20,face='bold'),
    legend.title=ggplot2::element_text(size=18,face='bold'),
    legend.text=ggplot2::element_text(size=16,face='bold'),
    plot.caption=ggplot2::element_text(size=15,face='bold',hjust=0,lineheight=1.15,margin=ggplot2::margin(t=14)),
    plot.caption.position='plot',plot.margin=ggplot2::margin(18,22,18,18))
}
mc19_p_mark <- function(p) {
  if(!is.finite(p)) 'NA' else if(p<.001) '***' else if(p<.01) '**' else if(p<.05) '*' else 'ns'
}
mc19_heatmap_forest <- function(results,cfg) {
  # A single row coordinate system guarantees alignment, including multiline rationales.
  par(mar=rep(0,4),family='sans',xpd=NA,font=2,font.axis=2,font.lab=2,lheight=.8)
  plot.new();plot.window(c(0,1),c(0,1),xaxs='i',yaxs='i')
  ids <- cfg$mc_ids; env <- cfg$environments; yy <- .925-(seq_along(ids)-.5)*.82/length(ids)
  display_ids <- sub('MC2s1','MC2S1',ids,fixed=TRUE)
  W <- grDevices::dev.size('in')[1]
  mx <- .20/W
  sx <- mx+max(strwidth(display_ids,units='inches',cex=.95)) / W+.25/W
  rationale <- mc19_selection_rationale()[ids]
  significant <- results$numerical_ok %in% TRUE & is.finite(results$beta_q_value) & results$beta_q_value<.05
  forest_env <- env[vapply(env,function(x) any(significant & results$environment==x),logical(1))]
  heatleft <- sx+max(strwidth(rationale,units='inches',cex=.85))/W+.20/W
  heatwidth <- 10/W
  forestleft <- heatleft+heatwidth+.18/W
  forestspan <- ((1-.25/W)-forestleft)/max(1L,length(forest_env))
  stopifnot(forestspan*W>2.5,max(strheight(rationale,units='inches',cex=.85))<.82*dev.size('in')[2]/length(ids))
  labels <- c('NO3-N','NO2-N','Ammonia-N','DIP')
  text(mx,.975,'MC',adj=0,font=2);text(sx,.975,'Selection rationale',adj=0,font=2)
  text(mx,yy,display_ids,adj=0,font=2,cex=.95)
  text(sx,yy,rationale,adj=0,cex=.85)
  lim <- max(abs(results$beta_per_TS_SD[results$numerical_ok %in% TRUE]),na.rm=TRUE)
  pal <- colorRampPalette(c('#083DCC','white','#D00000'))(201)
  visible <- results[results$environment %in% forest_env & results$numerical_ok %in% TRUE,]
  forestlim <- if(nrow(visible)) range(c(0,visible$CI95_low,visible$CI95_high),finite=TRUE) else c(-1,1)
  forestlim <- forestlim + c(-1,1)*diff(forestlim)*.04
  for(j in seq_along(env)) {
    rows <- results[results$environment==env[j],];rows <- rows[match(ids,rows$MC_ID),]
    stopifnot(!anyNA(rows$MC_ID))
    hx <- heatleft+(j-1)*heatwidth/4
    for(i in seq_along(ids)) {
      v <- rows$beta_per_TS_SD[i]; ok <- isTRUE(rows$numerical_ok[i]) && is.finite(v)
      rect(hx,yy[i]-.019,hx+heatwidth/4-.001,yy[i]+.019,col=if(ok) pal[1+round(200*(v+lim)/(2*lim))] else 'grey85',border='white')
      if(ok && is.finite(rows$beta_q_value[i]) && rows$beta_q_value[i]<.05) text(hx+heatwidth/8,yy[i],'*',cex=1.1,col=if(abs(v)>.65*lim) 'white' else 'black')
    }
    text(hx+heatwidth/8,.953,c('NO3','NO2','NH4','DIP')[j],cex=1.05)
    column <- match(env[j],forest_env)
    if(is.na(column)) next
    left <- forestleft+(column-1)*forestspan+.08/W; right <- left+forestspan-.35/W
    xx <- function(x) left+(x-forestlim[1])/diff(forestlim)*(right-left)
    segments(xx(0),.105,xx(0),.94,lty=3,col='grey55')
    good <- rows$numerical_ok %in% TRUE & is.finite(rows$CI95_low) & is.finite(rows$CI95_high)
    segments(xx(rows$CI95_low[good]),yy[good],xx(rows$CI95_high[good]),yy[good],lwd=1.1)
    points(xx(rows$beta_per_TS_SD[good]),yy[good],pch=ifelse(rows$beta_q_value[good]<.05,16,1),cex=.9)
    ticks <- pretty(forestlim,n=3);ticks <- ticks[ticks>=forestlim[1]&ticks<=forestlim[2]]
    segments(left,.095,right,.095);segments(xx(ticks),.095,xx(ticks),.09)
    text(xx(ticks),.077,format(ticks,trim=TRUE),cex=.95)
    text((left+right)/2,.975,labels[j],font=2,cex=1.1)
  }
  if(!length(forest_env)) text((forestleft+1)/2,.5,'No environmental panel has BH q < 0.05.',cex=1)
  text(forestleft+.30/W,.032,'Beta per TS nutrient SD; bars: pointwise 95% Wald CI.\nForest panels: nutrients with any MC q < 0.05. Filled points / *: saved BH q < 0.05 (76 tests).',adj=c(0,.5),cex=.8)
  # Numeric heatmap key; same range and colours as the effect cells.
  for(k in 1:201) rect(heatleft+(k-1)/201*heatwidth,.058,heatleft+k/201*heatwidth,.069,col=pal[k],border=NA)
  text(heatleft+c(0,.5,1)*heatwidth,.038,formatC(c(-lim,0,lim),digits=2,format='f'),cex=.75)
}
mc19_export_heatmap_forest <- function(stem,results,cfg) {
  mc19_main_export(stem,function() mc19_heatmap_forest(results,cfg),width=32,height=18,pointsize=24)
}
mc19_np_plot <- function(environment,out) {
  fields <- c('Nitrate_mgN_L','Nitrite_mgN_L','Ammonia_mgN_L','DIP_mgP_L')
  labs <- c('NO3-N (mg N/L)','NO2-N (mg N/L)','Ammonia nitrogen (mg N/L)','DIP (mg P/L)')
  e <- environment[order(as.integer(sub('TS','',environment$TS))),]
  e <- e[e$TS %in% paste0('TS',1:7),];stopifnot(nrow(e)==7,!anyDuplicated(e$TS))
  d <- do.call(rbind,lapply(seq_along(fields),function(j) {
    raw <- as.character(e[[paste0(fields[j],'_raw')]])
    data.frame(TS=e$TS,Date=as.Date(e$Water_quality_date),Indicator=labs[j],Value=e[[fields[j]]],Censored=grepl('<',raw,fixed=TRUE))
  }))
  d$Indicator <- factor(d$Indicator,levels=labs)
  lod <- data.frame(Indicator=factor(labs[c(2,4)],levels=labs),LOD=.002,
    Date=as.Date(e$Water_quality_date[1]),Label_y=c(.0033,.0045),Label='LOD = 0.002')
  p <- ggplot2::ggplot(d,ggplot2::aes(Date,Value))+
    ggplot2::geom_hline(data=lod,ggplot2::aes(yintercept=LOD),inherit.aes=FALSE,linetype=2,colour='grey35',linewidth=.75)+
    ggplot2::geom_line()+
    ggplot2::geom_point(ggplot2::aes(shape=Censored),size=3)+ggplot2::scale_shape_manual(values=c('FALSE'=16,'TRUE'=1),labels=c('Observed','Below LOD; L/2'))+
    ggplot2::geom_label(data=lod,ggplot2::aes(Date,Label_y,label=Label),inherit.aes=FALSE,
      hjust=0,size=5.5,fontface='bold',linewidth=0,fill='white',label.padding=grid::unit(.10,'lines'))+
    ggplot2::facet_wrap(~Indicator,ncol=2,scales='free_y')+
    ggplot2::scale_x_date(breaks=as.Date(e$Water_quality_date),labels=paste(e$TS,e$Sampling_month,sep='\n'))+
    ggplot2::labs(x=NULL,y=NULL,shape=NULL,caption='Dashed lines: NO2-N LOD = 0.002 mg N/L; DIP LOD = 0.002 mg P/L. Open circles: below LOD, shown at L/2.\nOne water-quality record per TS; TS8 has no contemporaneous environment record.')+
    mc19_talk_theme()+ggplot2::theme(axis.text.x=ggplot2::element_text(angle=45,hjust=1,size=17,face='bold'),legend.position='bottom',panel.spacing=grid::unit(1.2,'lines'))
  mc19_talk_save(p,file.path(out,'NP_by_TS'),20,11.25)
  write.csv(d,file.path(out,'NP_by_TS_data.csv'),row.names=FALSE)
}
mc19_talk_3d_widget <- function(d,results,scope,stress,initial='Nitrite_mgN_L') {
  fields <- c('Nitrite_mgN_L','Ammonia_mgN_L')
  stopifnot(initial %in% fields,all(is.finite(as.matrix(d[c('NMDS1','NMDS2','NMDS3')]))))
  short <- c(Nitrite_mgN_L='NO2-N',Ammonia_mgN_L='Ammonia-N')
  xyz <- as.matrix(d[c('NMDS1','NMDS2','NMDS3')])
  lower <- apply(xyz,2,min);upper <- apply(xyz,2,max)
  span <- pmax(upper-lower,.001)
  lower <- lower-.06*span;upper <- upper+.06*span
  span <- upper-lower
  colour_scale <- Map(function(x,c) list(x,c),seq(0,1,length.out=128),grDevices::hcl.colors(128,'Viridis'))
  caption <- function(field) {
    z <- results[results$Factor==field,];stopifnot(nrow(z)==1L)
    top <- sprintf('%s; n=%d; stress=%.4f; PERMANOVA R2=%.3f',scope,nrow(d),stress,z$R2_raw)
    bottom <- sprintf('raw P=%.4f %s; BH q=%.4f %s%s',z$P_raw,mc19_p_mark(z$P_raw),
      z$P_BH_raw,mc19_p_mark(z$P_BH_raw),if(scope=='coral' && field=='Ammonia_mgN_L') '; comparison panel' else '')
    paste0('<b>',top,'<br>',bottom,'</b>')
  }
  w <- plotly::plot_ly(height=850)
  for(field in fields) {
    value <- d[[field]];stopifnot(all(is.finite(value)))
    hover <- paste0('Sample: ',d$Sample_ID,'<br>TS: ',d$TS,'<br>',short[field],': ',value,' mg N/L')
    w <- plotly::add_trace(w,type='scatter3d',mode='markers',
      x=d$NMDS1,y=d$NMDS2,z=d$NMDS3,ids=as.character(d$New_ID),text=hover,hoverinfo='text',
      name=unname(short[field]),visible=field==initial,showlegend=FALSE,
      marker=list(size=3.3,opacity=.65,color=value,colorscale=colour_scale,
        cmin=min(value),cmax=max(value),showscale=TRUE,
        colorbar=list(title=list(text=paste0('<b>',short[field],'<br>(mg N/L)</b>'),font=list(size=20)),
          tickfont=list(size=18),x=1.035,y=.38,len=.60)),inherit=FALSE)
  }
  # Disable automatic scene axes: Plotly relocates their rendered edges with the
  # camera. Draw the cuboid, three axes and ticks as fixed world-space traces.
  edge <- function(a,b) rbind(a,b,rep(NA_real_,3))
  frame <- matrix(numeric(),ncol=3)
  for(axis in 1:3) for(a in 0:1) for(b in 0:1) {
    others <- setdiff(1:3,axis);start <- lower
    start[others] <- c(lower[others[1]],lower[others[2]])+c(a,b)*span[others]
    end <- start;end[axis] <- upper[axis]
    frame <- rbind(frame,edge(start,end))
  }
  w <- plotly::add_trace(w,type='scatter3d',mode='lines',x=frame[,1],y=frame[,2],z=frame[,3],
    line=list(color='#B5BCC4',width=2),name='Fixed bounding box',hoverinfo='skip',showlegend=FALSE,inherit=FALSE)
  colours <- c('#0072B2','#D55E00','#009E73')
  anchor <- c(upper[1],upper[2],lower[3])
  for(axis in 1:3) {
    end <- anchor;end[axis] <- if(axis<3L) lower[axis] else upper[axis]
    ticks <- pretty(c(lower[axis],upper[axis]),n=3)
    ticks <- ticks[ticks>lower[axis]+.07*span[axis] & ticks<upper[axis]-.01*span[axis]]
    if(length(ticks)<2L) {
      ticks <- pretty(c(lower[axis],upper[axis]),n=5)
      ticks <- ticks[ticks>lower[axis]+.07*span[axis] & ticks<upper[axis]-.01*span[axis]]
    }
    if(length(ticks)>3L) ticks <- ticks[unique(round(seq(1,length(ticks),length.out=3)))]
    line <- edge(anchor,end);label_pos <- matrix(numeric(),ncol=3)
    offset <- if(axis==1L) 2L else 1L
    for(tick in ticks) {
      start <- anchor;start[axis] <- tick;end_tick <- start
      end_tick[offset] <- end_tick[offset]+.025*span[offset]
      line <- rbind(line,edge(start,end_tick))
      label <- start;label[offset] <- label[offset]+.075*span[offset]
      label_pos <- rbind(label_pos,label)
    }
    w <- plotly::add_trace(w,type='scatter3d',mode='lines',x=line[,1],y=line[,2],z=line[,3],
      line=list(color=colours[axis],width=6),name=paste0('Fixed axis NMDS',axis),
      hoverinfo='skip',showlegend=FALSE,inherit=FALSE)
    w <- plotly::add_trace(w,type='scatter3d',mode='text',x=label_pos[,1],y=label_pos[,2],z=label_pos[,3],
      text=format(ticks,trim=TRUE,digits=3),textposition='middle center',
      textfont=list(family='Arial Black, Arial, sans-serif',size=17,color=colours[axis]),
      name=paste0('Fixed ticks NMDS',axis),hoverinfo='skip',showlegend=FALSE,inherit=FALSE)
  }
  axis_style <- function(axis) list(visible=FALSE,showgrid=FALSE,showbackground=FALSE,
    showline=FALSE,showticklabels=FALSE,zeroline=FALSE,title=list(text=''),
    range=c(lower[axis]-.16*span[axis],upper[axis]+.16*span[axis]),autorange=FALSE)
  fixed_axes <- paste0('<b><span style="color:#0072B2">X = NMDS1</span> &nbsp; | &nbsp; ',
    '<span style="color:#D55E00">Y = NMDS2</span> &nbsp; | &nbsp; ',
    '<span style="color:#009E73">Z = NMDS3</span></b>')
  # Explicit trace lists keep the fixed frame/axes visible during field switches.
  n_traces <- 9L
  buttons <- lapply(seq_along(fields),function(i) {
    visible <- c(seq_along(fields)==i,rep(TRUE,n_traces-2L))
    layout_change <- setNames(list(caption(fields[i])),'annotations[0].text')
    list(label=unname(short[fields[i]]),method='update',args=list(list(visible=visible),layout_change))
  })
  plotly::layout(w,font=list(family='Arial Black, Arial, sans-serif',size=20),
    uirevision='MC19-fixed-coordinates',margin=list(l=25,r=225,b=85,t=130),
    scene=list(aspectmode='data',uirevision='MC19-fixed-camera',
      camera=list(eye=list(x=1.30,y=1.30,z=1.00),center=list(x=0,y=0,z=0)),
      xaxis=axis_style(1),yaxis=axis_style(2),zaxis=axis_style(3)),
    updatemenus=list(list(type='buttons',direction='down',active=match(initial,fields)-1L,
      x=1.025,y=.98,xanchor='left',yanchor='top',font=list(size=19),buttons=buttons)),
    annotations=list(
      list(text=caption(initial),x=0,y=1.08,xref='paper',yref='paper',xanchor='left',yanchor='bottom',
        showarrow=FALSE,font=list(size=21),align='left'),
      list(text=fixed_axes,x=.5,y=-.07,xref='paper',yref='paper',xanchor='center',yanchor='top',
        showarrow=FALSE,font=list(size=20)),
      list(text='<b>Color by</b>',x=1.025,y=1.04,xref='paper',yref='paper',xanchor='left',yanchor='bottom',
        showarrow=FALSE,font=list(size=20))))
}
mc19_save_widget <- function(widget,path) {
  dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
  path <- file.path(normalizePath(dirname(path),winslash='/',mustWork=TRUE),basename(path))
  stage <- tempfile('mc19-html-');dir.create(stage)
  on.exit({
    resolved <- normalizePath(stage,winslash='/',mustWork=FALSE)
    if(startsWith(resolved,paste0(normalizePath(tempdir(),winslash='/'),'/'))) unlink(stage,recursive=TRUE)
  },add=TRUE)
  target <- file.path(stage,basename(path))
  withCallingHandlers(htmlwidgets::saveWidget(widget,target,selfcontained=FALSE,libdir='lib'),
    warning=function(w) stop('HTML resource export failed: ',conditionMessage(w),call.=FALSE))
  assets <- list.files(file.path(stage,'lib'),recursive=TRUE,full.names=FALSE)
  destinations <- file.path(dirname(path),'lib',assets)
  if(any(nchar(c(path,destinations),type='chars')>245L))
    stop('HTML output path is too deep for reliable Windows export; choose a shorter output directory.')
  for(i in seq_along(assets)) {
    dir.create(dirname(destinations[i]),recursive=TRUE,showWarnings=FALSE)
    if(!file.copy(file.path(stage,'lib',assets[i]),destinations[i],overwrite=TRUE) ||
       unname(tools::md5sum(file.path(stage,'lib',assets[i])))!=unname(tools::md5sum(destinations[i]))) stop('HTML dependency copy/verification failed.')
  }
  html <- paste(readLines(target,warn=FALSE),collapse='\n')
  refs <- regmatches(html,gregexpr('(?:src|href)="[^"]+"',html,perl=TRUE))[[1]]
  refs <- sub('^[^=]+="','',sub('"$','',refs))
  local <- refs[!grepl('^([a-zA-Z]+:|//|#)',refs)]
  if(!all(file.exists(file.path(dirname(path),local)))) stop('HTML contains missing local dependencies.')
  if(!file.copy(target,path,overwrite=TRUE)) stop('Cannot publish verified HTML.')
  invisible(path)
}
mc19_talk_nmds <- function(context,out,dimensions=c(2L,3L),trymax=500L,cache_only=FALSE) {
  stopifnot(requireNamespace('vegan',quietly=TRUE))
  ids <- as.character(context$input$metadata$New_ID)
  results <- context$saved$results
  fields <- c('Nitrite_mgN_L','Ammonia_mgN_L')
  diag <- list()
  for(k in dimensions) {
    cachepath <- file.path(out,paste0('NMDS',k,'D.rds'))
    oldpath <- mc19_support_file(context$result_dir,paste0('plot_NMDS',k,'D.rds'))
    prior <- if(file.exists(cachepath)) readRDS(cachepath) else if(file.exists(oldpath)) readRDS(oldpath) else NULL
    if(!is.null(prior)) stopifnot(identical(prior$input_md5,context$input$input_md5),identical(prior$New_ID,ids),identical(prior$MC,mc_cols),identical(prior$k,as.integer(k)))
    if(cache_only && is.null(prior)) stop('Saved NMDS cache required for style-only rendering: ',cachepath)
    reused <- !is.null(prior) && (cache_only || !isTRUE(getOption('mc19.nmds.refit_cached',FALSE)) || isTRUE(prior$fit$converged>0))
    if(reused) fit <- prior$fit else {
      set.seed(20261031L+k)
      fit <- vegan::metaMDS(context$saved$distance,k=k,trymax=trymax,autotransform=FALSE,wascores=FALSE,trace=FALSE,
        previous.best=if(!is.null(prior)) prior$fit else NULL)
    }
    xy <- as.matrix(vegan::scores(fit,display='sites',choices=seq_len(k)))
    stopifnot(setequal(rownames(xy),ids),all(is.finite(xy)),is.finite(fit$stress))
    xy <- xy[match(ids,rownames(xy)),,drop=FALSE];colnames(xy) <- paste0('NMDS',seq_len(k))
    if(!cache_only) saveRDS(list(fit=fit,coordinates=xy,k=as.integer(k),New_ID=ids,MC=mc_cols,input_md5=context$input$input_md5),cachepath)
    d <- cbind(context$input$metadata,as.data.frame(xy))
    if(!cache_only) write.csv(d,file.path(out,paste0('NMDS',k,'D_coordinates.csv')),row.names=FALSE)
    diag[[as.character(k)]] <- data.frame(k=k,stress=fit$stress,tries=fit$tries,best_repeats=fit$converged,reused=reused)
    grDevices::pdf(file.path(out,paste0('NMDS',k,'D_Shepard.pdf')),width=12.8,height=7.2,pointsize=20)
    graphics::par(font=2,font.lab=2,font.axis=2,font.main=2,mar=c(5,5,3,2));vegan::stressplot(fit);grDevices::dev.off()
    if(!isTRUE(fit$converged>0)) warning('NMDS best solution not repeated: inspect diagnostics; k=',k)
    for(field in fields) {
      z <- results[results$Factor==field,];stopifnot(nrow(z)==1)
      d$Env <- d[[field]]
      caption <- sprintf('%s; n=%d; stress=%.4f; PERMANOVA R2=%.3f; raw P=%.4f %s; BH q=%.4f %s%s',
        analysis_scope,nrow(d),fit$stress,z$R2_raw,z$P_raw,mc19_p_mark(z$P_raw),z$P_BH_raw,mc19_p_mark(z$P_BH_raw),
        if(analysis_scope=='coral' && field=='Ammonia_mgN_L') '; comparison panel' else '')
      if(k==2L) {
        p <- ggplot2::ggplot(d,ggplot2::aes(NMDS1,NMDS2,colour=Env))+ggplot2::geom_point(size=2,alpha=.8)+
          ggplot2::scale_colour_viridis_c(limits=range(context$input$metadata[[field]]))+
          ggplot2::coord_equal()+ggplot2::labs(colour=unname(context$labels[field]),caption=paste(strwrap(caption,85),collapse='\n'))+
          mc19_talk_theme()
        mc19_talk_save(p,file.path(out,paste0('NMDS2D_',field)),16,9)
      } else {
        stopifnot(requireNamespace('plotly',quietly=TRUE),requireNamespace('htmlwidgets',quietly=TRUE))
        w <- mc19_talk_3d_widget(d,results,analysis_scope,fit$stress,field)
        mc19_save_widget(w,file.path(out,paste0('NMDS3D_',field,'.html')))
        for(pair in list(c(1,2),c(1,3),c(2,3))) {
          p <- ggplot2::ggplot(d,ggplot2::aes(x=.data[[paste0('NMDS',pair[1])]],y=.data[[paste0('NMDS',pair[2])]],colour=Env))+
            ggplot2::geom_point(size=2,alpha=.8)+ggplot2::scale_colour_viridis_c()+ggplot2::coord_equal()+
            ggplot2::labs(colour=unname(context$labels[field]),caption=paste(strwrap(paste('Projection of 3D NMDS;',caption),85),collapse='\n'))+mc19_talk_theme()
          mc19_talk_save(p,file.path(out,paste0('NMDS3D_',field,'_axes',paste(pair,collapse=''))),16,9)
        }
      }
    }
    if(k==3L) {
      combined <- mc19_talk_3d_widget(d,results,analysis_scope,fit$stress)
      mc19_save_widget(combined,file.path(out,'NMDS3D_nutrients.html'))
    }
  }
  if(!cache_only) write.csv(do.call(rbind,diag),file.path(out,'NMDS_diagnostics.csv'),row.names=FALSE)
}
