# Independent 2D NMDS reproducibility check. Never overwrites plots/caches/PPT.
(function(){
 library(vegan)
 project<-normalizePath('Research/Msc main project-time series',winslash='/')
 out<-file.path(project,'results/nmds_2d_stability_v1');dir.create(out,recursive=TRUE,showWarnings=FALSE)
 sources<-c(file.path(project,'presentations/time_series.pptx'),unlist(lapply(c('coral','ecosystem'),function(s)c(file.path(project,'results',paste0('MC19_PERMANOVA_TS7_env6_',s,'_run001/supporting_files/'),paste0(s,'_environment_analysis.rds')),file.path(project,'figures',if(s=='coral')'overview_coral_v1' else 'overview_all_hosts_v1',paste0('NMDS',2:3,'D.rds'))))))
 hashes<-tools::md5sum(sources)
 inputs<-lapply(c('coral','ecosystem'),function(s){
  a<-readRDS(file.path(project,'results',paste0('MC19_PERMANOVA_TS7_env6_',s,'_run001/supporting_files/'),paste0(s,'_environment_analysis.rds')))
  old<-readRDS(file.path(project,'figures',if(s=='coral')'overview_coral_v1' else 'overview_all_hosts_v1','NMDS2D.rds'))
  stopifnot(identical(attr(a$distance,'Labels'),old$New_ID),max(abs(as.numeric(a$distance)-as.numeric(vegdist(sqrt(a$input$community),method='bray'))))<1e-7)
  list(scope=s,d=a$distance,old=old)
 });names(inputs)<-c('coral','ecosystem')
 jobs<-expand.grid(scope=names(inputs),mode=c('fresh','warm'),seed=c(101L,203L,307L,409L,503L),stringsAsFactors=FALSE)
 worker<-function(j,inputs,out){
  library(vegan);x<-inputs[[j$scope]];key<-paste(j$scope,j$mode,j$seed,sep='_');path<-file.path(out,paste0(key,'.rds'))
  if(file.exists(path))return(readRDS(path))
  set.seed(j$seed);t<-proc.time()[3]
  args<-list(comm=x$d,k=2,try=50,trymax=50,maxit=1000,autotransform=FALSE,wascores=FALSE,trace=FALSE,parallel=1)
  if(j$mode=='warm')args$previous.best<-x$old$fit
  fit<-do.call(metaMDS,args)
  z<-list(key=key,scope=j$scope,mode=j$mode,seed=j$seed,fit=fit,elapsed=unname(proc.time()[3]-t));saveRDS(z,path);z
 }
 cl<-parallel::makePSOCKcluster(min(4L,parallel::detectCores()));on.exit(parallel::stopCluster(cl),add=TRUE)
 message('Running 20 batches, 50 starts each, 4 workers; maxit=1000. Output: ',out)
 fits<-parallel::parLapplyLB(cl,split(jobs,seq_len(nrow(jobs))),worker,inputs=inputs,out=out)
 rows<-list();pairrows<-list();pointrows<-list()
 pdf(file.path(out,'diagnostics.pdf'),width=12,height=6)
 for(s in names(inputs)){
  x<-inputs[[s]];fs<-Filter(function(z)z$scope==s,fits);best<-fs[[which.min(vapply(fs,function(z)z$fit$stress,numeric(1)))]]
  ref<-scores(best$fit,display='sites');old<-scores(x$old$fit,display='sites');saveRDS(best,file.path(out,paste0(s,'_candidate_2D.rds')))
  all<-c(list(list(key='original',mode='original',seed=NA,fit=x$old$fit,elapsed=NA)),fs)
  for(z in all){
   xy<-scores(z$fit,display='sites');p<-procrustes(ref,xy,symmetric=TRUE);r<-residuals(p)
   rows[[length(rows)+1]]<-data.frame(scope=s,key=z$key,mode=z$mode,seed=z$seed,stress=z$fit$stress,best_repeats=z$fit$converged,tries=z$fit$tries,iterations=z$fit$iters,stop=z$fit$icause,proc_r=sqrt(max(0,1-p$ss)),rmse=sqrt(mean(r^2)),max_residual=max(r),elapsed=z$elapsed)
   pointrows[[length(pointrows)+1]]<-data.frame(scope=s,key=z$key,ID=rownames(ref),residual=r)
  }
  for(i in seq_along(fs))for(j in seq_len(i-1L)){
   p<-procrustes(scores(fs[[i]]$fit,display='sites'),scores(fs[[j]]$fit,display='sites'),symmetric=TRUE)
   pairrows[[length(pairrows)+1]]<-data.frame(scope=s,a=fs[[i]]$key,b=fs[[j]]$key,proc_r=sqrt(max(0,1-p$ss)),rmse=sqrt(max(0,p$ss)/nrow(ref)),max_residual=max(residuals(p)))
  }
  par(mfrow=c(1,2));stressplot(x$old$fit,x$d,main=paste(s,'original 2D'));stressplot(best$fit,x$d,main=paste(s,'candidate 2D'))
  par(mfrow=c(1,2));plot(procrustes(ref,old,symmetric=TRUE),main=paste(s,'original vs candidate'));plot(procrustes(ref,ref,symmetric=TRUE)$X,type='n',main=paste(s,'aligned solutions'));for(z in fs){p<-procrustes(ref,scores(z$fit,display='sites'),symmetric=TRUE);points(p$Yrot,pch=16,cex=.25,col=adjustcolor('blue',.2))};points(procrustes(ref,ref,symmetric=TRUE)$X,pch=1,cex=.5)
 }
 dev.off()
 write.csv(do.call(rbind,rows),file.path(out,'run_summary.csv'),row.names=FALSE)
 write.csv(do.call(rbind,pairrows),file.path(out,'pairwise_alignment.csv'),row.names=FALSE)
 write.csv(do.call(rbind,pointrows),file.path(out,'point_residuals.csv'),row.names=FALSE)
 for(s in names(inputs)){
  d<-subset(do.call(rbind,pointrows),scope==s & grepl('fresh',key));d$relative_rms<-d$residual*sqrt(attr(inputs[[s]]$d,'Size'))
  m<-aggregate(relative_rms~ID,d,max);a<-readRDS(file.path(project,'results',paste0('MC19_PERMANOVA_TS7_env6_',s,'_run001/supporting_files/'),paste0(s,'_environment_analysis.rds')))
  z<-merge(m,a$input$metadata,by.x='ID',by.y='New_ID');z<-z[order(-z$relative_rms),];write.csv(z,file.path(out,paste0(s,'_point_instability.csv')),row.names=FALSE)
 }
 stopifnot(identical(hashes,tools::md5sum(sources)))
 writeLines(c('Completed 20 independent seed/mode batches; 1000 random starts total.','Original sources and PPT MD5 unchanged.',capture.output(sessionInfo())),file.path(out,'COMPLETE.txt'))
 print(do.call(rbind,rows));message('COMPLETE')
})()
