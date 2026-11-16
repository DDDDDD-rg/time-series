# Background worker: one scope per process; results/diagnostics only, no figures.
(function() {
  args <- commandArgs(trailingOnly=TRUE)
  stopifnot(length(args)==4L,args[1] %in% c("coral","ecosystem"))
  scope <- args[1];threads <- as.integer(args[2]);project <- args[3];logs <- args[4]
  stopifnot(is.finite(threads),threads>=1L,dir.exists(project),dir.exists(logs))
  status <- file.path(logs,paste0(scope,"_status.txt"))
  stamp <- function(text) writeLines(paste(format(Sys.time(),"%Y-%m-%d %H:%M:%S %Z"),text),status)
  options(mc19.project_dir=project,mc19.analysis_scope=scope,mc19.temperature.stage="results",
    mc19.glmm.result_dir="",mc19.result_dir="",mc19.two_part.run_dir="",
    mc19.talk.output_dir="",mc19.figure_dir=NULL,mc19.run_id=NULL,
    mc19.run_plots=FALSE,mc19.run_diagnostics=TRUE,glmmTMB.cores=threads)
  stamp("RUNNING: GLMM + diagnostics -> PERMANOVA -> selected two-part; no figures")
  cat("BACKGROUND START:",scope,"PID",Sys.getpid(),"model threads",threads,"\n");flush.console()
  tryCatch({
    source(file.path(project,"scripts","TS7_temperature.R"),encoding="UTF-8",local=new.env(parent=globalenv()))
    stamp("COMPLETE: all requested result stages returned")
    cat("BACKGROUND COMPLETE:",scope,"\n")
  },error=function(e) {
    stamp(paste("FAILED:",conditionMessage(e)))
    message("BACKGROUND FAILED: ",scope,": ",conditionMessage(e))
    quit(save="no",status=1L)
  },interrupt=function(e) {
    stamp("INTERRUPTED; partial outputs retained")
    quit(save="no",status=130L)
  })
})()
