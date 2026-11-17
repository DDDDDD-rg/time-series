# Read-only reuse audits and user-run incremental writers. No execution on source().
mc19_incremental_same <- function(a,b,tolerance=1e-12) {
  length(a)==length(b) && all((is.na(a)&is.na(b)) | (!is.na(a)&!is.na(b)&abs(a-b)<=tolerance))
}
mc19_incremental_assert <- function(ok, message) {
  if(!isTRUE(ok)) stop("Reuse audit failed: ",message,call.=FALSE)
}
mc19_incremental_equal_data <- function(a,b,fields) {
  mc19_incremental_assert(nrow(a)==nrow(b) && all(fields %in% names(a)) && all(fields %in% names(b)),"data shape/fields")
  for(f in fields) {
    if(f=="Water_quality_date") { a[[f]] <- as.character(as.Date(a[[f]]));b[[f]] <- as.character(as.Date(b[[f]])) }
    ok <- if(is.numeric(a[[f]]) && is.numeric(b[[f]])) mc19_incremental_same(a[[f]],b[[f]],1e-10) else
      identical(as.character(a[[f]]),as.character(b[[f]]))
    mc19_incremental_assert(ok,paste("sample-level",f))
  }
  invisible(TRUE)
}
mc19_incremental_permanova_plan <- function(project,scope,input) {
  source <- getOption("mc19.temperature.legacy_permanova",file.path(project,"results",
    if(scope=="coral")"permanova_coral_v1" else "permanova_all_hosts_v1"))
  file <- file.path(source,"supporting_files",paste0(scope,"_environment_analysis.rds"))
  old <- readRDS(file)
  fields <- mc19_profile()$fields
  mc19_incremental_assert(length(fields)==5L && !"DIN_mgN_L" %in% fields,"PERMANOVA five fields")
  mc19_incremental_assert(all(old$results$Scope==scope) && all(old$results$Status=="ok"),"PERMANOVA scope/status")
  mc19_incremental_equal_data(old$input$metadata,input$metadata,c("New_ID","Sample_ID","Host_type","Host_niche","TS",fields))
  mc19_incremental_assert(identical(colnames(old$input$community),colnames(input$community)) &&
    identical(dim(old$input$community),dim(input$community)) &&
    mc19_incremental_same(as.numeric(old$input$community),as.numeric(input$community)),"PERMANOVA composition")
  mc19_incremental_equal_data(old$input$batches,input$batches,c("TS","Sampling_month","Water_quality_date",fields))
  mc19_incremental_assert(identical(old$input$row_batch,input$row_batch) &&
    identical(as.character(attr(old$distance,"Labels")),as.character(input$metadata$New_ID)) &&
    identical(attr(old$distance,"method"),"bray"),"PERMANOVA distance and TS mapping")
  idx <- match(fields,old$results$Factor)
  mc19_incremental_assert(!anyNA(idx) && !anyDuplicated(old$results$Factor) &&
    nrow(old$TS_orders)==5040L && ncol(old$TS_orders)==7L &&
    identical(as.integer(old$TS_orders[1,]),1:7) && nrow(old$F_per_order)==5040L &&
    all(fields %in% colnames(old$F_per_order)),"saved exhaustive TS permutations")
  for(f in fields) {
    row <- old$results[old$results$Factor==f,];values <- old$F_per_order[,f]
    tol <- 1e-10*max(1,abs(values[1]))
    mc19_incremental_assert(nrow(row)==1L && row$N_TS==7L && row$N_samples==nrow(input$metadata) &&
      mc19_incremental_same(values[1],row$F_raw,1e-10) &&
      mc19_incremental_same(sum(values>=values[1]-tol)/length(values),row$P_raw,1e-12),paste("saved F/raw P",f))
  }
  list(source=source,file=file,md5=unname(tools::md5sum(file)),old=old,indices=idx)
}
mc19_incremental_permanova_run <- function(project,scope,input,input_file) {
  plan <- mc19_incremental_permanova_plan(project,scope,input)
  old <- plan$old;fields <- mc19_profile()$fields;version <- mc19_method_version(scope,"PERMANOVA")
  run <- mc19_new_run(project,"PERMANOVA",scope,version,input_file,getOption("mc19.run_id",NULL))
  support <- file.path(run,"supporting_files");dir.create(support)
  writeLines("running",file.path(support,"00_run_status.txt"))
  tab <- old$results[plan$indices,,drop=FALSE]
  tab$Original_P_BH_raw <- tab$P_BH_raw
  tab$P_BH_raw <- p.adjust(tab$P_raw,"BH",n=5L)
  tab$Significance <- ifelse(tab$P_BH_raw<=.001,"***",ifelse(tab$P_BH_raw<=.01,"**",ifelse(tab$P_BH_raw<=.05,"*","ns")))
  tab$Source_version <- old$version;tab$Analysis_version <- version;tab$Input_MD5 <- input$input_md5
  tab$Reused_from <- normalizePath(plan$file,winslash="/",mustWork=TRUE)
  tab$Adjustment <- "BH across five selected single-environment tests; raw P/F/R2 reused"
  main <- c("Scope","Indicator","Factor","N_samples","N_MC","N_TS","R2_raw","F_raw","P_raw","P_BH_raw","Significance")
  write_csv_excel(tab[main],file.path(run,"environment_results.csv"))
  write_csv_excel(tab,file.path(support,"environment_results_full.csv"))
  write_csv_excel(input$batches,file.path(support,"environment_TS_values.csv"))
  write_csv_excel(cbind(input$metadata,as.data.frame(input$community)),file.path(support,"environment_analysis_input.csv"))
  saveRDS(list(results=tab,input=input,distance=old$distance,TS_orders=old$TS_orders,
    F_per_order=old$F_per_order[,fields,drop=FALSE],version=version,scope=scope),
    file.path(support,paste0(scope,"_environment_analysis.rds")))
  file.copy(input_file,file.path(support,"00_input_used_MC19_TS8_env5.xlsx"))
  file.copy(mc19_temperature_source(project),file.path(support,"PM11_temperature_source.csv"))
  write_csv_excel(input$temperature_source,file.path(support,"temperature_source_matching.csv"))
  write_csv_excel(data.frame(path=plan$file,md5=plan$md5),file.path(support,"reuse_fingerprints.csv"))
  # Coordinates depend on the unchanged community, not the selected colour variables.
  for(k in 2:3) {
    candidates <- c(file.path(plan$source,"supporting_files",paste0("plot_NMDS",k,"D.rds")),
      list.files(file.path(project,"figures"),pattern=paste0("^NMDS",k,"D[.]rds$"),recursive=TRUE,full.names=TRUE))
    for(path in candidates[file.exists(candidates)]) {
      cache <- tryCatch(readRDS(path),error=function(e)NULL)
      if(is.null(cache) || !identical(cache$input_md5,old$input$input_md5) ||
        !identical(as.character(cache$New_ID),as.character(input$metadata$New_ID)) ||
        !identical(cache$MC,colnames(input$community)) || !identical(as.integer(cache$k),as.integer(k))) next
      cache$reused_from <- normalizePath(path,winslash="/",mustWork=TRUE)
      cache$version <- version;cache$input_md5 <- input$input_md5
      saveRDS(cache,file.path(support,paste0("plot_NMDS",k,"D.rds")));break
    }
  }
  mc19_incremental_assert(identical(unname(tools::md5sum(plan$file)),plan$md5),"PERMANOVA source changed")
  writeLines(c("Five tests reused; no adonis2, statistical permutations or NMDS executed.",
    "Four nutrients plus temperature; DIN, salinity and pH excluded; BH recalculated from five raw P values."),file.path(support,"00_summary.txt"))
  writeLines(c("complete",version),file.path(support,"00_run_status.txt"))
  mc19_publish_result(project,"PERMANOVA",scope,run)
  invisible(run)
}
mc19_incremental_bind <- function(rows) {
  fields <- unique(unlist(lapply(rows,names)))
  do.call(rbind,lapply(rows,function(r) { for(f in setdiff(fields,names(r))) r[[f]] <- NA; r[fields] }))
}
mc19_incremental_entry <- function(entry,d) {
  mc19_incremental_equal_data(entry$data,d,c("New_ID","Sample_ID","Host_type","Host_niche","TS",
    "MC_ID","environment","MC_reads","Other_reads","Ruegeria_reads_total","MC_proportion","E_raw","E_z"))
  mc19_incremental_assert(identical(entry$model,"BB_TS"),"model must have TS random intercept")
  if(!is.null(entry$fit)) {
    fit <- entry$fit
    mc19_incremental_assert(family(fit)$family=="betabinomial" && family(fit)$link=="logit","family/link")
    form <- paste(deparse(formula(fit)),collapse=" ")
    mc19_incremental_assert(gsub("\\s+","",form)=="cbind(MC_reads,Other_reads)~E_z+(1|TS)","single-indicator formula")
    frame <- model.frame(fit)
    mc19_incremental_assert(mc19_same(as.numeric(model.response(frame)),as.numeric(cbind(d$MC_reads,d$Other_reads))) &&
      mc19_same(frame$E_z,d$E_z) && identical(as.character(frame$TS),as.character(d$TS)),"saved model design/response")
    if(isTRUE(entry$result$numerical_ok)) {
      co <- summary(fit)$coefficients$cond["E_z",]
      mc19_incremental_assert(mc19_same(unname(co[c(1,2,4)]),
        c(entry$result$beta_per_TS_SD,entry$result$SE,entry$result$beta_P_value),1e-10),"saved beta/SE/raw P")
    }
  }
  invisible(TRUE)
}
mc19_incremental_glmm_source <- function(cfg,input,source,fields) {
  support <- file.path(source,"supporting_files")
  m <- readRDS(file.path(support,"wald_manifest.rds"))
  tab <- read.csv(file.path(support,"all_results.csv"),stringsAsFactors=FALSE)
  mc19_incremental_assert(isTRUE(m$fit_complete),"source GLMM incomplete")
  scope <- if(is.null(m$analysis_scope)) "coral" else m$analysis_scope
  mc19_incremental_assert(identical(scope,cfg$analysis_scope),"source scope")
  entries <- hashes <- list()
  for(f in fields) for(mc in cfg$mc_ids) {
    key <- mc19_wald_key(mc,f,"BB_TS")
    path <- unname(m$model_files[key])
    mc19_incremental_assert(length(path)==1L && !is.na(path) && file.exists(path),paste("missing model",key))
    hash <- mc19_hash(path)
    old <- m$model_hashes$md5[match(hash$path,normalizePath(m$model_hashes$path,winslash="/",mustWork=TRUE))]
    mc19_incremental_assert(length(old)==1L && !is.na(old) && identical(hash$md5,old),paste("model hash",key))
    entry <- readRDS(path);d <- input$blocks[[paste(mc,f,sep="__")]]
    mc19_incremental_entry(entry,d)
    row <- tab[tab$MC_ID==mc & tab$environment==f & tab$model=="BB_TS",]
    mc19_incremental_assert(nrow(row)==1L && mc19_same(row$beta_P_value,entry$result$beta_P_value) &&
      mc19_same(row$beta_per_TS_SD,entry$result$beta_per_TS_SD),paste("source table",key))
    entry$source_model_file <- normalizePath(path,winslash="/",mustWork=TRUE)
    entry$source_diagnostic <- file.path(support,"diagnostics",paste0(key,".rds"))
    entries[[key]] <- entry;hashes[[key]] <- hash
  }
  list(entries=entries,hashes=do.call(rbind,hashes),source=source,version=m$version)
}
mc19_incremental_glmm_plan <- function(cfg,input=mc19_wald_load_input(cfg)) {
  fields <- setdiff(cfg$environments,"SST_C")
  mc19_incremental_assert(length(fields)==4L && !"DIN_mgN_L" %in% cfg$environments,"four nutrients plus temperature")
  legacy <- getOption("mc19.temperature.legacy_glmm",file.path(cfg$project_dir,"results",
    if(cfg$analysis_scope=="coral") "glmm_coral_v1" else "glmm_all_hosts_v1"))
  old <- mc19_incremental_glmm_source(cfg,input,legacy,fields)
  mc19_incremental_assert(length(old$entries)==76L,"76 legacy nutrient models")
  temperature <- NULL;notices <- character()
  # Reuse already saved temperature fits too, if the full source audit passes.
  if(isTRUE(getOption("mc19.temperature.reuse_saved",TRUE))) {
    candidates <- list.dirs(file.path(cfg$project_dir,"results"),recursive=FALSE,full.names=TRUE)
    candidates <- candidates[grepl(paste0("^MC19_GLMM.*_",cfg$analysis_scope,"_run[0-9]+$"),basename(candidates))]
    candidates <- candidates[order(file.info(candidates)$mtime,decreasing=TRUE)]
    for(path in candidates) {
      attempt <- tryCatch({
        m <- readRDS(file.path(path,"supporting_files","wald_manifest.rds"))
        allowed <- c(paste0("MC19_env6_TS7_temperature_",cfg$analysis_scope,"_Wald_v1"),cfg$version)
        if(!m$version %in% allowed) stop("other workflow")
        mc19_incremental_glmm_source(cfg,input,path,"SST_C")
      },error=function(e) { notices <<- c(notices,paste(basename(path),conditionMessage(e)));NULL })
      if(!is.null(attempt)) { temperature <- attempt;break }
    }
  }
  list(input=input,nutrients=old,temperature=temperature,
    new_temperature_models=if(is.null(temperature))19L else 0L,notices=notices)
}
mc19_incremental_copy_diagnostic <- function(entry,new_model,support) {
  path <- entry$source_diagnostic
  if(is.null(path) || !file.exists(path)) return(invisible(FALSE))
  diagnosis <- tryCatch(readRDS(path),error=function(e)NULL)
  if(is.null(diagnosis) || is.null(diagnosis$model_md5) ||
      !identical(diagnosis$model_md5,unname(tools::md5sum(entry$source_model_file)))) return(invisible(FALSE))
  if(!all(c("summaries","residual_table","tests","TS_residuals","n_sim","seed") %in% names(diagnosis))) return(invisible(FALSE))
  diagnosis$reused_from <- normalizePath(path,winslash="/",mustWork=TRUE)
  diagnosis$model_md5 <- unname(tools::md5sum(new_model))
  dir.create(file.path(support,"diagnostics"),showWarnings=FALSE)
  saveRDS(diagnosis,file.path(support,"diagnostics",basename(path)))
  invisible(TRUE)
}
mc19_incremental_glmm_run <- function(cfg) {
  plan <- mc19_incremental_glmm_plan(cfg);input <- plan$input
  scripts <- file.path(cfg$project_dir,"scripts","GLMM",cfg$scripts)
  script_hashes <- mc19_hash(scripts)
  reused <- c(plan$nutrients$entries,if(!is.null(plan$temperature))plan$temperature$entries)
  source_hashes <- rbind(plan$nutrients$hashes,if(!is.null(plan$temperature))plan$temperature$hashes)
  run <- mc19_new_run(cfg$project_dir,"GLMM",cfg$analysis_scope,cfg$version,cfg$workbook,getOption("mc19.run_id",NULL))
  support <- file.path(run,"supporting_files");dir.create(file.path(support,"models"),recursive=TRUE)
  dir.create(file.path(support,"scripts_used"))
  on.exit(mc19_write_text(capture.output(sessionInfo()),file.path(support,"sessionInfo.txt")),add=TRUE)
  for(p in scripts) if(!file.copy(p,file.path(support,"scripts_used",basename(p)))) stop("Script snapshot failed.")
  for(n in c("scaling","environment_TS","audit","coverage","inclusion","fingerprints"))
    mc19_write_csv(input[[n]],file.path(support,paste0("input_",n,".csv")))
  mc19_write_csv(do.call(rbind,input$blocks),file.path(support,"analysis_input.csv"))
  mc19_write_csv(source_hashes,file.path(support,"reused_source_fingerprints.csv"))
  mc19_write_csv(script_hashes,file.path(support,"script_fingerprints.csv"))
  file.copy(cfg$workbook,file.path(support,"input_workbook.xlsx"))
  file.copy(mc19_temperature_source(cfg$project_dir),file.path(support,"PM11_temperature_source.csv"))
  m <- list(version=cfg$version,cfg=cfg,run_dir=run,analysis_scope=cfg$analysis_scope,
    n_samples=input$scaling$n_samples[1],n_TS=7L,scaling=input$scaling,model_files=character(),
    fit_complete=FALSE,expected_models=95L,q_family_size=95L,input_hashes=input$fingerprints,
    new_fit_count=0L,reused_fit_count=0L,bootstrap=FALSE,permutations=FALSE,LOTO=FALSE)
  rows <- audit <- list()
  mc19_write_text("incremental_models_in_progress",file.path(support,"run_status.txt"))
  for(f in cfg$environments) for(mc in cfg$mc_ids) {
    key <- mc19_wald_key(mc,f,"BB_TS");d <- input$blocks[[paste(mc,f,sep="__")]]
    entry <- reused[[key]]
    if(is.null(entry)) {
      mc19_incremental_assert(f=="SST_C","attempt to refit a nutrient")
      message("New temperature model: ",mc)
      entry <- mc19_wald_fit_one(mc,d,cfg);m$new_fit_count <- m$new_fit_count+1L;origin <- "new_temperature_fit"
    } else {
      m$reused_fit_count <- m$reused_fit_count+1L;origin <- if(f=="SST_C")"reused_temperature_fit" else "reused_legacy_nutrient_fit"
      entry$data <- d
    }
    entry <- mc19_wald_annotate(entry,f,input$scaling[input$scaling$environment==f,],origin,cfg)
    path <- file.path(support,"models",paste0(key,".rds"));saveRDS(entry,path);m$model_files[key] <- path
    copied <- mc19_incremental_copy_diagnostic(entry,path,support)
    audit[[key]] <- data.frame(key=key,environment=f,origin=origin,source=entry$result$source_model_file,diagnostic_reused=isTRUE(copied))
    rows[[key]] <- entry$result
    mc19_wald_write_result(mc19_incremental_bind(rows),run)
    saveRDS(m,file.path(support,"wald_manifest.rds"))
  }
  results <- mc19_incremental_bind(rows)
  valid <- is.finite(results$beta_P_value) & results$numerical_ok %in% TRUE
  results$beta_q_value <- results$q_value <- NA_real_
  results$beta_q_value[valid] <- p.adjust(results$beta_P_value[valid],"BH",n=95L)
  results$q_value <- results$beta_q_value
  mc19_incremental_assert(nrow(results)==95L && m$new_fit_count==plan$new_temperature_models &&
    m$reused_fit_count+m$new_fit_count==95L,"incremental model counts")
  mc19_incremental_assert(identical(mc19_hash(source_hashes$path)$md5,source_hashes$md5) &&
    identical(mc19_hash(input$fingerprints$path)$md5,input$fingerprints$md5) &&
    identical(mc19_hash(scripts)$md5,script_hashes$md5),"source changed during run")
  mc19_wald_write_result(results,run);mc19_write_csv(do.call(rbind,audit),file.path(support,"reuse_audit.csv"))
  m$model_hashes <- mc19_hash(unname(m$model_files));m$fit_complete <- TRUE
  saveRDS(m,file.path(support,"wald_manifest.rds"));mc19_write_csv(m$model_hashes,file.path(support,"model_fingerprints.csv"))
  mc19_write_text(c("Four nutrients + temperature; DIN excluded. BH recomputed across all 95 raw coefficient P values.",
    paste("New temperature fits:",m$new_fit_count,"; reused fits:",m$reused_fit_count),
    "Old results/models are read-only; reused diagnostics retain their original simulations.",
    "No PERMANOVA or NMDS is run by this GLMM step."),file.path(run,"00_readme.txt"))
  mc19_write_text("Wald_fits_complete_diagnostics_pending",file.path(support,"run_status.txt"))
  mc19_publish_result(cfg$project_dir,"GLMM",cfg$analysis_scope,run)
  invisible(run)
}

mc19_incremental_tp_cache <- function(input,project,key,component) {
  source <- getOption("mc19.temperature.legacy_two_part",file.path(project,"results",
    if(input$scope=="coral")"two_part_coral_v1" else "two_part_all_hosts_v1"))
  support <- file.path(source,"supporting_files")
  path <- file.path(support,"models",paste0(key,"__",component,".rds"))
  if(!file.exists(path)) return(NULL)
  # Check that the original fitting procedure is unchanged; input selection is intentionally new.
  saved_script <- file.path(support,"mc19_two_part_helpers.R")
  expressions <- parse(saved_script)
  saved_body <- NULL
  for(x in expressions) if(is.call(x) && identical(x[[1]],as.name("<-")) &&
    identical(x[[2]],as.name("mc19_tp_fit"))) saved_body <- x[[3]][[3]]
  canonical <- function(x) paste(deparse(x,width.cutoff=500L),collapse="\n")
  mc19_incremental_assert(!is.null(saved_body) && identical(canonical(saved_body),canonical(body(mc19_tp_fit))),"two-part fitting method changed")
  original <- read.csv(file.path(support,"analysis_input.csv"),stringsAsFactors=FALSE)
  d <- input$data[input$data$key==key,,drop=FALSE]
  old <- original[original$key==key,,drop=FALSE]
  mc19_incremental_equal_data(old,d,c("New_ID","Sample_ID","Host_type","Host_niche","TS","MC_ID","environment",
    "MC_reads","Ruegeria_reads_total","MC_proportion","E_raw","E_z"))
  cache <- readRDS(path)
  selected <- if(component=="positive") d[d$MC_reads>0,,drop=FALSE] else d
  if(!is.null(cache$fit)) {
    mf <- model.frame(cache$fit)
    response <- if(component=="positive") selected$MC_proportion else as.integer(selected$MC_reads>0)
    mc19_incremental_assert(mc19_incremental_same(model.response(mf),response) &&
      mc19_incremental_same(mf$E_z,selected$E_z) &&
      identical(as.character(mf$TS),as.character(selected$TS)),"two-part response/design")
    if(isTRUE(cache$result$numerical_ok)) {
      co <- summary(cache$fit)$coefficients$cond["E_z",]
      mc19_incremental_assert(mc19_incremental_same(unname(co[c(1,2,4)]),
        c(cache$result$beta,cache$result$SE,cache$result$P),1e-10),"two-part raw estimates")
    }
  }
  cache$reused_from <- normalizePath(path,winslash="/",mustWork=TRUE)
  cache$source_md5 <- unname(tools::md5sum(path))
  cache
}
