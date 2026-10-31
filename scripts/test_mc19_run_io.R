# File-I/O regression checks only. No analysis entrypoint, model, NMDS or plot runs.
local({
  args <- commandArgs(trailingOnly=TRUE)
  project <- if (length(args)) normalizePath(args[1], winslash="/", mustWork=TRUE) else
    normalizePath(file.path(getwd(), "Research", "Msc main project-time series"), winslash="/", mustWork=TRUE)
  sys.source(file.path(project, "scripts", "mc19_run_io.R"), envir=environment())
  root <- tempfile("mc19-naming-")
  dir.create(root); dir.create(file.path(root, "results")); dir.create(file.path(root, "figures"))
  on.exit(unlink(root, recursive=TRUE), add=TRUE) # Explicitly created temporary fixture only.
  fail <- function(expr) stopifnot(inherits(tryCatch({force(expr); NULL}, error=identity), "error"))
  complete <- function(path, kind, scope) {
    s <- file.path(path, "supporting_files"); dir.create(s, recursive=TRUE, showWarnings=FALSE)
    if (kind == "GLMM") {
      saveRDS(list(analysis_scope=scope, version=if(scope=="coral") "MC19_env4_coral_TS_Wald_v2" else "MC19_env4_ecosystem_TS_Wald_v1",
        fit_complete=TRUE, model_files=rep("fixture",76)),file.path(s,"wald_manifest.rds"))
      writeLines("fixture",file.path(path,"00_result.csv")); writeLines("fixture",file.path(s,"all_results.csv"))
    } else if (kind == "PERMANOVA") {
      writeLines("complete",file.path(s,"00_run_status.txt"))
      write.csv(data.frame(Scope=rep(scope,7),Status="ok"),file.path(path,"environment_results.csv"),row.names=FALSE)
      write.csv(data.frame(Scope=rep(scope,7),Status="ok"),file.path(s,"environment_results_full.csv"),row.names=FALSE)
      saveRDS(list(),file.path(s,paste0(scope,"_environment_analysis.rds")))
    } else {
      saveRDS(list(scope=scope,keys=c("a","b")),file.path(s,"manifest.rds"))
      writeLines("fits_complete_some_require_review",file.path(path,"run_status.txt"))
      write.csv(data.frame(beta=1:4),file.path(s,"all_results.csv"),row.names=FALSE)
      writeLines("fixture",file.path(path,"00_result.csv"))
      saveRDS(data.frame(),file.path(s,"predictions.rds")); saveRDS(data.frame(),file.path(s,"observations.rds"))
    }
  }
  for (kind in c("GLMM","PERMANOVA","two_part")) for (scope in c("coral","ecosystem")) {
    one <- mc19_new_run(root,kind,scope,"method_v1")
    stopifnot(endsWith(one,"run001"))
    fail(mc19_publish_result(root,kind,scope,one))
    complete(one,kind,scope); mc19_publish_result(root,kind,scope,one)
    two <- mc19_new_run(root,kind,scope,"method_v1")
    stopifnot(endsWith(two,"run002"), identical(mc19_latest_result(root,kind,scope),one))
    fail(mc19_latest_result(root,kind,scope,two))
    fail(mc19_latest_result(root,kind,if(scope=="coral") "ecosystem" else "coral",one))
    complete(two,kind,scope)
    # A stale or missing pointer must not hide a newer complete directory.
    stopifnot(identical(mc19_latest_result(root,kind,scope),two))
    writeLines("missing/path",mc19_io_pointer(root,kind,scope))
    stopifnot(identical(mc19_latest_result(root,kind,scope),two))
    mc19_publish_result(root,kind,scope,two)
    alias <- mc19_io_aliases(root,kind,scope)
    stopifnot(all(vapply(alias,function(p) identical(readLines(p,warn=FALSE),two),logical(1))))
    # Moving an output to archive does not release its ID.
    dir.create(file.path(root,"archive"),showWarnings=FALSE)
    moved <- file.path(root,"archive",basename(two)); complete(moved,kind,scope)
    # Synthetic fixture only: simulate an archived run without Windows directory
    # rename permissions, while keeping the persistent reservation untouched.
    stopifnot(startsWith(normalizePath(two,winslash="/"),paste0(normalizePath(root,winslash="/"),"/")))
    unlink(two,recursive=TRUE); stopifnot(!dir.exists(two))
    three <- mc19_new_run(root,kind,scope,"method_v1")
    stopifnot(endsWith(three,"run003"))
    fail(mc19_new_run(root,kind,scope,requested="run002"))
    ten <- mc19_new_run(root,kind,scope,requested="run010"); complete(ten,kind,scope)
    stopifnot(identical(mc19_latest_result(root,kind,scope),ten))
  }
  n999 <- mc19_new_run(root,"GLMM","coral",requested="run999");complete(n999,"GLMM","coral")
  n1000 <- mc19_new_run(root,"GLMM","coral",requested="run1000");complete(n1000,"GLMM","coral")
  stopifnot(identical(mc19_latest_result(root,"GLMM","coral"),n1000))
  # Legacy timestamps are selected chronologically, and cannot cross analysis/scope.
  legacy_root <- file.path(root,"legacy");dir.create(legacy_root);dir.create(file.path(legacy_root,"results"))
  for (kind in c("GLMM","PERMANOVA","two_part")) for (scope in c("coral","ecosystem")) {
    prefix <- switch(kind,GLMM=paste0("GLMM_MC19_env4_TS_Wald_",scope),
      PERMANOVA=paste0("PERMANOVA_MC19_env7_",scope,"_samples"),two_part=paste0("GLMM_MC19_two_part_",scope))
    a <- file.path(legacy_root,"results",paste0(prefix,"_20261030_120000"));complete(a,kind,scope)
    b <- file.path(legacy_root,"results",paste0(prefix,"_20261031_120000"));complete(b,kind,scope)
    stopifnot(identical(mc19_latest_result(legacy_root,kind,scope),normalizePath(b,winslash="/")))
  }
  # Explicit complete archived paths and the old coral-only alias remain readable.
  archived <- file.path(root,"archive","MC19_GLMM_coral_run002")
  stopifnot(identical(mc19_latest_result(root,"GLMM","coral",archived),normalizePath(archived,winslash="/")))
  p1 <- mc19_plot_directory(root, archived,"GLMM","coral")
  p2 <- mc19_plot_directory(root, archived,"GLMM","coral")
  stopifnot(p1 != p2, endsWith(p1,"plot001"),endsWith(p2,"plot002"))
  mc19_publish_plot(root,"GLMM","coral",p1)
  fail(mc19_plot_directory(root, "another/source", "GLMM","coral",p1))
  talk <- mc19_new_run(root,"talk","coral","method_v1",archived)
  fail(mc19_latest_talk(root,"coral"))
  mc19_publish_plot(root,"talk","coral",talk)
  stopifnot(identical(mc19_latest_talk(root,"coral"),talk))
  # Verify existing production directories through read-only code only.
  production <- c()
  for (kind in c("GLMM","PERMANOVA","two_part")) for (scope in c("coral","ecosystem")) {
    p <- mc19_latest_result(project,kind,scope)
    stopifnot(dir.exists(p)); production <- c(production,p)
  }
  e <- new.env(parent=globalenv())
  for (f in c("glmm_MC19_pilot_helpers.R","glmm_MC19_wald_helpers.R"))
    sys.source(file.path(project,"scripts","GLMM",f),e)
  for (scope in c("coral","ecosystem")) {
    x <- e$mc19_wald_read_run(e$mc19_wald_config(project,scope))
    stopifnot(nrow(x$results)==76L)
  }
  # Actual downstream readers must still accept both existing scopes. Sourcing
  # these definition files does not invoke GLMM, PERMANOVA, NMDS or graphics.
  sys.source(file.path(project,"scripts","GLMM","mc19_two_part_helpers.R"),e)
  old_options <- options(mc19.analysis_scope="coral")
  on.exit(options(old_options),add=TRUE)
  for (scope in c("coral","ecosystem")) {
    expected <- if(scope=="coral") 8L else 14L
    stopifnot(nrow(e$mc19_tp_input(project,scope)$selected)==expected)
    options(mc19.analysis_scope=scope)
    for (f in c("permanova_MC19_data_settings.R","permanova_MC19_plot_helpers.R"))
      sys.source(file.path(project,"scripts","PERMANOVA",f),e)
    ctx <- e$mc19_plot_context()
    stopifnot(identical(ctx$input$scope,scope) || identical(ctx$saved$scope,scope) || all(ctx$saved$results$Scope==scope))
  }
  cat("PASS: numbering, archived IDs, completion filtering, stale pointers, scope isolation, legacy paths, plot revisions; six existing results and both saved GLMM readers verified read-only.\n")
})
