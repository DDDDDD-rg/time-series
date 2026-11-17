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
      saveRDS(list(analysis_scope=scope, version=paste0("MC19_env5_TS8_",scope,"_Wald_v1"),
        n_TS=8L, fit_complete=TRUE, model_files=rep("fixture",95)),file.path(s,"wald_manifest.rds"))
      writeLines("fixture",file.path(path,"00_result.csv")); writeLines("fixture",file.path(s,"all_results.csv"))
    } else if (kind == "PERMANOVA") {
      writeLines("complete",file.path(s,"00_run_status.txt"))
      write.csv(data.frame(Scope=rep(scope,5),Status="ok",N_TS=8L,Analysis_version=paste0("MC19_env5_TS8_",scope,"_individual_TSlabels_v1_vegan")),file.path(path,"environment_results.csv"),row.names=FALSE)
      write.csv(data.frame(Scope=rep(scope,5),Status="ok",N_TS=8L,Analysis_version=paste0("MC19_env5_TS8_",scope,"_individual_TSlabels_v1_vegan")),file.path(s,"environment_results_full.csv"),row.names=FALSE)
      saveRDS(list(),file.path(s,paste0(scope,"_environment_analysis.rds")))
    } else {
      parent <- file.path(root,paste0("parent_",scope))
      complete(parent,"GLMM",scope)
      saveRDS(list(scope=scope,source=parent,keys=c("a","b")),file.path(s,"manifest.rds"))
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
  # An old seven-TS manifest must not be accepted by current result discovery.
  old <- file.path(root,"results","MC19_GLMM_coral_run1001")
  complete(old,"GLMM","coral")
  mf <- file.path(old,"supporting_files","wald_manifest.rds")
  m <- readRDS(mf);m$n_TS <- 7L;m$version <- "MC19_env4_coral_TS_Wald_v2";saveRDS(m,mf)
  stopifnot(!mc19_result_complete(old,"GLMM","coral"),
    identical(mc19_latest_result(root,"GLMM","coral"),n1000))
  # Explicit archived paths remain available if they contain current 8-TS results.
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
  cat("PASS: numbering, reserved IDs, completion checks, stale pointers, scope isolation and exclusion of old 7TS results. No production output changed.\n")
})
