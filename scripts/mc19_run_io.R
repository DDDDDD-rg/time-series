# Shared naming and completed-run discovery. Definitions only; no fitting or plotting.
# Method versions remain in the analysis manifests, independent of run numbers.
mc19_io_scope <- function(scope) {
  stopifnot(length(scope) == 1L, scope %in% c("coral", "ecosystem"))
  scope
}
mc19_io_prefix <- function(kind, scope) {
  stopifnot(kind %in% c("GLMM", "PERMANOVA", "two_part", "talk"))
  paste0("MC19_", kind, "_", mc19_io_scope(scope))
}
mc19_io_time <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%OS6Z", tz="UTC")
mc19_io_record <- function(path) {
  f <- file.path(path, "run_record.csv")
  if (!file.exists(f)) return(NULL)
  tryCatch(utils::read.csv(f, stringsAsFactors=FALSE, check.names=FALSE), error=function(e) NULL)
}
mc19_io_update <- function(path, ...) {
  x <- mc19_io_record(path)
  if (is.null(x) || nrow(x) != 1L) stop("Missing run record: ", path)
  for (n in names(list(...))) x[[n]] <- list(...)[[n]]
  utils::write.csv(x, file.path(path, "run_record.csv"), row.names=FALSE, fileEncoding="UTF-8")
  invisible(path)
}
mc19_io_aliases <- function(project_dir, kind, scope, area="results") {
  names <- switch(kind,
    GLMM=c(paste0("00_latest_GLMM_MC19_env4_TS_Wald_", scope, ".txt"),
      if (scope == "coral") "00_latest_GLMM_MC19_env4_TS_Wald.txt"),
    PERMANOVA=paste0("00_latest_PERMANOVA_MC19_env7_", scope, ".txt"),
    two_part=paste0("00_latest_GLMM_MC19_two_part_", scope, ".txt"), character())
  file.path(project_dir, area, names)
}
mc19_io_pointer <- function(project_dir, kind, scope, area="results")
  file.path(project_dir, area, paste0("00_latest_", mc19_io_prefix(kind, scope), ".txt"))

mc19_io_reserve <- function(project_dir, prefix, suffix, area, requested=NULL) {
  # Persistent reservations survive output archiving. Atomic mkdir avoids two
  # sessions allocating the same ID; failed attempts leave a gap, never reuse.
  registry <- file.path(project_dir, ".run_registry")
  dir.create(registry, recursive=TRUE, showWarnings=FALSE)
  if (!dir.exists(registry)) stop("Cannot create run registry: ", registry)
  root <- file.path(project_dir, area)
  dir.create(root, recursive=TRUE, showWarnings=FALSE)
  pattern <- paste0("^", prefix, "_", suffix, "([0-9]{3,})$")
  used <- unique(c(list.files(registry, all.files=TRUE),
    list.files(file.path(project_dir, "results"), all.files=TRUE),
    list.files(file.path(project_dir, "figures"), all.files=TRUE)))
  ids <- as.numeric(sub(pattern, "\\1", used[grepl(pattern, used)]))
  next_id <- if (length(ids)) max(ids) + 1 else 1
  if (!is.null(requested) && nzchar(as.character(requested))) {
    value <- sub(paste0("^", suffix), "", as.character(requested))
    if (length(value) != 1L || !grepl("^[0-9]+$", value) || !is.finite(as.numeric(value)) || as.numeric(value) < next_id)
      stop("Run ID must be an unused increasing number (e.g. run001); next available: ", next_id)
    next_id <- as.numeric(value)
  }
  repeat {
    name <- paste0(prefix, "_", suffix, sprintf("%03.0f", next_id))
    path <- file.path(root, name)
    reservation <- file.path(registry, name)
    if (dir.exists(path) || file.exists(path)) {
      if (!is.null(requested)) stop("Run ID already reserved: ", name)
      next_id <- next_id + 1; next
    }
    if (!dir.create(reservation, showWarnings=FALSE)) {
      if (!dir.exists(reservation)) stop("Cannot reserve run ID: ", reservation)
      if (!is.null(requested)) stop("Run ID already reserved: ", name)
      next_id <- next_id + 1; next
    }
    writeLines(c(paste("created_at:", mc19_io_time()), paste("target:", path)),
      file.path(reservation, "reservation.txt"))
    if (!dir.create(path, showWarnings=FALSE)) stop("Cannot create reserved output: ", path)
    return(normalizePath(path, winslash="/", mustWork=TRUE))
  }
}
mc19_new_run <- function(project_dir, kind, scope, version="", source="", requested=NULL) {
  area <- if (kind == "talk") "figures" else "results"
  path <- mc19_io_reserve(project_dir, mc19_io_prefix(kind, scope), "run", area, requested)
  record <- data.frame(kind=kind, scope=scope, method_version=version, source=source,
    created_at=mc19_io_time(), finished_at="", state="running", stringsAsFactors=FALSE)
  utils::write.csv(record, file.path(path, "run_record.csv"), row.names=FALSE, fileEncoding="UTF-8")
  path
}
mc19_result_complete <- function(path, kind, scope) {
  # Only read small manifests/tables. Discovery never creates files or loads fits.
  tryCatch({
    s <- file.path(path, "supporting_files")
    if (!dir.exists(path)) return(FALSE)
    if (kind == "GLMM") {
      m <- readRDS(file.path(s, "wald_manifest.rds"))
      saved_scope <- if (is.null(m$analysis_scope)) "coral" else m$analysis_scope
      expected <- if (scope == "coral") "MC19_env4_coral_TS_Wald_v2" else "MC19_env4_ecosystem_TS_Wald_v1"
      return(identical(saved_scope, scope) && identical(m$version, expected) &&
        isTRUE(m$fit_complete) && length(m$model_files) == 76L &&
        file.exists(file.path(path, "00_result.csv")) && file.exists(file.path(s, "all_results.csv")))
    }
    if (kind == "PERMANOVA") {
      status <- readLines(file.path(s, "00_run_status.txt"), n=1L, warn=FALSE)
      t <- utils::read.csv(file.path(s, "environment_results_full.csv"), fileEncoding="UTF-8-BOM")
      return(identical(status, "complete") && nrow(t) == 7L &&
        all(c("Scope", "Status") %in% names(t)) && !anyNA(t$Scope) && !anyNA(t$Status) &&
        all(t$Scope == scope) && all(t$Status == "ok") &&
        file.exists(file.path(path, "environment_results.csv")) &&
        file.exists(file.path(s, paste0(scope, "_environment_analysis.rds"))))
    }
    if (kind == "two_part") {
      status <- readLines(file.path(path, "run_status.txt"), n=1L, warn=FALSE)
      m <- readRDS(file.path(s, "manifest.rds"))
      t <- utils::read.csv(file.path(s, "all_results.csv"))
      return(identical(m$scope, scope) && length(status) == 1L &&
        status %in% c("fits_complete_review_diagnostics", "fits_complete_some_require_review") &&
        nrow(t) == 2L * length(m$keys) && file.exists(file.path(path, "00_result.csv")) &&
        all(file.exists(file.path(s, c("predictions.rds", "observations.rds")))))
    }
    FALSE
  }, error=function(e) FALSE, warning=function(w) FALSE)
}
mc19_latest_result <- function(project_dir, kind, scope, explicit="") {
  mc19_io_prefix(kind, scope)
  if (length(explicit) != 1L || is.na(explicit)) stop("Supply one result directory.")
  if (nzchar(explicit) && explicit != "latest") {
    path <- normalizePath(explicit, winslash="/", mustWork=TRUE)
    if (!mc19_result_complete(path, kind, scope)) stop("Result is incomplete or belongs to a different analysis/scope: ", path)
    return(path)
  }
  pointers <- c(mc19_io_pointer(project_dir, kind, scope), mc19_io_aliases(project_dir, kind, scope))
  pointed <- unlist(lapply(pointers[file.exists(pointers)], function(p)
    trimws(readLines(p, n=1L, warn=FALSE, encoding="UTF-8"))), use.names=FALSE)
  legacy <- switch(kind, GLMM=paste0("GLMM_MC19_env4_TS_Wald_", scope),
    PERMANOVA=paste0("PERMANOVA_MC19_env7_", scope, "_samples"),
    two_part=paste0("GLMM_MC19_two_part_", scope))
  pattern <- paste0("^(", mc19_io_prefix(kind, scope), "_run[0-9]{3,}|", legacy, "_[0-9]{8}_[0-9]{6})$")
  discovered <- list.files(file.path(project_dir, "results"), pattern=pattern, full.names=TRUE)
  paths <- unique(c(pointed[nzchar(pointed)], discovered))
  paths <- paths[vapply(paths, mc19_result_complete, logical(1), kind=kind, scope=scope)]
  if (!length(paths)) stop("No completed ", kind, " result for ", scope, ". Supply an existing directory explicitly.")
  # New numbered runs follow legacy timestamps. Compare run numbers numerically,
  # not lexical names or directory mtimes (which change when plotting/archiving).
  score <- lapply(paths, function(p) {
    name <- basename(p); record <- mc19_io_record(p)
    is_new <- grepl(paste0("^", mc19_io_prefix(kind, scope), "_run[0-9]+$"), name)
    number <- if (is_new) as.numeric(sub(".*_run", "", name)) else 0
    stamp <- if (grepl("_[0-9]{8}_[0-9]{6}$", name)) sub(".*_([0-9]{8})_([0-9]{6})$", "\\1\\2", name) else ""
    if (!is.null(record) && identical(record$kind, kind) && identical(record$scope, scope)) stamp <- record$created_at
    c(new=as.integer(is_new), number=number, stamp=stamp)
  })
  scores <- do.call(rbind, score)
  index <- order(as.numeric(scores[,"new"]), as.numeric(scores[,"number"]), scores[,"stamp"], decreasing=TRUE)[1]
  normalizePath(paths[index], winslash="/", mustWork=TRUE)
}
mc19_publish_result <- function(project_dir, kind, scope, path) {
  if (!mc19_result_complete(path, kind, scope)) stop("Refusing to publish incomplete result: ", path)
  if (!is.null(mc19_io_record(path))) mc19_io_update(path, state="complete", finished_at=mc19_io_time())
  for (p in c(mc19_io_pointer(project_dir, kind, scope), mc19_io_aliases(project_dir, kind, scope)))
    writeLines(enc2utf8(normalizePath(path, winslash="/", mustWork=TRUE)), p, useBytes=TRUE)
  invisible(path)
}

mc19_plot_directory <- function(project_dir, source, kind, scope, requested=NULL) {
  if (!is.null(requested) && nzchar(requested)) {
    if (!dir.exists(requested)) dir.create(requested, recursive=TRUE)
    path <- normalizePath(requested, winslash="/", mustWork=TRUE)
  } else {
    prefix <- mc19_io_prefix(kind, scope)
    numbered <- grepl(paste0("^", prefix, "_run[0-9]+$"), basename(source))
    base <- if (numbered) basename(source) else paste0(prefix, "_legacy")
    path <- mc19_io_reserve(project_dir, base, "plot", "figures")
  }
  old <- mc19_io_record(path)
  if (!is.null(old) && (!identical(old$source, source) || !identical(old$scope, scope)))
    stop("Plot directory belongs to a different source/scope: ", path)
  record <- data.frame(kind=kind, scope=scope, method_version="", source=source,
    created_at=mc19_io_time(), finished_at="", state="running", stringsAsFactors=FALSE)
  utils::write.csv(record, file.path(path, "run_record.csv"), row.names=FALSE, fileEncoding="UTF-8")
  path
}
mc19_publish_plot <- function(project_dir, kind, scope, path) {
  mc19_io_update(path, state="complete", finished_at=mc19_io_time())
  for (p in c(mc19_io_pointer(project_dir, kind, scope, "figures"), mc19_io_aliases(project_dir, kind, scope, "figures")))
    writeLines(enc2utf8(normalizePath(path, winslash="/", mustWork=TRUE)), p, useBytes=TRUE)
  invisible(path)
}
mc19_latest_talk <- function(project_dir, scope) {
  prefix <- mc19_io_prefix("talk", scope)
  paths <- list.files(file.path(project_dir, "figures"),
    pattern=paste0("^(", prefix, "_run[0-9]+|MC19_talk_", scope, "_[0-9]{8}_[0-9]{6})$"), full.names=TRUE)
  good <- vapply(paths, function(p) {
    m <- mc19_io_record(p)
    if (!is.null(m)) return(identical(m$kind, "talk") && identical(m$scope, scope) && identical(m$state, "complete"))
    file.exists(file.path(p, "README.txt"))
  }, logical(1))
  paths <- paths[good]
  if (!length(paths)) stop("No completed talk figures for ", scope)
  numbered <- grepl("_run[0-9]+$", basename(paths))
  numbers <- numeric(length(paths))
  numbers[numbered] <- as.numeric(sub(".*_run", "", basename(paths[numbered])))
  normalizePath(paths[order(numbered, numbers, basename(paths), decreasing=TRUE)[1]], winslash="/", mustWork=TRUE)
}
