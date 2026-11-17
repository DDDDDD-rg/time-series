# Resolve archived manifest paths in memory; preserve saved RDS and fingerprints.
mc19_relocate_saved_paths <- function(x) {
  project <- getOption("mc19.project_dir", "")
  if (!nzchar(project)) return(x)
  relocate <- function(v) {
    v <- gsub("\\\\", "/", v)
    v <- gsub("results/joint_glmm_TS7_coral_run001", "archive/superseded_results/joint_env4_coral_v1", v, fixed=TRUE)
    v <- gsub("results/joint_glmm_TS7_ecosystem_run001", "archive/superseded_results/joint_env4_ecosystem_v1", v, fixed=TRUE)
    v <- gsub("results/MC19_GLMM_TS7_temp_coral_run001", "archive/superseded_results/glmm_temp_coral_v1", v, fixed=TRUE)
    v <- gsub("results/MC19_GLMM_TS7_temperature_coral_run001", "archive/superseded_results/glmm_temp_early_coral_v1", v, fixed=TRUE)
    v <- gsub("results/MC19_GLMM_TS7_temp_ecosystem_run001", "archive/superseded_results/glmm_temp_ecosystem_v1", v, fixed=TRUE)
    v <- gsub("results/MC19_PERMANOVA_TS7_temp_coral_run001", "archive/superseded_results/permanova_temp_coral_v1", v, fixed=TRUE)
    v <- gsub("results/MC19_PERMANOVA_TS7_temp_ecosystem_run001", "archive/superseded_results/permanova_temp_ecosystem_v1", v, fixed=TRUE)
    v <- gsub("results/MC19_two_part_TS7_temp_coral_run001", "archive/superseded_results/two_part_temp_coral_v1", v, fixed=TRUE)
    v <- gsub("results/MC19_two_part_TS7_temp_ecosystem_run001", "archive/superseded_results/two_part_temp_ecosystem_v1", v, fixed=TRUE)
    v
  }
  if (is.character(x)) { n <- names(x); x <- relocate(x); if(!is.null(n)) names(x) <- relocate(n) }
  else if (is.list(x)) { x[] <- lapply(x, mc19_relocate_saved_paths); if(!is.null(names(x))) names(x) <- relocate(names(x)) }
  x
}

# Current analysis: definitions only. Source after glmm_MC19_pilot_helpers.R.
# 19 MC x 5 nutrients, TS random intercept only; one BH family of 95.

mc19_wald_config <- function(project_dir, analysis_scope="coral") {
  cfg <- mc19_pilot_config(project_dir)
  sys.source(file.path(project_dir,"scripts","mc19_temperature_incremental.R"),envir=environment(mc19_wald_config))
  stopifnot(length(analysis_scope) == 1L, analysis_scope %in% c("coral", "ecosystem"))
  cfg$analysis_scope <- analysis_scope
  cfg$expected_samples <- if(mc19_profile()$temperature) c(coral=98L,ecosystem=310L)[[analysis_scope]] else NULL
  cfg$version <- mc19_method_version(analysis_scope,"GLMM")
  cfg$mc_ids <- cfg$all_mc_ids
  cfg$environments <- cfg$nutrient_fields
  cfg$model_ids <- "BB_TS"
  cfg$short_names <- c(Nitrate_mgN_L = "NO3", Nitrite_mgN_L = "NO2",
    Ammonia_mgN_L = "NH4", DIN_mgN_L = "DIN", DIP_mgP_L = "DIP", SST_C="Temp")
  cfg$axis_labels <- c(Nitrate_mgN_L = "Nitrate (mg N/L)",
    Nitrite_mgN_L = "Nitrite (mg N/L)", Ammonia_mgN_L = "Ammonia (mg N/L)",
    DIN_mgN_L = "DIN (mg N/L)", DIP_mgP_L = "DIP (mg P/L)", SST_C="Temperature (deg C)")
  cfg$n_models <- cfg$q_family_size <- 19L*length(cfg$environments)
  cfg$scripts <- c("GLMM.R", "GLMM_reads_diagnostics_full.R", "plot_GLMM_MC19.R",
    "glmm_MC19_wald_helpers.R", "glmm_MC19_wald_plot_helpers.R", "glmm_MC19_pilot_helpers.R", "../mc19_run_io.R", "../mc19_input.R", "../mc19_profile.R", "../TS7_temperature.R", "../mc19_temperature_incremental.R")
  cfg
}

mc19_wald_key <- function(mc, environment, model) paste(mc, environment, model, sep = "__")

mc19_wald_load_input <- function(cfg) {
  # The established reader verifies IDs, five nutrient fields, censoring/L/2,
  # all 19 proportions and the denominator from ALL MC rows in the raw table.
  reader_cfg <- mc19_pilot_config(cfg$project_dir)
  reader_cfg$workbook <- cfg$workbook
  reader_cfg$counts_file <- cfg$counts_file
  reader_cfg$analysis_scope <- cfg$analysis_scope
  reader_cfg$expected_samples <- cfg$expected_samples
  reader_cfg$complete_fields <- cfg$environments
  verified <- mc19_load_input(reader_cfg)
  template <- verified$data[verified$data$MC_ID == "MC10", , drop = FALSE]
  template <- template[, setdiff(names(template),
    c("MC_ID", "MC_reads", "Other_reads", "MC_proportion", "E_raw", "E_z")), drop = FALSE]
  if ((!is.null(cfg$expected_samples) && nrow(template) != cfg$expected_samples) ||
      nlevels(template$TS) != cfg$expected_TS || anyDuplicated(template$New_ID)) stop("样本或TS配置不符。")
  raw <- utils::read.delim(cfg$counts_file, check.names = FALSE, stringsAsFactors = FALSE,
    quote = "", comment.char = "")
  counts <- as.matrix(raw[, template$New_ID, drop = FALSE])
  if (!mc19_same(colSums(counts), template$Ruegeria_reads_total)) stop("两次读取的reads分母不一致。")
  environment_TS <- verified$ts_values
  environment_TS$E_z <- NULL
  scaling <- list(); blocks <- list(); audit <- list()
  for (environment in cfg$environments) {
    values <- environment_TS[[environment]]
    if (any(!is.finite(values)) || any(values < 0) || length(values) != nlevels(template$TS) ||
        any(!is.finite(template[[environment]]))) stop("All selected environments must cover the same samples: ", environment)
    center <- mean(values)
    scale <- stats::sd(values)
    if (!is.finite(scale) || scale <= 0) stop("Environment has no estimable variation: ", environment)
    scaling[[environment]] <- data.frame(environment = environment, center = center, scale = scale,
      E_low = min(values), E_high = max(values),
      censored_TS = sum(environment_TS[[paste0(environment, "_censored")]]),
      n_TS = nlevels(template$TS), n_samples = nrow(template), stringsAsFactors = FALSE)
    for (mc in cfg$mc_ids) {
      d <- template
      d$MC_ID <- mc; d$environment <- environment
      d$E_raw <- d[[environment]]
      d$E_z <- (d$E_raw - center) / scale
      d$MC_reads <- unname(counts[match(mc, raw$MC), ])
      d$Other_reads <- d$Ruegeria_reads_total - d$MC_reads
      d$MC_proportion <- d$MC_reads / d$Ruegeria_reads_total
      blocks[[paste(mc, environment, sep = "__")]] <- d
      if (environment == cfg$environments[1]) audit[[mc]] <- data.frame(MC_ID = mc,
        n_samples = nrow(d), zero_count = sum(d$MC_reads == 0),
        positive_TS = length(unique(d$TS[d$MC_reads > 0])),
        median_proportion = stats::median(d$MC_proportion), max_proportion = max(d$MC_proportion))
    }
  }
  if (length(blocks) != cfg$n_models || !identical(mc19_hash(verified$fingerprints$path)$md5,
      verified$fingerprints$md5)) stop("任务数或读取中的输入指纹变化。")
  list(blocks = blocks, scaling = do.call(rbind, scaling), environment_TS = environment_TS,
    audit = do.call(rbind, audit), coverage = verified$coverage, inclusion = verified$exclusion,
    fingerprints = verified$fingerprints, all_MC_count = verified$all_MC_count)
}

mc19_wald_fit_one <- function(mc, d, cfg) {
  row <- mc19_empty_result(mc, "BB_TS", d)
  formula <- stats::as.formula("cbind(MC_reads, Other_reads) ~ E_z + (1 | TS)", env = baseenv())
  begin <- proc.time()[[3]]
  attempt <- mc19_capture(function() glmmTMB::glmmTMB(formula = formula, data = d,
    family = glmmTMB::betabinomial(link = "logit"), ziformula = ~0, dispformula = ~1,
    REML = FALSE, na.action = stats::na.fail,
    control = glmmTMB::glmmTMBControl(parallel = getOption("glmmTMB.cores",1L),
      optCtrl = list(iter.max = 10000L, eval.max = 10000L))))
  fit <- if (nzchar(attempt$error)) NULL else attempt$value
  row$warnings <- paste(attempt$warnings, collapse = " | ")
  row$error <- attempt$error
  if (!is.null(fit)) {
    extracted <- mc19_capture(function() mc19_extract_result(fit, row, d))
    if (nzchar(extracted$error)) {
      row$fit_status <- "result_extraction_failed"; row$error <- extracted$error
    } else row <- extracted$value
    row$warnings <- paste(unique(c(attempt$warnings, extracted$warnings)), collapse = " | ")
  }
  row$elapsed_seconds <- proc.time()[[3]] - begin
  list(MC_ID = mc, model = "BB_TS", data = d, fit = fit, result = row,
    formula = paste(deparse(formula), collapse = " "), version = cfg$version)
}

mc19_wald_annotate <- function(entry, environment, scale_row, origin, cfg) {
  r <- entry$result
  r$environment <- environment
  r$analysis_scope <- cfg$analysis_scope
  r$P_value <- r$q_value <- r$beta_P_value <- r$beta_q_value <- NA_real_
  r$z_value <- r$CI95_low <- r$CI95_high <- r$OR95_low <- r$OR95_high <- NA_real_
  r$P_method <- "glmmTMB summary conditional-coefficient Wald z test"
  r$CI_method <- "Wald normal approximation, beta +/- qnorm(0.975)*SE"
  r$q_family_size <- cfg$q_family_size
  r$inference_status <- "unavailable_numerical_issue"
  r$diagnostic_status <- "not_run"
  r$fit_origin <- origin
  r$source_model_file <- if (!is.null(entry$source_model_file)) entry$source_model_file else ""
  r$E_center <- scale_row$center; r$E_SD <- scale_row$scale
  r$E_low <- scale_row$E_low; r$E_high <- scale_row$E_high
  r$censored_TS <- scale_row$censored_TS
  if (isTRUE(r$numerical_ok) && !is.null(entry$fit)) {
    extracted <- mc19_capture(function() {
      tab <- summary(entry$fit)$coefficients$cond
      required <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
      if (!"E_z" %in% rownames(tab) || !all(required %in% colnames(tab))) stop("缺少E_z的Wald系数行。")
      values <- unname(tab["E_z", required])
      if (any(!is.finite(values)) || values[2] <= 0 || values[4] < 0 || values[4] > 1 ||
          !mc19_same(values[1], r$beta_per_TS_SD) || !mc19_same(values[2], r$SE)) stop("Wald系数与原估计不一致。")
      values
    })
    if (!nzchar(extracted$error)) {
      values <- extracted$value
      r$beta_per_TS_SD <- values[1]; r$SE <- values[2]; r$z_value <- values[3]
      r$beta_P_value <- r$P_value <- values[4] # P is supplied by glmmTMB, not simulated.
      r$CI95_low <- values[1] - stats::qnorm(.975) * values[2]
      r$CI95_high <- values[1] + stats::qnorm(.975) * values[2]
      r$Wald95_low_pilot <- r$CI95_low; r$Wald95_high_pilot <- r$CI95_high
      r$OR95_low <- exp(r$CI95_low); r$OR95_high <- exp(r$CI95_high)
      r$inference_status <- "Wald_approximation_available"
    } else {
      r$inference_status <- "Wald_extraction_failed"
      r$error <- paste(r$error, extracted$error)
    }
    r$warnings <- paste(unique(c(r$warnings, extracted$warnings)), collapse = " | ")
  }
  entry$environment <- environment
  entry$version <- cfg$version
  entry$result <- r
  entry
}

mc19_wald_result_columns <- function() {
  c("MC_ID", "environment", "analysis_scope", "beta_per_TS_SD", "CI95_low", "CI95_high",
    "beta_P_value", "beta_q_value", "marginal_change_pp", "odds_ratio_per_TS_SD",
    "SE", "TS_SD_logit", "numerical_ok", "diagnostic_status", "n_samples", "n_TS",
    "model", "z_value", "censored_TS", "fit_status", "inference_status", "q_family_size", "fit_origin")
}

mc19_wald_write_result <- function(results, run_dir) {
  columns <- mc19_wald_result_columns()
  detail_columns <- c(columns, setdiff(names(results), columns))
  mc19_write_csv(results[, detail_columns, drop = FALSE], file.path(run_dir, "supporting_files", "all_results.csv"))
  mc19_write_csv(results[, columns, drop = FALSE], file.path(run_dir, "00_result.csv"))
}

mc19_wald_run <- function(cfg) {
  if(mc19_profile()$incremental) return(mc19_incremental_glmm_run(cfg))
  if (!identical(cfg$model_ids, "BB_TS") || cfg$n_models != 19L*length(cfg$environments) || cfg$q_family_size != cfg$n_models)
    stop("Model count differs from the selected workflow.")
  mc19_require(c("readxl", "glmmTMB"))
  input <- mc19_wald_load_input(cfg)
  scripts <- file.path(cfg$project_dir, "scripts", "GLMM", cfg$scripts)
  script_hashes <- mc19_hash(scripts)
  run_dir <- mc19_new_run(cfg$project_dir, "GLMM", cfg$analysis_scope, cfg$version, cfg$workbook,
    requested=getOption("mc19.run_id", NULL))
  support <- file.path(run_dir, "supporting_files")
  dir.create(file.path(support, "models"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(support, "scripts_used"), showWarnings = FALSE)
  run_dir <- normalizePath(run_dir, winslash = "/", mustWork = TRUE)
  on.exit(mc19_write_text(capture.output(utils::sessionInfo()), file.path(support, "sessionInfo.txt")), add = TRUE)
  for (path in scripts) if (!file.copy(path, file.path(support, "scripts_used", basename(path))))
    stop("源码快照失败：", path)
  for (name in c("scaling", "environment_TS", "audit", "coverage", "inclusion", "fingerprints"))
    mc19_write_csv(input[[name]], file.path(support, paste0("input_", name, ".csv")))
  mc19_write_csv(do.call(rbind, input$blocks), file.path(support, "analysis_input.csv"))
  mc19_write_csv(script_hashes, file.path(support, "script_fingerprints.csv"))
  manifest <- list(version = cfg$version, cfg = cfg, run_dir = run_dir,
    analysis_scope = cfg$analysis_scope, n_samples = input$scaling$n_samples[1], n_TS = input$scaling$n_TS[1],
    scaling = input$scaling, model_files = character(), fit_complete = FALSE,
    expected_models = cfg$n_models, q_family_size = cfg$q_family_size, input_hashes = input$fingerprints,
    new_fit_count = 0L,
    bootstrap = FALSE, permutations = FALSE, LOTO = FALSE)
  saveRDS(manifest, file.path(support, "wald_manifest.rds"))
  mc19_write_text(c(sprintf("19 MC x %d environments; %d single-indicator TS models; BH n=%d.",length(cfg$environments),cfg$n_models,cfg$q_family_size),
    sprintf("Scope=%s; samples=%d; TS=%d; workflow=%s.",cfg$analysis_scope,manifest$n_samples,manifest$n_TS,mc19_profile()$id),
    "cbind(MC_reads, Other_reads) ~ E_z + (1 | TS); Beta-binomial/logit; no Host/Repeat adjustment or sample averaging.",
    "Each environment is fitted separately, including temperature when selected. This is not the joint GLMM.",
    "Denominator: all original Ruegeria MC rows including Unassigned. Z scores weight each TS once.",
    "Wald P values from glmmTMB; BH across the full planned family, including unavailable rows.",
    "Censored nutrients use L/2; DIN sums the three components. Temperature: observed PM11 surface-water degrees C.",
    "All selected models are newly fitted. No bootstrap, label permutation or LOTO in this GLMM workflow.",
    "Diagnostic simulations check distributional fit and do not supply the coefficient P values."),file.path(run_dir,"00_readme.txt"))
  if(mc19_profile()$temperature) file.copy(mc19_temperature_source(cfg$project_dir),file.path(support,"PM11_temperature_source.csv"))
  mc19_write_text(paste0("fitting_",cfg$n_models,"_TS_models"),file.path(support,"run_status.txt"))
  rows <- list()
  for (environment in cfg$environments) for (mc in cfg$mc_ids) for (model in cfg$model_ids) {
    key <- mc19_wald_key(mc, environment, model)
    d <- input$blocks[[paste(mc, environment, sep = "__")]]
    scale_row <- input$scaling[input$scaling$environment == environment, , drop = FALSE]
    origin <- "new_fit"
    cat("模型 ", length(rows) + 1L, "/",cfg$n_models,"：", key, "\n", sep = ""); flush.console()
    entry <- mc19_wald_fit_one(mc, d, cfg)
    manifest$new_fit_count <- manifest$new_fit_count + 1L
    entry <- mc19_wald_annotate(entry, environment, scale_row, origin, cfg)
    path <- file.path(support, "models", paste0(key, ".rds"))
    saveRDS(entry, path); manifest$model_files[key] <- path
    rows[[key]] <- entry$result
    results <- do.call(rbind, rows); rownames(results) <- NULL
    mc19_wald_write_result(results, run_dir)
    saveRDS(manifest, file.path(support, "wald_manifest.rds"))
  }
  if (nrow(results) != cfg$n_models || !identical(mc19_hash(input$fingerprints$path)$md5, input$fingerprints$md5) ||
      !identical(mc19_hash(scripts)$md5, script_hashes$md5))
    stop("任务数、输入、历史模型或源码变化，本批不发布latest。")
  valid <- is.finite(results$beta_P_value) & results$numerical_ok %in% TRUE
  results$beta_q_value[valid] <- stats::p.adjust(results$beta_P_value[valid], method = "BH", n = cfg$q_family_size)
  results$q_value <- results$beta_q_value
  mc19_wald_write_result(results, run_dir)
  if (manifest$new_fit_count != cfg$n_models) stop("新拟合数量不一致。")
  manifest$model_hashes <- mc19_hash(unname(manifest$model_files))
  mc19_write_csv(manifest$model_hashes, file.path(support, "model_fingerprints.csv"))
  manifest$fit_complete <- TRUE; manifest$finished_at <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
  saveRDS(manifest, file.path(support, "wald_manifest.rds"))
  mc19_write_text("Wald_fits_complete_diagnostics_pending", file.path(support, "run_status.txt"))
  mc19_publish_result(cfg$project_dir, "GLMM", cfg$analysis_scope, run_dir)
  mc19_write_text(c("# Latest GLMM", "", sprintf("Current: 19 MC x %d environments; %d TS; BH n=%d.",length(cfg$environments),cfg$expected_TS,cfg$q_family_size), "",
    paste0("- [主结果](", basename(run_dir), "/00_result.csv)"),
    paste0("- [完整配套](", basename(run_dir), "/supporting_files/)"), "",
    paste(cfg$analysis_scope, "original samples; single environment; no Host/Repeat adjustment or averaging."),
    "历史172226为152项两结构Wald，153110为bootstrap；原结果保留。本批使用TS-only独立目录和指针。"),
    file.path(cfg$results_root, paste0("00_LATEST_GLMM",mc19_profile()$output_tag,"_", cfg$analysis_scope, ".md")))
  invisible(run_dir)
}

mc19_wald_read_run <- function(cfg, run_dir = "") {
  force(cfg) # Config loads shared I/O before R resolves its function names.
  run_dir <- mc19_latest_result(cfg$project_dir, "GLMM", cfg$analysis_scope, run_dir)
  run_dir <- normalizePath(run_dir, winslash = "/", mustWork = TRUE)
  support <- file.path(run_dir, "supporting_files")
  manifest <- mc19_relocate_saved_paths(readRDS(file.path(support, "wald_manifest.rds")))
  results <- utils::read.csv(file.path(support, "all_results.csv"), stringsAsFactors = FALSE)
  primary <- utils::read.csv(file.path(run_dir, "00_result.csv"), stringsAsFactors = FALSE)
  saved_scope <- if (is.null(manifest$analysis_scope)) "coral" else manifest$analysis_scope
  if (!identical(saved_scope, cfg$analysis_scope)) stop("保存结果与所选分析范围不同。")
  expected_n <- if (!is.null(cfg$expected_samples)) cfg$expected_samples else manifest$n_samples
  expected_ts <- cfg$expected_TS
  if (!identical(as.integer(manifest$n_TS), expected_ts)) stop("Result TS count differs from the selected workflow.")
  if (length(expected_n) != 1L || !is.finite(expected_n) || expected_n < 3L) stop("缺少样本数记录。")
  if (is.null(results$analysis_scope)) results$analysis_scope <- saved_scope
  if (anyNA(results$analysis_scope) || any(results$analysis_scope != saved_scope)) stop("结果范围字段不符。")
  plan <- expand.grid(MC_ID = cfg$mc_ids, environment = cfg$environments, model = cfg$model_ids,
    stringsAsFactors = FALSE)
  keys <- mc19_wald_key(results$MC_ID, results$environment, results$model)
  expected <- mc19_wald_key(plan$MC_ID, plan$environment, plan$model)
  if (!identical(manifest$version, cfg$version) || !isTRUE(manifest$fit_complete) ||
      length(manifest$model_files) != cfg$n_models || nrow(results) != cfg$n_models || anyDuplicated(keys) ||
      !setequal(keys, expected) || !setequal(names(manifest$model_files), expected) ||
      any(results$n_samples != expected_n | results$n_TS != expected_ts | results$q_family_size != cfg$q_family_size))
    stop("Incomplete or incompatible single-environment Wald run.")
  primary_keys <- mc19_wald_key(primary$MC_ID, primary$environment, primary$model)
  if (nrow(primary) != cfg$n_models || anyDuplicated(primary_keys) || !setequal(primary_keys, keys))
    stop("主表与完整明细的模型键不一致。")
  primary <- primary[match(keys, primary_keys), , drop = FALSE]
  # The full-precision detail is authoritative; accept display rounding in the main CSV.
  for (field in c("beta_per_TS_SD", "SE", "CI95_low", "CI95_high", "beta_P_value", "beta_q_value")) {
    a <- primary[[field]]; b <- results[[field]]
    same <- (is.na(a) & is.na(b)) | (is.finite(a) & is.finite(b) & abs(a - b) <= 1e-7)
    if (!all(same)) stop("主表与完整明细不同步：", field)
  }
  if (is.null(manifest$model_hashes) ||
      !identical(mc19_hash(manifest$model_hashes$path)$md5, manifest$model_hashes$md5))
    stop("已保存模型发生变化，请先核对来源。")
  list(run_dir = run_dir, support = support, manifest = manifest, results = results)
}

mc19_wald_diagnostics <- function(cfg, run_dir = "") {
  mc19_require(c("glmmTMB", "DHARMa"))
  if (length(cfg$n_sim) != 1L || !is.finite(cfg$n_sim) || cfg$n_sim < 250 ||
      cfg$n_sim != as.integer(cfg$n_sim)) stop("诊断n_sim须为至少250的整数。")
  run <- mc19_wald_read_run(cfg, run_dir)
  destination <- file.path(run$support, "diagnostics")
  dir.create(destination, showWarnings = FALSE)
  rows <- list()
  for (i in seq_len(nrow(run$results))) {
    r <- run$results[i, ]; key <- mc19_wald_key(r$MC_ID, r$environment, r$model)
    entry <- readRDS(run$manifest$model_files[[key]])
    row <- data.frame(MC_ID = r$MC_ID, environment = r$environment, model = r$model,
      status = "skipped_numerical_issue", uniformity_p_diagnostic = NA_real_,
      dispersion_p_diagnostic = NA_real_, zero_inflation_p_diagnostic = NA_real_,
      summary_flags = NA_integer_, n_sim = cfg$n_sim, seed = cfg$seed + i, warnings = "", error = "")
    if (isTRUE(r$numerical_ok) && !is.null(entry$fit)) {
      cat("诊断 ", i, "/",cfg$n_models,"：", key, "（只模拟，不重新拟合）\n", sep = ""); flush.console()
      cached_path <- file.path(destination,paste0(key,".rds"))
      cached <- if(mc19_profile()$incremental && file.exists(cached_path)) tryCatch(readRDS(cached_path),error=function(e)NULL) else NULL
      reusable <- !is.null(cached) && identical(cached$model_md5,unname(tools::md5sum(run$manifest$model_files[[key]]))) &&
        identical(as.integer(cached$n_sim),as.integer(cfg$n_sim))
      if(reusable) {
        message("Reusing saved diagnostic: ",key)
        attempt <- list(value=cached,warnings=character(),error="")
        row$n_sim <- cached$n_sim;row$seed <- cached$seed
      } else attempt <- mc19_capture(function() mc19_diagnostic_one(entry, cfg, cfg$seed + i))
      row$warnings <- paste(attempt$warnings, collapse = " | "); row$error <- attempt$error
      if (!nzchar(attempt$error)) {
        diagnosis <- attempt$value; row$status <- diagnosis$status
        diagnosis$model_md5 <- unname(tools::md5sum(run$manifest$model_files[[key]]))
        row$uniformity_p_diagnostic <- diagnosis$tests$p_diagnostic_only[1]
        row$dispersion_p_diagnostic <- diagnosis$tests$p_diagnostic_only[2]
        row$zero_inflation_p_diagnostic <- diagnosis$tests$p_diagnostic_only[3]
        row$summary_flags <- sum(diagnosis$summaries$outside_simulated95, na.rm = TRUE)
        saveRDS(diagnosis, file.path(destination, paste0(key, ".rds")))
        for (name in c("summaries", "residual_table", "tests", "TS_residuals"))
          mc19_write_csv(diagnosis[[name]], file.path(destination, paste0(key, "__", name, ".csv")))
      } else row$status <- "diagnostic_failed"
    }
    rows[[key]] <- row
    mc19_write_csv(do.call(rbind, rows), file.path(destination, "diagnostic_overview.csv"))
  }
  diagnostics <- do.call(rbind, rows); rownames(diagnostics) <- NULL
  run$results$diagnostic_status <- diagnostics$status
  mc19_wald_write_result(run$results, run$run_dir) # Beta P/q are retained exactly.
  mc19_write_text(c("2000次（或用户指定次数）拟合后模拟，仅检查分布，不重拟合。",
    "诊断P与beta的Wald P分开。review_flags不是自动拒绝规则，不要求所有P>.05。",
    "不据少数极端点、最大值或LOTO翻转淘汰；中心/零值/系统性偏差结合研究目标解释。"),
    file.path(destination, "README.txt"))
  mc19_write_text("Wald_analysis_complete_diagnostics_review_required", file.path(run$support, "run_status.txt"))
  invisible(diagnostics)
}
