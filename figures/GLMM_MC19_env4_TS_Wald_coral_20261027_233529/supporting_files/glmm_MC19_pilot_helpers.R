# MC19 coral pilot: definitions only. Sourcing this file never fits a model.
# Four MCs x nitrate x {BB_GLM, BB_TS}; all original biological samples retained.

mc19_pilot_config <- function(project_dir) {
  list(
    version = "MC19_coral_nitrate_pilot_v1",
    project_dir = normalizePath(project_dir, winslash = "/", mustWork = TRUE),
    workbook = file.path(project_dir, "data", "PERMANOVA_input_MC19_env7_20261026.xlsx"),
    counts_file = file.path(project_dir, "data", "mc", "MC_abundance_absolute.tsv"),
    mc_ids = c("MC10", "MC6", "MC95", "MC376"),
    all_mc_ids = c("MC41", "MC95", "MC6", "MC507", "MC21", "MC451", "MC98",
      "MC2s1", "MC43", "MC25", "MC77", "MC376", "MC35", "MC472", "MC357",
      "MC80", "MC0", "MC10", "MC15"),
    environment = "Nitrate_mgN_L",
    nutrient_fields = c("Nitrate_mgN_L", "Nitrite_mgN_L", "Ammonia_mgN_L", "DIP_mgP_L"),
    model_ids = c("BB_GLM", "BB_TS"),
    analysis_scope = "coral", complete_fields = "Nitrate_mgN_L",
    expected_samples = 98L, expected_TS = 7L,
    n_sim = 2000L, seed = 20261027L,
    results_root = file.path(project_dir, "results"),
    figures_root = file.path(project_dir, "figures"),
    latest_file = file.path(project_dir, "results", "00_latest_GLMM_MC19_pilot.txt"),
    scripts = c("GLMM.R", "GLMM_reads_diagnostics_full.R", "plot_GLMM_MC19.R",
      "glmm_MC19_pilot_helpers.R")
  )
}

mc19_require <- function(packages) {
  absent <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(absent)) stop("缺少R包：", paste(absent, collapse = ", "),
    "。请先安装；本脚本不会自动安装。", call. = FALSE)
}

mc19_write_csv <- function(x, path) {
  utils::write.csv(x, path, row.names = FALSE, fileEncoding = "UTF-8", na = "")
}

mc19_write_text <- function(x, path) writeLines(enc2utf8(x), path, useBytes = TRUE)

mc19_hash <- function(paths) {
  if (any(!file.exists(paths))) stop("文件不存在：", paste(paths[!file.exists(paths)], collapse = ", "))
  data.frame(path = normalizePath(paths, winslash = "/", mustWork = TRUE),
    md5 = unname(tools::md5sum(paths)), stringsAsFactors = FALSE)
}

mc19_same <- function(a, b, tolerance = 1e-12) {
  length(a) == length(b) && all((is.na(a) & is.na(b)) |
    (!is.na(a) & !is.na(b) & abs(a - b) <= tolerance))
}

mc19_capture <- function(fun) {
  notices <- character()
  value <- tryCatch(withCallingHandlers(fun(), warning = function(w) {
    notices <<- c(notices, conditionMessage(w))
    invokeRestart("muffleWarning")
  }), error = function(e) e)
  list(value = value, warnings = unique(notices),
    error = if (inherits(value, "error")) conditionMessage(value) else "")
}

mc19_load_input <- function(cfg) {
  mc19_require("readxl")
  if (!identical(cfg$mc_ids, c("MC10", "MC6", "MC95", "MC376")) ||
      !identical(cfg$environment, "Nitrate_mgN_L") ||
      !identical(cfg$model_ids, c("BB_GLM", "BB_TS"))) {
    stop("本入口只执行4MC×硝酸盐×两种结构的8个先导模型。")
  }
  fingerprints <- mc19_hash(c(cfg$workbook, cfg$counts_file))
  sheets <- lapply(c("Metadata", "Community", "Environment_TS"), function(sheet)
    as.data.frame(readxl::read_excel(cfg$workbook, sheet = sheet)))
  names(sheets) <- c("metadata", "community", "environment")
  meta <- sheets$metadata; community <- sheets$community; env <- sheets$environment
  required <- c("New_ID", "Sample_ID", "Host_type", "Host_niche", "TS", cfg$nutrient_fields)
  if (!all(required %in% names(meta)) ||
      !all(c("New_ID", cfg$all_mc_ids) %in% names(community)) ||
      !all(c("TS", "Sampling_month", cfg$nutrient_fields,
        paste0(cfg$nutrient_fields, "_raw")) %in% names(env))) stop("输入工作表缺少约定字段。")
  for (name in c("metadata", "community")) {
    sheets[[name]]$New_ID <- trimws(as.character(sheets[[name]]$New_ID))
    ids <- sheets[[name]]$New_ID
    if (anyNA(ids) || any(ids == "") || anyDuplicated(ids)) stop(name, "的New_ID缺失或重复。")
  }
  meta <- sheets$metadata; community <- sheets$community
  meta$TS <- as.character(meta$TS); env$TS <- as.character(env$TS)
  if (anyNA(meta$TS) || anyNA(env$TS) || anyDuplicated(env$TS) ||
      !all(meta$TS %in% paste0("TS", 1:8)) ||
      !setequal(meta$New_ID, community$New_ID)) stop("样本ID或TS对应关系不完整。")
  community <- community[match(meta$New_ID, community$New_ID), , drop = FALSE]
  env_row <- match(meta$TS, env$TS)
  if (anyNA(env_row)) stop("Metadata的TS无法匹配Environment_TS。")
  for (field in cfg$nutrient_fields) {
    if (!is.numeric(meta[[field]]) || !is.numeric(env[[field]]) ||
        !mc19_same(meta[[field]], env[[field]][env_row])) stop(field, "的两个工作表数值不一致。")
    raw <- trimws(as.character(env[[paste0(field, "_raw")]]))
    censored <- !is.na(raw) & grepl("^<", raw)
    parsed <- suppressWarnings(as.numeric(sub("^<\\s*", "", raw)))
    parsed[censored] <- parsed[censored] / 2
    if (!mc19_same(parsed, env[[field]])) stop(field, "数值与原始值/L/2规则不一致。")
    env[[paste0(field, "_censored")]] <- censored
  }
  is_coral <- !is.na(meta$Host_niche) & meta$Host_niche == "coral"
  scope <- if (is.null(cfg$analysis_scope)) "coral" else cfg$analysis_scope
  if (length(scope) != 1L || !scope %in% c("coral", "ecosystem")) stop("Invalid analysis scope.")
  niche <- tolower(trimws(as.character(meta$Host_niche)))
  host <- tolower(trimws(as.character(meta$Host_type)))
  biological <- !is.na(niche) & nzchar(niche) & !niche %in% c("env", "seawater", "sediment") &
    !is.na(host) & !host %in% c("seawater", "sediment")
  in_scope <- if (scope == "coral") biological & is_coral else biological
  fields <- if (is.null(cfg$complete_fields)) cfg$environment else cfg$complete_fields
  complete <- apply(as.matrix(meta[fields]), 1L, function(x) all(is.finite(x)))
  keep <- in_scope & complete
  exclusion <- meta[, c("New_ID", "Sample_ID", "Host_type", "Host_niche", "TS"), drop = FALSE]
  exclusion$included <- keep
  exclusion$reason <- ifelse(!in_scope, paste0("outside_", scope, "_scope"),
    ifelse(!complete, "missing_required_environment", "included"))
  selected <- meta[keep, required, drop = FALSE]
  if (nrow(selected) < 3L || (!is.null(cfg$expected_samples) && nrow(selected) != cfg$expected_samples) ||
      length(unique(selected$TS)) != cfg$expected_TS) {
    stop("所选范围的样本数或有效TS与配置不符：", scope)
  }
  if (any(selected[[cfg$environment]] < 0)) stop("硝酸盐不能为负。")
  ts_values <- env[env$TS %in% selected$TS, , drop = FALSE]
  ts_values <- ts_values[order(match(ts_values$TS, paste0("TS", 1:8))), , drop = FALSE]
  # 每个TS只参与一次标准化参数计算；不平均MC或样本数据。
  center <- mean(ts_values[[cfg$environment]])
  scale <- stats::sd(ts_values[[cfg$environment]])
  if (!is.finite(scale) || scale <= 0) stop("硝酸盐在有效TS间没有可估计的变化。")
  ts_values$E_z <- (ts_values[[cfg$environment]] - center) / scale
  selected$TS <- factor(selected$TS, levels = ts_values$TS)
  selected$E_raw <- selected[[cfg$environment]]
  selected$E_z <- ts_values$E_z[match(as.character(selected$TS), ts_values$TS)]
  raw_counts <- utils::read.delim(cfg$counts_file, check.names = FALSE,
    stringsAsFactors = FALSE, quote = "", comment.char = "")
  if (!"MC" %in% names(raw_counts) || anyNA(raw_counts$MC) || anyDuplicated(raw_counts$MC) ||
      !all(cfg$all_mc_ids %in% raw_counts$MC) ||
      !all(selected$New_ID %in% names(raw_counts))) stop("原始计数表的MC或样本ID不完整。")
  count_matrix <- as.matrix(raw_counts[, selected$New_ID, drop = FALSE])
  if (!is.numeric(count_matrix) || anyNA(count_matrix) || any(!is.finite(count_matrix)) ||
      any(count_matrix < 0) || any(count_matrix != floor(count_matrix))) stop("原始reads须为非负整数。")
  totals <- colSums(count_matrix) # 必须加总全部MC；不是仅19MC或4MC。
  if (any(totals <= 0)) stop("Ruegeria总reads必须为正。")
  check_counts <- t(count_matrix[match(cfg$all_mc_ids, raw_counts$MC), , drop = FALSE])
  check_proportions <- sweep(check_counts, 1L, totals, "/")
  workbook_proportions <- as.matrix(community[match(selected$New_ID, community$New_ID),
    cfg$all_mc_ids, drop = FALSE])
  if (!is.numeric(workbook_proportions) || anyNA(workbook_proportions) ||
      any(abs(check_proportions - workbook_proportions) > 1e-10)) {
    stop("原始reads/全部Ruegeria reads与当前19MC比例不一致。")
  }
  selected$Ruegeria_reads_total <- unname(totals)
  blocks <- lapply(cfg$mc_ids, function(mc) {
    d <- selected
    d$MC_ID <- mc
    d$MC_reads <- unname(count_matrix[match(mc, raw_counts$MC), ])
    d$Other_reads <- d$Ruegeria_reads_total - d$MC_reads
    d$MC_proportion <- d$MC_reads / d$Ruegeria_reads_total
    d
  })
  long <- do.call(rbind, blocks); rownames(long) <- NULL
  audit <- do.call(rbind, lapply(blocks, function(d) data.frame(
    MC_ID = d$MC_ID[1], n_samples = nrow(d), n_TS = nlevels(d$TS),
    zero_count = sum(d$MC_reads == 0), positive_count = sum(d$MC_reads > 0),
    positive_TS_count = length(unique(d$TS[d$MC_reads > 0])),
    median_proportion = stats::median(d$MC_proportion), max_proportion = max(d$MC_proportion))))
  coverage <- data.frame(TS = env$TS,
    scope_samples = as.integer(table(factor(meta$TS[in_scope], levels = env$TS))),
    used_samples = as.integer(table(factor(as.character(selected$TS), levels = env$TS))))
  coverage$excluded_missing_environment <- coverage$scope_samples - coverage$used_samples
  if (scope == "coral") { # Historical pilot readers retain their old coverage columns.
    coverage$coral_samples <- coverage$scope_samples
    coverage$excluded_missing_nitrate <- coverage$excluded_missing_environment
  }
  list(data = long, ts_values = ts_values, audit = audit, coverage = coverage,
    exclusion = exclusion, center = center, scale = scale, fingerprints = fingerprints,
    all_MC_count = nrow(raw_counts), proportion_checks = length(check_proportions))
}

mc19_empty_result <- function(mc, model, d) {
  data.frame(MC_ID = mc, environment = "Nitrate_mgN_L", model = model,
    n_samples = nrow(d), n_TS = nlevels(d$TS), zero_count = sum(d$MC_reads == 0),
    fit_status = "fit_failed", convergence_code = NA_integer_, pdHess = FALSE,
    max_abs_gradient = NA_real_, gradient_review_flag = FALSE,
    large_coefficient_flag = FALSE, numerical_ok = FALSE,
    beta_per_TS_SD = NA_real_, SE = NA_real_,
    Wald95_low_pilot = NA_real_, Wald95_high_pilot = NA_real_,
    odds_ratio_per_TS_SD = NA_real_, phi = NA_real_, TS_SD_logit = NA_real_,
    TS_variance_boundary = NA, logLik = NA_real_, AIC = NA_real_,
    fixed_part_change_pp = NA_real_, marginal_change_pp = NA_real_,
    P_value = NA_real_, q_value = NA_real_, inference_status = "pilot_not_calibrated",
    diagnostic_status = "not_run", warnings = "", error = "", elapsed_seconds = NA_real_,
    stringsAsFactors = FALSE)
}

mc19_marginal_probability <- function(eta, sd_ts) {
  # 对随机截距分布作数值积分；与b_TS=0的固定部分预测分别保存。
  if (!is.finite(sd_ts) || sd_ts < 0) return(rep(NA_real_, length(eta)))
  if (sd_ts == 0) return(stats::plogis(eta))
  vapply(eta, function(x) stats::integrate(function(z)
    stats::plogis(x + sd_ts * z) * stats::dnorm(z), -Inf, Inf,
    rel.tol = 1e-7)$value, numeric(1))
}

mc19_extract_result <- function(fit, row, d) {
  row$convergence_code <- fit$fit$convergence
  row$pdHess <- isTRUE(fit$sdr$pdHess)
  gradient <- fit$sdr$gradient.fixed
  if (length(gradient) && all(is.finite(gradient))) row$max_abs_gradient <- max(abs(gradient))
  row$gradient_review_flag <- is.finite(row$max_abs_gradient) && row$max_abs_gradient > .01
  row$logLik <- as.numeric(stats::logLik(fit)); row$AIC <- stats::AIC(fit)
  beta <- glmmTMB::fixef(fit)$cond
  row$large_coefficient_flag <- any(abs(beta) > 10, na.rm = TRUE)
  row$beta_per_TS_SD <- unname(beta["E_z"])
  row$phi <- exp(unname(glmmTMB::fixef(fit)$disp["(Intercept)"]))
  row$TS_SD_logit <- if (row$model == "BB_TS")
    unname(attr(glmmTMB::VarCorr(fit)$cond$TS, "stddev")[1]) else 0
  row$TS_variance_boundary <- row$model == "BB_TS" &&
    is.finite(row$TS_SD_logit) && row$TS_SD_logit < 1e-4
  covariance <- tryCatch(as.matrix(stats::vcov(fit)$cond), error = function(e) NULL)
  if (!is.null(covariance) && "E_z" %in% rownames(covariance) &&
      is.finite(covariance["E_z", "E_z"]) && covariance["E_z", "E_z"] > 0) {
    row$SE <- sqrt(covariance["E_z", "E_z"])
  }
  row$numerical_ok <- isTRUE(row$convergence_code == 0L) && row$pdHess &&
    all(is.finite(c(beta, row$SE, row$logLik, row$phi, row$TS_SD_logit))) &&
    row$SE > 0 && row$phi > 0 && row$TS_SD_logit >= 0
  row$fit_status <- if (row$numerical_ok) "numerically_usable" else "needs_numerical_review"
  if (row$numerical_ok && (row$gradient_review_flag || row$large_coefficient_flag ||
      isTRUE(row$TS_variance_boundary))) row$fit_status <- "numerical_boundary_or_scale_review"
  if (row$numerical_ok) {
    row$Wald95_low_pilot <- row$beta_per_TS_SD - 1.96 * row$SE
    row$Wald95_high_pilot <- row$beta_per_TS_SD + 1.96 * row$SE
    row$odds_ratio_per_TS_SD <- exp(row$beta_per_TS_SD)
    endpoints <- unname(beta["(Intercept)"]) + row$beta_per_TS_SD * range(d$E_z)
    row$fixed_part_change_pp <- 100 * diff(stats::plogis(endpoints))
    row$marginal_change_pp <- 100 * diff(mc19_marginal_probability(endpoints, row$TS_SD_logit))
  }
  row
}

mc19_fit_one <- function(mc, model, d, cfg) {
  row <- mc19_empty_result(mc, model, d)
  formula <- stats::as.formula(if (model == "BB_TS")
    "cbind(MC_reads, Other_reads) ~ E_z + (1 | TS)" else
    "cbind(MC_reads, Other_reads) ~ E_z", env = baseenv())
  begin <- proc.time()[[3]]
  attempt <- mc19_capture(function() glmmTMB::glmmTMB(formula = formula, data = d,
    family = glmmTMB::betabinomial(link = "logit"), ziformula = ~0, dispformula = ~1,
    REML = FALSE, na.action = stats::na.fail,
    control = glmmTMB::glmmTMBControl(parallel = 1L,
      optCtrl = list(iter.max = 10000L, eval.max = 10000L))))
  fit <- if (nzchar(attempt$error)) NULL else attempt$value
  row$warnings <- paste(attempt$warnings, collapse = " | ")
  row$error <- attempt$error
  if (!is.null(fit)) {
    extracted <- mc19_capture(function() mc19_extract_result(fit, row, d))
    row$warnings <- paste(unique(c(attempt$warnings, extracted$warnings)), collapse = " | ")
    if (nzchar(extracted$error)) {
      row$fit_status <- "result_extraction_failed"; row$error <- extracted$error
    } else {
      row <- extracted$value
      row$warnings <- paste(unique(c(attempt$warnings, extracted$warnings)), collapse = " | ")
    }
  }
  row$elapsed_seconds <- proc.time()[[3]] - begin
  list(MC_ID = mc, model = model, data = d, fit = fit, result = row,
    formula = paste(deparse(formula), collapse = " "), version = cfg$version)
}

mc19_write_comparison <- function(results, support) {
  comparisons <- lapply(unique(results$MC_ID), function(mc) {
    a <- results[results$MC_ID == mc & results$model == "BB_GLM", , drop = FALSE]
    b <- results[results$MC_ID == mc & results$model == "BB_TS", , drop = FALSE]
    data.frame(MC_ID = mc, GLM_numerical_ok = a$numerical_ok, TS_numerical_ok = b$numerical_ok,
      beta_GLM = a$beta_per_TS_SD, beta_TS = b$beta_per_TS_SD,
      SE_GLM = a$SE, SE_TS = b$SE, TS_SD_logit = b$TS_SD_logit,
      marginal_change_pp_GLM = a$marginal_change_pp, marginal_change_pp_TS = b$marginal_change_pp,
      delta_AIC_TS_minus_GLM = if (a$numerical_ok && b$numerical_ok) b$AIC - a$AIC else NA_real_,
      interpretation = "Structure comparison only; no model selected by P or AIC")
  })
  mc19_write_csv(do.call(rbind, comparisons), file.path(support, "structure_comparison.csv"))
}

mc19_run_fit <- function(cfg) {
  mc19_require(c("readxl", "glmmTMB"))
  input <- mc19_load_input(cfg)
  script_paths <- file.path(cfg$project_dir, "scripts", "GLMM", cfg$scripts)
  script_hashes <- mc19_hash(script_paths)
  tag <- paste0("GLMM_MC19_nitrate_pilot_", format(Sys.time(), "%Y%m%d_%H%M%S"))
  run_dir <- file.path(cfg$results_root, tag)
  if (dir.exists(run_dir)) stop("运行目录已存在，请稍后重新运行。")
  support <- file.path(run_dir, "supporting_files")
  dir.create(file.path(support, "models"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(support, "scripts_used"), showWarnings = FALSE)
  run_dir <- normalizePath(run_dir, winslash = "/", mustWork = TRUE)
  on.exit(mc19_write_text(capture.output(utils::sessionInfo()),
    file.path(support, "sessionInfo.txt")), add = TRUE)
  for (path in script_paths) {
    if (!file.copy(path, file.path(support, "scripts_used", basename(path)), overwrite = FALSE)) {
      stop("无法保存脚本快照：", path)
    }
  }
  mc19_write_csv(input$data, file.path(support, "analysis_input.csv"))
  mc19_write_csv(input$ts_values, file.path(support, "environment_TS_values.csv"))
  mc19_write_csv(input$audit, file.path(support, "MC_data_audit.csv"))
  mc19_write_csv(input$coverage, file.path(support, "TS_sample_coverage.csv"))
  mc19_write_csv(input$exclusion, file.path(support, "sample_inclusion.csv"))
  mc19_write_csv(input$fingerprints, file.path(support, "input_fingerprints.csv"))
  mc19_write_csv(script_hashes, file.path(support, "script_fingerprints.csv"))
  manifest <- list(version = cfg$version, cfg = cfg, run_dir = run_dir,
    input = input, model_files = character(), fit_complete = FALSE,
    center = input$center, scale = input$scale,
    model_count = 8L, final_inference = FALSE, refitting_bootstrap = FALSE,
    label_permutations = FALSE, leave_one_TS_out = FALSE)
  saveRDS(manifest, file.path(support, "run_manifest.rds"))
  mc19_write_text(c("8个先导模型；P/q尚未校准，留空。",
    "每个MC分别拟合BB_GLM和BB_TS；Beta-binomial/logit，恒定离散度，无零膨胀。",
    "98个coral原始样本；全部Ruegeria reads作分母；不按Repeat_group或TS平均。",
    "E_z按7个TS的当期硝酸盐数值标准化，每个TS等权计算中心和标准差。",
    "Wald95区间仅作先导数值参考，不是已完成TS校准的正式区间。",
    "BB_TS斜率是条件logit斜率；跨结构另比较积分后的边际预测占比变化。",
    "诊断只模拟，不重新拟合；本入口不执行逐TS留一、置换或bootstrap。",
    "先看00_result.csv；模型、输入、诊断、运行状态及快照在supporting_files/。"),
    file.path(run_dir, "00_readme.txt"))
  rows <- list(); files <- character()
  for (mc in cfg$mc_ids) for (model in cfg$model_ids) {
    d <- input$data[input$data$MC_ID == mc, , drop = FALSE]
    cat("拟合 ", length(rows) + 1L, "/8：", mc, " / ", model, "\n", sep = ""); flush.console()
    entry <- mc19_fit_one(mc, model, d, cfg)
    key <- paste(mc, model, sep = "__")
    path <- file.path(support, "models", paste0(key, ".rds"))
    saveRDS(entry, path); files[key] <- path; rows[[key]] <- entry$result
    results <- do.call(rbind, rows); rownames(results) <- NULL
    mc19_write_csv(results, file.path(run_dir, "00_result.csv"))
    manifest$model_files <- files
    saveRDS(manifest, file.path(support, "run_manifest.rds"))
  }
  final_hashes <- mc19_hash(input$fingerprints$path)
  if (!identical(final_hashes$md5, input$fingerprints$md5) ||
      !identical(mc19_hash(script_paths)$md5, script_hashes$md5)) {
    stop("运行中输入或脚本变化；此批次不发布为latest。")
  }
  manifest$fit_complete <- TRUE
  manifest$fit_finished <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
  saveRDS(manifest, file.path(support, "run_manifest.rds"))
  mc19_write_comparison(results, support)
  status <- if (all(results$numerical_ok)) "fit_complete_diagnostics_pending" else
    "fit_complete_with_issues_diagnostics_pending"
  mc19_write_text(status, file.path(support, "run_status.txt"))
  # 专用新指针：不会读取或改写旧MC8机器指针。
  mc19_write_text(run_dir, cfg$latest_file)
  mc19_write_text(c("# Latest GLMM", "", "当前MC19先导运行：", "",
    paste0("- [8模型结果](", basename(run_dir), "/00_result.csv)"), "",
    "P/q尚未完成TS层级校准；先核查数值和诊断。", "",
    "历史MC8结果：GLMM_Z_TS_SD_Q_MC8_20261013_185720/。"),
    file.path(cfg$results_root, "00_LATEST_GLMM.md"))
  cat("8个拟合已结束。结果：", run_dir, "\n")
  invisible(run_dir)
}

mc19_read_run <- function(cfg, run_dir = "") {
  if (!nzchar(run_dir)) {
    if (!file.exists(cfg$latest_file)) stop("尚无MC19先导运行。请先由用户运行GLMM.R。")
    run_dir <- trimws(readLines(cfg$latest_file, warn = FALSE, encoding = "UTF-8")[1])
  }
  if (is.na(run_dir) || !nzchar(run_dir)) stop("MC19运行目录为空。")
  run_dir <- normalizePath(run_dir, winslash = "/", mustWork = TRUE)
  support <- file.path(run_dir, "supporting_files")
  manifest <- readRDS(file.path(support, "run_manifest.rds"))
  if (!identical(manifest$version, cfg$version) || !isTRUE(manifest$fit_complete) ||
      length(manifest$model_files) != 8L) stop("不是完整的当前MC19八模型先导运行。")
  results <- utils::read.csv(file.path(run_dir, "00_result.csv"), stringsAsFactors = FALSE)
  expected <- as.vector(outer(cfg$mc_ids, cfg$model_ids, paste, sep = "__"))
  keys <- paste(results$MC_ID, results$model, sep = "__")
  if (nrow(results) != 8L || anyDuplicated(keys) || !setequal(keys, expected) ||
      !setequal(names(manifest$model_files), expected)) stop("结果表与8模型清单不一致。")
  list(run_dir = run_dir, support = support, manifest = manifest, results = results)
}

mc19_simulated_counts <- function(simulations, d) {
  n <- nrow(d)
  columns <- lapply(simulations, function(x) {
    if (is.matrix(x) || is.data.frame(x)) {
      x <- as.matrix(x)
      if (nrow(x) != n || ncol(x) != 2L ||
          any(abs(rowSums(x) - d$Ruegeria_reads_total) > 1e-8)) stop("模拟两列计数分母不一致。")
      x <- x[, 1]
    }
    if (!is.numeric(x) || length(x) != n || any(!is.finite(x)) ||
        any(x < 0 | x > d$Ruegeria_reads_total) || any(abs(x - round(x)) > 1e-8)) {
      stop("模拟响应不是有效整数计数；未将比例擅自转换为计数。")
    }
    as.numeric(x)
  })
  result <- do.call(cbind, columns)
  if (!is.matrix(result) || nrow(result) != n) stop("模拟矩阵维度不正确。")
  result
}

mc19_stat_vector <- function(y, N) {
  p <- y / N; positive <- p[y > 0]
  c(zero_fraction = mean(y == 0), median_proportion = stats::median(p),
    q95_proportion = unname(stats::quantile(p, .95)), max_proportion = max(p),
    median_positive = if (length(positive)) stats::median(positive) else NA_real_)
}

mc19_simulation_summaries <- function(d, sims) {
  groups <- c(list(All = seq_len(nrow(d))), split(seq_len(nrow(d)), d$TS))
  do.call(rbind, lapply(names(groups), function(group) {
    ix <- groups[[group]]
    observed <- mc19_stat_vector(d$MC_reads[ix], d$Ruegeria_reads_total[ix])
    simulated <- vapply(seq_len(ncol(sims)), function(j)
      mc19_stat_vector(sims[ix, j], d$Ruegeria_reads_total[ix]), numeric(length(observed)))
    do.call(rbind, lapply(seq_along(observed), function(k) {
      values <- simulated[k, ]; values <- values[is.finite(values)]
      bounds <- if (length(values)) stats::quantile(values, c(.025, .5, .975), names = FALSE) else
        rep(NA_real_, 3)
      data.frame(group = group, statistic = names(observed)[k], observed = unname(observed[k]),
        simulated_low = bounds[1], simulated_median = bounds[2], simulated_high = bounds[3],
        outside_simulated95 = if (is.finite(observed[k]) && all(is.finite(bounds)))
          observed[k] < bounds[1] || observed[k] > bounds[3] else NA,
        valid_simulations = length(values))
    }))
  }))
}

mc19_diagnostic_one <- function(entry, cfg, seed) {
  d <- entry$data
  # glmmTMB::simulate会重抽TS随机截距；不是固定已拟合TS的条件模拟。
  simulated <- mc19_simulated_counts(stats::simulate(entry$fit, nsim = cfg$n_sim, seed = seed), d)
  residuals <- DHARMa::createDHARMa(simulatedResponse = simulated,
    observedResponse = d$MC_reads, fittedPredictedResponse = rowMeans(simulated),
    integerResponse = TRUE, seed = seed)
  test_functions <- list(uniformity = function() DHARMa::testUniformity(residuals, plot = FALSE),
    dispersion = function() DHARMa::testDispersion(residuals, plot = FALSE),
    zero_inflation = function() DHARMa::testZeroInflation(residuals, plot = FALSE))
  tests <- do.call(rbind, lapply(names(test_functions), function(name) {
    answer <- mc19_capture(test_functions[[name]])
    data.frame(test = name, p_diagnostic_only = if (!nzchar(answer$error))
      answer$value$p.value else NA_real_, error = answer$error,
      warnings = paste(answer$warnings, collapse = " | "))
  }))
  summaries <- mc19_simulation_summaries(d, simulated)
  fixed_p <- as.numeric(stats::predict(entry$fit, type = "response", re.form = NA))
  fitted_p <- as.numeric(stats::predict(entry$fit, type = "response", re.form = NULL))
  residual_table <- data.frame(New_ID = d$New_ID, TS = as.character(d$TS),
    E_raw = d$E_raw, E_z = d$E_z, observed_proportion = d$MC_proportion,
    fixed_part_proportion = fixed_p, conditional_proportion = fitted_p,
    simulation_mean_proportion = rowMeans(simulated) / d$Ruegeria_reads_total,
    scaled_residual = residuals$scaledResiduals)
  mismatch <- any(summaries$outside_simulated95 %in% TRUE)
  flagged_tests <- any(tests$p_diagnostic_only < .05, na.rm = TRUE)
  ts_residuals <- split(residuals$scaledResiduals, d$TS)
  ts_summary <- data.frame(TS = names(ts_residuals),
    n_samples = lengths(ts_residuals),
    residual_median = vapply(ts_residuals, stats::median, numeric(1)),
    residual_mean = vapply(ts_residuals, mean, numeric(1)))
  list(status = if (mismatch || flagged_tests || any(nzchar(tests$error))) "review_flags" else
      "completed_no_automatic_flags", tests = tests, summaries = summaries,
    residual_table = residual_table, TS_residuals = ts_summary, dharma = residuals,
    simulated_counts = simulated,
    method = "Unconditional glmmTMB simulation; TS random intercept redrawn; no refits",
    seed = seed, n_sim = cfg$n_sim)
}

mc19_run_diagnostics <- function(cfg, run_dir = "") {
  mc19_require(c("glmmTMB", "DHARMa"))
  if (length(cfg$n_sim) != 1L || !is.finite(cfg$n_sim) || cfg$n_sim < 250L ||
      cfg$n_sim != as.integer(cfg$n_sim)) stop("n_sim应为至少250的整数，默认2000。")
  run <- mc19_read_run(cfg, run_dir)
  destination <- file.path(run$support, "diagnostics")
  dir.create(destination, showWarnings = FALSE)
  rows <- list()
  for (i in seq_len(nrow(run$results))) {
    r <- run$results[i, ]; key <- paste(r$MC_ID, r$model, sep = "__")
    entry <- readRDS(run$manifest$model_files[[key]])
    row <- data.frame(MC_ID = r$MC_ID, model = r$model, status = "skipped_numerical_issue",
      uniformity_p_diagnostic = NA_real_, dispersion_p_diagnostic = NA_real_,
      zero_inflation_p_diagnostic = NA_real_, summary_flags = NA_integer_,
      n_sim = cfg$n_sim, seed = cfg$seed + i, warnings = "", error = "")
    if (isTRUE(r$numerical_ok) && !is.null(entry$fit)) {
      cat("模拟诊断 ", i, "/8：", key, "（不重新拟合）\n", sep = ""); flush.console()
      attempt <- mc19_capture(function() mc19_diagnostic_one(entry, cfg, cfg$seed + i))
      row$warnings <- paste(attempt$warnings, collapse = " | "); row$error <- attempt$error
      if (!nzchar(attempt$error)) {
        diagnosis <- attempt$value; row$status <- diagnosis$status
        diagnosis$model_md5 <- unname(tools::md5sum(run$manifest$model_files[[key]]))
        row$uniformity_p_diagnostic <- diagnosis$tests$p_diagnostic_only[1]
        row$dispersion_p_diagnostic <- diagnosis$tests$p_diagnostic_only[2]
        row$zero_inflation_p_diagnostic <- diagnosis$tests$p_diagnostic_only[3]
        row$summary_flags <- sum(diagnosis$summaries$outside_simulated95, na.rm = TRUE)
        saveRDS(diagnosis, file.path(destination, paste0(key, ".rds")))
        mc19_write_csv(diagnosis$summaries, file.path(destination, paste0(key, "__simulation_summary.csv")))
        mc19_write_csv(diagnosis$residual_table, file.path(destination, paste0(key, "__residuals.csv")))
        mc19_write_csv(diagnosis$tests, file.path(destination, paste0(key, "__tests.csv")))
        mc19_write_csv(diagnosis$TS_residuals, file.path(destination, paste0(key, "__TS_residuals.csv")))
      } else row$status <- "diagnostic_failed"
    }
    rows[[key]] <- row
  }
  diagnostics <- do.call(rbind, rows); rownames(diagnostics) <- NULL
  run$results$diagnostic_status <- diagnostics$status[match(
    paste(run$results$MC_ID, run$results$model), paste(diagnostics$MC_ID, diagnostics$model))]
  mc19_write_csv(run$results, file.path(run$run_dir, "00_result.csv"))
  mc19_write_csv(diagnostics, file.path(destination, "diagnostic_overview.csv"))
  mc19_write_text(c("模拟诊断完成；未新增模型拟合，未置换，未bootstrap。",
    "模拟重抽TS随机截距。诊断P值用于失配筛查，不是硝酸盐关联P值。",
    "同TS相关性会影响全局残差检验的校准；诊断P仅作筛查，结合分TS模拟图判断。",
    "多个分组/统计量触发标记时人工查看图，不以所有P>0.05作为通关条件。",
    "本阶段不进行自动模型选择或零膨胀追加。"), file.path(destination, "README.txt"))
  file.copy(file.path(cfg$project_dir, "scripts", "GLMM", "GLMM_reads_diagnostics_full.R"),
    file.path(destination, "diagnostic_entry_used.R"), overwrite = TRUE)
  mc19_write_text(capture.output(utils::sessionInfo()), file.path(destination, "sessionInfo.txt"))
  mc19_write_text(if (all(diagnostics$status == "completed_no_automatic_flags"))
    "pilot_complete_review_required" else "pilot_complete_with_review_flags",
    file.path(run$support, "run_status.txt"))
  invisible(diagnostics)
}

mc19_pdf <- function(path, fun, width = 11, height = 8.5) {
  grDevices::pdf(path, width = width, height = height, family = "sans", onefile = TRUE,
    useDingbats = FALSE)
  on.exit(grDevices::dev.off(), add = TRUE)
  fun()
}

mc19_note_page <- function(title, text) {
  graphics::par(mfrow = c(1, 1), mar = c(3, 3, 4, 2))
  graphics::plot.new(); graphics::title(main = title)
  graphics::text(.03, .92, labels = paste(text, collapse = "\n"), adj = c(0, 1), cex = .9)
}

mc19_predict_grid <- function(entry, manifest) {
  if (!isTRUE(entry$result$numerical_ok)) return(NULL)
  d <- entry$data; grid <- d[rep(1L, 160L), , drop = FALSE]
  grid$E_raw <- seq(min(d$E_raw), max(d$E_raw), length.out = nrow(grid))
  grid$E_z <- (grid$E_raw - manifest$center) / manifest$scale
  prediction <- stats::predict(entry$fit, newdata = grid, type = "link", re.form = NA, se.fit = TRUE)
  eta <- as.numeric(prediction$fit); se <- as.numeric(prediction$se.fit)
  data.frame(MC_ID = entry$MC_ID, model = entry$model, E_raw = grid$E_raw,
    E_z = grid$E_z, fixed_part = stats::plogis(eta),
    Wald95_low_pilot = stats::plogis(eta - 1.96 * se),
    Wald95_high_pilot = stats::plogis(eta + 1.96 * se),
    marginal = mc19_marginal_probability(eta, entry$result$TS_SD_logit))
}

mc19_plot_effects <- function(run, full_data_inference = FALSE) {
  r <- run$results
  good <- r$numerical_ok & is.finite(r$beta_per_TS_SD) & is.finite(r$Wald95_low_pilot) &
    is.finite(r$Wald95_high_pilot)
  if (!any(good)) return(mc19_note_page("No numerically usable fits",
    "Inspect 00_result.csv and supporting_files/models before interpreting effects."))
  graphics::par(mfrow = c(1, 1), mar = c(5, 6, 4, 2))
  y <- rev(seq_len(nrow(r)))
  limits <- range(c(0, r$Wald95_low_pilot[good], r$Wald95_high_pilot[good]))
  if (diff(limits) == 0) limits <- limits + c(-1, 1)
  graphics::plot(NA, xlim = limits, ylim = c(.4, nrow(r) + .6), yaxt = "n",
    xlab = "Nitrate coefficient per SD across 7 TS (logit scale)", ylab = "",
    main = if (full_data_inference) "Full-data effects: BB_GLM and BB_TS" else
      "Pilot effects: BB_GLM and BB_TS")
  graphics::abline(v = 0, lty = 3, col = "grey50")
  graphics::axis(2, at = y, labels = paste(r$MC_ID, r$model), las = 1)
  graphics::segments(r$Wald95_low_pilot[good], y[good], r$Wald95_high_pilot[good], y[good])
  graphics::points(r$beta_per_TS_SD[good], y[good], pch = ifelse(r$model[good] == "BB_TS", 17, 16))
  if (any(!good)) graphics::text(rep(limits[1], sum(!good)), y[!good], "Unavailable fit", adj = 0)
  graphics::mtext(if (full_data_inference)
      "Approximate Wald 95% intervals; bootstrap P/q reported separately; GLMM slopes are conditional." else
      "Approximate Wald 95% intervals; no calibrated P/q; GLMM slopes are conditional.",
    side = 1, line = 3.6, cex = .8)
}

mc19_plot_curves <- function(mc, entries, run, grids) {
  d <- entries[[1]]$data; ts <- levels(d$TS); pchs <- c(1, 2, 0, 5, 6, 16, 17)
  graphics::par(mfrow = c(1, 2), mar = c(5, 4, 4, 1), oma = c(2, 0, 0, 0))
  for (entry in entries) {
    grid <- grids[[paste(mc, entry$model, sep = "__")]]
    upper <- if (is.null(grid)) max(d$MC_proportion) else
      max(c(d$MC_proportion, grid$Wald95_high_pilot), na.rm = TRUE)
    graphics::plot(d$E_raw, 100 * d$MC_proportion, pch = pchs[as.integer(d$TS)],
      xlab = "Nitrate (mg N/L)", ylab = "MC proportion (%)", ylim = c(0, max(1, 105 * upper)),
      main = paste(mc, entry$model))
    if (!is.null(grid)) {
      graphics::polygon(c(grid$E_raw, rev(grid$E_raw)),
        100 * c(grid$Wald95_low_pilot, rev(grid$Wald95_high_pilot)),
        col = grDevices::adjustcolor("grey70", alpha.f = .4), border = NA)
      graphics::lines(grid$E_raw, 100 * grid$fixed_part, lwd = 2)
      graphics::lines(grid$E_raw, 100 * grid$marginal, lwd = 2, lty = 2)
      graphics::points(d$E_raw, 100 * d$MC_proportion, pch = pchs[as.integer(d$TS)])
      if (entry$model == "BB_TS") {
        conditional <- as.numeric(stats::predict(entry$fit, type = "response", re.form = NULL))
        first <- !duplicated(d$TS)
        graphics::points(d$E_raw[first], 100 * conditional[first], pch = 4, cex = 1.5, lwd = 2)
      }
    } else graphics::mtext("Fit requires numerical review", side = 3, cex = .8)
    graphics::legend("topleft", legend = ts, pch = pchs[seq_along(ts)], bty = "n", cex = .7, ncol = 2)
    graphics::legend("topright", legend = c("Fixed part (b_TS = 0)", "Marginal over TS",
      "Wald 95%: fixed part"), lty = c(1, 2, NA), pch = c(NA, NA, 15),
      col = c("black", "black", "grey70"), bty = "n", cex = .65)
  }
  graphics::mtext("All points are original samples. Crosses: fitted TS means. No TS-calibrated intervals yet.",
    outer = TRUE, side = 1, line = .4, cex = .8)
}

mc19_plot_TS_offsets <- function(entries, cfg) {
  graphics::par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
  for (mc in cfg$mc_ids) {
    entry <- entries[[paste(mc, "BB_TS", sep = "__")]]
    if (!isTRUE(entry$result$numerical_ok)) {
      graphics::plot.new(); graphics::title(main = paste(mc, "TS offsets unavailable")); next
    }
    effects <- glmmTMB::ranef(entry$fit)$cond$TS
    graphics::plot(seq_len(nrow(effects)), effects[, "(Intercept)"], xaxt = "n", pch = 16,
      xlab = "Sampling time", ylab = "Fitted TS intercept (logit)", main = mc)
    graphics::axis(1, at = seq_len(nrow(effects)), labels = rownames(effects))
    graphics::abline(h = 0, lty = 3)
  }
}

mc19_plot_distribution <- function(entry, diagnosis) {
  d <- entry$data
  graphics::par(mfrow = c(2, 2), mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0))
  DHARMa::plotQQunif(diagnosis$dharma, testUniformity = FALSE, testOutliers = FALSE,
    testDispersion = FALSE)
  graphics::mtext(paste(entry$MC_ID, entry$model), side = 3, outer = TRUE, line = .4)
  graphics::plot(d$E_raw, diagnosis$dharma$scaledResiduals,
    pch = c(1, 2, 0, 5, 6, 16, 17)[as.integer(d$TS)],
    xlab = "Nitrate (mg N/L)", ylab = "Simulation-scaled residual")
  graphics::abline(h = .5, lty = 3)
  graphics::boxplot(split(diagnosis$dharma$scaledResiduals, d$TS),
    xlab = "TS", ylab = "Simulation-scaled residual", col = "white")
  graphics::abline(h = .5, lty = 3)
  graphics::plot(d$Ruegeria_reads_total, diagnosis$dharma$scaledResiduals, log = "x", pch = 16,
    xlab = "Total Ruegeria reads (log scale)", ylab = "Simulation-scaled residual")
  graphics::abline(h = .5, lty = 3)
}

mc19_plot_simulation_summary <- function(entry, diagnosis) {
  stats_to_plot <- c("zero_fraction", "median_proportion", "q95_proportion", "median_positive")
  graphics::par(mfrow = c(2, 2), mar = c(5, 4, 3, 1), oma = c(0, 0, 0, 0))
  for (statistic in stats_to_plot) {
    s <- diagnosis$summaries[diagnosis$summaries$statistic == statistic &
      diagnosis$summaries$group != "All", , drop = FALSE]
    finite <- is.finite(s$observed) & is.finite(s$simulated_low) & is.finite(s$simulated_high)
    if (!any(finite)) {
      graphics::plot.new(); graphics::title(main = paste(statistic, "unavailable")); next
    }
    limits <- range(c(s$observed[finite], s$simulated_low[finite], s$simulated_high[finite]))
    if (diff(limits) == 0) limits <- limits + c(-.01, .01)
    x <- seq_len(nrow(s))
    graphics::plot(x, s$observed, xaxt = "n", ylim = limits, xlab = "TS", ylab = statistic,
      main = paste(entry$MC_ID, entry$model), pch = 16)
    graphics::axis(1, at = x, labels = s$group, las = 2)
    graphics::segments(x[finite], s$simulated_low[finite], x[finite], s$simulated_high[finite],
      lwd = 4, col = "grey65")
    graphics::points(x[finite], s$simulated_median[finite], pch = 1)
    graphics::points(x[finite], s$observed[finite], pch = 16)
    graphics::legend("topleft", c("Observed", "Simulated median", "Simulated 95%"),
      pch = c(16, 1, NA), lty = c(NA, NA, 1), col = c("black", "black", "grey65"),
      bty = "n", cex = .65)
  }
}

mc19_run_plots <- function(cfg, run_dir = "") {
  mc19_require("glmmTMB")
  run <- mc19_read_run(cfg, run_dir)
  entries <- lapply(run$manifest$model_files, readRDS)
  figure_dir <- file.path(cfg$figures_root, basename(run$run_dir))
  support <- file.path(figure_dir, "supporting_files")
  dir.create(support, recursive = TRUE, showWarnings = FALSE)
  grids <- list(); plot_errors <- list()
  for (key in names(entries)) {
    answer <- mc19_capture(function() mc19_predict_grid(entries[[key]], run$manifest))
    if (!nzchar(answer$error)) grids[key] <- list(answer$value)
    if (nzchar(answer$error) || length(answer$warnings)) plot_errors[[key]] <- data.frame(
      model = key, error = answer$error, warnings = paste(answer$warnings, collapse = " | "))
  }
  available <- grids[!vapply(grids, is.null, logical(1))]
  if (length(available)) mc19_write_csv(do.call(rbind, available), file.path(support, "prediction_curves.csv"))
  # 原页直接写入两份多页PDF，避免大量零碎单图。
  mc19_pdf(file.path(figure_dir, "00_1_effects_and_curves.pdf"), function() {
    mc19_note_page("MC19 nitrate pilot", c("Four MCs x two structures = eight models.",
      "Beta-binomial / logit; all coral samples retained; no sample averaging.",
      "P and q are intentionally blank: TS-level calibration is pending.",
      "Shaded bands are approximate fixed-part Wald intervals, not prediction intervals.",
      "Solid: b_TS=0; dashed: integrated over the fitted TS-intercept distribution.",
      "Conditional GLMM coefficients and marginal GLM coefficients have different interpretations.",
      "No structure is automatically selected by significance or AIC."))
    mc19_plot_effects(run)
    for (mc in cfg$mc_ids) mc19_plot_curves(mc,
      entries[paste(mc, cfg$model_ids, sep = "__")], run, grids)
    mc19_plot_TS_offsets(entries, cfg)
  })
  diagnostic_path <- file.path(run$support, "diagnostics", "diagnostic_overview.csv")
  diagnostic_overview <- if (file.exists(diagnostic_path))
    utils::read.csv(diagnostic_path, stringsAsFactors = FALSE) else NULL
  if (!is.null(diagnostic_overview)) mc19_require("DHARMa")
  mc19_pdf(file.path(figure_dir, "00_2_diagnostics.pdf"), function() {
    mc19_note_page("Pilot diagnostic guide", c(
      "Unconditional simulation: glmmTMB redraws the TS random intercept.",
      "Diagnostic P values describe fit checks, not nitrate associations.",
      "Observed zeros, positive values and tails are compared with model simulations.",
      "Flags require review; passing every residual test is not an acceptance rule.",
      "Skipped or failed fits are kept visible. No refitting is performed here."))
    for (key in names(entries)) {
      path <- file.path(run$support, "diagnostics", paste0(key, ".rds"))
      current <- if (!is.null(diagnostic_overview)) diagnostic_overview[
        paste(diagnostic_overview$MC_ID, diagnostic_overview$model, sep = "__") == key, , drop = FALSE] else
        NULL
      diagnosis <- if (file.exists(path)) readRDS(path) else NULL
      usable_diagnostic <- !is.null(current) && nrow(current) == 1L &&
        current$status %in% c("completed_no_automatic_flags", "review_flags") &&
        !is.null(diagnosis) && identical(diagnosis$n_sim, as.integer(current$n_sim)) &&
        identical(diagnosis$seed, as.integer(current$seed)) &&
        identical(diagnosis$model_md5, unname(tools::md5sum(run$manifest$model_files[[key]])))
      if (usable_diagnostic && isTRUE(entries[[key]]$result$numerical_ok)) {
        mc19_plot_distribution(entries[[key]], diagnosis)
        mc19_plot_simulation_summary(entries[[key]], diagnosis)
      } else mc19_note_page(key, c("Diagnostics unavailable or skipped.",
        "Inspect the result table and diagnostic overview before interpretation."))
    }
  })
  if (length(plot_errors)) mc19_write_csv(do.call(rbind, plot_errors), file.path(support, "plot_notices.csv"))
  file.copy(file.path(run$run_dir, "00_result.csv"), file.path(support, "result_snapshot.csv"), overwrite = TRUE)
  file.copy(file.path(cfg$project_dir, "scripts", "GLMM", "plot_GLMM_MC19.R"),
    file.path(support, "plot_entry_used.R"), overwrite = TRUE)
  mc19_write_text(c(paste("Source results:", run$run_dir),
    "Replotting reads saved fits and diagnostics; it does not fit or simulate models.",
    "The b_TS=0 curve is distinct from the numerically integrated marginal curve."),
    file.path(support, "README.txt"))
  mc19_write_text(capture.output(utils::sessionInfo()), file.path(support, "sessionInfo.txt"))
  mc19_write_text(figure_dir, file.path(cfg$figures_root, "00_latest_GLMM_MC19_pilot.txt"))
  cat("合并图目录：", figure_dir, "\n")
  invisible(figure_dir)
}
