# Definitions only. Full-data, model-based parametric-bootstrap LRT.
# Requires pilot and LOTO helpers for saved-model validation; never calls LOTO fitting.

mc19_inference_version <- function() "MC19_coral_nitrate_full_inference_v1"

mc19_inference_latest <- function(cfg) file.path(cfg$results_root,
  "00_latest_GLMM_MC19_full_inference.txt")

mc19_inference_write_result <- function(result, run_dir) {
  primary <- c("MC_ID", "environment", "model", "n_samples", "n_TS",
    "beta_per_TS_SD", "SE", "CI95_low_approx", "CI95_high_approx",
    "odds_ratio_per_TS_SD", "marginal_change_pp", "TS_SD_logit",
    "beta_P_value", "beta_q_value", "P_interpretation", "q_family_size",
    "bootstrap_B_planned", "bootstrap_valid", "bootstrap_failed",
    "inference_status", "fit_status", "diagnostic_status", "environment_context")
  mc19_write_csv(result, file.path(run_dir, "supporting_files", "inference_details.csv"))
  mc19_write_csv(result[, primary, drop = FALSE], file.path(run_dir, "00_result.csv"))
}

mc19_inference_scripts <- function() c("GLMM.R", "plot_GLMM_MC19_full_data.R",
  "glmm_MC19_inference_helpers.R", "glmm_MC19_pilot_helpers.R", "glmm_MC19_loto_helpers.R")

mc19_inference_fit <- function(d, model, null = FALSE, start = NULL) {
  rhs <- if (null) "1" else "E_z"
  if (model == "BB_TS") rhs <- paste(rhs, "+ (1 | TS)")
  f <- stats::as.formula(paste("cbind(MC_reads, Other_reads) ~", rhs), env = baseenv())
  captured <- mc19_capture(function() glmmTMB::glmmTMB(formula = f, data = d,
    family = glmmTMB::betabinomial(link = "logit"), ziformula = ~0, dispformula = ~1,
    REML = FALSE, na.action = stats::na.fail, start = start,
    control = glmmTMB::glmmTMBControl(parallel = 1L,
      optCtrl = list(iter.max = 10000L, eval.max = 10000L))))
  fit <- if (nzchar(captured$error)) NULL else captured$value
  answer <- list(fit = fit, ok = FALSE, logLik = NA_real_, pdHess = FALSE,
    convergence = NA_integer_, TS_SD = NA_real_, boundary = FALSE,
    warnings = paste(captured$warnings, collapse = " | "), error = captured$error)
  if (is.null(fit)) return(answer)
  inspected <- mc19_capture(function() {
    ll <- as.numeric(stats::logLik(fit))
    sd_ts <- if (model == "BB_TS")
      unname(attr(glmmTMB::VarCorr(fit)$cond$TS, "stddev")[1]) else 0
    parameters <- unlist(glmmTMB::fixef(fit), use.names = FALSE)
    list(logLik = ll, TS_SD = sd_ts, pdHess = isTRUE(fit$sdr$pdHess),
      convergence = fit$fit$convergence, boundary = model == "BB_TS" && sd_ts < 1e-4,
      ok = isTRUE(fit$fit$convergence == 0L) && isTRUE(fit$sdr$pdHess) &&
        all(is.finite(c(ll, parameters, sd_ts))) && sd_ts >= 0)
  })
  answer$warnings <- paste(unique(c(captured$warnings, inspected$warnings)), collapse = " | ")
  if (nzchar(inspected$error)) answer$error <- inspected$error else
    answer[names(inspected$value)] <- inspected$value
  answer
}

mc19_inference_lrt <- function(full_ll, null_ll, tolerance = 1e-6) {
  value <- 2 * (full_ll - null_ll)
  if (length(value) != 1L || !is.finite(value) || value < -tolerance) return(NA_real_)
  max(0, value) # Only roundoff-sized negative values are clamped.
}

mc19_inference_start <- function(null_fit, model) {
  start <- list(beta = c(unname(glmmTMB::fixef(null_fit)$cond["(Intercept)"]), 0),
    betadisp = unname(glmmTMB::fixef(null_fit)$disp))
  if (model == "BB_TS") start$theta <- unname(null_fit$fit$par[
    names(null_fit$fit$par) == "theta"])
  start
}

mc19_inference_one <- function(entry, B, seed, work_dir) {
  # Fixed B, no significance-driven stopping, no replacement of failed draws.
  d <- entry$data; model <- entry$model
  checkpoint <- file.path(work_dir, "checkpoint.rds")
  if (file.exists(checkpoint)) {
    state <- readRDS(checkpoint)
    if (!identical(state$B, B) || !identical(state$seed, seed) ||
        !identical(state$key, paste(entry$MC_ID, model, sep = "__"))) stop("Bootstrap检查点不匹配。")
  } else {
    state <- list(key = paste(entry$MC_ID, model, sep = "__"), B = B, seed = seed,
      null = mc19_inference_fit(d, model, null = TRUE), draws = NULL,
      observed_LRT = NA_real_, status = "null_fit_failed", simulations = NULL,
      simulation_warnings = "", simulation_error = "")
    if (isTRUE(state$null$ok)) {
      state$observed_LRT <- mc19_inference_lrt(
        as.numeric(stats::logLik(entry$fit)), state$null$logLik)
      state$status <- if (is.finite(state$observed_LRT)) "running" else "observed_LRT_invalid"
      if (state$status == "running") {
        simulated <- mc19_capture(function() mc19_simulated_counts(
          stats::simulate(state$null$fit, nsim = B, seed = seed), d))
        state$simulation_warnings <- paste(simulated$warnings, collapse = " | ")
        state$simulation_error <- simulated$error
        if (nzchar(simulated$error)) state$status <- "simulation_failed" else
          state$simulations <- simulated$value
      }
    }
    saveRDS(state, checkpoint)
  }
  if (!state$status %in% c("running", "complete")) return(state)
  n_done <- if (is.null(state$draws)) 0L else nrow(state$draws)
  started <- proc.time()[[3]]
  if (n_done < B) for (b in seq.int(n_done + 1L, B)) {
    simulated_data <- d
    simulated_data$MC_reads <- state$simulations[, b]
    simulated_data$Other_reads <- d$Ruegeria_reads_total - simulated_data$MC_reads
    simulated_data$MC_proportion <- simulated_data$MC_reads / d$Ruegeria_reads_total
    null <- mc19_inference_fit(simulated_data, model, null = TRUE)
    full <- mc19_inference_fit(simulated_data, model, start = if (isTRUE(null$ok))
      mc19_inference_start(null$fit, model) else NULL)
    statistic <- if (isTRUE(null$ok) && isTRUE(full$ok))
      mc19_inference_lrt(full$logLik, null$logLik) else NA_real_
    row <- data.frame(draw = b, LRT = statistic, valid = is.finite(statistic),
      null_ok = null$ok, full_ok = full$ok, null_logLik = null$logLik, full_logLik = full$logLik,
      null_convergence = null$convergence, full_convergence = full$convergence,
      null_pdHess = null$pdHess, full_pdHess = full$pdHess,
      null_TS_SD = null$TS_SD, full_TS_SD = full$TS_SD,
      null_boundary = null$boundary, full_boundary = full$boundary,
      null_warnings = null$warnings, full_warnings = full$warnings,
      null_error = null$error, full_error = full$error, stringsAsFactors = FALSE)
    state$draws <- if (is.null(state$draws)) row else rbind(state$draws, row)
    if (b %% 25L == 0L || b == B) {
      saveRDS(state, checkpoint)
      mc19_write_csv(state$draws, file.path(work_dir, "bootstrap_draws.csv"))
      elapsed <- proc.time()[[3]] - started
      cat(state$key, ": ", b, "/", B, "; 本段耗时 ", round(elapsed, 1),
        " s; 本模型预计剩余 ", round(elapsed / (b - n_done) * (B - b), 1), " s\n", sep = "")
      flush.console()
    }
  }
  state$status <- "complete"
  saveRDS(state, checkpoint)
  state
}

mc19_inference_summary <- function(row, state, B, seed) {
  row$P_value <- row$q_value <- NA_real_
  row$beta_P_value <- row$beta_q_value <- NA_real_
  row$beta_test_H0 <- "beta_Nitrate_mgN_L = 0 (two-sided association, LRT)"
  row$CI95_low_approx <- row$Wald95_low_pilot
  row$CI95_high_approx <- row$Wald95_high_pilot
  row$CI_method <- "Full-data Wald approximation; not inverted bootstrap LRT"
  row$inference_status <- state$status
  row$P_method <- "Full-data model-based parametric-bootstrap LRT, beta_nitrate=0"
  row$P_interpretation <- "unavailable"
  row$LRT <- state$observed_LRT
  row$bootstrap_B_planned <- B; row$bootstrap_seed <- seed
  row$bootstrap_B_finished <- if (is.null(state$draws)) 0L else nrow(state$draws)
  row$bootstrap_valid <- row$bootstrap_failed <- row$bootstrap_exceedances <- 0L
  row$P_numerical_lower <- row$P_numerical_upper <- NA_real_
  row$MC95_lower <- row$MC95_upper <- NA_real_
  row$null_fit_status <- if (isTRUE(state$null$ok)) "numerically_usable" else "needs_review"
  row$null_TS_boundary <- isTRUE(state$null$boundary)
  row$null_warnings <- state$null$warnings; row$null_error <- state$null$error
  row$simulation_warnings <- state$simulation_warnings
  row$simulation_error <- state$simulation_error
  row$q_family <- "This stage: 4 MC x nitrate x 2 structures; BH n=8"
  row$q_family_size <- 8L
  row$main_analysis_scope <- "All 98 coral samples, TS1-TS7 including TS4 and TS6"
  row$environment_context <- "High nitrate concentrated in TS4 (0.098) and TS6 (0.110 mg N/L)"
  row$LOTO_role <- "Supplementary influence analysis; never a model-selection or P-value rule"
  if (state$status == "complete" && nrow(state$draws) == B) {
    valid <- state$draws$valid %in% TRUE
    k <- sum(state$draws$LRT[valid] >= state$observed_LRT - 1e-10)
    failed <- sum(!valid)
    row$bootstrap_valid <- sum(valid); row$bootstrap_failed <- failed
    row$bootstrap_exceedances <- k
    row$P_numerical_lower <- (1 + k) / (B + 1)
    row$P_numerical_upper <- (1 + k + failed) / (B + 1)
    # Failed draws remain in the denominator and all count as exceedances.
    # These bounds concern unknown numerical outcomes, not sampling uncertainty.
    row$P_value <- row$P_numerical_upper
    row$beta_P_value <- row$P_value
    row$P_interpretation <- if (failed) "Conservative upper bound for failed refits" else
      "Parametric-bootstrap P, (1+exceedances)/(1+B)"
    row$inference_status <- if (failed) "complete_with_failed_refits_upper_bound" else
      "complete_model_based_bootstrap"
    # Clopper-Pearson simulation uncertainty; combine worst cases if draws failed.
    row$MC95_lower <- stats::binom.test(k, B)$conf.int[1]
    row$MC95_upper <- stats::binom.test(k + failed, B)$conf.int[2]
  }
  row
}

mc19_inference_run <- function(cfg, source_results_dir, loto_results_dir,
    B = 999L, seed = 20261028L, resume_results_dir = "") {
  mc19_require("glmmTMB")
  if (length(B) != 1L || !is.finite(B) || B < 999 || B != as.integer(B) ||
      length(seed) != 1L || !is.finite(seed) || seed != as.integer(seed) ||
      seed < 1 || seed > .Machine$integer.max - 8000L) stop("B须为至少999的整数，seed须为有效正整数。")
  B <- as.integer(B); seed <- as.integer(seed)
  source <- mc19_read_run(cfg, source_results_dir)
  prepared <- mc19_loto_prepare(cfg, source) # Validation only: no LOTO fits are run.
  loto <- mc19_loto_read_run(cfg, loto_results_dir)
  if (!identical(normalizePath(loto$manifest$source_run_dir, winslash = "/", mustWork = TRUE),
      source$run_dir)) stop("LOTO来源与完整模型来源不同。")
  if (!identical(mc19_hash(loto$manifest$source_fingerprints$path)$md5,
      loto$manifest$source_fingerprints$md5)) stop("LOTO记录的完整模型来源已变化。")
  loto_keys <- paste(loto$summary$MC_ID, loto$summary$model, sep = "__")
  full_keys <- paste(source$results$MC_ID, source$results$model, sep = "__")
  match_loto <- match(full_keys, loto_keys)
  if (anyNA(match_loto) || anyDuplicated(loto_keys) ||
      !mc19_same(source$results$beta_per_TS_SD, loto$summary$full_beta[match_loto]) ||
      !mc19_same(source$results$SE, loto$summary$full_SE[match_loto])) stop("LOTO完整估计与原模型不符。")
  source_paths <- c(unname(source$manifest$model_files),
    file.path(source$run_dir, "00_result.csv"), file.path(source$support, "run_manifest.rds"),
    file.path(loto$run_dir, "00_result.csv"), file.path(loto$support, "loto_manifest.rds"))
  source_hashes <- mc19_hash(source_paths)
  scripts <- file.path(cfg$project_dir, "scripts", mc19_inference_scripts())
  script_hashes <- mc19_hash(scripts)
  # Entry settings may change to point to a resume directory; helper code may not.
  computational_scripts <- scripts[grepl("helpers[.]R$", scripts)]
  helper_hashes <- mc19_hash(computational_scripts)
  resumed <- nzchar(resume_results_dir)
  if (resumed) {
    run_dir <- normalizePath(resume_results_dir, winslash = "/", mustWork = TRUE)
    support <- file.path(run_dir, "supporting_files")
    manifest <- readRDS(file.path(support, "inference_manifest.rds"))
    if (!identical(manifest$version, mc19_inference_version()) ||
        !identical(manifest$B, B) || !identical(manifest$seed, seed) ||
        !identical(manifest$source_run_dir, source$run_dir) ||
        !identical(manifest$source_hashes, source_hashes) ||
        !identical(manifest$helper_hashes, helper_hashes) ||
        !identical(manifest$glmmTMB_version, as.character(utils::packageVersion("glmmTMB"))) ||
        !identical(manifest$R_version, R.version.string)) stop("恢复时输入、辅助源码、R/包版本或参数发生变化。")
    mc19_write_csv(script_hashes, file.path(support,
      paste0("resume_scripts_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")))
  } else {
    tag <- paste0("GLMM_MC19_nitrate_full_data_", format(Sys.time(), "%Y%m%d_%H%M%S"))
    run_dir <- file.path(cfg$results_root, tag)
    if (dir.exists(run_dir)) stop("推断运行目录已存在。")
    support <- file.path(run_dir, "supporting_files")
    dir.create(file.path(support, "bootstrap"), recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(support, "scripts_used"), showWarnings = FALSE)
    run_dir <- normalizePath(run_dir, winslash = "/", mustWork = TRUE)
    for (path in scripts) if (!file.copy(path, file.path(support, "scripts_used", basename(path))))
      stop("无法保存脚本快照：", path)
    manifest <- list(version = mc19_inference_version(), run_dir = run_dir, cfg = cfg,
      B = B, seed = seed, source_run_dir = source$run_dir, loto_run_dir = loto$run_dir,
      source_hashes = source_hashes, helper_hashes = helper_hashes,
      source_model_files = source$manifest$model_files, model_keys = names(prepared$entries),
      glmmTMB_version = as.character(utils::packageVersion("glmmTMB")), R_version = R.version.string,
      completed_keys = character(), analysis_complete = FALSE, model_selection = FALSE,
      n_full_samples = 98L, n_TS = 7L, null_models = 8L,
      planned_bootstrap_refits = 16L * B, q_family_size = 8L,
      reused_full_models = 8L, new_observed_full_models = 0L,
      center = prepared$center, scale = prepared$scale)
    mc19_write_csv(source_hashes, file.path(support, "source_fingerprints.csv"))
    mc19_write_csv(script_hashes, file.path(support, "script_fingerprints.csv"))
    mc19_write_csv(source$results, file.path(support, "source_full_results.csv"))
    mc19_write_csv(loto$summary, file.path(support, "LOTO_supplement_only.csv"))
    mc19_write_csv(do.call(rbind, prepared$data_by_mc), file.path(support, "frozen_full_input.csv"))
    mc19_write_csv(source$manifest$input$ts_values, file.path(support, "environment_TS_values.csv"))
    mc19_write_text(c("完整数据硝酸盐推断：4MC × 两结构，主估计沿用原98样本/7TS模型。",
      paste("原完整结果：", source$run_dir), paste("LOTO补充结果：", loto$run_dir),
      "TS4、TS6全部保留；高硝酸盐信息集中于两次采样。无水温、Host、Repeat或样本平均。",
      "两结构分别比较含硝酸盐的完整模型与去掉硝酸盐的零模型；其余结构不变。",
      paste0("每结构/MC固定", B, "次零模型参数模拟，每次用全部98行重拟合零/完整模型。"),
      "BB_TS模拟重抽7个TS截距；BB_GLM模拟依赖样本独立假设。均固定实际TS、硝酸盐与reads分母。",
      "检验目标是硝酸盐斜率beta=0（双侧关联）；不是检验TS随机项是否存在或诊断残差。",
      "beta_P_value/P_value=(1+LRT不小于实际值的次数)/(B+1)；beta_q_value/q_value为对应BH调整。",
      "模型参数bootstrap依赖所选分布/随机项假设，不能分离未测量季节混杂。",
      "失败模拟不替换、不删除：P_numerical_lower/upper界定未知数值结果，P_value使用保守上界。",
      "MC95_lower/upper为模拟次数引入的二项区间，与效应CI及数值失败上下界不同。",
      "q_value=本轮预定8项P的BH调整；失败模型仍保留在n=8中。不是未来19MC×4营养盐的最终q。",
      "CI95_*_approx仍是完整模型Wald近似区间，不是bootstrap反演区间；不据此反推bootstrap P。",
      "LOTO仅在supporting_files提供补充；不按翻转次数、LOTO显著性、诊断P或AIC筛除模型。",
      "主表00_result.csv优先展示beta/SE/近似CI/幅度/beta_P_value/beta_q_value及状态。",
      "supporting_files/inference_details.csv保留全部精度/失败/数值界限；其余draws/检查点和源码同放配套目录。",
      "中断后设置GLMM.R的resume_results_dir并保持B/seed及helper不变；最多重算最近未保存的24次draw。",
      "模型比较方法参考：https://hojsgaard.github.io/pbkrtest/reference/pb-modcomp.html （不调用不支持glmmTMB的PBmodcomp）。",
      "glmmTMB模拟：https://glmmtmb.github.io/glmmTMB/reference/simulate.glmmTMB.html",
      "BH：https://stat.ethz.ch/R-manual/R-devel/library/stats/html/p.adjust.html"),
      file.path(run_dir, "00_readme.txt"))
    saveRDS(manifest, file.path(support, "inference_manifest.rds"))
  }
  on.exit(mc19_write_text(capture.output(utils::sessionInfo()),
    file.path(support, "sessionInfo.txt")), add = TRUE)
  manifest$analysis_complete <- FALSE
  saveRDS(manifest, file.path(support, "inference_manifest.rds"))
  mc19_write_text("running_full_data_bootstrap", file.path(support, "run_status.txt"))
  cat("完整数据推断目录：", run_dir, "\n计划最多8个零模型 + ", 16L * B,
    " 次bootstrap重拟合；不重跑原8模型或56次LOTO。\n", sep = "")
  rows <- list()
  for (i in seq_along(prepared$entries)) {
    key <- names(prepared$entries)[i]; entry <- prepared$entries[[i]]
    work <- file.path(support, "bootstrap", key)
    dir.create(work, recursive = TRUE, showWarnings = FALSE)
    state <- mc19_inference_one(entry, B, seed + 1000L * i, work)
    baseline <- source$results[paste(source$results$MC_ID, source$results$model, sep = "__") == key,
      , drop = FALSE]
    rows[[key]] <- mc19_inference_summary(baseline, state, B, seed + 1000L * i)
    # Q-values are assigned only after the entire predeclared family has finished.
    result <- do.call(rbind, rows); rownames(result) <- NULL
    mc19_inference_write_result(result, run_dir)
    manifest$completed_keys <- names(rows)
    saveRDS(manifest, file.path(support, "inference_manifest.rds"))
  }
  if (nrow(result) != 8L || !identical(mc19_hash(source_paths)$md5, source_hashes$md5) ||
      !identical(mc19_hash(scripts)$md5, script_hashes$md5)) stop("任务数/原结果/源码在运行中变化，不发布latest。")
  available <- is.finite(result$P_value)
  if (any(available)) result$q_value[available] <- stats::p.adjust(
    result$P_value[available], method = "BH", n = 8L)
  result$beta_q_value <- result$q_value
  mc19_inference_write_result(result, run_dir)
  mc19_write_comparison(result, support)
  manifest$analysis_complete <- TRUE
  manifest$finished_at <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
  saveRDS(manifest, file.path(support, "inference_manifest.rds"))
  mc19_write_text(if (all(available) && all(result$bootstrap_failed == 0L))
    "full_data_inference_complete_model_review_pending" else
    "full_data_inference_complete_with_numerical_issues", file.path(support, "run_status.txt"))
  mc19_write_text(run_dir, mc19_inference_latest(cfg))
  mc19_write_text(c("# Latest GLMM", "", "当前MC19完整数据硝酸盐推断：", "",
    paste0("- [完整数据beta P/q主表](", basename(run_dir), "/00_result.csv)"),
    paste0("- [全部推断明细](", basename(run_dir), "/supporting_files/inference_details.csv)"),
    paste0("- [原8模型及诊断](", basename(source$run_dir), "/00_result.csv)"),
    paste0("- [LOTO补充敏感性](", basename(loto$run_dir), "/00_result.csv)"), "",
    "主分析全98个coral，含TS4/TS6，两结构并列；高硝酸盐信息集中于两TS。",
    "beta P为完整数据零模型参数bootstrap LRT；q为本阶段8项BH。",
    "LOTO只作补充，不筛模型、不决定显著性；数值失败界限及模型诊断需一并查看。"),
    file.path(cfg$results_root, "00_LATEST_GLMM.md"))
  invisible(run_dir)
}

mc19_inference_read <- function(cfg, run_dir = "") {
  if (!nzchar(run_dir)) run_dir <- trimws(readLines(mc19_inference_latest(cfg),
    warn = FALSE, encoding = "UTF-8")[1])
  run_dir <- normalizePath(run_dir, winslash = "/", mustWork = TRUE)
  support <- file.path(run_dir, "supporting_files")
  manifest <- readRDS(file.path(support, "inference_manifest.rds"))
  if (!identical(manifest$version, mc19_inference_version()) ||
      !isTRUE(manifest$analysis_complete)) stop("完整数据推断尚未结束。")
  results <- utils::read.csv(file.path(support, "inference_details.csv"), stringsAsFactors = FALSE)
  primary <- utils::read.csv(file.path(run_dir, "00_result.csv"), stringsAsFactors = FALSE)
  expected <- as.vector(outer(cfg$mc_ids, cfg$model_ids, paste, sep = "__"))
  keys <- paste(results$MC_ID, results$model, sep = "__")
  if (nrow(results) != 8L || anyDuplicated(keys) || !setequal(keys, expected) ||
      any(results$n_samples != 98L | results$n_TS != 7L)) stop("完整结果范围不符。")
  if (nrow(primary) != 8L || !identical(paste(primary$MC_ID, primary$model, sep = "__"), keys) ||
      !mc19_same(primary$beta_per_TS_SD, results$beta_per_TS_SD) ||
      !mc19_same(primary$SE, results$SE) ||
      !mc19_same(primary$CI95_low_approx, results$CI95_low_approx) ||
      !mc19_same(primary$CI95_high_approx, results$CI95_high_approx) ||
      !mc19_same(primary$beta_P_value, results$beta_P_value) ||
      !mc19_same(primary$beta_q_value, results$beta_q_value)) stop("主表与完整明细不同步。")
  list(run_dir = run_dir, support = support, manifest = manifest, results = results)
}

mc19_inference_plots <- function(cfg, run_dir = "") {
  mc19_require("glmmTMB")
  run <- mc19_inference_read(cfg, run_dir)
  source <- mc19_read_run(cfg, run$manifest$source_run_dir)
  if (!identical(mc19_hash(run$manifest$source_hashes$path)$md5,
      run$manifest$source_hashes$md5)) stop("原完整/LOTO结果变化，请先检查来源指纹。")
  entries <- lapply(source$manifest$model_files, readRDS)
  figure_dir <- file.path(cfg$figures_root, basename(run$run_dir))
  support <- file.path(figure_dir, "supporting_files")
  dir.create(support, recursive = TRUE, showWarnings = FALSE)
  grids <- lapply(entries, mc19_predict_grid, manifest = source$manifest)
  plot_run <- source; plot_run$results <- run$results
  mc19_pdf(file.path(figure_dir, "00_full_data_inference.pdf"), function() {
    mc19_note_page("Full-data nitrate analysis: all 98 coral samples", c(
      "Four MCs, both structures; TS1-TS7 including all samples from TS4 and TS6.",
      "High nitrate is concentrated in TS4 (0.098) and TS6 (0.110 mg N/L).",
      "Original full-data effects and approximate Wald intervals are retained.",
      "P: model-based null parametric-bootstrap LRT; q: BH across this stage's eight tests.",
      "A failed bootstrap draw counts as an exceedance; P is then a conservative numerical upper bound.",
      "LOTO describes influence only; it never filters models or tests significance.",
      "No SST/Host/Repeat adjustment, no biological-sample averaging.",
      "Random-intercept assumptions and unmeasured seasonal confounding still affect interpretation."))
    mc19_plot_effects(plot_run, full_data_inference = TRUE)
    for (mc in cfg$mc_ids) mc19_plot_curves(mc,
      entries[paste(mc, cfg$model_ids, sep = "__")], source, grids)
    r <- run$results
    display <- c("MC / structure       beta       SE     beta P   beta BH q  valid/failed draws",
      vapply(seq_len(nrow(r)), function(i) sprintf("%-6s %-7s %+.4f   %.4f   %.4f   %.4f       %d/%d",
        r$MC_ID[i], r$model[i], r$beta_per_TS_SD[i], r$SE[i], r$P_value[i], r$q_value[i],
        r$bootstrap_valid[i], r$bootstrap_failed[i]), character(1)),
      "", "Intervals on preceding pages are approximate Wald intervals, not bootstrap-test inversion.",
      "Monte Carlo intervals and numerical P bounds are retained in 00_result.csv.",
      "Both structures are reported; no minimum-P or LOTO-based structure selection.")
    mc19_note_page("Full-data effects and model-based inference", display)
  })
  mc19_write_csv(do.call(rbind, grids), file.path(support, "prediction_curves.csv"))
  file.copy(file.path(run$run_dir, "00_result.csv"), file.path(support, "result_snapshot.csv"),
    overwrite = TRUE)
  for (script in mc19_inference_scripts()) file.copy(file.path(cfg$project_dir, "scripts", script),
    file.path(support, script), overwrite = TRUE)
  mc19_write_text(c(paste("Source complete-data inference:", run$run_dir),
    paste("Existing diagnostics:", file.path(cfg$figures_root, basename(source$run_dir),
      "00_2_diagnostics.pdf")),
    paste("Existing LOTO supplement:", file.path(cfg$figures_root, basename(run$manifest$loto_run_dir),
      "00_LOTO_stability.pdf")),
    "This plot entry reads saved fits and tables only; no fitting, simulation or LOTO refits."),
    file.path(support, "README.txt"))
  mc19_write_text(figure_dir, file.path(cfg$figures_root, "00_latest_GLMM_MC19_full_inference.txt"))
  cat("完整数据合并图目录：", figure_dir, "\n", sep = "")
  invisible(figure_dir)
}
