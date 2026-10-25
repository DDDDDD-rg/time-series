# 当前入口：已有结果的三个结构精简复核；由用户运行。
# 默认只读取D0、Z_TS_SD_Q、SD_M；不重新拟合、不生成新模拟、不做CV。
# 生成新的效应/校准比较目录；输入、模型及历史结果保持只读。
# 无效拟合列入待核查清单；兼容的上一轮有效结果可作为明确标注的来源。
# 图中使用英语；先查看00_next_steps.html/txt和00_pending_checks.csv。
# 目标8个轻量复核进程，按CPU、任务数和可用内存限制。
# GLMM_PARALLEL_WORKERS可覆盖进程数；MC_selection可限定需要复核的MC。
# 下方保留历史拟合函数及其配置；文件底部只自动调用run_saved_MC_review()。

project_dir <- "C:/Users/xzy/Desktop/CUHK_Study/GNBF6010_Research_Projects/Research/Msc main project-time series"
input_file <- file.path(project_dir, "data", "GLMM_input_MC8_20261011.xlsx")
input_sheet <- "Input_long"
MC_selection <- NULL  # NULL 表示输入表中的全部 MC。
host_reference <- "Acropora sp."
n_sim <- 2000L
simulation_seed <- 123L
quadrature_nodes <- 31L
quadrature_check_nodes <- 61L
maximum_fit_nodes <- 121L
alpha_starts <- c(-6, -2, 2, 6)
reuse_baseline_cache <- TRUE
baseline_cache_dir <- file.path(project_dir, "results", "GLMM_all_MC_20261011_191248")
output_root <- file.path(project_dir, "results")
review_source_dirs <- file.path(output_root,c(
  "GLMM_structure_suite_MC8_20261012_235153",
  "GLMM_structure_suite_MC8_20261012_185348"))
review_models <- c("D0","Z_TS_SD_Q","SD_M")
review_parallel_workers <- 8L # 复核任务较轻；仍按CPU、可用内存和MC数限制。
review_worker_memory_gb <- 0.4
review_min_positive_sim <- 50L
review_save_plots <- TRUE
# 以下并行参数仅用于保留的历史拟合入口；当前复核入口使用review_parallel_workers。
parallel_model_workers <- 12L # GLMM_PARALLEL_WORKERS可覆盖；各进程内部单线程。
parallel_start_workers <- 6L  # 仅供下方历史入口；当前入口不会再建立起点子进程池。
resume_run_dir <- Sys.getenv("GLMM_ZI_RESUME_DIR", unset = "")  # 一般保持空；中断后可恢复。

# 下列为保留的历史D三档惩罚配置；当前结构筛查的设置在suite_version附近。
# Q = 总负对数似然 + lambda * gamma_reads^2；lambda = 1/(2*s^2)。
# x为中心化log10 reads，未按SD缩放；一个单位对应reads增加10倍。
zero_penalty_prior_sds <- c(middle = 2, stronger = 1, weaker = 4)
zero_penalty_reference_sd <- 2
zero_penalty_seed_dir <- file.path(project_dir, "results", "GLMM_joint_reads_MC8_20261012_001527")
zero_penalty_resume_dir <- Sys.getenv("GLMM_ZERO_PENALTY_RESUME_DIR", unset = "")

align_observation_effects <- function(used_data, effects, covariance) {
  n <- nrow(used_data)
  if (!"Observation_ID" %in% names(used_data) || anyNA(used_data$Observation_ID) ||
      length(unique(used_data$Observation_ID)) != n || is.null(effects)) {
    stop("需要每行唯一且完整的 Observation_ID。")
  }
  ri <- match(as.character(used_data$Observation_ID), rownames(effects))
  sd_values <- attr(covariance, "stddev")
  if (is.null(sd_values)) sd_values <- sqrt(diag(covariance))
  if (is.null(names(sd_values))) names(sd_values) <- colnames(covariance)
  if (anyNA(ri) || is.null(names(sd_values)) || any(!is.finite(sd_values)) ||
      any(sd_values < 0)) stop("无法对齐 Observation 随机效应或标准差。")
  if (ncol(effects) == 1L && identical(colnames(effects), "(Intercept)")) {
    term <- rep("(Intercept)", n)
    group_variable <- "common"
    group_level <- rep("All", n)
    structure <- "Common Observation SD"
  } else {
    matching <- list()
    for (column in intersect(c("reads_band", "TS"), names(used_data))) {
      if (anyNA(used_data[[column]])) next
      design <- stats::model.matrix(stats::reformulate(column, intercept = FALSE),
        data = used_data)
      if (ncol(effects) == ncol(design) &&
          setequal(colnames(effects), colnames(design)) &&
          all(rowSums(design != 0) == 1L) && all(design %in% c(0, 1))) {
        matching[[column]] <- design
      }
    }
    if (length(matching) != 1L) {
      stop("Observation 项必须为共同截距或一个分类变量的 diag(0+group|Observation_ID)。")
    }
    group_variable <- names(matching)[1]
    design <- matching[[1]]
    off_diagonal <- covariance
    diag(off_diagonal) <- 0
    if (any(abs(off_diagonal) > 1e-10)) stop("分组 OLRE 模拟要求对角协方差结构。")
    term <- colnames(design)[max.col(design, ties.method = "first")]
    group_level <- as.character(used_data[[group_variable]])
    structure <- paste0(group_variable, "-specific Observation SD")
  }
  ci <- match(term, colnames(effects))
  si <- match(term, names(sd_values))
  if (anyNA(ci) || anyNA(si)) stop("Observation 项的设计列、偏移和标准差名称不一致。")
  offsets <- as.numeric(as.matrix(effects)[cbind(ri, ci)])
  row_sd <- as.numeric(sd_values[si])
  if (any(!is.finite(offsets))) stop("Observation 偏移中有无效值。")
  parameters <- unique(data.frame(group_variable = group_variable,
    group_level = group_level, term = term, observation_sd_logit = row_sd,
    observation_variance_logit = row_sd^2, structure = structure))
  list(offsets = offsets, row_sd = row_sd, term = term,
    parameters = parameters, structure = structure)
}


make_olre_residuals <- function(fit_model, n_sim = 2000L, seed = 123L) {
  if (length(n_sim) != 1L || !is.finite(n_sim) || n_sim < 250L ||
      n_sim != floor(n_sim)) stop("n_sim 应为至少 250 的整数；正式分析建议 2000。")
  if (stats::family(fit_model)$family != "binomial" ||
      stats::family(fit_model)$link != "logit") stop("本模拟用于 binomial(logit) + OLRE。")
  used_data <- stats::model.frame(fit_model)
  counts <- stats::model.response(used_data)
  if (!is.matrix(counts) || ncol(counts) != 2L || any(!is.finite(counts)) ||
      any(counts < 0) || any(counts != floor(counts))) stop("响应须为两列非负整数。")
  total_reads <- rowSums(counts)
  if (any(total_reads <= 0)) stop("Ruegeria 总 reads 必须为正。")
  observation <- align_observation_effects(used_data,
    glmmTMB::ranef(fit_model)$cond$Observation_ID,
    glmmTMB::VarCorr(fit_model)$cond$Observation_ID)

  # 在拟合行预测：移除本样本已拟合的 Observation 项，保留 TS / Repeat。
  eta_all <- as.numeric(stats::predict(fit_model, type = "link", re.form = NULL))
  eta_base <- eta_all - observation$offsets
  if (length(eta_base) != nrow(counts) || any(!is.finite(eta_base))) {
    stop("无法获得有效的样本线性预测值。")
  }
  set.seed(seed)
  simulated_MC <- replicate(n_sim, {
    # 使用每行所属波动组的 SD；共同模型每行相同，reads 分组模型每组不同。
    new_offsets <- stats::rnorm(nrow(counts), mean = 0, sd = observation$row_sd)
    probabilities <- stats::plogis(eta_base + new_offsets)
    stats::rbinom(nrow(counts), size = total_reads, prob = probabilities)
  })
  expected_counts <- rowMeans(simulated_MC)
  res_check <- DHARMa::createDHARMa(simulatedResponse = simulated_MC,
    observedResponse = counts[, 1], fittedPredictedResponse = expected_counts,
    integerResponse = TRUE, seed = seed)
  res_check$fittedModel <- fit_model
  list(res_check = res_check, total_reads = total_reads,
    expected_prop = expected_counts / total_reads,
    observation_sd = observation$row_sd, observation_parameters = observation$parameters,
    observation_structure = observation$structure,
    simulation_method = paste("Conditional on fitted TS/Repeat; redraw independent",
      "Observation effects with each row's model SD, then binomial counts"))
}

# p=0 可以是有限次模拟的结果；负数、>1、NA / Inf 均不能作为有效 p 值。

valid_p <- function(p) is.finite(p) & p >= 0 & p <= 1


adjust_valid_p <- function(p, method = "BH", n = length(p)) {
  answer <- rep(NA_real_, length(p))
  keep <- valid_p(p)
  if (any(keep)) answer[keep] <- stats::p.adjust(p[keep], method = method, n = n)
  answer
}


quantile_fit_state <- function(fit) {
  converged <- !is.null(fit) && isTRUE(fit$converged)
  outer_status <- if (is.null(fit) || is.null(fit$outer.info$conv)) NA_character_
    else paste(fit$outer.info$conv, collapse = ";")
  list(converged = converged, outer_status = outer_status,
    fit_ok = converged && !is.na(outer_status) &&
      identical(tolower(trimws(outer_status)), "full convergence"))
}

# 使用 DHARMa testQuantiles 的 qgam 方法与 BH 合并方式，但先逐项验证。
# 同一 k 规则用于两个模型；失败结果保留为 NA，不自动选择更容易通过的 k。

run_quantile_checks <- function(residuals, predictor, diagnostic,
                                quantiles = c(0.25, 0.5, 0.75)) {
  if (length(residuals) != length(predictor) || any(!is.finite(residuals)) ||
      any(!is.finite(predictor)) || length(unique(predictor)) < 4L) {
    stop("分位数检查需要对齐的有限数据及至少四个不同的自变量值。")
  }
  dimSmooth <- min(length(unique(predictor)), 10L)
  s <- mgcv::s
  fits <- vector("list", length(quantiles))
  rows <- vector("list", length(quantiles))
  curves <- vector("list", length(quantiles))
  messages <- vector("list", length(quantiles))
  for (i in seq_along(quantiles)) {
    tau <- quantiles[i]
    fit_warnings <- character(0)
    dat <- data.frame(res = residuals - tau, pred = predictor)
    fit <- tryCatch(withCallingHandlers({
      invisible(utils::capture.output(qfit <- qgam::qgam(
        res ~ s(pred, k = dimSmooth), data = dat, qu = tau)))
      qfit
    }, warning = function(w) {
      fit_warnings <<- unique(c(fit_warnings, conditionMessage(w)))
      invokeRestart("muffleWarning")
    }), error = function(e) e)
    optimizer_used <- "qgam/newton"
    initial_outer_status <- if (inherits(fit, "error")) "fit_error" else quantile_fit_state(fit)$outer_status
    # qgam 校准结束后，若最后的平滑拟合 step failed，则保留同一 ELF family、
    # 同一学习率、同一 k 和数据，仅用 BFGS 重新完成最终 REML 拟合。
    if (!inherits(fit, "error") && !quantile_fit_state(fit)$fit_ok) {
      retry <- tryCatch(withCallingHandlers(
        mgcv::gam(res ~ s(pred, k = dimSmooth), data = dat, family = fit$family,
          method = "REML", optimizer = c("outer", "bfgs"),
          control = mgcv::gam.control(epsilon = 1e-6, maxit = 500),
          in.out = list(sp = fit$sp * 1.1, scale = fit$sig2)),
        warning = function(w) {
          fit_warnings <<- unique(c(fit_warnings, paste0("REML retry: ", conditionMessage(w))))
          invokeRestart("muffleWarning")
        }), error = function(e) e)
      if (!inherits(retry, "error") && quantile_fit_state(retry)$fit_ok) {
        fit <- retry
        optimizer_used <- "qgam ELF retained; final REML/BFGS retry"
      }
    }
    row <- data.frame(diagnostic = diagnostic, quantile = tau, k = dimSmooth,
      converged = FALSE, outer_status = NA_character_, valid_for_inference = FALSE,
      curve_valid = FALSE, status = "fit_error", reason = "",
      p_intercept_raw = NA_real_, p_trend_raw = NA_real_,
      p_intercept = NA_real_, p_trend = NA_real_,
      p_combined_reference_and_trend = NA_real_, optimizer_used = optimizer_used,
      initial_outer_status = initial_outer_status, warnings = paste(fit_warnings, collapse = " | "))
    if (inherits(fit, "error")) {
      row$reason <- conditionMessage(fit)
    } else {
      fits[i] <- list(fit)
      state <- quantile_fit_state(fit)
      row$converged <- state$converged
      row$outer_status <- state$outer_status
      sm <- tryCatch(base::summary(fit), error = function(e) e)
      if (inherits(sm, "error") || is.null(sm$p.table) || is.null(sm$s.table) ||
          nrow(sm$p.table) < 1L || nrow(sm$s.table) != 1L) {
        row$status <- "summary_error"
        row$reason <- if (inherits(sm, "error")) conditionMessage(sm) else "Missing coefficient tests"
      } else {
        row$p_intercept_raw <- as.numeric(sm$p.table[1, ncol(sm$p.table)])
        row$p_trend_raw <- as.numeric(sm$s.table[1, ncol(sm$s.table)])
        if (!state$fit_ok) {
          row$status <- "not_converged"
          row$reason <- paste("converged:", state$converged, "; outer:", state$outer_status)
        } else if (!all(valid_p(c(row$p_intercept_raw, row$p_trend_raw)))) {
          row$status <- "invalid_p"
          row$reason <- "A component p value is missing, non-finite or outside [0,1]"
        } else {
          row$valid_for_inference <- TRUE
          row$status <- "ok"
          row$p_intercept <- row$p_intercept_raw
          row$p_trend <- row$p_trend_raw
          row$p_combined_reference_and_trend <- min(stats::p.adjust(
            c(row$p_intercept, row$p_trend), method = "BH"))
          pred <- tryCatch(stats::predict(fit, se.fit = TRUE), error = function(e) e)
          if (inherits(pred, "error") || length(pred$fit) != length(predictor) ||
              length(pred$se.fit) != length(predictor) ||
              any(!is.finite(pred$fit)) || any(!is.finite(pred$se.fit)) || any(pred$se.fit < 0)) {
            row$status <- "prediction_error"
            row$reason <- if (inherits(pred, "error")) conditionMessage(pred) else "Invalid quantile predictions"
          } else {
            row$curve_valid <- TRUE
            curves[[i]] <- data.frame(diagnostic = diagnostic, predictor = predictor,
              quantile = factor(tau, levels = quantiles),
              fitted_quantile = as.numeric(pred$fit) + tau,
              lower = as.numeric(pred$fit) + tau - 1.96 * as.numeric(pred$se.fit),
              upper = as.numeric(pred$fit) + tau + 1.96 * as.numeric(pred$se.fit))
          }
        }
      }
    }
    rows[[i]] <- row
    messages[[i]] <- paste(diagnostic, tau, row$status, row$reason, row$warnings, sep = ": ")
  }
  tab <- do.call(rbind, rows)
  curve_table <- do.call(rbind, curves)
  if (is.null(curve_table)) curve_table <- data.frame(diagnostic = character(),
    predictor = numeric(), quantile = factor(levels = quantiles),
    fitted_quantile = numeric(), lower = numeric(), upper = numeric())
  all_valid <- all(tab$valid_for_inference) &&
    all(valid_p(tab$p_combined_reference_and_trend))
  combined <- if (all_valid) min(stats::p.adjust(
    tab$p_combined_reference_and_trend, method = "BH")) else NA_real_
  list(p.value = combined, pvals = tab$p_combined_reference_and_trend,
    qgamFits = fits, quantile_table = tab, curve_table = curve_table,
    diagnostic = diagnostic, status = if (all(tab$status == "ok")) "ok" else "partial",
    messages = unlist(messages[tab$status != "ok"], use.names = FALSE))
}


extract_component_table <- function(coefficient_matrix) {
  if (is.null(coefficient_matrix) || !nrow(coefficient_matrix)) return(NULL)
  cm <- coefficient_matrix
  tab <- data.frame(term = rownames(cm), estimate = cm[, 1],
    standard_error = cm[, 2], z_value = cm[, 3], p_raw = cm[, 4], row.names = NULL)
  tab$CI_low <- tab$estimate - 1.96 * tab$standard_error
  tab$CI_high <- tab$estimate + 1.96 * tab$standard_error
  tab$exp_estimate <- exp(tab$estimate)
  tab$exp_CI_low <- exp(tab$CI_low)
  tab$exp_CI_high <- exp(tab$CI_high)
  tab
}

# Levene / Brown-Forsythe 检验：对各组到组中位数的绝对距离做单因素 F 检验。
# 单独计算统计量，避免 DHARMa 0.5.0 testCategorical(plot=FALSE) 仍调用 title()。

median_levene_test <- function(residuals, grouping) {
  grouping <- droplevels(factor(grouping))
  if (length(residuals) != length(grouping) || anyNA(grouping) ||
      any(!is.finite(residuals)) || nlevels(grouping) < 2L ||
      length(residuals) <= nlevels(grouping)) stop("组间方差检验没有足够的有效数据。")
  centered <- abs(residuals - stats::ave(residuals, grouping, FUN = stats::median))
  tab <- stats::anova(stats::lm(centered ~ grouping))
  p <- as.numeric(tab[1, "Pr(>F)"])
  if (!valid_p(p)) stop("组间方差检验产生无效 p 值。")
  list(p.value = p, statistic = as.numeric(tab[1, "F value"]),
    df_group = as.numeric(tab[1, "Df"]), df_residual = as.numeric(tab[2, "Df"]), table = tab)
}


diagnose_residual_groups <- function(res_check, data_test) {
  groupings <- list(TS = data_test$TS, Host = data_test$Host_type,
                   DO = factor(data_test$DO_c), Reads_band = data_test$reads_band)
  residuals_used <- as.numeric(res_check$scaledResiduals)
  if (length(residuals_used) != nrow(data_test)) stop("分组诊断的残差与样本没有对齐。")
  group_tables <- homogeneity_tables <- objects <- list()
  errors <- character(0)
  for (name in names(groupings)) {
    grouping <- droplevels(factor(groupings[[name]]))
    tab <- do.call(rbind, lapply(levels(grouping), function(g) {
      x <- residuals_used[grouping == g]
      ks <- tryCatch(stats::ks.test(x, "punif", exact = FALSE), error = function(e) e)
      p <- if (inherits(ks, "error")) NA_real_ else ks$p.value
      ok <- valid_p(p)
      reason <- if (inherits(ks, "error")) conditionMessage(ks) else if (!ok) "Invalid KS p value" else ""
      data.frame(grouping = name, group = g, n = length(x),
        median_residual = stats::median(x), residual_IQR = stats::IQR(x),
        p_raw = if (ok) p else NA_real_, status = if (ok) "ok" else "failed", reason = reason)
    }))
    tab$p_BH_within_grouping <- adjust_valid_p(tab$p_raw, n = nrow(tab))
    group_tables[[name]] <- tab
    levene <- tryCatch(median_levene_test(residuals_used, grouping), error = function(e) e)
    failed <- inherits(levene, "error")
    if (failed) errors <- c(errors, paste(name, conditionMessage(levene), sep = ": "))
    else objects[[name]] <- levene
    homogeneity_tables[[name]] <- data.frame(grouping = name,
      method = "Median-centered Levene (Brown-Forsythe)",
      p_value = if (failed) NA_real_ else levene$p.value,
      F_value = if (failed) NA_real_ else levene$statistic,
      df_group = if (failed) NA_real_ else levene$df_group,
      df_residual = if (failed) NA_real_ else levene$df_residual,
      status = if (failed) "failed" else "ok",
      reason = if (failed) conditionMessage(levene) else "")
  }
  group_table <- do.call(rbind, group_tables)
  homogeneity_table <- do.call(rbind, homogeneity_tables)
  list(group_table = group_table, homogeneity_table = homogeneity_table,
    categorical_objects = objects, errors = errors,
    status = if (all(group_table$status == "ok") &&
      all(homogeneity_table$status == "ok")) "ok" else "partial")
}


extract_quantile_fit_status <- function(check, diagnostic = check$diagnostic) {
  check$quantile_table[, c("diagnostic", "quantile", "k", "converged", "outer_status",
    "valid_for_inference", "curve_valid", "status", "reason", "optimizer_used", "initial_outer_status", "warnings"), drop = FALSE]
}


diagnose_reads_relationship <- function(fit_model, res_check, plot = TRUE) {
  needed_packages <- c("glmmTMB", "DHARMa", "qgam", "ggplot2")
  missing_packages <- needed_packages[
    !vapply(needed_packages, requireNamespace, logical(1), quietly = TRUE)
  ]
  if (length(missing_packages)) {
    stop("缺少 R 包：", paste(missing_packages, collapse = ", "))
  }
  if (!(inherits(fit_model, "glmmTMB") || inherits(fit_model, "mc_olre_fit")) ||
      !inherits(res_check, "DHARMa")) {
    stop("需要 glmmTMB / mc_olre_fit 模型及对应 DHARMa 对象。")
  }
  if (stats::family(fit_model)$family != "binomial" ||
      stats::family(fit_model)$link != "logit") {
    stop("本脚本用于 binomial(logit) + OLRE 模型。")
  }

  # 只使用模型实际拟合的行，避免缺测删行后发生错位。
  used_data <- stats::model.frame(fit_model)
  counts <- stats::model.response(used_data)
  if (!is.matrix(counts) || ncol(counts) != 2L ||
      any(!is.finite(counts)) || any(counts < 0) ||
      any(counts != floor(counts))) {
    stop("响应必须是两列非负整数计数：MC_reads 和 Other_reads。")
  }
  total_reads <- rowSums(counts)
  residuals_used <- as.numeric(res_check$scaledResiduals)
  if (nrow(counts) < 40L || any(total_reads <= 0) ||
      length(unique(total_reads)) < 8L) {
    stop("需要至少 40 行、正的总 reads，以及至少 8 个不同的总 reads 值。")
  }
  if (length(residuals_used) != nrow(counts) ||
      any(!is.finite(residuals_used)) ||
      any(residuals_used < 0 | residuals_used > 1) ||
      !isTRUE(all.equal(as.numeric(res_check$observedResponse),
                        as.numeric(counts[, 1]), check.attributes = FALSE))) {
    stop("res_check 与 fit_model 的观测行不匹配，请先重新生成正确的残差对象。")
  }
  if (!is.null(res_check$fittedModel) &&
      (!identical(stats::model.frame(res_check$fittedModel), used_data) ||
       !identical(res_check$fittedModel$fit$par, fit_model$fit$par))) {
    stop("res_check 与当前模型不匹配；修改或重新拟合模型后，需要重新生成残差对象。")
  }
  required_columns <- c("Host_type", "TS")
  if (!all(required_columns %in% names(used_data))) {
    stop("模型数据缺少 Host_type 或 TS 列。")
  }

  log_reads <- log10(total_reads)
  absolute_deviation <- 2 * abs(residuals_used - 0.5)
  quantiles_checked <- c(0.25, 0.5, 0.75)

  location_check <- run_quantile_checks(residuals_used, log_reads, "Residual_location")
  deviation_check <- run_quantile_checks(absolute_deviation, log_reads, "Absolute_deviation")
  quantile_table <- rbind(location_check$quantile_table, deviation_check$quantile_table)
  main_table <- quantile_table[quantile_table$quantile == 0.5,
    c("diagnostic", "p_trend", "valid_for_inference", "status", "reason"), drop = FALSE]
  rownames(main_table) <- NULL
  main_table$p_BH <- adjust_valid_p(main_table$p_trend, n = 2L)

  # 分组仅用于描述和画图，不将连续 reads 改成模型里的类别变量。
  breaks <- unique(as.numeric(stats::quantile(log_reads,
                                              probs = seq(0, 1, 0.25))))
  reads_group <- cut(log_reads, breaks = breaks, include.lowest = TRUE,
                     labels = paste0("Q", seq_len(length(breaks) - 1L)))
  row_table <- data.frame(
    fitted_row = rownames(used_data),
    MC_reads = counts[, 1], Other_reads = counts[, 2],
    Ruegeria_reads_total = total_reads, log10_Ruegeria_reads = log_reads,
    MC_relative_abundance = counts[, 1] / total_reads,
    Host_type = used_data$Host_type, TS = used_data$TS,
    residual = residuals_used, absolute_deviation = absolute_deviation,
    reads_group = reads_group
  )
  group_table <- do.call(rbind, lapply(levels(reads_group), function(g) {
    keep <- reads_group == g
    data.frame(
      reads_group = g, n = sum(keep),
      min_reads = min(total_reads[keep]),
      max_reads = max(total_reads[keep]),
      median_residual = stats::median(residuals_used[keep]),
      residual_IQR = stats::IQR(residuals_used[keep]),
      median_absolute_deviation = stats::median(absolute_deviation[keep])
    )
  }))

  curve_table <- rbind(location_check$curve_table, deviation_check$curve_table)
  names(curve_table)[names(curve_table) == "predictor"] <- "log10_Ruegeria_reads"

  plot_rows <- rbind(
    data.frame(log10_Ruegeria_reads = log_reads,
               value = residuals_used, diagnostic = "Residual_location"),
    data.frame(log10_Ruegeria_reads = log_reads,
               value = absolute_deviation, diagnostic = "Absolute_deviation")
  )
  reference_lines <- expand.grid(
    diagnostic = c("Residual_location", "Absolute_deviation"),
    reference = quantiles_checked
  )
  residual_plot <- ggplot2::ggplot(
    plot_rows, ggplot2::aes(x = log10_Ruegeria_reads, y = value)
  ) +
    ggplot2::geom_point(alpha = 0.35, size = 1.3) +
    ggplot2::geom_hline(data = reference_lines,
                        ggplot2::aes(yintercept = reference),
                        linetype = "dashed", colour = "grey55") +
    ggplot2::geom_ribbon(data = curve_table, inherit.aes = FALSE,
      ggplot2::aes(x = log10_Ruegeria_reads, ymin = lower, ymax = upper,
                   fill = quantile), alpha = 0.12) +
    ggplot2::geom_line(data = curve_table, inherit.aes = FALSE,
      ggplot2::aes(x = log10_Ruegeria_reads, y = fitted_quantile,
                   colour = quantile), linewidth = 0.8) +
    ggplot2::facet_wrap(~ diagnostic, ncol = 1) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::labs(x = "Ruegeria reads per sample (log10)",
                   y = "Residual / 2 * abs(residual - 0.5)",
                   colour = "Quantile", fill = "Quantile") +
    ggplot2::theme_bw()
  raw_plot <- ggplot2::ggplot(row_table,
    ggplot2::aes(x = log10_Ruegeria_reads,
                 y = 100 * MC_relative_abundance, colour = TS)) +
    ggplot2::geom_point(alpha = 0.7, size = 1.5) +
    ggplot2::facet_wrap(~ Host_type, ncol = 4) +
    ggplot2::labs(x = "Ruegeria reads per sample (log10)",
                   y = "MC relative abundance within Ruegeria (%)") +
    ggplot2::theme_bw()

  location_flag <- main_table$p_BH[1] < 0.05
  deviation_flag <- main_table$p_BH[2] < 0.05
  interpretation <- if (anyNA(c(location_flag, deviation_flag))) {
    "部分趋势检验无效或未完成，请查看分位数状态和原因；有效项目单独保留。"
  } else if (location_flag && deviation_flag) {
    "残差中心与绝对偏离均随 reads 改变。先检查中心关系；改善后再检查是否仍有波动大小的趋势。"
  } else if (location_flag) {
    "检测到残差中心随 reads 改变；优先检查预期比例与 reads 的关系。"
  } else if (deviation_flag) {
    "检测到绝对偏离随 reads 改变，支持进一步检查随 reads 变化的变异结构；结合分位数曲线及分组 IQR 判断。"
  } else {
    "两项主要趋势检验经 BH 校正后均未检出显著关系；结合分位数曲线查看局部变化。"
  }

  cat("\n主要结果（p_trend 检验是否随 reads 变化）：\n")
  print(main_table, row.names = FALSE)
  cat("\n分位数细节（p_combined 同时涉及参考水平和趋势）：\n")
  print(quantile_table, row.names = FALSE)
  cat("\n按 reads 四分位分组的描述摘要：\n")
  print(group_table, row.names = FALSE)
  cat("\n", interpretation, "\n", sep = "")
  cat("绝对偏离同时受中心偏移和分布形状影响，不能单独当作随机效应方差的检验。\n")
  if (isTRUE(plot)) {
    print(raw_plot)
    print(residual_plot)
  }
  invisible(list(main_table = main_table, quantile_table = quantile_table,
                 group_table = group_table, row_table = row_table,
                 curve_table = curve_table, interpretation = interpretation,
                 status = if (location_check$status == "ok" && deviation_check$status == "ok") "ok" else "partial",
                  location_check = location_check,
                 deviation_check = deviation_check,
                 plots = list(raw_relationship = raw_plot,
                              residual_relationship = residual_plot),
                 package_versions = vapply(needed_packages,
                   function(p) as.character(utils::packageVersion(p)), character(1))))
}



write_result_table <- function(x, output_dir, filename) {
  utils::write.csv(x, file = file.path(output_dir, filename),
                   row.names = FALSE, fileEncoding = "UTF-8")
}


make_reads_band_setup <- function(total_reads) {
  if (any(!is.finite(total_reads)) || any(total_reads <= 0)) stop("reads 分组需要正的有限分母。")
  limits <- as.numeric(stats::quantile(total_reads, probs = c(0, 1/3, 2/3, 1),
    type = 7L, names = FALSE))
  if (any(diff(limits) <= 0)) stop("reads 分位数分界重复，无法固定三个波动组。")
  labels <- c("Low", "Middle", "High")
  band <- cut(total_reads, breaks = limits, labels = labels,
    include.lowest = TRUE, right = TRUE, ordered_result = FALSE)
  if (anyNA(band) || any(table(band) < 10L)) stop("reads 分组缺失或有组少于 10 个样本。")
  table <- do.call(rbind, lapply(seq_along(labels), function(i) {
    values <- total_reads[band == labels[i]]
    data.frame(reads_band = labels[i], lower_boundary = limits[i],
      upper_boundary = limits[i + 1L], lower_inclusive = i == 1L,
      upper_inclusive = TRUE, n = length(values), min_observed_reads = min(values),
      max_observed_reads = max(values), median_reads = stats::median(values),
      rule = "Type-7 tertiles of model-row Ruegeria total reads; no response-based selection")
  }))
  list(band = band, breaks = limits, labels = labels, table = table)
}


summarise_dispersion_contributions <- function(data_test, expected_prop) {
  residual_count <- data_test$MC_reads - data_test$Ruegeria_reads_total * expected_prop
  squared_centered <- (residual_count - mean(residual_count))^2
  total <- sum(squared_centered)
  if (!is.finite(total) || total <= 0) stop("无法计算计数残差方差贡献。")
  columns <- intersect(c(".input_row", "Sample_ID", "New_ID", "Host_type", "TS",
    "Repeat_group_ID", "Observation_ID", "reads_band", "MC_reads", "Other_reads",
    "Ruegeria_reads_total"), names(data_test))
  samples <- data_test[, columns, drop = FALSE]
  samples$observed_MC_proportion <- data_test$MC_reads / data_test$Ruegeria_reads_total
  samples$expected_MC_proportion <- expected_prop
  samples$count_residual <- residual_count
  samples$squared_centered_count_residual <- squared_centered
  samples$variance_contribution_fraction <- squared_centered / total
  samples <- samples[order(-samples$variance_contribution_fraction), , drop = FALSE]
  samples$cumulative_contribution_fraction <- cumsum(samples$variance_contribution_fraction)
  bands <- do.call(rbind, lapply(levels(data_test$reads_band), function(g) {
    keep <- data_test$reads_band == g
    data.frame(reads_band = g, n = sum(keep),
      count_variance_contribution_fraction = sum(squared_centered[keep]) / total,
      mean_count_residual = mean(residual_count[keep]),
      median_MC_proportion = stats::median(data_test$MC_reads[keep] /
        data_test$Ruegeria_reads_total[keep]))
  }))
  summary <- data.frame(top1_fraction = samples$variance_contribution_fraction[1],
    top5_fraction = sum(head(samples$variance_contribution_fraction, 5L)),
    top10_fraction = sum(head(samples$variance_contribution_fraction, 10L)))
  list(samples = samples, bands = bands, summary = summary)
}

# 对新 reads 网格沿用拟合时的样条节点。

set_reads_basis_columns <- function(data, centered_reads, mean_basis, prefix = "mean_reads_s") {
  values <- stats::predict(mean_basis, newx = centered_reads)
  for (i in seq_len(ncol(values))) data[[paste0(prefix, i)]] <- values[, i]
  data
}


prepare_mc_data <- function(raw_data, MC_test) {
  required_packages <- c("glmmTMB", "DHARMa", "qgam", "ggplot2")
  missing_packages <- required_packages[
    !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
  ]
  if (length(missing_packages)) {
    stop("缺少 R 包：", paste(missing_packages, collapse = ", "))
  }
  required_columns <- c("MC_ID", "MC_reads", "Other_reads", "DO_c",
                         "Host_type", "Repeat_group_ID", "TS")
  if (!all(required_columns %in% names(raw_data))) {
    stop("输入表缺少：", paste(setdiff(required_columns, names(raw_data)),
                                 collapse = ", "))
  }
  if (length(MC_test) != 1L || is.na(MC_test) || !nzchar(MC_test)) {
    stop("MC_test 应为一个 MC 名称，例如 'MC10'。")
  }
  keep_MC <- !is.na(raw_data$MC_ID) & as.character(raw_data$MC_ID) == MC_test
  if (!any(keep_MC)) stop("输入表中没有找到 ", MC_test)
  selected <- as.data.frame(raw_data[keep_MC, , drop = FALSE])
  selected$.input_row <- which(keep_MC) + 1L
  numeric_columns <- c("MC_reads", "Other_reads", "DO_c")
  if ("Ruegeria_reads_total" %in% names(selected)) {
    numeric_columns <- c(numeric_columns, "Ruegeria_reads_total")
  }
  for (column in numeric_columns) {
    if (!is.numeric(selected[[column]])) {
      stop(column, " 必须为数值列；请检查 Excel 是否混入文本。")
    }
    if (any(!is.finite(selected[[column]]) & !is.na(selected[[column]]))) {
      stop(column, " 中存在无穷大等无效数值。")
    }
  }
  for (column in c("MC_reads", "Other_reads")) {
    values <- selected[[column]]
    if (any(values < 0 | values != floor(values), na.rm = TRUE)) {
      stop(column, " 必须是非负整数。零读数会保留。")
    }
  }
  for (column in c("Host_type", "Repeat_group_ID", "TS")) {
    selected[[column]] <- as.character(selected[[column]])
    blank <- !is.na(selected[[column]]) & !nzchar(trimws(selected[[column]]))
    selected[[column]][blank] <- NA_character_
  }
  total_reads <- selected$MC_reads + selected$Other_reads
  if ("Ruegeria_reads_total" %in% names(selected)) {
    comparable <- !is.na(total_reads) & !is.na(selected$Ruegeria_reads_total)
    if (any(total_reads[comparable] != selected$Ruegeria_reads_total[comparable])) {
      stop("MC_reads + Other_reads 与 Ruegeria_reads_total 不一致。")
    }
  }

  # 明确记录删行原因；MC_reads = 0 不构成删行原因。
  reason <- vapply(seq_len(nrow(selected)), function(i) {
    missing <- required_columns[vapply(selected[i, required_columns, drop = FALSE],
                                       function(x) is.na(x[1]), logical(1))]
    parts <- if (length(missing)) paste0("Missing: ", paste(missing, collapse = ";"))
             else character(0)
    if (!is.na(total_reads[i]) && total_reads[i] <= 0) {
      parts <- c(parts, "Ruegeria total reads = 0")
    }
    paste(parts, collapse = " | ")
  }, character(1))
  audit_columns <- intersect(c(".input_row", "New_ID", "Sample_ID", "MC_ID",
                                "Host_type", "TS", "Repeat_group_ID",
                                "MC_reads", "Other_reads", "DO_c"), names(selected))
  sample_audit <- selected[, audit_columns, drop = FALSE]
  sample_audit$included <- reason == ""
  sample_audit$exclusion_reason <- reason
  data_test <- selected[reason == "", , drop = FALSE]
  if (nrow(data_test) < 40L) stop("完整且有正分母的样本少于 40 行。")
  data_test <- data_test[order(data_test$Repeat_group_ID, data_test$.input_row), ,
                           drop = FALSE]
  data_test$Host_type <- factor(data_test$Host_type)
  if ("Acropora sp." %in% levels(data_test$Host_type)) {
    data_test$Host_type <- stats::relevel(data_test$Host_type, ref = "Acropora sp.")
  }
  data_test$Repeat_group_ID <- factor(data_test$Repeat_group_ID)
  data_test$TS <- factor(data_test$TS)
  data_test$Observation_ID <- factor(seq_len(nrow(data_test)))
  rownames(data_test) <- as.character(data_test$.input_row)
  if (any(vapply(data_test[c("Host_type", "Repeat_group_ID", "TS")],
                 nlevels, integer(1)) < 2L)) {
    stop("当前公式需要至少两个宿主、Repeat group 和 TS 水平。")
  }


  # log10 reads 的中心由同一批完整样本计算，两个模型共用。
  data_test$Ruegeria_reads_total <- data_test$MC_reads + data_test$Other_reads
  data_test$log10_Ruegeria_reads <- log10(data_test$Ruegeria_reads_total)
  log_reads_center <- mean(data_test$log10_Ruegeria_reads)
  data_test$log_reads_c <- data_test$log10_Ruegeria_reads - log_reads_center
  if (length(unique(data_test$log_reads_c)) < 8L) {
    stop("候选模型比较需要至少 8 个不同的正 Ruegeria 总 reads 值。")
  }

  # 固定均值曲线 df=3，只比较观察层随机效应的波动分组。
  mean_basis <- splines::ns(data_test$log_reads_c, df = 3L, intercept = FALSE)
  data_test <- set_reads_basis_columns(data_test, data_test$log_reads_c, mean_basis)
  # reads_band 只用于诊断汇总，不进入本轮模型的波动结构。
  band_setup <- make_reads_band_setup(data_test$Ruegeria_reads_total)
  data_test$reads_band <- band_setup$band
  list(data_test = data_test, sample_audit = sample_audit, mean_basis = mean_basis,
    log_reads_center = log_reads_center, n_source = nrow(selected))
}

.normal_quadrature_cache <- new.env(parent = emptyenv())
normal_quadrature <- function(n = 31L) {
  # 每个进程缓存与数据/参数无关的节点，避免各起点重复做同一个特征分解。
  key <- as.character(n)
  if (exists(key, envir = .normal_quadrature_cache, inherits = FALSE))
    return(get(key, envir = .normal_quadrature_cache, inherits = FALSE))
  J <- matrix(0, n, n)
  J[cbind(seq_len(n - 1L), 2:n)] <- sqrt(seq_len(n - 1L))
  J <- J + t(J)
  e <- eigen(J, symmetric = TRUE)
  o <- order(e$values)
  value <- list(nodes = e$values[o], weights = e$vectors[1, o]^2)
  assign(key, value, envir = .normal_quadrature_cache)
  value
}

make_olre_objective <- function(d, X, quadrature, laplace_observation = FALSE, zi = FALSE) {
  if (zi && laplace_observation) stop("ZI comparisons require Observation quadrature.")
  y <- as.numeric(d$MC_reads)
  N <- as.numeric(d$MC_reads + d$Other_reads)
  ir <- as.integer(d$Repeat_group_ID)
  it <- as.integer(d$TS)
  observed_logit <- qlogis((y + 0.5) / (N + 1))
  function(p) {
    sr <- exp(p$log_sd_repeat)
    st <- exp(p$log_sd_ts)
    so <- exp(p$log_sd_observation)
    eta <- as.vector(X %*% p$beta) + p$b_repeat[ir] + p$b_ts[it]
    nll <- -sum(RTMB::dnorm(p$b_repeat, 0, sr, log = TRUE)) -
      sum(RTMB::dnorm(p$b_ts, 0, st, log = TRUE))
    if (laplace_observation) {
      nll <- nll - sum(RTMB::dnorm(p$u_observation, 0, so, log = TRUE)) -
        sum(RTMB::dbinom_robust(y, N, eta + p$u_observation, log = TRUE))
    } else {
      delta <- p$alpha / sqrt(1 + p$alpha^2)
      omega <- so / sqrt(1 - 2 * delta^2 / pi)
      xi <- -omega * delta * sqrt(2 / pi)
      eta0 <- eta + xi
      # Centre at the FULL binomial * skew-normal kernel mode, including zeros.
      zmode <- 10 * tanh(((observed_logit - eta0) / omega) / 10)
      for (iteration in 1:30) {
        probability <- RTMB::plogis(eta0 + omega * zmode)
        az <- p$alpha * zmode
        mills <- exp(RTMB::dnorm(az, log = TRUE) - RTMB::pnorm(az, log.p = TRUE))
        curvature <- 1 + omega^2 * N * probability * (1 - probability) +
          p$alpha^2 * mills * (az + mills)
        step <- (omega * (y - N * probability) - zmode + p$alpha * mills) / curvature
        zmode <- zmode + 3 * tanh(step / 3)
      }
      probability <- RTMB::plogis(eta0 + omega * zmode)
      az <- p$alpha * zmode
      mills <- exp(RTMB::dnorm(az, log = TRUE) - RTMB::pnorm(az, log.p = TRUE))
      zsd <- 1 / sqrt(1 + omega^2 * N * probability * (1 - probability) +
        p$alpha^2 * mills * (az + mills))
      for (j in seq_along(quadrature$nodes)) {
        z <- zmode + zsd * quadrature$nodes[j]
        term <- log(quadrature$weights[j]) +
          RTMB::dbinom_robust(y, N, eta0 + omega * z, log = TRUE) +
          RTMB::dnorm(z, 0, 1, log = TRUE) + log(2) +
          RTMB::pnorm(p$alpha * z, log.p = TRUE) -
          RTMB::dnorm(z, zmode, zsd, log = TRUE)
        if (j == 1L) ll <- term else ll <- RTMB::logspace_add(ll, term)
      }
      # Integrate the count component first, then mix in the point mass at zero.
      # pi is independent of the Observation effect. The count branch still allows zero.
      if (isTRUE(zi)) {
        log_pi <- -RTMB::logspace_add(0, -p$gamma_zi)
        log_one_minus_pi <- -RTMB::logspace_add(0, p$gamma_zi)
        positive <- which(y > 0)
        zero <- which(y == 0)
        ll[positive] <- log_one_minus_pi + ll[positive]
        ll[zero] <- RTMB::logspace_add(log_pi, log_one_minus_pi + ll[zero])
      }
      nll <- nll - sum(ll)
    }
    nll
  }
}

initial_olre_parameters <- function(fit, d, laplace_observation = FALSE, alpha = 0) {
  vc <- glmmTMB::VarCorr(fit)$cond
  re <- glmmTMB::ranef(fit)$cond
  answer <- list(beta = unname(glmmTMB::fixef(fit)$cond),
    log_sd_repeat = log(unname(attr(vc$Repeat_group_ID, "stddev"))),
    log_sd_ts = log(unname(attr(vc$TS, "stddev"))),
    log_sd_observation = log(unname(attr(vc$Observation_ID, "stddev"))),
    alpha = alpha,
    b_repeat = re$Repeat_group_ID[levels(d$Repeat_group_ID), "(Intercept)"],
    b_ts = re$TS[levels(d$TS), "(Intercept)"])
  if (laplace_observation) answer$u_observation <-
    re$Observation_ID[as.character(d$Observation_ID), "(Intercept)"]
  answer
}

fit_rtmb_olre <- function(d, X, initial, quadrature_n = 31L,
                          skew = FALSE, laplace_observation = FALSE, zi = FALSE,
                          control = list(iter.max = 600, eval.max = 1200, rel.tol = 1e-9, x.tol = 1e-8)) {
  quadrature <- normal_quadrature(quadrature_n)
  fun <- make_olre_objective(d, X, quadrature, laplace_observation, zi = zi)
  mapping <- if (!skew) list(alpha = factor(NA)) else list()
  random <- c("b_repeat", "b_ts", if (laplace_observation) "u_observation")
  obj <- RTMB::MakeADFun(fun, parameters = initial, random = random,
    map = mapping, silent = TRUE)
  lower <- rep(-Inf, length(obj$par)); upper <- rep(Inf, length(obj$par))
  names(lower) <- names(upper) <- names(obj$par)
  lower[grepl("^log_sd", names(obj$par))] <- -8
  upper[grepl("^log_sd", names(obj$par))] <- 3
  lower[names(obj$par) == "alpha"] <- -20
  upper[names(obj$par) == "alpha"] <- 20
  if (isTRUE(zi)) {
    lower[names(obj$par) == "gamma_zi"] <- -18
    upper[names(obj$par) == "gamma_zi"] <- stats::qlogis(.99)
  }
  opt <- nlminb(obj$par, obj$fn, obj$gr, lower = lower, upper = upper, control = control)
  initial_convergence <- opt$convergence
  optimizer_path <- "nlminb"
  if (opt$convergence != 0L || max(abs(obj$gr(opt$par))) > .005) {
    retry <- optim(opt$par, obj$fn, obj$gr, method = "L-BFGS-B", lower = lower,
      upper = upper, control = list(maxit = 800, factr = 1e6, pgtol = 1e-6))
    if (is.finite(retry$value) && retry$value <= opt$objective + 1e-6 &&
        retry$convergence == 0L && max(abs(obj$gr(retry$par))) < .01) {
      opt <- list(par = retry$par, objective = retry$value,
        convergence = retry$convergence, message = retry$message)
      optimizer_path <- "nlminb -> L-BFGS-B"
    }
  }
  val <- obj$fn(opt$par)
  gradient <- obj$gr(opt$par)
  estimates <- obj$env$parList(opt$par)
  sdr <- tryCatch(RTMB::sdreport(obj, par.fixed = opt$par), error = function(e) e)
  list(objective = obj, opt = opt, parameters = estimates, value = val,
    max_gradient = max(abs(gradient)), sdr = sdr,
    pdHess = !inherits(sdr, "error") && isTRUE(sdr$pdHess),
    quadrature_n = quadrature_n, laplace_observation = laplace_observation,
    skew = skew, zi = zi, k = length(opt$par), optimizer_path = optimizer_path,
    initial_convergence = initial_convergence)
}

# ---- RTMB 模型对象与结果读取方法 ----
model.frame.mc_olre_fit <- function(formula, ...) formula$frame
family.mc_olre_fit <- function(object, ...) stats::binomial(link = "logit")
logLik.mc_olre_fit <- function(object, ...) structure(-object$nll,
  df = object$k, nobs = nrow(object$data), class = "logLik")
AIC.mc_olre_fit <- function(object, ..., k = 2) {
  if (isTRUE(object$penalty_requested)) return(NA_real_)
  NextMethod("AIC")
}
BIC.mc_olre_fit <- function(object, ...) {
  if (isTRUE(object$penalty_requested)) return(NA_real_)
  NextMethod("BIC")
}
nobs.mc_olre_fit <- function(object, ...) nrow(object$data)
summary.mc_olre_fit <- function(object, ...) {
  beta <- object$parameters$beta
  se <- sqrt(diag(object$cov_fixed))[seq_along(beta)]
  cm <- cbind(Estimate = beta, `Std. Error` = se, `z value` = beta / se,
    `Pr(>|z|)` = 2 * stats::pnorm(-abs(beta / se)))
  if (isTRUE(object$penalty_requested)) cm[, -1L] <- NA_real_
  rownames(cm) <- colnames(object$X)
  answer <- list(model = object$model_id, family = "binomial (logit)",
    formula = object$formula, coefficients = list(cond = cm),
    alpha = object$parameters$alpha,
    alpha_se = if (isTRUE(object$penalty_requested)) NA_real_ else object$alpha_se,
    observation_SD = exp(object$parameters$log_sd_observation),
    Repeat_SD = exp(object$parameters$log_sd_repeat),
    TS_SD = exp(object$parameters$log_sd_ts),
    AIC = stats::AIC(object), BIC = stats::BIC(object),
    convergence = object$fit$convergence, pdHess = object$pdHess,
    max_gradient = object$max_gradient, n = nrow(object$data),
    engine = object$engine, quadrature_n = object$quadrature_n,
    zero_inflation_probability = if (isTRUE(object$zi)) stats::plogis(object$parameters$gamma_zi) else 0,
    penalty_requested = isTRUE(object$penalty_requested), raw_logLik = -object$nll,
    penalty_value = object$penalty_value, objective_value = object$objective_value,
    lambda_zi_reads = object$joint_spec$lambda_zi_reads)
  class(answer) <- "summary.mc_olre_fit"
  answer
}
print.summary.mc_olre_fit <- function(x, ...) {
  cat("Model:", x$model, "\nFamily:", x$family, "\nEngine:", x$engine, "\n")
  print(x$formula)
  print(data.frame(n = x$n, AIC = x$AIC, BIC = x$BIC,
    convergence = x$convergence, pdHess = x$pdHess, max_gradient = x$max_gradient))
  if (isTRUE(x$penalty_requested)) {
    cat("Penalized fit: raw logLik at the penalized estimate =", x$raw_logLik,
      "; penalty =", x$penalty_value, "; total objective =", x$objective_value,
      "; lambda =", x$lambda_zi_reads, "\n")
    cat("Point estimates only; ordinary AIC/BIC, Wald intervals and p values are left blank.\n")
  }
  cat("Conditional fixed effects (logit scale):\n")
  printCoefmat(x$coefficients$cond)
  cat("Random-effect SDs (logit): Repeat =", x$Repeat_SD,
    "; TS =", x$TS_SD, "; Observation =", x$observation_SD, "\n")
  cat("Observation shape alpha =", x$alpha, "; SE =", x$alpha_se, "\n")
  cat("Extra-zero probability pi =", x$zero_inflation_probability, "\n")
  invisible(x)
}

wrap_rtmb_fit <- function(result, reference_fit, d, X, model_id) {
  alpha_index <- which(names(result$opt$par) == "alpha")
  if (inherits(result$sdr, "error")) {
    cov_fixed <- matrix(NA_real_, result$k, result$k)
  } else cov_fixed <- result$sdr$cov.fixed
  answer <- list(model_id = model_id, frame = stats::model.frame(reference_fit),
    data = d, X = X, formula = stats::formula(reference_fit),
    parameters = result$parameters, cov_fixed = cov_fixed,
    alpha_se = if (length(alpha_index)) sqrt(diag(cov_fixed))[alpha_index] else NA_real_,
    fit = list(par = result$opt$par, convergence = result$opt$convergence,
      message = result$opt$message), pdHess = result$pdHess,
    max_gradient = result$max_gradient, nll = as.numeric(result$value), k = result$k,
    engine = "RTMB: adaptive Gaussian quadrature for Observation; Laplace for Repeat/TS",
    quadrature_n = result$quadrature_n, zi = isTRUE(result$zi))
  p <- answer$parameters
  answer$eta_base <- as.vector(X %*% p$beta) +
    p$b_repeat[as.integer(d$Repeat_group_ID)] + p$b_ts[as.integer(d$TS)]
  class(answer) <- "mc_olre_fit"
  answer
}

check_quadrature_precision <- function(result, d, X, check_n) {
  fun <- make_olre_objective(d, X, normal_quadrature(check_n), zi = isTRUE(result$zi))
  obj <- RTMB::MakeADFun(fun, result$parameters,
    random = c("b_repeat", "b_ts"),
    map = if (result$skew) list() else list(alpha = factor(NA)), silent = TRUE)
  val <- obj$fn(result$opt$par)
  gr <- obj$gr(result$opt$par)
  data.frame(fitted_nodes = result$quadrature_n, check_nodes = check_n,
    fitted_nll = as.numeric(result$value), check_nll = as.numeric(val),
    absolute_nll_difference = abs(val - result$value),
    check_max_gradient = max(abs(gr)),
    passed = is.finite(val) && abs(val - result$value) < 0.01 && max(abs(gr)) < 0.02)
}

# 另用 R 的自适应积分检查若干观测的条件似然，不依赖 RTMB 自动微分。
check_conditional_integrals <- function(fit, rows = NULL, nodes = fit$quadrature_n) {
  d <- fit$data; y <- d$MC_reads; N <- y + d$Other_reads
  if (is.null(rows)) rows <- unique(c(order(N)[1:3], order(N, decreasing = TRUE)[1:3],
    order(y / N, decreasing = TRUE)[1:5], order(y / N)[1:3]))
  p <- fit$parameters; alpha <- p$alpha
  delta <- alpha / sqrt(1 + alpha^2)
  omega <- exp(p$log_sd_observation) / sqrt(1 - 2 * delta^2 / pi)
  xi <- -omega * delta * sqrt(2 / pi)
  gh <- normal_quadrature(nodes)
  do.call(rbind, lapply(rows, function(i) {
    eta <- fit$eta_base[i] + xi
    z <- 10 * tanh(((qlogis((y[i] + .5) / (N[i] + 1)) - eta) / omega) / 10)
    for (j in 1:30) {
      pr <- plogis(eta + omega * z)
      az <- alpha * z
      mills <- exp(stats::dnorm(az, log = TRUE) - stats::pnorm(az, log.p = TRUE))
      step <- (omega * (y[i] - N[i] * pr) - z + alpha * mills) /
        (1 + omega^2 * N[i] * pr * (1-pr) + alpha^2 * mills * (az + mills))
      z <- z + 3 * tanh(step / 3)
    }
    pr <- plogis(eta + omega * z)
    az <- alpha * z
    mills <- exp(stats::dnorm(az, log = TRUE) - stats::pnorm(az, log.p = TRUE))
    zs <- 1 / sqrt(1 + omega^2 * N[i] * pr * (1-pr) + alpha^2 * mills * (az + mills))
    log_kernel <- function(a) stats::dbinom(y[i], N[i], plogis(eta + omega * a), log = TRUE) +
      stats::dnorm(a, log = TRUE) + log(2) + stats::pnorm(alpha * a, log.p = TRUE)
    zz <- z + zs * gh$nodes
    terms <- log(gh$weights) + log_kernel(zz) - stats::dnorm(zz, z, zs, log = TRUE)
    m <- max(terms); log_gh <- m + log(sum(exp(terms - m)))
    anchor <- log_kernel(z)
    integral <- stats::integrate(function(a) exp(log_kernel(z + zs * a) - anchor) * zs,
      lower = -Inf, upper = Inf, rel.tol = 1e-9, subdivisions = 300L)
    log_direct <- anchor + log(integral$value)
    count_log_gh <- log_gh; count_log_direct <- log_direct
    if (isTRUE(fit$zi)) {
      log_pi <- stats::plogis(p$gamma_zi, log.p = TRUE)
      log_one_minus_pi <- stats::plogis(p$gamma_zi, lower.tail = FALSE, log.p = TRUE)
      mix_zero <- function(log_count) {
        v <- c(log_pi, log_one_minus_pi + log_count)
        m <- max(v); m + log(sum(exp(v - m)))
      }
      log_gh <- if (y[i] == 0) mix_zero(log_gh) else log_one_minus_pi + log_gh
      log_direct <- if (y[i] == 0) mix_zero(log_direct) else log_one_minus_pi + log_direct
    }
    data.frame(fitted_row = rownames(d)[i], Sample_ID = d$Sample_ID[i], TS = d$TS[i],
      N = N[i], y = y[i], log_likelihood_quadrature = log_gh,
      log_likelihood_direct_integral = log_direct,
      absolute_difference = abs(log_gh - log_direct),
      count_log_likelihood_quadrature = count_log_gh,
      count_log_likelihood_direct_integral = count_log_direct,
      count_absolute_difference = abs(count_log_gh - count_log_direct),
      passed = is.finite(log_direct) && abs(log_gh - log_direct) < 1e-4 &&
        abs(count_log_gh - count_log_direct) < 1e-4)
  }))
}

simulate_centered_olre <- function(fit, n_sim, seed) {
  if (inherits(fit, "glmmTMB")) return(make_olre_residuals(fit, n_sim, seed))
  d <- fit$data; N <- d$MC_reads + d$Other_reads
  alpha <- fit$parameters$alpha
  delta <- alpha / sqrt(1 + alpha^2)
  sd_obs <- exp(fit$parameters$log_sd_observation)
  omega <- sd_obs / sqrt(1 - 2 * delta^2 / pi)
  xi <- -omega * delta * sqrt(2 / pi)
  set.seed(seed)
  simulated <- replicate(n_sim, {
    if (identical(as.numeric(alpha), 0)) u <- stats::rnorm(nrow(d), sd = sd_obs)
    else {
      z0 <- stats::rnorm(nrow(d)); z1 <- stats::rnorm(nrow(d))
      u <- xi + omega * (delta * abs(z0) + sqrt(1 - delta^2) * z1)
    }
    counts <- stats::rbinom(nrow(d), size = N, prob = stats::plogis(fit$eta_base + u))
    if (isTRUE(fit$zi)) {
      extra_zero <- stats::rbinom(nrow(d), size = 1L, prob = stats::plogis(fit$parameters$gamma_zi)) == 1L
      counts[extra_zero] <- 0L
    }
    counts
  })
  expected <- rowMeans(simulated)
  res <- DHARMa::createDHARMa(simulatedResponse = simulated,
    observedResponse = d$MC_reads, fittedPredictedResponse = expected,
    integerResponse = TRUE, seed = seed)
  res$fittedModel <- fit
  list(res_check = res, total_reads = N, expected_prop = expected / N,
    observation_sd = rep(sd_obs, nrow(d)),
    simulation_method = "Conditional on fitted Repeat/TS; redraw centered Observation, binomial counts, and the extra-zero indicator when present")
}

check_distribution_tails <- function(res, N) {
  summaries <- list(
    Q95_MC_proportion = function(x) unname(stats::quantile(x / N, .95)),
    Q99_MC_proportion = function(x) unname(stats::quantile(x / N, .99)),
    Max_MC_proportion = function(x) max(x / N),
    N_MC_proportion_above_25pct = function(x) sum(x / N > .25),
    N_MC_proportion_above_50pct = function(x) sum(x / N > .50),
    N_MC_proportion_below_1pct = function(x) sum(x / N < .01))
  tab <- do.call(rbind, lapply(names(summaries), function(name) {
    fun <- summaries[[name]]
    test <- DHARMa::testGeneric(res, summary = fun, plot = FALSE)
    sim <- apply(res$simulatedResponse, 2L, fun)
    data.frame(check = name, observed = fun(res$observedResponse),
      simulated_mean = mean(sim), simulated_95_low = unname(stats::quantile(sim, .025)),
      simulated_95_high = unname(stats::quantile(sim, .975)), p_raw = test$p.value)
  }))
  tab$p_BH_within_six <- adjust_valid_p(tab$p_raw, n = 6L)
  tab
}

diagnose_skew_comparison_model <- function(id, fit, d, n_sim, seed, output_dir) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  warnings_seen <- character()
  capture_warning <- function(w) {
    warnings_seen <<- unique(c(warnings_seen, conditionMessage(w)))
    invokeRestart("muffleWarning")
  }
  answer <- withCallingHandlers({
    simulation <- simulate_centered_olre(fit, n_sim, seed)
    res <- simulation$res_check; N <- simulation$total_reads
    tests <- list(
      Uniformity = DHARMa::testUniformity(res, plot = FALSE),
      Dispersion = DHARMa::testDispersion(res, plot = FALSE),
      Mean_relative_abundance = DHARMa::testGeneric(res,
        summary = function(x) mean(x / N), plot = FALSE),
      Q90_relative_abundance = DHARMa::testGeneric(res,
        summary = function(x) unname(stats::quantile(x / N, .9)), plot = FALSE))
    set.seed(seed)
    tests$Outliers <- DHARMa::testOutliers(res, type = "bootstrap", nBoot = 200, plot = FALSE)
    reads <- diagnose_reads_relationship(fit, res, plot = FALSE)
    pred <- run_quantile_checks(as.numeric(res$scaledResiduals),
      simulation$expected_prop, "Predicted_proportion")
    grouped <- diagnose_residual_groups(res, d)
    overall <- do.call(rbind, lapply(names(tests), function(name) data.frame(
      check = name, p_value = tests[[name]]$p.value, status = "ok")))
    overall <- rbind(overall,
      data.frame(check = "Quantiles_vs_predicted_proportion", p_value = pred$p.value, status = pred$status),
      data.frame(check = "Quantiles_vs_Ruegeria_reads", p_value = reads$location_check$p.value,
        status = reads$location_check$status))
    overall <- overall[match(c("Uniformity", "Dispersion", "Outliers",
      "Quantiles_vs_predicted_proportion", "Quantiles_vs_Ruegeria_reads",
      "Mean_relative_abundance", "Q90_relative_abundance"), overall$check), ]
    status <- rbind(extract_quantile_fit_status(reads$location_check),
      extract_quantile_fit_status(reads$deviation_check), extract_quantile_fit_status(pred))
    tails <- check_distribution_tails(res, N)
    rows <- cbind(d, expected_MC_proportion = simulation$expected_prop,
      DHARMa_residual = res$scaledResiduals,
      absolute_deviation = 2 * abs(res$scaledResiduals - .5))
    contributions <- summarise_dispersion_contributions(d, simulation$expected_prop)
    write_result_table(overall, output_dir, "04_overall_diagnostics.csv")
    write_result_table(tails, output_dir, "05_distribution_tail_checks.csv")
    write_result_table(reads$main_table, output_dir, "06_reads_trend_tests_BH.csv")
    write_result_table(reads$quantile_table, output_dir, "07_reads_quantile_details.csv")
    write_result_table(reads$group_table, output_dir, "08_reads_group_summary.csv")
    write_result_table(rows, output_dir, "09_fitted_rows_and_residuals.csv")
    write_result_table(grouped$group_table, output_dir, "11_group_residual_tests.csv")
    write_result_table(grouped$homogeneity_table, output_dir, "11_group_homogeneity.csv")
    write_result_table(status, output_dir, "12_quantile_fit_status.csv")
    write_result_table(rbind(reads$quantile_table, pred$quantile_table), output_dir, "12_all_quantile_details.csv")
    write_result_table(contributions$samples, output_dir, "17_dispersion_contributions_by_sample.csv")
    writeLines(enc2utf8(capture.output(base::summary(fit))),
      file.path(output_dir, "01_model_summary.txt"))
    cm <- base::summary(fit)$coefficients$cond
    write_result_table(extract_component_table(cm), output_dir, "02_conditional_coefficients.csv")
    residual_plot <- reads$plots$residual_relationship + ggplot2::labs(title = id)
    ggplot2::ggsave(file.path(output_dir, "10_residual_location_and_deviation.png"),
      residual_plot, width = 8, height = 9, dpi = 300, bg = "white")
    ggplot2::ggsave(file.path(output_dir, "10_residual_location_and_deviation.pdf"),
      residual_plot, width = 8, height = 9)
    # 独立绘图设备，避免 RStudio 的 invalid graphics state 影响保存。
    grDevices::png(file.path(output_dir, "10_DHARMa_overview.png"), width = 1800, height = 900, res = 180)
    tryCatch(plot(res), finally = grDevices::dev.off())
    complete <- all(overall$status == "ok") && all(status$status == "ok") && grouped$status == "ok"
    list(fit_model = fit, res_check = res, overall_table = overall,
      tail_table = tails, reads_diagnostics = reads, grouped_diagnostics = grouped,
      quantile_fit_status = status, fitted_rows = rows, prediction_quantiles = pred,
      dispersion_ratio = unname(tests$Dispersion$statistic),
      diagnostics_complete = complete, diagnostic_status = if (complete) "ok" else "partial",
      simulation_method = simulation$simulation_method, plots = reads$plots)
  }, warning = capture_warning)
  answer$warnings <- warnings_seen
  writeLines(enc2utf8(if (length(warnings_seen)) warnings_seen else "No warnings."),
    file.path(output_dir, "warnings.txt"))
  saveRDS(answer, file.path(output_dir, "model_analysis_objects.rds"))
  answer
}

model_key_metrics <- function(id, fit, diagnostics) {
  sm <- base::summary(fit); cm <- sm$coefficients$cond
  if (inherits(fit, "glmmTMB")) {
    sd <- vapply(glmmTMB::VarCorr(fit)$cond, function(x) unname(attr(x, "stddev")), numeric(1))
    alpha <- 0; alpha_se <- NA_real_
    converged <- fit$fit$convergence; pdHess <- fit$sdr$pdHess; maxgr <- NA_real_
    engine <- "glmmTMB Laplace (historical M2 reference)"; precision <- NA_integer_
  } else {
    sd <- c(Repeat_group_ID = sm$Repeat_SD, TS = sm$TS_SD, Observation_ID = sm$observation_SD)
    alpha <- sm$alpha; alpha_se <- sm$alpha_se
    converged <- sm$convergence; pdHess <- sm$pdHess; maxgr <- sm$max_gradient
    engine <- sm$engine; precision <- fit$quadrature_n
  }
  getp <- function(name) diagnostics$overall_table$p_value[match(name, diagnostics$overall_table$check)]
  data.frame(model = id, n = nrow(stats::model.frame(fit)), k = attr(stats::logLik(fit), "df"),
    logLik = as.numeric(stats::logLik(fit)), AIC = stats::AIC(fit), BIC = stats::BIC(fit),
    convergence = converged, pdHess = pdHess, max_gradient = maxgr,
    alpha = alpha, alpha_SE = alpha_se,
    Repeat_SD = unname(sd["Repeat_group_ID"]), TS_SD = unname(sd["TS"]),
    Observation_SD = unname(sd["Observation_ID"]),
    DO_beta = cm["DO_c", 1], DO_SE = cm["DO_c", 2], DO_p_raw = cm["DO_c", 4],
    DO_CI_low = cm["DO_c", 1] - 1.96 * cm["DO_c", 2],
    DO_CI_high = cm["DO_c", 1] + 1.96 * cm["DO_c", 2],
    Uniformity_p = getp("Uniformity"), Dispersion_ratio = diagnostics$dispersion_ratio,
    Dispersion_p = getp("Dispersion"), Outliers_p = getp("Outliers"),
    Predicted_quantiles_p = getp("Quantiles_vs_predicted_proportion"),
    Reads_quantiles_p = getp("Quantiles_vs_Ruegeria_reads"),
    Mean_proportion_p = getp("Mean_relative_abundance"), Q90_proportion_p = getp("Q90_relative_abundance"),
    diagnostic_status = diagnostics$diagnostic_status, quadrature_nodes = precision,
    engine = engine)
}


# ---- 统一设计、批量拟合和一张总表 ----

formula.mc_olre_fit <- function(x, ...) x$formula

write_excel_csv <- function(data, path) {
  con <- file(path, open = "wt", encoding = "UTF-8")
  on.exit(close(con))
  cat("\ufeff", file = con)
  utils::write.table(data, con, sep = ",", row.names = FALSE,
    col.names = TRUE, na = "", qmethod = "double")
  invisible(path)
}

natural_MC_order <- function(ids) {
  number <- suppressWarnings(as.numeric(sub("^[^0-9]*([0-9]+).*$", "\\1", ids)))
  ids[order(number, ids, na.last = TRUE)]
}

build_shared_design <- function(raw_data, MC_ids, host_reference) {
  required <- c("New_ID", "Sample_ID", "MC_ID", "MC_reads", "Other_reads", "DO_c",
    "Host_type", "Repeat_group_ID", "TS")
  if (!all(required %in% names(raw_data))) stop("缺少列：", paste(setdiff(required, names(raw_data)), collapse = ", "))
  raw <- as.data.frame(raw_data[as.character(raw_data$MC_ID) %in% MC_ids, , drop = FALSE])
  for (column in c("MC_reads", "Other_reads", "DO_c")) {
    if (!is.numeric(raw[[column]])) stop(column, " 必须为数值列。")
    if (any(!is.finite(raw[[column]]) & !is.na(raw[[column]]))) stop(column, " 中有无穷大。")
  }
  for (column in c("MC_reads", "Other_reads")) {
    x <- raw[[column]]
    if (any(x < 0 | x != floor(x), na.rm = TRUE)) stop(column, " 必须是非负整数。")
  }
  for (column in c("New_ID", "Sample_ID", "Host_type", "Repeat_group_ID", "TS")) {
    raw[[column]] <- as.character(raw[[column]])
    raw[[column]][!is.na(raw[[column]]) & !nzchar(trimws(raw[[column]]))] <- NA_character_
  }
  if (anyNA(raw$New_ID)) stop("New_ID 存在空值，无法统一各 MC 的设计。")
  if (anyDuplicated(raw[c("New_ID", "MC_ID")])) stop("同一 New_ID、MC_ID 出现多行；请核对输入。")
  total <- raw$MC_reads + raw$Other_reads
  if ("Ruegeria_reads_total" %in% names(raw)) {
    if (!is.numeric(raw$Ruegeria_reads_total)) stop("Ruegeria_reads_total 必须为数值。")
    comparable <- !is.na(total) & !is.na(raw$Ruegeria_reads_total)
    if (any(total[comparable] != raw$Ruegeria_reads_total[comparable])) stop("MC_reads + Other_reads 与 Ruegeria_reads_total 不一致。")
    missing <- is.na(total)
    total[missing] <- raw$Ruegeria_reads_total[missing]
  }
  raw$.shared_N <- total
  split_rows <- split(seq_len(nrow(raw)), raw$New_ID)
  fields <- c("Sample_ID", "DO_c", "Host_type", "Repeat_group_ID", "TS", ".shared_N")
  inconsistent <- vapply(split_rows, function(i) any(vapply(fields,
    function(column) length(unique(raw[[column]][i])) > 1L, logical(1))), logical(1))
  if (any(inconsistent)) stop("不同 MC 中相同 New_ID 的标签、协变量或分母不一致：",
    paste(head(names(inconsistent)[inconsistent], 8L), collapse = ", "))
  shared <- raw[!duplicated(raw$New_ID), c("New_ID", fields), drop = FALSE]
  keep <- stats::complete.cases(shared) & is.finite(shared$.shared_N) & shared$.shared_N > 0
  shared <- shared[keep, , drop = FALSE]
  if (nrow(shared) < 40L) stop("共同设计中有效样本少于 40 个。")
  levels_host <- sort(unique(shared$Host_type))
  if (!host_reference %in% levels_host) stop("没有找到统一宿主参照：", host_reference)
  levels_host <- c(host_reference, setdiff(levels_host, host_reference))
  log_center <- mean(log10(shared$.shared_N))
  x <- log10(shared$.shared_N) - log_center
  if (length(unique(x)) < 8L) stop("共同设计的正 reads 数不足 8 种。")
  mean_basis <- splines::ns(x, df = 3L, intercept = FALSE)
  shared$Ruegeria_reads_total <- shared$.shared_N
  shared$.shared_N <- NULL
  shared$log_reads_c <- x
  shared <- set_reads_basis_columns(shared, x, mean_basis)
  list(samples = shared, log_reads_center = log_center, mean_basis = mean_basis,
    host_levels = levels_host, host_reference = host_reference, MC_ids = MC_ids)
}

prepare_mc_shared <- function(raw_data, MC_id, shared) {
  prepared <- prepare_mc_data(raw_data, MC_id)
  d <- prepared$data_test
  if (anyDuplicated(d$New_ID)) stop("当前 MC 的 New_ID 不唯一。")
  if (anyNA(match(d$New_ID, shared$samples$New_ID))) stop("模型记录未匹配到共同设计。")
  if (!shared$host_reference %in% d$Host_type) stop("当前 MC 缺少统一宿主参照。")
  d$Host_type <- droplevels(factor(as.character(d$Host_type), levels = shared$host_levels))
  d$log_reads_c <- log10(d$Ruegeria_reads_total) - shared$log_reads_center
  d <- set_reads_basis_columns(d, d$log_reads_c, shared$mean_basis)
  X <- stats::model.matrix(~ DO_c + Host_type + mean_reads_s1 + mean_reads_s2 + mean_reads_s3, d)
  if (qr(X)$rank < ncol(X)) stop("固定效应设计矩阵秩不足，无法估计统一公式。")
  prepared$data_test <- d
  prepared$mean_basis <- shared$mean_basis
  prepared$log_reads_center <- shared$log_reads_center
  prepared$X <- X
  prepared
}

make_model_stub <- function(d) {
  fit_formula <- stats::as.formula(paste("cbind(MC_reads, Other_reads) ~ DO_c + Host_type +",
    "mean_reads_s1 + mean_reads_s2 + mean_reads_s3 +",
    "(1 | Repeat_group_ID) + (1 | TS) + (1 | Observation_ID)"), env = baseenv())
  frame_formula <- stats::as.formula(paste("cbind(MC_reads, Other_reads) ~ DO_c + Host_type +",
    "mean_reads_s1 + mean_reads_s2 + mean_reads_s3 + Repeat_group_ID + TS + Observation_ID"), env = baseenv())
  frame <- stats::model.frame(frame_formula, data = d, na.action = stats::na.fail)
  structure(list(frame = frame, formula = fit_formula), class = "mc_olre_fit")
}

make_initial_parameters <- function(d, X, stub) {
  # 正态 glmmTMB 只提供起点；失败时仍尝试 RTMB 主模型，不自动更换主模型。
  reference <- tryCatch(glmmTMB::glmmTMB(stats::formula(stub), data = d,
    family = stats::binomial(link = "logit"), REML = FALSE), error = function(e) NULL)
  initial <- NULL
  if (!is.null(reference)) initial <- tryCatch(initial_olre_parameters(reference, d), error = function(e) NULL)
  if (is.null(initial) || any(!is.finite(initial$beta))) {
    initial <- list(beta = rep(0, ncol(X)), log_sd_repeat = log(.25),
      log_sd_ts = log(.25), log_sd_observation = log(.75), alpha = 0,
      b_repeat = rep(0, nlevels(d$Repeat_group_ID)), b_ts = rep(0, nlevels(d$TS)))
    initial$beta[1] <- stats::qlogis((sum(d$MC_reads) + .5) / (sum(d$Ruegeria_reads_total) + 1))
  }
  for (name in c("log_sd_repeat", "log_sd_ts", "log_sd_observation")) {
    if (!is.finite(initial[[name]])) initial[[name]] <- -6
    initial[[name]] <- max(-7.5, min(2.5, initial[[name]]))
  }
  initial$alpha <- 0
  initial
}

basic_fit_valid <- function(fit) {
  !is.null(fit) && is.finite(fit$value) && isTRUE(fit$opt$convergence == 0L) &&
    isTRUE(fit$pdHess) && is.finite(fit$max_gradient) && fit$max_gradient < .02 &&
    !inherits(fit$sdr, "error") && is.matrix(fit$sdr$cov.fixed) &&
    all(is.finite(fit$sdr$cov.fixed)) && all(diag(fit$sdr$cov.fixed) > 0)
}

select_best_likelihood <- function(trials, tolerance = .001) {
  finite <- which(vapply(trials, function(x) !is.null(x) && is.finite(x$value), logical(1)))
  if (!length(finite)) return(list(fit = NULL, best = NA_integer_, reproducible = FALSE))
  values <- vapply(trials[finite], function(x) as.numeric(x$value), numeric(1))
  best <- finite[which.min(values)]
  # 即使更好的解触及边界，也保留并标记，不能用较差的内部局部解冒充成功。
  near <- finite[abs(values - min(values)) < tolerance]
  list(fit = trials[[best]], best = best, near_best = near,
    reproducible = sum(vapply(trials[near], basic_fit_valid, logical(1))) >= 2L)
}

fit_shared_candidates <- function(d, X, stub, nodes, check_nodes, starts, output_dir) {
  attempts <- list(); precision_rows <- list(); warnings_seen <- character()
  record_warning <- function(w) { warnings_seen <<- unique(c(warnings_seen, conditionMessage(w))); invokeRestart("muffleWarning") }
  initial <- withCallingHandlers(make_initial_parameters(d, X, stub), warning = record_warning)
  run_fit <- function(p, skew) {
    tryCatch(withCallingHandlers(fit_rtmb_olre(d, X, p, node_n, skew = skew),
      warning = record_warning), error = function(e) { warnings_seen <<- unique(c(warnings_seen, conditionMessage(e))); NULL })
  }
  check_precision <- function(fit) {
    if (is.null(fit)) return(data.frame(fitted_nodes = node_n, check_nodes = check_n,
      fitted_nll = NA_real_, check_nll = NA_real_, absolute_nll_difference = NA_real_,
      check_max_gradient = NA_real_, passed = FALSE))
    tryCatch(check_quadrature_precision(fit, d, X, check_n), error = function(e) {
      warnings_seen <<- unique(c(warnings_seen, conditionMessage(e)))
      data.frame(fitted_nodes = node_n, check_nodes = check_n, fitted_nll = fit$value,
        check_nll = NA_real_, absolute_nll_difference = NA_real_, check_max_gradient = NA_real_, passed = FALSE)
    })
  }
  node_n <- as.integer(nodes)
  for (stage in 1:2) {
    cat("  正态参照：", node_n, " 个积分节点。\n", sep = "")
    normal <- run_fit(initial, FALSE)
    trials <- vector("list", length(starts))
    for (j in seq_along(starts)) {
      cat("  偏正态起点 ", j, "/", length(starts), "：alpha = ", starts[j], "\n", sep = "")
      p <- if (!is.null(normal) && all(is.finite(normal$parameters$beta))) normal$parameters else initial
      p$alpha <- starts[j]
      trial <- run_fit(p, TRUE)
      trials[j] <- list(trial)
      attempts[[length(attempts) + 1L]] <- data.frame(stage = stage, nodes = node_n,
        alpha_start = starts[j], nll = if (is.null(trial)) NA_real_ else as.numeric(trial$value),
        alpha = if (is.null(trial)) NA_real_ else trial$parameters$alpha,
        convergence = if (is.null(trial)) NA_integer_ else trial$opt$convergence,
        pdHess = !is.null(trial) && isTRUE(trial$pdHess),
        max_gradient = if (is.null(trial)) NA_real_ else trial$max_gradient,
        basic_valid = basic_fit_valid(trial))
    }
    selected <- select_best_likelihood(trials)
    skew <- selected$fit
    check_n <- max(as.integer(check_nodes), 2L * node_n - 1L)
    p_normal <- check_precision(normal); p_skew <- check_precision(skew)
    precision_rows[[length(precision_rows) + 1L]] <- cbind(model = "M2_NORMAL_AGHQ", stage = stage, p_normal)
    precision_rows[[length(precision_rows) + 1L]] <- cbind(model = "M2_SKEW_AGHQ", stage = stage, p_skew)
    need_more_nodes <- !isTRUE(p_skew$passed) || (!is.null(normal) && !isTRUE(p_normal$passed))
    if (stage == 1L && need_more_nodes) {
      cat("  提高积分节点数，并在统一节点下重新拟合两个候选。\n")
      node_n <- check_n
      if (!is.null(normal)) { initial <- normal$parameters; initial$alpha <- 0 }
      next
    }
    break
  }
  attempt_table <- do.call(rbind, attempts)
  attempt_table$near_best_final_stage <- attempt_table$stage == stage & is.finite(attempt_table$nll) &
    !is.null(skew) & abs(attempt_table$nll - if (is.null(skew)) Inf else skew$value) < .001
  write_result_table(attempt_table, output_dir, "01_multistart_fit_attempts.csv")
  write_result_table(do.call(rbind, precision_rows), output_dir, "01_quadrature_precision_checks.csv")
  wrapped <- list(normal = if (is.null(normal)) NULL else wrap_rtmb_fit(normal, stub, d, X, "M2_NORMAL_AGHQ"),
    primary = if (is.null(skew)) NULL else wrap_rtmb_fit(skew, stub, d, X, "M2_SKEW_AGHQ"))
  integral_tables <- list()
  for (name in names(wrapped)) if (!is.null(wrapped[[name]])) {
    integral_tables[[name]] <- tryCatch(cbind(model = name,
      check_conditional_integrals(wrapped[[name]])), error = function(e) {
        warnings_seen <<- unique(c(warnings_seen, conditionMessage(e)))
        data.frame(model = name, passed = FALSE)
      })
  }
  bind_flexible <- function(tables) {
    if (!length(tables)) return(data.frame())
    columns <- unique(unlist(lapply(tables, names)))
    do.call(rbind, lapply(tables, function(x) { for (name in setdiff(columns, names(x))) x[[name]] <- NA; x[, columns, drop = FALSE] }))
  }
  integral_table <- bind_flexible(integral_tables)
  write_result_table(integral_table, output_dir, "01_independent_conditional_integrals.csv")
  direct_ok <- function(name) !is.null(integral_tables[[name]]) && all(integral_tables[[name]]$passed %in% TRUE)
  alpha_interior <- !is.null(skew) && is.finite(skew$parameters$alpha) && abs(skew$parameters$alpha) < 19.99
  observation_interior <- !is.null(skew) && is.finite(skew$parameters$log_sd_observation) &&
    skew$parameters$log_sd_observation > -7.99 && skew$parameters$log_sd_observation < 2.99
  primary_valid <- basic_fit_valid(skew) && isTRUE(p_skew$passed) && direct_ok("primary") &&
    alpha_interior && observation_interior && isTRUE(selected$reproducible)
  normal_valid <- basic_fit_valid(normal) && isTRUE(p_normal$passed) && direct_ok("normal")
  reasons <- c(if (is.null(skew)) "所有偏正态起点均失败",
    if (!is.null(skew) && !basic_fit_valid(skew)) "优化或 Hessian/协方差检查未通过",
    if (!is.null(skew) && !alpha_interior) "alpha 触及边界",
    if (!is.null(skew) && !observation_interior) "Observation SD 触及边界",
    if (!is.null(skew) && !isTRUE(selected$reproducible)) "最优解未由至少两个有效起点复现",
    if (!isTRUE(p_skew$passed)) "积分节点核对未通过",
    if (!direct_ok("primary")) "独立条件积分核对未通过")
  writeLines(enc2utf8(if (length(warnings_seen)) warnings_seen else "No warnings."), file.path(output_dir, "01_fitting_warnings.txt"))
  list(primary = wrapped$primary, normal = wrapped$normal,
    primary_valid = primary_valid, normal_valid = normal_valid,
    reasons = reasons, alpha_interior = alpha_interior, reproducible = selected$reproducible,
    quadrature_ok = isTRUE(p_skew$passed), direct_integrals_ok = direct_ok("primary"),
    warnings = warnings_seen, attempts = attempt_table)
}

wald_block <- function(fit, pattern) {
  if (isTRUE(fit$penalty_requested)) return(list(statistic = NA_real_, df = NA_integer_, p = NA_real_))
  index <- grep(pattern, colnames(fit$X))
  if (!length(index)) return(list(statistic = NA_real_, df = 0L, p = NA_real_))
  beta <- fit$parameters$beta[index]
  V <- fit$cov_fixed[index, index, drop = FALSE]
  value <- tryCatch(as.numeric(crossprod(beta, solve(V, beta))), error = function(e) NA_real_)
  if (!is.finite(value) || value < 0) return(list(statistic = NA_real_, df = length(index), p = NA_real_))
  list(statistic = value, df = length(index), p = stats::pchisq(value, length(index), lower.tail = FALSE))
}

diagnose_primary <- function(fit, d, n_sim, seed, output_dir) {
  answer <- diagnose_skew_comparison_model("M2_SKEW_AGHQ", fit, d, n_sim, seed, output_dir)
  res <- answer$res_check
  N <- d$Ruegeria_reads_total
  zero_test <- DHARMa::testGeneric(res, summary = function(x) sum(x == 0), plot = FALSE)
  one_test <- DHARMa::testGeneric(res, summary = function(x) sum(x == N), plot = FALSE)
  answer$overall_table <- rbind(answer$overall_table,
    data.frame(check = c("Zero_MC_count", "All_reads_MC_count"),
      p_value = c(zero_test$p.value, one_test$p.value), status = "ok"))
  answer$zero_count_check <- zero_test
  answer$all_reads_count_check <- one_test
  write_result_table(answer$overall_table, output_dir, "04_overall_diagnostics.csv")
  saveRDS(answer, file.path(output_dir, "model_analysis_objects.rds"))
  answer
}

summary_row <- function(MC_id, n_source) {
  data.frame(MC_ID = MC_id, status = "未完成", n_source = as.integer(n_source),
    n_used = NA_integer_, n_excluded = NA_integer_, n_repeat_groups = NA_integer_, n_TS = NA_integer_,
    n_zero_MC = NA_integer_, n_all_reads_MC = NA_integer_, observed_mean_proportion = NA_real_,
    host_reference = host_reference, DO_beta = NA_real_, DO_SE = NA_real_,
    DO_CI_low = NA_real_, DO_CI_high = NA_real_, DO_p_raw = NA_real_, DO_p_BH = NA_real_,
    DO_OR = NA_real_, DO_OR_CI_low = NA_real_, DO_OR_CI_high = NA_real_,
    Host_Wald_chisq = NA_real_, Host_Wald_df = NA_integer_, Host_p_raw = NA_real_, Host_p_BH = NA_real_,
    Reads_Wald_chisq = NA_real_, Reads_Wald_df = NA_integer_, Reads_p_raw = NA_real_, Reads_p_BH = NA_real_,
    alpha = NA_real_, alpha_SE = NA_real_, Repeat_SD = NA_real_, TS_SD = NA_real_, Observation_SD = NA_real_,
    AIC = NA_real_, BIC = NA_real_, logLik = NA_real_, k = NA_integer_,
    normal_AIC = NA_real_, delta_AIC_skew_minus_normal = NA_real_, normal_DO_beta = NA_real_,
    convergence = NA_integer_, pdHess = FALSE, max_gradient = NA_real_, numerical_valid = FALSE,
    multistart_reproducible = FALSE, quadrature_ok = FALSE, direct_integrals_ok = FALSE,
    diagnostic_status = "未运行", Uniformity_p = NA_real_, Dispersion_p = NA_real_, Dispersion_ratio = NA_real_,
    Outliers_p = NA_real_, Predicted_quantiles_p = NA_real_, Reads_quantiles_p = NA_real_,
    Mean_proportion_p = NA_real_, Q90_proportion_p = NA_real_, Zero_count_p = NA_real_, All_reads_count_p = NA_real_,
    Tail_min_BH_p = NA_real_, Reads_location_trend_BH_p = NA_real_, Reads_deviation_trend_BH_p = NA_real_,
    diagnostic_flags = "", notes = "", stringsAsFactors = FALSE)
}

fill_fit_summary <- function(row, candidates) {
  fit <- candidates$primary
  row$numerical_valid <- isTRUE(candidates$primary_valid)
  row$multistart_reproducible <- isTRUE(candidates$reproducible)
  row$quadrature_ok <- isTRUE(candidates$quadrature_ok)
  row$direct_integrals_ok <- isTRUE(candidates$direct_integrals_ok)
  row$notes <- paste(candidates$reasons, collapse = "; ")
  if (!is.null(fit)) {
    row$convergence <- fit$fit$convergence; row$pdHess <- fit$pdHess
    row$max_gradient <- fit$max_gradient; row$alpha <- fit$parameters$alpha
    row$alpha_SE <- fit$alpha_se
    row$Repeat_SD <- exp(fit$parameters$log_sd_repeat)
    row$TS_SD <- exp(fit$parameters$log_sd_ts)
    row$Observation_SD <- exp(fit$parameters$log_sd_observation)
    row$AIC <- stats::AIC(fit); row$BIC <- stats::BIC(fit)
    row$logLik <- as.numeric(stats::logLik(fit)); row$k <- fit$k
  }
  if (isTRUE(candidates$normal_valid)) {
    row$normal_AIC <- stats::AIC(candidates$normal)
    row$normal_DO_beta <- candidates$normal$parameters$beta[match("DO_c", colnames(candidates$normal$X))]
  }
  if (!row$numerical_valid) { row$status <- "数值拟合需复核"; return(row) }
  cm <- base::summary(fit)$coefficients$cond
  row$DO_beta <- cm["DO_c", 1]; row$DO_SE <- cm["DO_c", 2]; row$DO_p_raw <- cm["DO_c", 4]
  row$DO_CI_low <- row$DO_beta - 1.96 * row$DO_SE
  row$DO_CI_high <- row$DO_beta + 1.96 * row$DO_SE
  row$DO_OR <- exp(row$DO_beta)
  row$DO_OR_CI_low <- exp(row$DO_CI_low); row$DO_OR_CI_high <- exp(row$DO_CI_high)
  host <- wald_block(fit, "^Host_type"); reads <- wald_block(fit, "^mean_reads_s")
  row$Host_Wald_chisq <- host$statistic; row$Host_Wald_df <- host$df; row$Host_p_raw <- host$p
  row$Reads_Wald_chisq <- reads$statistic; row$Reads_Wald_df <- reads$df; row$Reads_p_raw <- reads$p
  if (isTRUE(candidates$normal_valid)) row$delta_AIC_skew_minus_normal <- row$AIC - row$normal_AIC
  row$status <- "数值有效；等待诊断"
  row
}

fill_diagnostic_summary <- function(row, diagnostics) {
  getp <- function(name) {
    value <- diagnostics$overall_table$p_value[match(name, diagnostics$overall_table$check)]
    if (length(value) != 1L) NA_real_ else value
  }
  mapping <- c(Uniformity_p = "Uniformity", Dispersion_p = "Dispersion", Outliers_p = "Outliers",
    Predicted_quantiles_p = "Quantiles_vs_predicted_proportion", Reads_quantiles_p = "Quantiles_vs_Ruegeria_reads",
    Mean_proportion_p = "Mean_relative_abundance", Q90_proportion_p = "Q90_relative_abundance",
    Zero_count_p = "Zero_MC_count", All_reads_count_p = "All_reads_MC_count")
  for (name in names(mapping)) row[[name]] <- getp(mapping[[name]])
  row$Dispersion_ratio <- diagnostics$dispersion_ratio
  row$diagnostic_status <- diagnostics$diagnostic_status
  if (any(is.finite(diagnostics$tail_table$p_BH_within_six))) row$Tail_min_BH_p <- min(diagnostics$tail_table$p_BH_within_six, na.rm = TRUE)
  trends <- diagnostics$reads_diagnostics$main_table
  row$Reads_location_trend_BH_p <- trends$p_BH[match("Residual_location", trends$diagnostic)]
  row$Reads_deviation_trend_BH_p <- trends$p_BH[match("Absolute_deviation", trends$diagnostic)]
  values <- unlist(row[, c(names(mapping), "Tail_min_BH_p", "Reads_location_trend_BH_p", "Reads_deviation_trend_BH_p")])
  flags <- names(values)[is.finite(values) & values < .05]
  row$diagnostic_flags <- paste(flags, collapse = "; ")
  row$status <- if (!isTRUE(diagnostics$diagnostics_complete)) "数值有效；诊断未完成" else
    if (length(flags)) "数值有效；残差需复核" else "数值有效；未检出主要残差偏差"
  row
}

add_cross_MC_adjustment <- function(table) {
  for (term in c("DO", "Host", "Reads")) {
    p <- table[[paste0(term, "_p_raw")]]
    p[!table$numerical_valid] <- NA_real_
    table[[paste0(term, "_p_BH")]] <- adjust_valid_p(p, method = "BH", n = nrow(table))
  }
  table
}

fixed_coefficient_table <- function(fit, MC_id) {
  cm <- base::summary(fit)$coefficients$cond
  data.frame(MC_ID = MC_id, term = rownames(cm), estimate = cm[, 1], SE = cm[, 2],
    CI_low = cm[, 1] - 1.96 * cm[, 2], CI_high = cm[, 1] + 1.96 * cm[, 2],
    z = cm[, 3], p_raw = cm[, 4], row.names = NULL)
}

make_reads_curve <- function(fit, shared, MC_id) {
  x <- seq(min(shared$samples$log_reads_c), max(shared$samples$log_reads_c), length.out = 150L)
  B <- stats::predict(shared$mean_basis, newx = x)
  B0 <- as.numeric(stats::predict(shared$mean_basis, newx = 0))
  index <- match(paste0("mean_reads_s", 1:3), colnames(fit$X))
  contrast <- sweep(B, 2L, B0)
  value <- as.vector(contrast %*% fit$parameters$beta[index])
  V <- fit$cov_fixed[index, index, drop = FALSE]
  SE <- sqrt(pmax(0, rowSums((contrast %*% V) * contrast)))
  data.frame(MC_ID = MC_id, Ruegeria_reads_total = 10^(x + shared$log_reads_center),
    log10_reads = x + shared$log_reads_center, reads_contribution_logit = value,
    CI_low = value - 1.96 * SE, CI_high = value + 1.96 * SE)
}

html_escape <- function(x) {
  x <- gsub("&", "&amp;", as.character(x), fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub('"', "&quot;", x, fixed = TRUE)
  x
}

fmt_number <- function(x, digits = 3L) {
  ifelse(is.finite(x), formatC(x, format = "f", digits = digits), "—")
}

fmt_p <- function(x) {
  vapply(x, function(p) if (!is.finite(p)) "—" else if (p == 0) "0*" else
    if (p < .001) formatC(p, format = "e", digits = 2) else formatC(p, format = "f", digits = 3), character(1))
}

write_summary_html <- function(table, run_dir) {
  view <- data.frame(MC = table$MC_ID, 状态 = table$status, 纳入记录 = table$n_used,
    排除数 = table$n_excluded, 零值数 = table$n_zero_MC,
    DO_beta及95CI = ifelse(is.finite(table$DO_beta), paste0(fmt_number(table$DO_beta), " [",
      fmt_number(table$DO_CI_low), ", ", fmt_number(table$DO_CI_high), "]"), "—"),
    DO_p = fmt_p(table$DO_p_raw), DO_BH = fmt_p(table$DO_p_BH),
    DO_OR及95CI = ifelse(is.finite(table$DO_OR), paste0(fmt_number(table$DO_OR), " [",
      fmt_number(table$DO_OR_CI_low), ", ", fmt_number(table$DO_OR_CI_high), "]"), "—"),
    Host_p = fmt_p(table$Host_p_raw), Host_BH = fmt_p(table$Host_p_BH),
    Reads_p = fmt_p(table$Reads_p_raw), Reads_BH = fmt_p(table$Reads_p_BH),
    AIC = fmt_number(table$AIC, 2L), BIC = fmt_number(table$BIC, 2L),
    ΔAIC对正态 = fmt_number(table$delta_AIC_skew_minus_normal, 2L),
    alpha = fmt_number(table$alpha), Repeat_SD = fmt_number(table$Repeat_SD),
    TS_SD = fmt_number(table$TS_SD), Observation_SD = fmt_number(table$Observation_SD),
    Uniformity_p = fmt_p(table$Uniformity_p), Dispersion_p = fmt_p(table$Dispersion_p),
    Predicted_quantiles_p = fmt_p(table$Predicted_quantiles_p), Reads_quantiles_p = fmt_p(table$Reads_quantiles_p),
    其他诊断提示 = table$diagnostic_flags, 说明 = table$notes, check.names = FALSE)
  group_titles <- c("MC", "数据与状态", "DO 关联", "Host 和 reads 整体关联", "拟合与随机效应", "残差诊断")
  group_spans <- c(1L, 4L, 4L, 4L, 7L, 6L)
  header1 <- paste0('<tr>', paste0('<th colspan="', group_spans, '">', group_titles, '</th>', collapse = ''), '</tr>')
  header2 <- paste0('<tr>', paste0('<th>', html_escape(names(view)), '</th>', collapse = ''), '</tr>')
  rows <- vapply(seq_len(nrow(view)), function(i) {
    cells <- vapply(seq_along(view), function(j) {
      value <- view[[j]][i]
      if (is.na(value)) value <- "—"
      significant <- names(view)[j] %in% c("DO_BH", "Host_BH", "Reads_BH") &&
        is.finite(table[[c(DO_BH = "DO_p_BH", Host_BH = "Host_p_BH", Reads_BH = "Reads_p_BH")[[names(view)[j]]]]][i]) &&
        table[[c(DO_BH = "DO_p_BH", Host_BH = "Host_p_BH", Reads_BH = "Reads_p_BH")[[names(view)[j]]]]][i] < .05
      paste0('<td', if (significant) ' class="sig"' else '', '>', html_escape(value), '</td>')
    }, character(1))
    paste0('<tr class="', if (table$numerical_valid[i]) 'valid' else 'review', '">', paste(cells, collapse = ''), '</tr>')
  }, character(1))
  text <- c('<!doctype html><html lang="zh-CN"><meta charset="utf-8"><title>所有 MC 的 GLMM 结果</title>',
    '<style>body{font:15px/1.5 Arial,"Microsoft YaHei",sans-serif;margin:24px;color:#172b3a;background:#f6f8fa}h1{font-size:25px;margin-bottom:6px}.intro{max-width:1100px}.scroll{overflow:auto;max-height:76vh;border:1px solid #c9d4dd;background:white}table{border-collapse:separate;border-spacing:0;min-width:2900px}th,td{padding:10px 12px;border-bottom:1px solid #e0e6eb;border-right:1px solid #eef2f5;white-space:nowrap;text-align:right}th{background:#193c55;color:white}thead tr:first-child th{position:sticky;top:0;z-index:4;height:26px;text-align:center}thead tr:nth-child(2) th{position:sticky;top:47px;z-index:4}tbody td:first-child,thead tr:nth-child(2) th:first-child{position:sticky;left:0;z-index:3;min-width:72px;text-align:left;background:#edf3f7;font-weight:bold;color:#172b3a}tbody tr:nth-child(even) td{background:#f5f8fa}tbody tr.review td{background:#fff2df}.sig{font-weight:bold;color:#005845;background:#e0f2e9!important}td:nth-child(2),td:nth-last-child(1),td:nth-last-child(2){white-space:normal;min-width:220px;text-align:left}.foot{max-width:1100px;color:#405363}code{background:#e6edf3;padding:2px 5px}</style>',
    '<h1>所有 MC 的 GLMM 结果</h1>',
    paste0('<p class="intro">统一主模型：二项 logit + DO + Host + reads 自然样条 df=3 + Repeat/TS 正态随机截距 + 中心化偏正态样本偏移。共 ', nrow(table), ' 个 MC；数值核对通过 ', sum(table$numerical_valid), ' 个。每个 MC 一行，可横向滚动查看。记录按现有 New_ID 匹配，Sample_ID 标签可重复。</p>'),
    '<div class="scroll"><table><thead>', header1, header2, '</thead><tbody>', rows, '</tbody></table></div>',
    '<div class="foot"><p>DO β 是所有宿主共用的条件 logit 斜率。OR=exp(β)，表示 DO 增加 1 个原始单位时的 odds 比值。95% CI 为逐项 Wald 区间，未做同时区间校正。</p>',
    '<p>DO、Host 整体检验、reads 整体检验分别跨全部计划 MC 做 BH 校正。Host/reads 为近似联合 Wald 检验；Host 检验不等于两两比较。ΔAIC=偏正态−正态，仅在同一 MC、相同数据和积分框架内比较。</p>',
    '<p>偏正态始终是主模型。数值无效时，主模型推断列留空，alpha/AIC 等仅用于定位失败。残差诊断 p 值用于定位模型偏差；未检出偏差不等于已证明拟合正确。0* 表示有限模拟中未观察到同样极端的统计量。</p>',
    '<p>完整数值见 <a href="00_all_MC_summary.csv">00_all_MC_summary.csv</a>；逐列定义见 <a href="00_column_dictionary.csv">00_column_dictionary.csv</a>。本轮为完整数据主分析，未运行逐批次剔除检查。</p></div></html>')
  writeLines(enc2utf8(text), file.path(run_dir, "00_all_MC_summary.html"))
}

make_column_dictionary <- function() {
  columns <- names(summary_row("MC", 0L))
  description <- setNames(rep("", length(columns)), columns)
  definitions <- c(MC_ID = "输入表中的 MC 名称，每个 MC 独立拟合同一主模型。",
    status = "同时说明数值和残差诊断状态；未检出主要偏差不等于证明模型正确。",
    n_source = "输入表中这个 MC 的原始行数。", n_used = "完整且分母为正的原始记录数，每个 New_ID 一行；Sample_ID 标签可能重复。",
    n_excluded = "因缺失或分母为零排除的行数。", n_repeat_groups = "实际纳入的 Repeat_group 数。",
    n_TS = "实际纳入的 TS 数。", n_zero_MC = "MC_reads=0 的样本数，这些样本会保留。",
    n_all_reads_MC = "Other_reads=0 的样本数，即 MC 在 Ruegeria 内比例为 100%。",
    observed_mean_proportion = "纳入样本 MC_reads/Ruegeria_reads_total 的算术平均，范围 0-1。",
    host_reference = "统一的宿主参考水平。", DO_beta = "中心化 DO 的条件 logit 系数，宿主之间共用一条 DO 斜率。",
    DO_SE = "DO 系数标准误。", DO_CI_low = "DO 系数 95% Wald CI 下限。", DO_CI_high = "DO 系数 95% Wald CI 上限。",
    DO_p_raw = "DO 系数双侧 Wald 原始 p 值。", DO_p_BH = "所有计划 MC 的 DO 检验组成一个 BH 校正家族。",
    DO_OR = "exp(DO_beta)，DO 每增加 1 原始单位时的条件 odds 比值。",
    DO_OR_CI_low = "DO odds ratio 的 95% CI 下限。", DO_OR_CI_high = "DO odds ratio 的 95% CI 上限。",
    Host_Wald_chisq = "所有宿主对照系数共同等于 0 的近似联合 Wald 统计量。",
    Host_Wald_df = "宿主整体检验自由度。", Host_p_raw = "宿主整体关联检验原始 p，不是两两比较。",
    Host_p_BH = "宿主整体检验在全部计划 MC 之间单独做 BH 校正。",
    Reads_Wald_chisq = "三项 reads 样条系数共同等于 0 的近似联合 Wald 统计量。",
    Reads_Wald_df = "reads 整体检验自由度，通常为 3。", Reads_p_raw = "reads 曲线整体关联的原始 p。",
    Reads_p_BH = "reads 整体检验在全部计划 MC 之间单独做 BH 校正。",
    alpha = "中心化偏正态样本偏移的形状参数：正为右偏、负为左偏、0 为正态。",
    alpha_SE = "形状参数标准误，只辅助数值判断；不进行 alpha 的普通 Wald 显著性检验。",
    Repeat_SD = "Repeat_group 随机截距 SD，logit 尺度。", TS_SD = "TS 随机截距 SD，logit 尺度。",
    Observation_SD = "中心化偏正态样本偏移的实际 SD，logit 尺度；同一 MC 共用。",
    AIC = "偏正态主模型 AIC；数值无效时只作排错信息。", BIC = "偏正态主模型 BIC。",
    logLik = "偏正态主模型最大化对数似然。", k = "估计参数数目。",
    normal_AIC = "数值有效的正态对照 AIC。", delta_AIC_skew_minus_normal = "主模型 AIC 减同一 MC 的正态对照 AIC；负值偏向偏正态。",
    normal_DO_beta = "数值有效的正态对照 DO 系数，仅用于敏感性对照。",
    convergence = "优化器收敛码，0 为正常结束。", pdHess = "Hessian 是否正定。",
    max_gradient = "最优点的最大绝对梯度。", numerical_valid = "优化、Hessian、协方差、多起点、参数边界及两类积分检查的综合结果。",
    multistart_reproducible = "是否至少两个数值有效起点复现最优负对数似然（差小于 0.001）。",
    quadrature_ok = "提高积分节点数后，负对数似然差小于 0.01 且梯度小于 0.02。",
    direct_integrals_ok = "代表性样本的 AGHQ 条件似然与独立直接积分差小于 0.0001。",
    diagnostic_status = "ok 表示诊断计算完成，partial 表示部分未完成；不代表所有诊断均不显著。",
    Uniformity_p = "模拟残差整体是否偏离均匀分布。", Dispersion_p = "整体残差离散程度是否超出模拟预期。",
    Dispersion_ratio = "观测与模拟离散程度之比。", Outliers_p = "模拟残差异常值检验（bootstrap）。",
    Predicted_quantiles_p = "残差 25%、50%、75% 分位数是否偏离参考水平或随预测比例变化。",
    Reads_quantiles_p = "残差分位数是否偏离参考水平或随 log10 Ruegeria reads 变化。",
    Mean_proportion_p = "观测相对丰度均值与模型模拟是否吻合。", Q90_proportion_p = "观测相对丰度第 90 百分位与模型模拟是否吻合。",
    Zero_count_p = "MC_reads=0 的样本数是否超出模型模拟预期。", All_reads_count_p = "MC_reads=N 的样本数是否超出模型模拟预期。",
    Tail_min_BH_p = "六项高低比例端诊断在同一 MC 内 BH 校正后的最小 p。",
    Reads_location_trend_BH_p = "残差中位数随 reads 的平滑趋势 p，与绝对偏离趋势组成两项 BH 校正。",
    Reads_deviation_trend_BH_p = "2*abs(残差-0.5) 中位数随 reads 的趋势；兼受中心偏移和分布形状影响。",
    diagnostic_flags = "原始诊断 p<0.05 或指定组内校正 p<0.05 的项目，供定位问题。",
    notes = "数值失败、诊断失败或输入异常的具体原因。")
  description[names(definitions)] <- definitions
  stopifnot(all(nzchar(description)))
  data.frame(column = columns, description = unname(description), row.names = NULL)
}

write_final_outputs <- function(table, run_dir, completed, planned) {
  table <- add_cross_MC_adjustment(table)
  write_excel_csv(table, file.path(run_dir, "00_all_MC_summary.csv"))
  write_summary_html(table, run_dir)
  write_excel_csv(make_column_dictionary(), file.path(run_dir, "00_column_dictionary.csv"))
  text <- capture.output({
    cat("============================================================\n所有 MC 的统一 GLMM 主分析\n")
    cat("已完成：", completed, "/", planned, "；数值有效：", sum(table$numerical_valid), "。\n", sep = "")
    cat("主模型：二项 logit；DO + Host + reads 自然样条 df=3；Repeat/TS 正态随机截距；Observation 中心化偏正态。\n")
    cat("1. 数据、拟合和诊断状态\n")
    print(table[, c("MC_ID", "n_source", "n_used", "n_excluded", "n_zero_MC", "numerical_valid", "status")], row.names = FALSE)
    cat("\n2. DO 主要结果（BH 跨所有计划 MC；95% CI 为逐项 Wald 区间）\n")
    print(table[, c("MC_ID", "DO_beta", "DO_CI_low", "DO_CI_high", "DO_OR", "DO_p_raw", "DO_p_BH")], row.names = FALSE)
    cat("\n3. 宿主与 reads 整体关联（近似联合 Wald；两个检验家族分别 BH）\n")
    print(table[, c("MC_ID", "Host_p_raw", "Host_p_BH", "Reads_p_raw", "Reads_p_BH")], row.names = FALSE)
    cat("\n4. 偏斜、波动及正态对照\n")
    print(table[, c("MC_ID", "alpha", "Repeat_SD", "TS_SD", "Observation_SD", "AIC", "BIC", "delta_AIC_skew_minus_normal")], row.names = FALSE)
    cat("\n5. 需要复核的项目\n")
    review <- table[!table$numerical_valid | table$status != "数值有效；未检出主要残差偏差", c("MC_ID", "diagnostic_flags", "notes"), drop = FALSE]
    if (nrow(review)) print(review, row.names = FALSE) else cat("本轮未检出主要数值或残差问题。\n")
    cat("\n优先查看：00_all_MC_summary.html；完整数值：00_all_MC_summary.csv。\n")
    cat("正态是对照，主模型始终为偏正态。数值无效的推断留空并保留失败行。\n")
    cat("诊断 ok 表示计算完成；诊断 p 值用于定位偏差。所有零值保留，TS8 若缺失 DO 会按相同规则排除。\n")
    cat("本轮未做逐 TS 剔除敏感性分析；没有对被剔除 TS 计算预测评分。\n")
    cat("结果目录：", normalizePath(run_dir, winslash = "/", mustWork = FALSE), "\n", sep = "")
  })
  writeLines(enc2utf8(text), file.path(run_dir, "00_key_results.txt"))
  invisible(list(table = table, text = text))
}


# ---- 零膨胀候选与全部 MC 的比较流程 ----

zi_worker_fit <- function(p, d, design, nodes) {
  messages <- character()
  fit <- tryCatch(withCallingHandlers(
    fit_rtmb_olre(d, design, p, nodes, skew = TRUE, zi = TRUE),
    warning = function(w) {
      messages <<- unique(c(messages, conditionMessage(w))); invokeRestart("muffleWarning")
    }), error = function(e) e)
  if (inherits(fit, "error")) return(list(fit = NULL, warnings = messages, error = conditionMessage(fit)))
  # AD tapes are process-local. Only numerical estimates/covariance cross workers.
  fit$objective <- NULL
  if (!inherits(fit$sdr, "error")) fit$sdr <- list(pdHess = fit$sdr$pdHess, cov.fixed = fit$sdr$cov.fixed)
  list(fit = fit, warnings = messages, error = "")
}

bind_flexible_tables <- function(tables) {
  tables <- Filter(function(x) is.data.frame(x) && nrow(x), tables)
  if (!length(tables)) return(data.frame())
  columns <- unique(unlist(lapply(tables, names)))
  do.call(rbind, lapply(tables, function(x) {
    for (nm in setdiff(columns, names(x))) x[[nm]] <- NA
    x[, columns, drop = FALSE]
  }))
}

cached_fit_result <- function(fit, d, X) {
  fun <- make_olre_objective(d, X, normal_quadrature(fit$quadrature_n))
  obj <- RTMB::MakeADFun(fun, fit$parameters, random = c("b_repeat", "b_ts"), silent = TRUE)
  value <- obj$fn(obj$par)
  gradient <- obj$gr(obj$par)
  if (abs(value - fit$nll) > .01) stop("缓存似然与当前实现不一致。")
  list(objective = obj, opt = list(par = obj$par, objective = value,
    convergence = fit$fit$convergence, message = "Verified cached fit"),
    parameters = obj$env$parList(obj$par), value = value,
    max_gradient = max(abs(gradient)), pdHess = fit$pdHess,
    sdr = list(pdHess = fit$pdHess, cov.fixed = fit$cov_fixed),
    quadrature_n = fit$quadrature_n, skew = TRUE, zi = FALSE, k = fit$k)
}

load_verified_baseline <- function(mc, d, X, output_dir) {
  if (!isTRUE(reuse_baseline_cache)) return(NULL)
  model_path <- file.path(baseline_cache_dir, mc, "fitted_models.rds")
  if (!file.exists(model_path)) return(NULL)
  saved <- readRDS(model_path)
  candidates <- saved$candidates
  fit <- candidates$primary
  if (is.null(fit) || !inherits(fit, "mc_olre_fit")) return(NULL)
  columns <- c("New_ID", "Sample_ID", "MC_ID", "MC_reads", "Other_reads", "DO_c",
    "Host_type", "Repeat_group_ID", "TS", "Ruegeria_reads_total", "log_reads_c",
    "mean_reads_s1", "mean_reads_s2", "mean_reads_s3")
  same_data <- nrow(fit$data) == nrow(d) && all(columns %in% names(fit$data)) &&
    all(vapply(columns, function(nm) {
      if (is.numeric(d[[nm]])) isTRUE(all.equal(as.numeric(d[[nm]]), as.numeric(fit$data[[nm]]), tolerance = 1e-12))
      else identical(as.character(d[[nm]]), as.character(fit$data[[nm]]))
    }, logical(1)))
  if (!same_data || !identical(colnames(X), colnames(fit$X)) ||
      !isTRUE(all.equal(unname(X), unname(fit$X), tolerance = 1e-12))) {
    warning("缓存样本或设计不一致，重新拟合基线：", mc)
    return(NULL)
  }
  numerical <- tryCatch(cached_fit_result(fit, d, X), error = function(e) e)
  precision <- tryCatch({
    if (inherits(numerical, "error")) stop(conditionMessage(numerical))
    check_quadrature_precision(numerical, d, X, 2L * fit$quadrature_n - 1L)
  }, error = function(e) data.frame(passed = FALSE, reason = conditionMessage(e)))
  direct <- tryCatch(check_conditional_integrals(fit), error = function(e)
    data.frame(passed = FALSE, reason = conditionMessage(e)))
  covariance_ok <- is.matrix(fit$cov_fixed) && all(is.finite(fit$cov_fixed)) &&
    all(diag(fit$cov_fixed) > 0)
  candidates$primary_valid <- isTRUE(candidates$primary_valid) &&
    isTRUE(precision$passed) && all(direct$passed) && covariance_ok &&
    !inherits(numerical, "error") && numerical$max_gradient < .02
  candidates$quadrature_ok <- isTRUE(precision$passed)
  candidates$direct_integrals_ok <- all(direct$passed)
  if (!candidates$primary_valid && !length(candidates$reasons)) candidates$reasons <- "缓存复核未通过"
  if (inherits(numerical, "error")) candidates$reasons <- unique(c(candidates$reasons, conditionMessage(numerical)))
  write_result_table(precision, output_dir, "01_baseline_precision_recheck.csv")
  write_result_table(direct, output_dir, "01_baseline_direct_integral_recheck.csv")
  writeLines(normalizePath(model_path, winslash = "/"), file.path(output_dir, "01_baseline_cache_source.txt"))
  diag_path <- file.path(baseline_cache_dir, mc, "primary_diagnostics", "model_analysis_objects.rds")
  list(candidates = candidates, source = "verified_previous_run", numerical = numerical,
    diagnostics = if (candidates$primary_valid && file.exists(diag_path)) readRDS(diag_path) else NULL,
    diagnostic_source = if (file.exists(diag_path)) diag_path else "")
}

check_pi0_likelihood_boundary <- function(zi_fit, baseline, output_dir) {
  base_fit <- baseline$primary
  pi_hat <- if (isTRUE(zi_fit$zi)) stats::plogis(zi_fit$parameters$gamma_zi) else 0
  if (!isTRUE(baseline$primary_valid) || is.null(base_fit) || !is.finite(pi_hat) || pi_hat >= 1e-6)
    return(list(passed = FALSE, checked = FALSE))
  # At pi=0 the count model is exactly the baseline. Check the direction into pi>0,
  # in addition to the multi-start likelihood, rather than trusting a flat logit Hessian.
  small_pi <- 1e-5
  p <- base_fit$parameters; p$gamma_zi <- stats::qlogis(small_pi)
  check <- tryCatch({
    obj <- RTMB::MakeADFun(make_olre_objective(base_fit$data, base_fit$X,
      normal_quadrature(max(base_fit$quadrature_n, zi_fit$quadrature_n)), zi = TRUE),
      p, random = c("b_repeat", "b_ts"), silent = TRUE)
    value <- obj$fn(obj$par)
    score <- obj$gr(obj$par)[names(obj$par) == "gamma_zi"] / (small_pi * (1-small_pi))
    data.frame(pi_hat = pi_hat, positive_pi_nll = zi_fit$nll, exact_pi0_nll = base_fit$nll,
      small_pi = small_pi, nll_at_small_pi = value, nll_score_at_small_pi = unname(score),
      passed = is.finite(score) && score >= 0 && zi_fit$nll >= base_fit$nll - .001 &&
        abs(zi_fit$nll - base_fit$nll) < .001)
  }, error = function(e) data.frame(pi_hat = pi_hat, passed = FALSE, reason = conditionMessage(e)))
  write_result_table(check, output_dir, "01_ZI_pi0_likelihood_boundary_check.csv")
  list(passed = isTRUE(check$passed), checked = TRUE, table = check)
}

pi0_boundary_candidate <- function(baseline, mode, attempts = data.frame(), warnings = character()) {
  fit <- baseline$primary
  if (is.null(fit)) stop("pi=0 边界缺少基线拟合。")
  fit$model_id <- "M2_SKEW_ZI_PI0_BOUNDARY"; fit$zi <- FALSE
  list(primary = fit, primary_valid = isTRUE(baseline$primary_valid), normal = NULL, normal_valid = FALSE,
    reproducible = baseline$reproducible, quadrature_ok = baseline$quadrature_ok,
    direct_integrals_ok = baseline$direct_integrals_ok,
    reasons = if (isTRUE(baseline$primary_valid)) character() else baseline$reasons,
    zi_mode = mode, pi = 0, attempts = attempts, warnings = warnings, boundary_zero = TRUE)
}

append_parameter_notes <- function(row, candidates) {
  fit <- candidates$primary
  if (is.null(fit)) return(row)
  notes <- Filter(nzchar, row$notes)
  if (isTRUE(candidates$boundary_zero)) notes <- c(notes,
    "额外零值概率为 pi=0 边界；该候选按相同的基线模型报告")
  if (is.matrix(fit$cov_fixed)) {
    bad <- which(!is.finite(diag(fit$cov_fixed)) | diag(fit$cov_fixed) <= 0)
    if (length(bad)) notes <- c(notes, paste0("标准误估计无效的参数：",
      paste(names(fit$fit$par)[bad], collapse = ", ")))
  }
  d <- fit$data
  for (host in levels(d$Host_type)) {
    keep <- d$Host_type == host
    if (any(keep) && (all(d$MC_reads[keep] == 0) || all(d$Other_reads[keep] == 0)))
      notes <- c(notes, paste0(host, " 的 ", sum(keep), " 条记录全部为",
        if (all(d$MC_reads[keep] == 0)) "零值" else "比例 1",
        "；该宿主系数处于边界，估计不稳定"))
  }
  row$notes <- paste(unique(notes), collapse = "; ")
  row
}

fit_zero_inflated_candidate <- function(d, X, stub, baseline, output_dir) {
  # With no observed zeros, pi=0 is the exact ML boundary, not an arbitrary cutoff.
  if (!any(d$MC_reads == 0)) {
    return(pi0_boundary_candidate(baseline$candidates, "exact_pi0_no_observed_zeros"))
  }
  warnings_seen <- character(); errors_seen <- character(); attempts <- list(); precisions <- list()
  record_warning <- function(w) {
    warnings_seen <<- unique(c(warnings_seen, conditionMessage(w))); invokeRestart("muffleWarning")
  }
  seed_normal <- baseline$candidates$normal
  seed_skew <- baseline$candidates$primary
  initial <- if (!is.null(seed_normal)) seed_normal$parameters else
    withCallingHandlers(make_initial_parameters(d, X, stub), warning = record_warning)
  skew_initial <- if (!is.null(seed_skew)) seed_skew$parameters else initial
  zero_fraction <- mean(d$MC_reads == 0)
  seed_alpha <- if (is.finite(skew_initial$alpha)) max(-8, min(8, skew_initial$alpha)) else -2
  starts <- data.frame(alpha = c(alpha_starts, seed_alpha, seed_alpha),
    pi = pmax(1e-4, pmin(.8, zero_fraction * c(.5, .5, .5, .5, .2, .8))),
    seed = c(rep("normal", length(alpha_starts)), "skew", "skew"))
  node_n <- max(quadrature_nodes, if (is.null(seed_skew)) quadrature_nodes else seed_skew$quadrature_n)
  run_fit <- function(p) tryCatch(withCallingHandlers(
    fit_rtmb_olre(d, X, p, node_n, skew = TRUE, zi = TRUE), warning = record_warning),
    error = function(e) { errors_seen <<- unique(c(errors_seen, conditionMessage(e))); NULL })
  for (stage in seq_len(3L)) {
    trials <- vector("list", nrow(starts))
    parameters <- lapply(seq_len(nrow(starts)), function(j) {
      p <- if (starts$seed[j] == "skew") skew_initial else initial
      p$alpha <- starts$alpha[j]; p$gamma_zi <- stats::qlogis(starts$pi[j]); p
    })
    pool <- getOption("GLMM_ZI_worker_pool", NULL)
    parallel_trials <- NULL
    if (!is.null(pool)) {
      cat("    ", length(parameters), " 个独立起点使用 ", length(pool), " 个本地进程；节点 ", node_n, "。\n", sep = "")
      parallel_trials <- parallel::parLapply(pool, parameters, zi_worker_fit, d = d, design = X, nodes = node_n)
    }
    for (j in seq_len(nrow(starts))) {
      cat("    ZI 起点 ", j, "/", nrow(starts), "；节点 ", node_n,
        "；alpha=", signif(starts$alpha[j], 3), "；pi=", signif(starts$pi[j], 3), "\n", sep = "")
      if (is.null(parallel_trials)) trial <- run_fit(parameters[[j]]) else {
        trial <- parallel_trials[[j]]$fit
        warnings_seen <- unique(c(warnings_seen, parallel_trials[[j]]$warnings))
        if (nzchar(parallel_trials[[j]]$error)) errors_seen <- unique(c(errors_seen, parallel_trials[[j]]$error))
      }
      trials[j] <- list(trial)
      attempts[[length(attempts) + 1L]] <- data.frame(stage = stage, nodes = node_n,
        alpha_start = starts$alpha[j], pi_start = starts$pi[j],
        nll = if (is.null(trial)) NA_real_ else as.numeric(trial$value),
        alpha = if (is.null(trial)) NA_real_ else trial$parameters$alpha,
        pi = if (is.null(trial)) NA_real_ else stats::plogis(trial$parameters$gamma_zi),
        convergence = if (is.null(trial)) NA_integer_ else trial$opt$convergence,
        pdHess = !is.null(trial) && isTRUE(trial$pdHess),
        max_gradient = if (is.null(trial)) NA_real_ else trial$max_gradient,
        basic_valid = basic_fit_valid(trial))
    }
    selected <- select_best_likelihood(trials)
    best <- selected$fit
    if (is.null(best)) break
    check_n <- max(quadrature_check_nodes, 2L * node_n - 1L)
    precision <- tryCatch(check_quadrature_precision(best, d, X, check_n), error = function(e) {
      errors_seen <<- unique(c(errors_seen, conditionMessage(e)))
      data.frame(fitted_nodes = node_n, check_nodes = check_n, fitted_nll = best$value,
        check_nll = NA_real_, absolute_nll_difference = NA_real_, check_max_gradient = NA_real_, passed = FALSE)
    })
    wrapped <- wrap_rtmb_fit(best, stub, d, X, "M2_SKEW_ZI_AGHQ")
    direct <- tryCatch(check_conditional_integrals(wrapped), error = function(e) {
      errors_seen <<- unique(c(errors_seen, conditionMessage(e)))
      data.frame(passed = FALSE)
    })
    precisions[[length(precisions) + 1L]] <- cbind(stage = stage, precision)
    needs_more <- !isTRUE(precision$passed) || !all(direct$passed %in% TRUE)
    if (!needs_more || node_n >= maximum_fit_nodes) break
    node_n <- min(maximum_fit_nodes, check_n)
    cat("    增加节点至 ", node_n, "，重新核对积分精度。\n", sep = "")
    skew_initial <- best$parameters
    skew_initial$gamma_zi <- NULL
  }
  attempt_table <- bind_flexible_tables(attempts)
  write_result_table(attempt_table, output_dir, "01_ZI_multistart_fit_attempts.csv")
  write_result_table(bind_flexible_tables(precisions), output_dir, "01_ZI_quadrature_precision_checks.csv")
  writeLines(enc2utf8(unique(c(warnings_seen, errors_seen))), file.path(output_dir, "01_ZI_fitting_messages.txt"))
  if (is.null(best)) return(list(primary = NULL, primary_valid = FALSE, normal = NULL,
    normal_valid = FALSE, reproducible = FALSE, quadrature_ok = FALSE, direct_integrals_ok = FALSE,
    reasons = c("所有零膨胀起点均失败", errors_seen), zi_mode = "failed", pi = NA_real_,
    attempts = attempt_table, warnings = warnings_seen, boundary_zero = FALSE))
  write_result_table(direct, output_dir, "01_ZI_independent_conditional_integrals.csv")
  boundary <- check_pi0_likelihood_boundary(wrapped, baseline$candidates, output_dir)
  if (isTRUE(boundary$passed)) {
    cat("    新增概率收敛到 pi=0；核对似然方向后精确退回基线。\n")
    answer <- pi0_boundary_candidate(baseline$candidates, "pi0_verified_likelihood_boundary", attempt_table, warnings_seen)
    answer$boundary_check <- boundary$table
    answer$interior_attempt_before_boundary_check <- wrapped
    return(answer)
  }
  shape_ok <- is.finite(best$parameters$alpha) && abs(best$parameters$alpha) < 19.99
  obs_sd_ok <- best$parameters$log_sd_observation > -7.99 && best$parameters$log_sd_observation < 2.99
  pi_ok <- stats::plogis(best$parameters$gamma_zi) >= 1e-6 &&
    best$parameters$gamma_zi < stats::qlogis(.99) - .01
  valid <- basic_fit_valid(best) && shape_ok && obs_sd_ok && pi_ok &&
    isTRUE(selected$reproducible) && isTRUE(precision$passed) && all(direct$passed %in% TRUE)
  reasons <- c(if (!basic_fit_valid(best)) "优化、Hessian 或协方差检查未通过",
    if (!shape_ok) "alpha 触及边界", if (!obs_sd_ok) "Observation SD 触及边界",
    if (!pi_ok) "额外零值概率触及数值边界",
    if (!isTRUE(selected$reproducible)) "最优解未由至少两个有效起点复现",
    if (!isTRUE(precision$passed)) "积分节点核对未通过",
    if (!all(direct$passed %in% TRUE)) "独立积分核对未通过")
  list(primary = wrapped, primary_valid = valid, normal = NULL, normal_valid = FALSE,
    reproducible = isTRUE(selected$reproducible), quadrature_ok = isTRUE(precision$passed),
    direct_integrals_ok = all(direct$passed %in% TRUE), reasons = reasons,
    zi_mode = "estimated_common_pi", pi = stats::plogis(best$parameters$gamma_zi),
    attempts = attempt_table, warnings = warnings_seen, boundary_zero = FALSE)
}

zero_count_diagnostics <- function(diagnostics, d, output_dir) {
  res <- diagnostics$res_check
  sim <- colSums(res$simulatedResponse == 0)
  test <- DHARMa::testGeneric(res, summary = function(x) sum(x == 0), plot = FALSE)
  overall <- data.frame(observed = sum(d$MC_reads == 0), simulated_mean = mean(sim),
    simulated_95_low = unname(stats::quantile(sim, .025)),
    simulated_95_high = unname(stats::quantile(sim, .975)),
    ratioObsSim = unname(test$statistic), p_raw = test$p.value)
  groups <- list(TS = factor(d$TS), Host = factor(d$Host_type), Reads_quartile =
    factor(paste0("Q", cut(rank(d$Ruegeria_reads_total, ties.method = "average"),
      breaks = stats::quantile(rank(d$Ruegeria_reads_total, ties.method = "average"), seq(0, 1, .25)),
      include.lowest = TRUE, labels = FALSE))))
  rows <- list()
  for (grouping in names(groups)) {
    group <- groups[[grouping]]
    these <- lapply(levels(group), function(g) {
      keep <- group == g
      simulated <- colSums(res$simulatedResponse[keep, , drop = FALSE] == 0)
      t <- DHARMa::testGeneric(res, summary = function(x) sum(x[keep] == 0), plot = FALSE)
      data.frame(grouping = grouping, group = g, n = sum(keep), observed = sum(d$MC_reads[keep] == 0),
        simulated_mean = mean(simulated), simulated_95_low = unname(stats::quantile(simulated, .025)),
        simulated_95_high = unname(stats::quantile(simulated, .975)), p_raw = t$p.value)
    })
    these <- do.call(rbind, these)
    these$p_BH_within_grouping <- adjust_valid_p(these$p_raw)
    rows[[grouping]] <- these
  }
  grouped <- do.call(rbind, rows)
  diagnostics$zero_count_table <- overall
  diagnostics$zero_group_table <- grouped
  diagnostics$zero_count_check <- test
  N <- d$Ruegeria_reads_total
  one_test <- DHARMa::testGeneric(res, summary = function(x) sum(x == N), plot = FALSE)
  keep <- !diagnostics$overall_table$check %in% c("Zero_MC_count", "All_reads_MC_count")
  diagnostics$overall_table <- rbind(diagnostics$overall_table[keep, , drop = FALSE],
    data.frame(check = c("Zero_MC_count", "All_reads_MC_count"),
      p_value = c(test$p.value, one_test$p.value), status = "ok"))
  write_result_table(overall, output_dir, "13_zero_count_summary.csv")
  write_result_table(grouped, output_dir, "14_zero_count_by_group.csv")
  write_result_table(diagnostics$overall_table, output_dir, "04_overall_diagnostics.csv")
  diagnostics
}

new_comparison_row <- function(mc, model, n_source) {
  row <- summary_row(mc, n_source)
  row$model <- model; row$extra_zero_pi <- if (model == "Baseline") 0 else NA_real_
  row$zi_mode <- if (model == "Baseline") "no_ZI" else "not_fitted"
  row$baseline_source <- ""
  row$expected_zero_mean <- row$expected_zero_95_low <- row$expected_zero_95_high <- NA_real_
  row$Zero_TS_min_BH_p <- row$Zero_Host_min_BH_p <- row$Zero_Reads_min_BH_p <- NA_real_
  row$top1_count_variance_fraction <- NA_real_
  row
}

fill_zero_diagnostic_row <- function(row, diagnostics) {
  row <- fill_diagnostic_summary(row, diagnostics)
  z <- diagnostics$zero_count_table
  row$expected_zero_mean <- z$simulated_mean
  row$expected_zero_95_low <- z$simulated_95_low
  row$expected_zero_95_high <- z$simulated_95_high
  mapping <- c(TS = "Zero_TS_min_BH_p", Host = "Zero_Host_min_BH_p", Reads_quartile = "Zero_Reads_min_BH_p")
  for (grouping in names(mapping)) {
    nm <- mapping[[grouping]]
    vals <- diagnostics$zero_group_table$p_BH_within_grouping[diagnostics$zero_group_table$grouping == grouping]
    row[[nm]] <- if (any(is.finite(vals))) min(vals, na.rm = TRUE) else NA_real_
  }
  extra <- c(Zero_TS_min_BH_p = row$Zero_TS_min_BH_p, Zero_Host_min_BH_p = row$Zero_Host_min_BH_p,
    Zero_Reads_min_BH_p = row$Zero_Reads_min_BH_p)
  flags <- names(extra)[is.finite(extra) & extra < .05]
  row$diagnostic_flags <- paste(Filter(nzchar, c(row$diagnostic_flags, flags)), collapse = "; ")
  if (length(flags)) row$status <- "数值有效；残差需复核"
  contributions <- summarise_dispersion_contributions(diagnostics$fitted_rows,
    diagnostics$fitted_rows$expected_MC_proportion)
  row$top1_count_variance_fraction <- contributions$summary$top1_fraction
  row
}

adjust_comparison_p <- function(table, MC_ids) {
  for (model in c("Baseline", "ZI")) {
    index <- which(table$model == model)
    for (term in c("DO", "Host", "Reads")) {
      p <- table[[paste0(term, "_p_raw")]][index]
      p[!table$numerical_valid[index]] <- NA_real_
      table[[paste0(term, "_p_BH")]][index] <- adjust_valid_p(p, n = length(MC_ids))
    }
  }
  table
}

make_comparison_table <- function(table, MC_ids) {
  rows <- lapply(MC_ids, function(mc) {
    b <- table[table$MC_ID == mc & table$model == "Baseline", , drop = FALSE]
    z <- table[table$MC_ID == mc & table$model == "ZI", , drop = FALSE]
    valid_pair <- isTRUE(b$numerical_valid) && isTRUE(z$numerical_valid)
    delta <- if (valid_pair) z$AIC - b$AIC else NA_real_
    decision <- if (!z$numerical_valid) "ZI 数值需复核" else if (!b$numerical_valid) "ZI 有效；基线无效，AIC 比较留空" else
      if (grepl("pi0", z$zi_mode)) "pi=0；与基线相同" else
      if (delta < -2) "ZI 的 AIC 改善" else if (delta > 2) "基线的 AIC 较低" else "AIC 接近"
    data.frame(MC_ID = mc, n_used = z$n_used, observed_zero = z$n_zero_MC,
      baseline_valid = b$numerical_valid, ZI_valid = z$numerical_valid, extra_zero_pi = z$extra_zero_pi,
      baseline_AIC = b$AIC, ZI_AIC = z$AIC, delta_AIC_ZI_minus_baseline = delta,
      baseline_expected_zero = b$expected_zero_mean, ZI_expected_zero = z$expected_zero_mean,
      baseline_zero_p = b$Zero_count_p, ZI_zero_p = z$Zero_count_p,
      baseline_DO_beta = b$DO_beta, baseline_DO_p_BH = b$DO_p_BH,
      ZI_DO_beta = z$DO_beta, ZI_DO_CI_low = z$DO_CI_low, ZI_DO_CI_high = z$DO_CI_high,
      ZI_DO_p_raw = z$DO_p_raw, ZI_DO_p_BH = z$DO_p_BH,
      ZI_Host_p_BH = z$Host_p_BH, ZI_Reads_p_BH = z$Reads_p_BH,
      ZI_Uniformity_p = z$Uniformity_p, ZI_Dispersion_p = z$Dispersion_p,
      ZI_Dispersion_ratio = z$Dispersion_ratio, ZI_Predicted_quantiles_p = z$Predicted_quantiles_p,
      ZI_Reads_quantiles_p = z$Reads_quantiles_p,
      ZI_Reads_location_trend_BH_p = z$Reads_location_trend_BH_p,
      ZI_Reads_deviation_trend_BH_p = z$Reads_deviation_trend_BH_p,
      ZI_Zero_TS_min_BH_p = z$Zero_TS_min_BH_p, ZI_Zero_Host_min_BH_p = z$Zero_Host_min_BH_p,
      ZI_Zero_Reads_min_BH_p = z$Zero_Reads_min_BH_p, ZI_Tail_min_BH_p = z$Tail_min_BH_p,
      ZI_alpha = z$alpha, ZI_Repeat_SD = z$Repeat_SD, ZI_TS_SD = z$TS_SD,
      ZI_Observation_SD = z$Observation_SD, ZI_max_gradient = z$max_gradient,
      baseline_status = b$status, ZI_status = z$status, AIC_comparison = decision,
      baseline_flags = b$diagnostic_flags, ZI_flags = z$diagnostic_flags,
      ZI_notes = z$notes, zi_mode = z$zi_mode, baseline_source = b$baseline_source,
      stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

write_ZI_html <- function(comparison, table, run_dir) {
  view <- data.frame(MC = comparison$MC_ID, 零值数 = comparison$observed_zero,
    基线数值有效 = comparison$baseline_valid, ZI数值有效 = comparison$ZI_valid,
    额外零值概率 = fmt_number(comparison$extra_zero_pi),
    基线AIC = fmt_number(comparison$baseline_AIC, 2), ZI_AIC = fmt_number(comparison$ZI_AIC, 2),
    ΔAIC = fmt_number(comparison$delta_AIC_ZI_minus_baseline, 2),
    基线预测零值 = fmt_number(comparison$baseline_expected_zero, 1),
    ZI预测零值 = fmt_number(comparison$ZI_expected_zero, 1),
    ZI零值p = fmt_p(comparison$ZI_zero_p),
    基线DO_BH = fmt_p(comparison$baseline_DO_p_BH),
    ZI_DO_beta = fmt_number(comparison$ZI_DO_beta), ZI_DO_BH = fmt_p(comparison$ZI_DO_p_BH),
    ZI状态 = comparison$ZI_status, AIC比较 = comparison$AIC_comparison,
    ZI残差提示 = comparison$ZI_flags, ZI说明 = comparison$ZI_notes, check.names = FALSE)
  heads <- paste0("<tr>", paste0("<th>", html_escape(names(view)), "</th>", collapse = ""), "</tr>")
  rows <- vapply(seq_len(nrow(view)), function(i) paste0("<tr>", paste0("<td>",
    vapply(view, function(v) html_escape(if (is.na(v[i])) "—" else as.character(v[i])), character(1)),
    "</td>", collapse = ""), "</tr>"), character(1))
  html <- c('<!doctype html><html lang="zh-CN"><meta charset="utf-8"><title>全部 MC 零膨胀比较</title>',
    '<style>body{font:15px/1.5 Arial,"Microsoft YaHei",sans-serif;margin:24px;color:#20303e}table{border-collapse:collapse}th,td{padding:9px 12px;border:1px solid #d7e0e8;white-space:nowrap;text-align:right}th{background:#193c55;color:white;position:sticky;top:0}td:first-child{font-weight:bold;position:sticky;left:0;background:#edf3f7}tr:nth-child(even){background:#f5f8fa}.scroll{overflow:auto;max-height:75vh}td:nth-last-child(-n+4){white-space:normal;min-width:200px;text-align:left}</style>',
    '<h1>全部 MC：共同零膨胀概率候选</h1>',
    '<p>基线与 ZI 均使用中心化偏正态 Observation 效应。ZI 仅新增 logit(pi)=gamma0；其余比例方程和数据一致。ΔAIC=ZI−基线，负值偏向 ZI。每个 MC 独立估计参数。</p>',
    paste0('<p>ZI 数值核对通过 ', sum(comparison$ZI_valid), '/', nrow(comparison), '。pi=0 边界候选按同一个基线模型报告；相应 AIC 不另计零值参数。诊断 ok 表示计算完成。</p>'),
    '<div class="scroll"><table><thead>', heads, '</thead><tbody>', rows, '</tbody></table></div>',
    '<p>DO、Host、reads 的 BH 分别跨全部计划 MC 校正，基线和 ZI 分开。DO β 属于计数部分；总体预期比例还乘以 (1−pi)。无效拟合的推断留空。AIC 只在同一 MC、数值有效的候选间比较。</p>',
    '<p>完整比较：<a href="00_all_MC_ZI_comparison.csv">00_all_MC_ZI_comparison.csv</a>；两候选全部数值：<a href="01_all_model_metrics.csv">01_all_model_metrics.csv</a>。</p>',
    '<p>本轮不做逐 TS 剔除，也不计算被剔除批次的预测评分。pi=0 在边界，不报告普通卡方 LRT 或零值部分 Wald p。</p></html>')
  writeLines(enc2utf8(html), file.path(run_dir, "00_all_MC_ZI_comparison.html"))
}

write_ZI_outputs <- function(table, MC_ids, run_dir, completed) {
  table <- adjust_comparison_p(table, MC_ids)
  comparison <- make_comparison_table(table, MC_ids)
  write_excel_csv(table, file.path(run_dir, "01_all_model_metrics.csv"))
  write_excel_csv(comparison, file.path(run_dir, "00_all_MC_ZI_comparison.csv"))
  write_ZI_html(comparison, table, run_dir)
  txt <- capture.output({
    cat("============================================================\n全部 MC：偏正态基线 vs 共同零膨胀概率候选\n")
    cat("已完成：", completed, "/", length(MC_ids), "；ZI 数值有效：", sum(comparison$ZI_valid), "。\n", sep = "")
    cat("1. 数值状态、额外零值概率和 AIC\n")
    print(comparison[, c("MC_ID", "baseline_valid", "ZI_valid", "extra_zero_pi", "delta_AIC_ZI_minus_baseline", "AIC_comparison")], row.names = FALSE)
    cat("\n2. 零值预测（模型模拟 ", n_sim, " 次；实际零值全部保留）\n", sep = "")
    print(comparison[, c("MC_ID", "observed_zero", "baseline_expected_zero", "ZI_expected_zero", "baseline_zero_p", "ZI_zero_p")], row.names = FALSE)
    cat("\n3. ZI 候选 DO 结果（BH 跨全部计划 MC；95% CI 为逐项 Wald）\n")
    print(comparison[, c("MC_ID", "ZI_DO_beta", "ZI_DO_CI_low", "ZI_DO_CI_high", "ZI_DO_p_raw", "ZI_DO_p_BH")], row.names = FALSE)
    cat("\n4. ZI 候选宿主、reads 整体关联\n")
    print(comparison[, c("MC_ID", "ZI_Host_p_BH", "ZI_Reads_p_BH")], row.names = FALSE)
    cat("\n5. ZI 剩余问题\n")
    print(comparison[, c("MC_ID", "ZI_status", "ZI_flags", "ZI_notes")], row.names = FALSE)
    cat("\n无零值 MC 采用精确 pi=0 边界；极小 pi 候选另核查边界似然方向。\n")
    cat("确认 pi=0 后按同一个基线模型报告，AIC 不另计零值参数。\n")
    cat("没有报告普通卡方 LRT / 零值部分 Wald p。本轮未逐 TS 剔除。\n")
    cat("基线缓存仅在所有记录、协变量、分母和样条设计核对一致后复用，并重新核查积分。\n")
    cat("优先阅读：00_all_MC_ZI_comparison.html；00_all_MC_ZI_comparison.csv；00_key_results.txt。\n")
    cat("结果目录：", normalizePath(run_dir, winslash = "/", mustWork = FALSE), "\n", sep = "")
  })
  writeLines(enc2utf8(txt), file.path(run_dir, "00_key_results.txt"))
  invisible(list(table = table, comparison = comparison, text = txt))
}

run_all_MC_zero_inflation <- function() {
  previous_options <- options(contrasts = c("contr.treatment", "contr.poly"))
  on.exit(options(previous_options), add = TRUE)
  packages <- c("readxl", "glmmTMB", "RTMB", "DHARMa", "qgam", "ggplot2")
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("缺少包：", paste(missing, collapse = ", "))
  workers <- max(1L, min(as.integer(parallel_start_workers), parallel::detectCores()))
  if (workers > 1L) {
    cat("启动 ", workers, " 个本地计算进程，用于独立起点。\n", sep = "")
    pool <- parallel::makePSOCKcluster(workers)
    on.exit(parallel::stopCluster(pool), add = TRUE)
    worker_source <- getOption("GLMM_ZI_worker_source", file.path(project_dir, "scripts", "GLMM_reads_diagnostics_full.R"))
    parallel::clusterCall(pool, function(path) {
      options(GLMM_ZI_definitions_only = TRUE)
      Sys.setenv(OMP_NUM_THREADS = 1, OPENBLAS_NUM_THREADS = 1, MKL_NUM_THREADS = 1)
      source(path, encoding = "UTF-8")
      requireNamespace("RTMB", quietly = TRUE)
      TRUE
    }, worker_source)
    options(GLMM_ZI_worker_pool = pool)
    on.exit(options(GLMM_ZI_worker_pool = NULL), add = TRUE)
  }
  raw <- readxl::read_excel(input_file, sheet = input_sheet)
  MC_ids <- natural_MC_order(unique(as.character(raw$MC_ID[!is.na(raw$MC_ID)])))
  if (!is.null(MC_selection)) {
    if (any(!MC_selection %in% MC_ids)) stop("指定的 MC 不存在。")
    MC_ids <- natural_MC_order(unique(MC_selection))
  }
  shared <- build_shared_design(raw, MC_ids, host_reference)
  resuming <- nzchar(resume_run_dir)
  run_dir <- if (resuming) resume_run_dir else
    file.path(output_root, paste0("GLMM_all_MC_ZI_", format(Sys.time(), "%Y%m%d_%H%M%S")))
  if (!resuming && dir.exists(run_dir)) run_dir <- paste0(run_dir, "_", Sys.getpid())
  dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(shared, file.path(run_dir, "01_shared_model_design.rds"))
  write_excel_csv(shared$samples, file.path(run_dir, "01_shared_sample_design.csv"))
  writeLines(capture.output(utils::sessionInfo()), file.path(run_dir, "sessionInfo.txt"))
  settings <- list(input_file = input_file, input_md5 = unname(tools::md5sum(input_file)), MC_ids = MC_ids,
    n_sim = n_sim, seed = simulation_seed, nodes = quadrature_nodes, maximum_fit_nodes = maximum_fit_nodes,
    mean_df = 3L, zi_formula = "~1", baseline_cache_dir = baseline_cache_dir,
    cache_used_only_after_record_and_design_checks = TRUE, batch_sensitivity = FALSE)
  if (resuming) {
    previous <- readRDS(file.path(run_dir, "01_run_settings.rds"))
    for (nm in c("input_md5", "MC_ids", "n_sim", "seed", "mean_df", "zi_formula")) {
      if (!identical(previous[[nm]], settings[[nm]])) stop("恢复运行设置不一致：", nm)
    }
  }
  saveRDS(settings, file.path(run_dir, "01_run_settings.rds"))
  rows <- list(); model_files <- list(); all_zeros <- list(); all_overall <- list(); coefficients <- list()
  # First confirm the no-zero boundary and the previously most problematic valid fit.
  # Output tables retain the natural MC order regardless of fitting order.
  fit_order <- c(intersect(c("MC10", "MC43"), MC_ids), setdiff(MC_ids, c("MC10", "MC43")))
  for (i in seq_along(fit_order)) {
    mc <- fit_order[i]
    previous_model <- file.path(run_dir, mc, "fitted_comparison.rds")
    if (resuming && file.exists(previous_model)) {
      saved <- readRDS(previous_model)
      if (!is.null(saved$metrics) && nrow(saved$metrics) == 2L) {
        rows[[mc]] <- saved$metrics
        model_files[[mc]] <- previous_model
        for (model in names(saved$diagnostic_files)) {
          diag <- readRDS(saved$diagnostic_files[[model]])
          if (is.null(diag$zero_group_table)) diag <- zero_count_diagnostics(diag, diag$fitted_rows,
            file.path(run_dir, mc, paste0(model, "_diagnostics")))
          all_zeros[[paste(mc, model)]] <- cbind(MC_ID = mc, model = model, diag$zero_group_table)
          all_overall[[paste(mc, model)]] <- cbind(MC_ID = mc, model = model, diag$overall_table)
          coefficients[[paste(mc, model)]] <- cbind(model = model,
            fixed_coefficient_table(saved$candidates[[model]]$primary, mc))
        }
        cat("恢复已完成的 ", mc, "。\n", sep = "")
        next
      }
    }
    cat("\n========== ", i, "/", length(MC_ids), "：", mc, " ==========\n", sep = "")
    mc_dir <- file.path(run_dir, mc); dir.create(mc_dir, recursive = TRUE, showWarnings = FALSE)
    n_source <- sum(as.character(raw$MC_ID) == mc, na.rm = TRUE)
    pair_rows <- list(Baseline = new_comparison_row(mc, "Baseline", n_source),
      ZI = new_comparison_row(mc, "ZI", n_source))
    captured <- tryCatch({
      prepared <- prepare_mc_shared(raw, mc, shared); d <- prepared$data_test; X <- prepared$X
      stub <- make_model_stub(d)
      write_result_table(prepared$sample_audit, mc_dir, "01_sample_audit.csv")
      write_result_table(d, mc_dir, "01_model_input_rows.csv")
      for (model in names(pair_rows)) {
        r <- pair_rows[[model]]; r$n_used <- nrow(d); r$n_excluded <- n_source - nrow(d)
        r$n_repeat_groups <- nlevels(d$Repeat_group_ID); r$n_TS <- nlevels(d$TS)
        r$n_zero_MC <- sum(d$MC_reads == 0); r$n_all_reads_MC <- sum(d$Other_reads == 0)
        r$observed_mean_proportion <- mean(d$MC_reads / d$Ruegeria_reads_total)
        pair_rows[[model]] <- r
      }
      cat("  纳入 ", nrow(d), "；零值 ", sum(d$MC_reads == 0), "。\n", sep = "")
      baseline <- load_verified_baseline(mc, d, X, mc_dir)
      if (is.null(baseline)) {
        cat("  重新拟合偏正态基线。\n")
        base_dir <- file.path(mc_dir, "baseline_fit"); dir.create(base_dir)
        candidates <- fit_shared_candidates(d, X, stub, quadrature_nodes,
          quadrature_check_nodes, alpha_starts, base_dir)
        baseline <- list(candidates = candidates, source = "refitted", diagnostics = NULL,
          diagnostic_source = "")
      } else cat("  基线样本及设计一致；缓存似然与积分复核完成。\n")
      zi <- fit_zero_inflated_candidate(d, X, stub, baseline, mc_dir)
      candidates <- list(Baseline = baseline$candidates, ZI = zi)
      diag_paths <- list()
      saveRDS(list(MC_ID = mc, candidates = candidates, shared_design = shared,
        baseline_source = baseline$source), file.path(mc_dir, "fitted_comparison.rds"))
      for (model in names(candidates)) {
        r <- pair_rows[[model]]
        r <- append_parameter_notes(fill_fit_summary(r, candidates[[model]]), candidates[[model]])
        r$baseline_source <- baseline$source
        r$extra_zero_pi <- if (model == "ZI") zi$pi else 0
        r$zi_mode <- if (model == "ZI") zi$zi_mode else "no_ZI"
        if (r$numerical_valid) {
          diagnostic_dir <- file.path(mc_dir, paste0(model, "_diagnostics"))
          dir.create(diagnostic_dir, recursive = TRUE, showWarnings = FALSE)
          diag <- if (model == "Baseline") baseline$diagnostics else NULL
          if (model == "ZI" && isTRUE(zi$boundary_zero)) diag <- baseline$diagnostics
          if (is.null(diag)) {
            cat("  ", model, "：残差模拟及诊断（", n_sim, " 次）。\n", sep = "")
            diag <- tryCatch(diagnose_skew_comparison_model(paste0(mc, "_", model),
              candidates[[model]]$primary, d, n_sim, simulation_seed + i * 100L,
              diagnostic_dir), error = function(e) e)
          }
          if (inherits(diag, "error")) {
            r$status <- "数值有效；诊断失败"; r$diagnostic_status <- "failed"
            r$notes <- paste(Filter(nzchar, c(r$notes, conditionMessage(diag))), collapse = "; ")
          } else {
            diag <- zero_count_diagnostics(diag, d, diagnostic_dir)
            r <- fill_zero_diagnostic_row(r, diag)
            # Baseline simulation matrices are already saved in the verified previous run.
            if (model == "Baseline" && nzchar(baseline$diagnostic_source)) {
              writeLines(normalizePath(baseline$diagnostic_source, winslash = "/"),
                file.path(diagnostic_dir, "00_cached_diagnostic_source.txt"))
              diag_paths[[model]] <- baseline$diagnostic_source
            } else {
              saveRDS(diag, file.path(diagnostic_dir, "model_analysis_objects.rds"))
              diag_paths[[model]] <- file.path(diagnostic_dir, "model_analysis_objects.rds")
            }
            write_result_table(diag$tail_table, diagnostic_dir, "05_distribution_tail_checks.csv")
            all_zeros[[paste(mc, model)]] <- cbind(MC_ID = mc, model = model, diag$zero_group_table)
            all_overall[[paste(mc, model)]] <- cbind(MC_ID = mc, model = model, diag$overall_table)
            coefficients[[paste(mc, model)]] <- cbind(model = model,
              fixed_coefficient_table(candidates[[model]]$primary, mc))
          }
        }
        pair_rows[[model]] <- r
      }
      saveRDS(list(MC_ID = mc, candidates = candidates, metrics = do.call(rbind, pair_rows),
        diagnostic_files = diag_paths, baseline_source = baseline$source),
        file.path(mc_dir, "fitted_comparison.rds"))
      model_files[[mc]] <- file.path(mc_dir, "fitted_comparison.rds")
      TRUE
    }, error = function(e) {
      writeLines(enc2utf8(conditionMessage(e)), file.path(mc_dir, "00_error.txt"))
      for (model in names(pair_rows)) {
        r <- pair_rows[[model]]
        if (!r$numerical_valid) r$status <- "处理失败"
        r$notes <- paste(Filter(nzchar, c(r$notes, conditionMessage(e))), collapse = "; ")
        pair_rows[[model]] <<- r
      }
      cat("  错误已保留：", conditionMessage(e), "\n", sep = "")
      FALSE
    })
    rows[[mc]] <- do.call(rbind, pair_rows)
    table <- do.call(rbind, rows)
    # Include not-yet-run MC rows so BH always uses the same planned family.
    pending <- setdiff(MC_ids, names(rows))
    for (name in pending) table <- rbind(table,
      new_comparison_row(name, "Baseline", sum(as.character(raw$MC_ID) == name, na.rm = TRUE)),
      new_comparison_row(name, "ZI", sum(as.character(raw$MC_ID) == name, na.rm = TRUE)))
    final <- write_ZI_outputs(table, MC_ids, run_dir, i)
    saveRDS(list(summary_table = final$table, comparison = final$comparison, model_files = model_files,
      MC_ids = MC_ids, shared_design = shared, output_dir = run_dir, completed = i, settings = settings),
      file.path(run_dir, "GLMM_all_MC_ZI_results.rds"))
    print(final$comparison[final$comparison$MC_ID == mc,
      c("MC_ID", "ZI_valid", "extra_zero_pi", "delta_AIC_ZI_minus_baseline", "ZI_status")], row.names = FALSE)
    gc(FALSE)
  }
  final <- write_ZI_outputs(do.call(rbind, rows), MC_ids, run_dir, length(rows))
  write_excel_csv(bind_flexible_tables(all_zeros), file.path(run_dir, "02_all_zero_group_diagnostics.csv"))
  write_excel_csv(bind_flexible_tables(all_overall), file.path(run_dir, "03_all_overall_diagnostics.csv"))
  write_excel_csv(bind_flexible_tables(coefficients), file.path(run_dir, "04_all_fixed_coefficients.csv"))
  writeLines(enc2utf8(c(
    "All MCs: centered skew-normal binomial mixed model, baseline versus one common extra-zero probability per MC.",
    "Mean: DO_c + Host + natural spline df=3 of centered log10 Ruegeria reads; Gaussian Repeat/TS intercepts.",
    "ZI: logit(pi)=gamma0. Count Observation effects are integrated before mixing in a point mass at zero.",
    "P(Y=0)=pi+(1-pi)*P_count(Y=0); P(Y=y>0)=(1-pi)*P_count(Y=y).",
    "With no observed zeros, the exact ML estimate is pi=0 and the candidate equals the baseline with no extra parameter.",
    "For near-zero pi fits, the likelihood direction at small positive pi and the multi-start likelihood are checked before reducing to pi=0.",
    "Confirmed pi=0 candidates are reported as the same baseline model, with the baseline parameter count and AIC.",
    "All New_ID records and zero counts retained; missing covariates / non-positive denominators excluded and audited.",
    "Existing baseline fits reused only after exact record/design matching and fresh likelihood/integration checks.",
    "New fits use multiple alpha/pi starts, increased quadrature nodes if required, and independent integration checks.",
    "Inference blank for numerically invalid models. AIC compared only within a numerically valid same-MC pair.",
    "DO, Host and reads BH families each include all planned MCs, separately for baseline and ZI.",
    "DO is the count-component conditional logit slope; overall mean proportion=(1-pi)*E_count(p).",
    "Diagnostic simulations condition on fitted Repeat/TS modes, redraw Observation effects, counts and the extra-zero indicator.",
    "No ordinary chi-square LRT or ZI Wald test at pi=0. No leave-one-TS or held-out prediction checks this run.",
    "Zero-group tests use separate within-grouping BH families (TS, Host, reads quartile).",
    "https://journal.r-project.org/articles/RJ-2017-066/")), file.path(run_dir, "METHODS.txt"))
  if (!identical(unname(tools::md5sum(input_file)), settings$input_md5)) warning("输入文件在分析期间发生变化。")
  cat("\n", paste(final$text, collapse = "\n"), "\n", sep = "")
  invisible(list(summary_table = final$table, comparison = final$comparison,
    model_files = model_files, output_dir = run_dir, shared_design = shared, settings = settings))
}


# ===== Current experiment: a fixed four-candidate comparison for every MC =====
# A: constant extra-zero probability, constant Observation SD.
# B: logit(pi_i)=gamma0+gamma_reads*x_i, constant Observation SD.
# C: constant pi, log(SD_i)=log_sd_observation+sd_reads*x_i.
# D: both reads relationships; x_i is centered log10 total Ruegeria reads.
# Mean, df=3, data, denominator and skew-normal family are shared.
joint_cache_dir <- file.path(project_dir, "results", "GLMM_all_MC_ZI_20261011_232226")
joint_resume_dir <- Sys.getenv("GLMM_JOINT_RESUME_DIR", unset = "")
joint_specs <- list(
  A = list(id="A", zi=TRUE, zi_reads=FALSE, sd_reads=FALSE),
  B = list(id="B", zi=TRUE, zi_reads=TRUE, sd_reads=FALSE),
  C = list(id="C", zi=TRUE, zi_reads=FALSE, sd_reads=TRUE),
  D = list(id="D", zi=TRUE, zi_reads=TRUE, sd_reads=TRUE))

joint_zero_reads_separation <- function(d) {
  zero <- d$log_reads_c[d$MC_reads==0];positive <- d$log_reads_c[d$MC_reads>0]
  if(!length(zero)||!length(positive))return("")
  if(max(zero)<min(positive))return("零值的reads全部低于非零值；零值logit斜率趋于负无穷，无有限最大似然估计")
  if(min(zero)>max(positive))return("零值的reads全部高于非零值；零值logit斜率趋于正无穷，无有限最大似然估计")
  ""
}

# Fixed-coefficient penalty: independent of random effects, hence it can be
# subtracted exactly from the Laplace objective to recover the data-only nll.
joint_penalty_value <- function(parameters, spec) {
  lambda <- spec$lambda_zi_reads
  if (is.null(lambda)) lambda <- 0
  if (!isTRUE(spec$zi) || !isTRUE(spec$zi_reads) || lambda == 0) return(0)
  lambda * parameters$zi_reads^2
}

joint_objective <- function(d, X, nodes, spec) {
  y <- as.numeric(d$MC_reads); N <- as.numeric(d$Ruegeria_reads_total)
  x <- as.numeric(d$log_reads_c)
  ir <- as.integer(d$Repeat_group_ID); it <- as.integer(d$TS)
  observed_logit <- qlogis((y+.5)/(N+1))
  gh <- normal_quadrature(nodes)
  function(p) {
    eta <- as.vector(X %*% p$beta) + p$b_repeat[ir] + p$b_ts[it]
    nll <- 0
    if (!isTRUE(spec$zero_repeat)) nll <- nll - sum(RTMB::dnorm(p$b_repeat,0,exp(p$log_sd_repeat),log=TRUE))
    if (!isTRUE(spec$zero_ts)) nll <- nll - sum(RTMB::dnorm(p$b_ts,0,exp(p$log_sd_ts),log=TRUE))
    so <- exp(p$log_sd_observation + if (isTRUE(spec$sd_reads)) p$sd_reads*x else 0)
    delta <- p$alpha/sqrt(1+p$alpha^2)
    omega <- so/sqrt(1-2*delta^2/pi)
    xi <- -omega*delta*sqrt(2/pi)
    eta0 <- eta+xi
    zmode <- 10*tanh(((observed_logit-eta0)/omega)/10)
    for (iteration in 1:30) {
      probability <- RTMB::plogis(eta0+omega*zmode)
      az <- p$alpha*zmode
      mills <- exp(RTMB::dnorm(az,log=TRUE)-RTMB::pnorm(az,log.p=TRUE))
      curvature <- 1+omega^2*N*probability*(1-probability)+p$alpha^2*mills*(az+mills)
      step <- (omega*(y-N*probability)-zmode+p$alpha*mills)/curvature
      zmode <- zmode+3*tanh(step/3)
    }
    probability <- RTMB::plogis(eta0+omega*zmode); az <- p$alpha*zmode
    mills <- exp(RTMB::dnorm(az,log=TRUE)-RTMB::pnorm(az,log.p=TRUE))
    zsd <- 1/sqrt(1+omega^2*N*probability*(1-probability)+p$alpha^2*mills*(az+mills))
    for (j in seq_along(gh$nodes)) {
      z <- zmode+zsd*gh$nodes[j]
      term <- log(gh$weights[j])+RTMB::dbinom_robust(y,N,eta0+omega*z,log=TRUE)+
        RTMB::dnorm(z,0,1,log=TRUE)+log(2)+RTMB::pnorm(p$alpha*z,log.p=TRUE)-
        RTMB::dnorm(z,zmode,zsd,log=TRUE)
      if (j==1L) ll <- term else ll <- RTMB::logspace_add(ll,term)
    }
    if (isTRUE(spec$zi)) {
      g <- p$gamma_zi + if (isTRUE(spec$zi_reads)) p$zi_reads*x else rep(0,length(x))
      lp <- -RTMB::logspace_add(0,-g); lq <- -RTMB::logspace_add(0,g)
      positive <- which(y>0); zero <- which(y==0)
      ll[positive] <- lq[positive]+ll[positive]
      ll[zero] <- RTMB::logspace_add(lp[zero],lq[zero]+ll[zero])
    }
    nll-sum(ll)+joint_penalty_value(p,spec)
  }
}

joint_map <- function(p,spec) {
  m <- list()
  for (name in c("repeat","ts")) if (isTRUE(spec[[paste0("zero_",name)]])) {
    m[[paste0("log_sd_",name)]] <- factor(NA)
    m[[paste0("b_",name)]] <- factor(rep(NA,length(p[[paste0("b_",name)]])))
  }
  m
}
joint_random <- function(spec) c(if(!isTRUE(spec$zero_repeat)) "b_repeat", if(!isTRUE(spec$zero_ts)) "b_ts")
joint_parameters <- function(p,d,spec) {
  keep <- c("beta","log_sd_repeat","log_sd_ts","log_sd_observation","alpha","b_repeat","b_ts")
  p <- p[c(keep,intersect(c("gamma_zi","zi_reads","sd_reads"),names(p)))]
  if (isTRUE(spec$zi)) {
    if (is.null(p$gamma_zi)) p$gamma_zi <- qlogis(max(1e-4,min(.5,mean(d$MC_reads==0)/2)))
    if (isTRUE(spec$zi_reads) && is.null(p$zi_reads)) p$zi_reads <- 0
    if (!isTRUE(spec$zi_reads)) p$zi_reads <- NULL
  } else {p$gamma_zi <- NULL; p$zi_reads <- NULL}
  if (isTRUE(spec$sd_reads) && is.null(p$sd_reads)) p$sd_reads <- 0
  if (!isTRUE(spec$sd_reads)) p$sd_reads <- NULL
  for (name in c("repeat","ts")) if (isTRUE(spec[[paste0("zero_",name)]])) {
    p[[paste0("log_sd_",name)]] <- -Inf
    p[[paste0("b_",name)]][] <- 0
  }
  p
}
# Only these two sdreport fields are consumed by basic_fit_valid / joint_wrap.
# This uses optimHess and one Cholesky factorization/inverse,
# without a second AD report tape or unused random-effect standard errors.
joint_fixed_covariance <- function(obj, par.fixed) {
  hessian.fixed <- stats::optimHess(par.fixed, obj$fn, obj$gr)
  factor <- try(chol(hessian.fixed), silent = TRUE)
  pdHess <- !inherits(factor, "try-error")
  covariance <- if(pdHess)chol2inv(factor) else hessian.fixed * NaN
  dimnames(covariance) <- dimnames(hessian.fixed)
  list(pdHess = pdHess, cov.fixed = covariance)
}

joint_fit <- function(p,d,X,nodes,spec) {
  p <- joint_parameters(p,d,spec)
  obj <- RTMB::MakeADFun(joint_objective(d,X,nodes,spec),p,
    random=joint_random(spec),map=joint_map(p,spec),silent=TRUE)
  lower <- rep(-Inf,length(obj$par)); upper <- rep(Inf,length(obj$par))
  names(lower) <- names(upper) <- names(obj$par)
  lower[grepl("^log_sd",names(lower))] <- -8; upper[grepl("^log_sd",names(upper))] <- 3
  lower[names(lower)=="alpha"] <- -20; upper[names(upper)=="alpha"] <- 20
  lower[names(lower)=="gamma_zi"] <- -18; upper[names(upper)=="gamma_zi"] <- qlogis(.99)
  lower[names(lower)%in%c("zi_reads","sd_reads","zi_DO")] <- -10
  upper[names(upper)%in%c("zi_reads","sd_reads","zi_DO")] <- 10
  lower[names(lower)=="sd_reads2"] <- -4
  upper[names(upper)=="sd_reads2"] <- 4
  lower[names(lower)=="sd_mean"] <- -3
  upper[names(upper)=="sd_mean"] <- 3
  lower[names(lower)=="alpha_mean"] <- -5
  upper[names(upper)=="alpha_mean"] <- 5
  opt <- nlminb(obj$par,obj$fn,obj$gr,lower=lower,upper=upper,
    control=list(iter.max=600,eval.max=1200,rel.tol=1e-9,x.tol=1e-8))
  path <- "nlminb"
  if (opt$convergence!=0L || max(abs(obj$gr(opt$par)))>.005) {
    retry <- optim(opt$par,obj$fn,obj$gr,method="L-BFGS-B",lower=lower,upper=upper,
      control=list(maxit=800,factr=1e6,pgtol=1e-6))
    if (is.finite(retry$value) && retry$value<=opt$objective+1e-6 && retry$convergence==0L &&
        max(abs(obj$gr(retry$par)))<.01) {
      opt <- list(par=retry$par,objective=retry$value,convergence=retry$convergence,message=retry$message)
      path <- "nlminb -> L-BFGS-B"
    }
  }
  value <- obj$fn(opt$par); gradient <- obj$gr(opt$par)
  parameters <- obj$env$parList(opt$par)
  penalty_value <- as.numeric(joint_penalty_value(parameters,spec))
  # A failed optimizer/gradient cannot pass basic_fit_valid. Keep its objective
  # for the multistart audit, but skip its otherwise wasted numerical Hessian.
  sdr <- if(isTRUE(opt$convergence==0L) && is.finite(value) &&
    all(is.finite(gradient)) && max(abs(gradient))<.02)
      tryCatch(joint_fixed_covariance(obj,par.fixed=opt$par),error=function(e)e) else
      list(pdHess=FALSE,cov.fixed=matrix(NA_real_,length(opt$par),length(opt$par)))
  ans <- list(opt=opt,parameters=parameters,value=as.numeric(value),max_gradient=max(abs(gradient)),
    sdr=sdr,pdHess=!inherits(sdr,"error")&&isTRUE(sdr$pdHess),quadrature_n=nodes,
    skew=TRUE,zi=isTRUE(spec$zi),k=length(opt$par),optimizer_path=path,spec=spec,
    penalty_value=penalty_value,data_nll=as.numeric(value)-penalty_value)
  if(isTRUE(spec$sd_mean)||isTRUE(spec$alpha_mean))ans$abundance_X <- X
  # Only numerical estimates and covariance cross process boundaries.
  if(!inherits(sdr,"error")) ans$sdr <- list(pdHess=sdr$pdHess,cov.fixed=sdr$cov.fixed)
  ans
}
joint_worker <- function(p,d,design,nodes,spec) {
  started <- proc.time()[["elapsed"]]
  warnings <- character()
  f <- tryCatch(withCallingHandlers(joint_fit(p,d,design,nodes,spec),warning=function(w){
    warnings <<- unique(c(warnings,conditionMessage(w))); invokeRestart("muffleWarning")}),error=function(e)e)
  elapsed <- proc.time()[["elapsed"]] - started
  # The tape is local to joint_fit and is never returned. Release it before the
  # next serial start so eight workers do not retain several old tapes each.
  gc(FALSE)
  if(inherits(f,"error")) list(fit=NULL,error=conditionMessage(f),warnings=warnings,elapsed_seconds=elapsed)
  else list(fit=f,error="",warnings=warnings,elapsed_seconds=elapsed)
}
joint_wrap <- function(result,stub,d,X) {
  fit <- wrap_rtmb_fit(result,stub,d,X,paste0("JOINT_",result$spec$id))
  fit$joint_spec <- result$spec
  fit$penalty_requested <- !is.null(result$spec$lambda_zi_reads) && result$spec$lambda_zi_reads > 0
  fit$penalty_value <- as.numeric(joint_penalty_value(result$parameters,result$spec))
  fit$objective_value <- as.numeric(result$value)
  fit$nll <- fit$objective_value-fit$penalty_value
  fit
}

joint_total_objective <- function(fit) {
  if (!is.null(fit$objective_value)) fit$objective_value else fit$nll
}
joint_row_values <- function(fit,x=fit$data$log_reads_c) {
  p <- fit$parameters
  spec <- fit$joint_spec
  if(is.null(spec)) spec <- list(zi=isTRUE(fit$zi),zi_reads=FALSE,sd_reads=FALSE)
  log_sd <- rep(p$log_sd_observation,length(x))
  if(isTRUE(spec$sd_reads)) log_sd <- log_sd+p$sd_reads*x
  pi <- if(isTRUE(spec$zi)) plogis(p$gamma_zi+if(isTRUE(spec$zi_reads))p$zi_reads*x else rep(0,length(x))) else rep(0,length(x))
  list(sd=exp(log_sd),pi=pi)
}
joint_precision <- function(result,d,X,check_n) {
  p <- result$parameters; spec <- result$spec
  obj <- RTMB::MakeADFun(joint_objective(d,X,check_n,spec),p,
    random=joint_random(spec),map=joint_map(p,spec),silent=TRUE)
  val <- obj$fn(result$opt$par); gr <- obj$gr(result$opt$par)
  data.frame(fitted_nodes=result$quadrature_n,check_nodes=check_n,
    fitted_nll=result$value,check_nll=val,absolute_nll_difference=abs(val-result$value),
    check_max_gradient=max(abs(gr)),passed=is.finite(val)&&abs(val-result$value)<.01&&max(abs(gr))<.02)
}
joint_integrals <- function(fit,rows=NULL,nodes=fit$quadrature_n) {
  d <- fit$data; y <- d$MC_reads; N <- d$Ruegeria_reads_total
  if(is.null(rows)) rows <- unique(c(order(N)[1:3],order(N,decreasing=TRUE)[1:3],
    order(y/N)[1:4],order(y/N,decreasing=TRUE)[1:5],which(y==0)[seq_len(min(4,sum(y==0)))]))
  rows <- rows[is.finite(rows)]
  values <- joint_row_values(fit); alpha <- fit$parameters$alpha
  delta <- alpha/sqrt(1+alpha^2); gh <- normal_quadrature(nodes)
  do.call(rbind,lapply(rows,function(i){
    omega <- values$sd[i]/sqrt(1-2*delta^2/pi)
    xi <- -omega*delta*sqrt(2/pi); eta <- fit$eta_base[i]+xi
    z <- 10*tanh(((qlogis((y[i]+.5)/(N[i]+1))-eta)/omega)/10)
    for(j in 1:30){
      pr <- plogis(eta+omega*z); az <- alpha*z
      mills <- exp(dnorm(az,log=TRUE)-pnorm(az,log.p=TRUE))
      step <- (omega*(y[i]-N[i]*pr)-z+alpha*mills)/(1+omega^2*N[i]*pr*(1-pr)+alpha^2*mills*(az+mills))
      z <- z+3*tanh(step/3)
    }
    pr <- plogis(eta+omega*z); az <- alpha*z
    mills <- exp(dnorm(az,log=TRUE)-pnorm(az,log.p=TRUE))
    zs <- 1/sqrt(1+omega^2*N[i]*pr*(1-pr)+alpha^2*mills*(az+mills))
    kernel <- function(a) dbinom(y[i],N[i],plogis(eta+omega*a),log=TRUE)+
      dnorm(a,log=TRUE)+log(2)+pnorm(alpha*a,log.p=TRUE)
    zz <- z+zs*gh$nodes
    terms <- log(gh$weights)+kernel(zz)-dnorm(zz,z,zs,log=TRUE)
    m <- max(terms); count_gh <- m+log(sum(exp(terms-m)))
    anchor <- kernel(z)
    integral <- integrate(function(a)exp(kernel(z+zs*a)-anchor)*zs,-Inf,Inf,rel.tol=1e-9,subdivisions=400L)
    count_direct <- anchor+log(integral$value)
    mix <- function(v){
      if(!isTRUE(fit$zi)) return(v)
      log_g <- log(values$pi[i]); log_q <- log1p(-values$pi[i])
      if(y[i]>0) return(log_q+v)
      a <- c(log_g,log_q+v); m <- max(a); m+log(sum(exp(a-m)))
    }
    ll_gh <- mix(count_gh); ll_direct <- mix(count_direct)
    data.frame(New_ID=d$New_ID[i],N=N[i],y=y[i],Observation_SD=values$sd[i],pi=values$pi[i],
      log_likelihood_quadrature=ll_gh,log_likelihood_direct=ll_direct,
      absolute_difference=abs(ll_gh-ll_direct),count_absolute_difference=abs(count_gh-count_direct),
      passed=is.finite(ll_direct)&&abs(ll_gh-ll_direct)<1e-4&&abs(count_gh-count_direct)<1e-4)
  }))
}

legacy_simulate_centered_olre <- simulate_centered_olre
simulate_centered_olre <- function(fit,n_sim,seed) {
  if(is.null(fit$joint_spec)) return(legacy_simulate_centered_olre(fit,n_sim,seed))
  d <- fit$data; N <- d$Ruegeria_reads_total; rv <- joint_row_values(fit)
  alpha <- fit$parameters$alpha; delta <- alpha/sqrt(1+alpha^2)
  omega <- rv$sd/sqrt(1-2*delta^2/pi); xi <- -omega*delta*sqrt(2/pi)
  set.seed(seed)
  sims <- replicate(n_sim,{
    u <- xi+omega*(delta*abs(rnorm(nrow(d)))+sqrt(1-delta^2)*rnorm(nrow(d)))
    counts <- rbinom(nrow(d),size=N,prob=plogis(fit$eta_base+u))
    if(isTRUE(fit$zi)) counts[rbinom(nrow(d),1,rv$pi)==1] <- 0L
    counts
  })
  expected <- rowMeans(sims)
  res <- DHARMa::createDHARMa(simulatedResponse=sims,observedResponse=d$MC_reads,
    fittedPredictedResponse=expected,integerResponse=TRUE,seed=seed)
  res$fittedModel <- fit
  list(res_check=res,total_reads=N,expected_prop=expected/N,observation_sd=rv$sd,
    simulation_method="Conditional Repeat/TS; redraw centered skew-normal Observation with row-specific SD, binomial counts and row-specific extra-zero indicator")
}

joint_interior <- function(result,d) {
  if(is.null(result)) return(FALSE)
  p <- result$parameters; spec <- result$spec
  ls <- p$log_sd_observation+if(isTRUE(spec$sd_reads))p$sd_reads*d$log_reads_c else 0
  all(is.finite(ls)) && max(abs(p$alpha))<19.99 &&
    p$log_sd_observation> -7.99 && p$log_sd_observation<2.99 &&
    (!isTRUE(spec$sd_reads)||abs(p$sd_reads)<9.99) &&
    (!isTRUE(spec$zi_reads)||abs(p$zi_reads)<9.99) &&
    (!isTRUE(spec$zi)||p$gamma_zi<qlogis(.99)-.01) &&
    (!isTRUE(spec$zi)||is.null(spec$lambda_zi_reads)||spec$lambda_zi_reads<=0||p$gamma_zi> -17.99)
}
joint_starts <- function(seed,d,spec) {
  a <- max(-8,min(8,seed$alpha))
  aa <- c(a,if(abs(a)>.1)-a else .5,-2,2,-6,6)
  zs <- c(-1,1,-.5,.5,-2,2); ss <- c(-.2,.2,-.5,.5,0,0)
  lapply(seq_along(aa),function(j){
    p <- joint_parameters(seed,d,spec); p$alpha <- aa[j]
    if(isTRUE(spec$zi_reads)) p$zi_reads <- zs[j]
    if(isTRUE(spec$sd_reads)) p$sd_reads <- ss[j]
    p
  })
}
joint_fit_pass <- function(seed,d,X,stub,spec,output_dir,pool=NULL,nodes=31L) {
  dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
  attempts <- list(); precision_rows <- list(); warnings <- character()
  best <- NULL; selected <- list(reproducible=FALSE)
  precision <- data.frame(passed=FALSE); direct <- data.frame(passed=FALSE)
  for(stage in 1:3){
    ps <- joint_starts(seed,d,spec)
    cat("    ",spec$id,": ",length(ps)," 起点；积分节点 ",nodes,
      if(isTRUE(spec$zero_repeat))"；Repeat方差=0" else "","\n",sep="")
    flush.console()
    args <- list(d=d,design=X,nodes=nodes,spec=spec)
    trials_out <- if(is.null(pool)) lapply(seq_along(ps),function(j) {
      cat("      起点 ",j,"/",length(ps)," 开始\n",sep=""); flush.console()
      result <- do.call(joint_worker,c(list(p=ps[[j]]),args))
      cat("      起点 ",j," 完成；",round(result$elapsed_seconds,1),"秒\n",sep=""); flush.console()
      result
    }) else
      parallel::parLapplyLB(pool,ps,joint_worker,d=d,design=X,nodes=nodes,spec=spec,chunk.size=1L)
    trials <- lapply(trials_out,`[[`,"fit")
    for(j in seq_along(trials)){
      f <- trials[[j]]
      warnings <- unique(c(warnings,trials_out[[j]]$warnings,trials_out[[j]]$error))
      attempts[[length(attempts)+1L]] <- data.frame(stage=stage,nodes=nodes,start=j,
        alpha_start=ps[[j]]$alpha,zi_reads_start=if(isTRUE(spec$zi_reads))ps[[j]]$zi_reads else 0,
        sd_reads_start=if(isTRUE(spec$sd_reads))ps[[j]]$sd_reads else 0,
        nll=if(is.null(f))NA_real_ else f$value,convergence=if(is.null(f))NA_integer_ else f$opt$convergence,
        raw_data_nll=if(is.null(f))NA_real_ else f$data_nll,
        objective_value=if(is.null(f))NA_real_ else f$value,
        penalty_value=if(is.null(f))NA_real_ else f$penalty_value,
        lambda_zi_reads=if(is.null(spec$lambda_zi_reads))0 else spec$lambda_zi_reads,
        pdHess=!is.null(f)&&isTRUE(f$pdHess),basic_valid=basic_fit_valid(f),
        max_gradient=if(is.null(f))NA_real_ else f$max_gradient,
        elapsed_seconds=trials_out[[j]]$elapsed_seconds)
    }
    selected <- select_best_likelihood(trials); best <- selected$fit
    if(is.null(best)) break
    check_started <- proc.time()[["elapsed"]]
    cat("    开始高节点精度及独立积分复核。\n"); flush.console()
    precision <- tryCatch(joint_precision(best,d,X,max(61L,2L*nodes-1L)),
      error=function(e)data.frame(passed=FALSE,reason=conditionMessage(e)))
    wrapped <- joint_wrap(best,stub,d,X)
    direct <- tryCatch(joint_integrals(wrapped),error=function(e)data.frame(passed=FALSE,reason=conditionMessage(e)))
    cat("    积分复核完成；",round(proc.time()[["elapsed"]]-check_started,1),"秒\n",sep=""); flush.console()
    precision_rows[[length(precision_rows)+1L]] <- cbind(stage=stage,precision)
    write_result_table(bind_flexible_tables(attempts),output_dir,"01_multistart_fit_attempts.csv")
    write_result_table(bind_flexible_tables(precision_rows),output_dir,"01_quadrature_precision_checks.csv")
    write_result_table(direct,output_dir,"01_independent_integrals.csv")
    saveRDS(list(best=best,selected=selected,precision=precision,direct=direct,
      attempts=bind_flexible_tables(attempts),spec=spec),file.path(output_dir,"01_latest_fit_stage.rds"))
    # A nonzero optimizer gradient is not evidence of inaccurate quadrature.
    # Refining nodes only helps when likelihood or independent integral agreement fails.
    integration_agrees <- nrow(precision)==1L &&
      "absolute_nll_difference"%in%names(precision) &&
      is.finite(precision$absolute_nll_difference) &&
      precision$absolute_nll_difference<.01 && all(direct$passed%in%TRUE)
    if(integration_agrees||nodes>=maximum_fit_nodes)break
    nodes <- min(maximum_fit_nodes,2L*nodes-1L); seed <- best$parameters
  }
  reasons <- c(if(is.null(best))"所有起点失败",
    if(!is.null(best)&&!basic_fit_valid(best))"优化、Hessian或协方差未通过",
    if(!isTRUE(selected$reproducible))"最优似然未由至少两个有效起点复现",
    if(!isTRUE(precision$passed))"积分精度复核未通过",
    if(!all(direct$passed%in%TRUE))"独立条件积分未通过",
    if(!is.null(best)&&!joint_interior(best,d))"参数触及边界")
  valid <- !is.null(best)&&basic_fit_valid(best)&&isTRUE(selected$reproducible)&&
    isTRUE(precision$passed)&&all(direct$passed%in%TRUE)&&joint_interior(best,d)
  writeLines(enc2utf8(Filter(nzchar,warnings)),file.path(output_dir,"01_fitting_messages.txt"))
  list(primary=if(is.null(best))NULL else joint_wrap(best,stub,d,X),primary_valid=valid,
    normal=NULL,normal_valid=FALSE,reproducible=isTRUE(selected$reproducible),
    quadrature_ok=isTRUE(precision$passed),direct_integrals_ok=all(direct$passed%in%TRUE),
    reasons=reasons,raw_fit=best,attempts=bind_flexible_tables(attempts),spec=spec)
}

joint_variance_boundaries <- function(candidate,d,X,stub,output_dir,pool) {
  f <- candidate$primary
  if(is.null(f))return(candidate)
  small <- c("repeat"=exp(f$parameters$log_sd_repeat)<.003,ts=exp(f$parameters$log_sd_ts)<.003)
  # A precisely zero variance is a genuine boundary of the same random-effect family.
  # Do not replace a better invalid likelihood with a materially worse valid solution.
  if(!any(small))return(candidate)
  spec <- candidate$spec
  if(small["repeat"])spec$zero_repeat <- TRUE
  if(small["ts"])spec$zero_ts <- TRUE
  boundary <- joint_fit_pass(f$parameters,d,X,stub,spec,file.path(output_dir,"variance_boundary"),pool,f$quadrature_n)
  b <- boundary$primary
  comparison <- data.frame(full_nll=f$nll,full_valid=candidate$primary_valid,
    boundary_nll=if(is.null(b))NA_real_ else b$nll,boundary_valid=boundary$primary_valid,
    zero_repeat=isTRUE(spec$zero_repeat),zero_ts=isTRUE(spec$zero_ts),chosen=FALSE)
  comparison$full_objective <- joint_total_objective(f)
  comparison$boundary_objective <- if(is.null(b))NA_real_ else joint_total_objective(b)
  if(!is.null(b)&&isTRUE(boundary$primary_valid)&&joint_total_objective(b)<=joint_total_objective(f)+.001){
    comparison$chosen <- TRUE
    boundary$reasons <- character(); candidate <- boundary
    candidate$boundary_note <- paste(c(if(isTRUE(spec$zero_repeat))"Repeat方差经似然核对取精确0",
      if(isTRUE(spec$zero_ts))"TS方差经似然核对取精确0"),collapse="; ")
  }
  write_result_table(comparison,output_dir,"01_variance_boundary_comparison.csv")
  candidate
}

joint_verify_cached <- function(mc,d,X,stub,output_dir) {
  path <- file.path(joint_cache_dir,mc,"fitted_comparison.rds")
  if(!file.exists(path))return(NULL)
  saved <- readRDS(path); candidate <- saved$candidates$ZI; f <- candidate$primary
  if(is.null(f))return(NULL)
  cols <- c("New_ID","Sample_ID","MC_ID","MC_reads","Other_reads","Ruegeria_reads_total",
    "DO_c","Host_type","Repeat_group_ID","TS","log_reads_c","mean_reads_s1","mean_reads_s2","mean_reads_s3")
  same <- nrow(f$data)==nrow(d)&&all(vapply(cols,function(nm){
    if(is.numeric(d[[nm]]))isTRUE(all.equal(as.numeric(d[[nm]]),as.numeric(f$data[[nm]]),tolerance=1e-12))
    else identical(as.character(d[[nm]]),as.character(f$data[[nm]]))},logical(1))) &&
    isTRUE(all.equal(unname(f$X),unname(X),tolerance=1e-12))&&identical(colnames(f$X),colnames(X))
  if(!same)stop("A缓存与当前数据/设计不一致：",mc)
  spec <- joint_specs$A; spec$zi <- isTRUE(f$zi)
  p <- joint_parameters(f$parameters,d,spec)
  obj <- RTMB::MakeADFun(joint_objective(d,X,f$quadrature_n,spec),p,random=joint_random(spec),silent=TRUE)
  val <- obj$fn(obj$par); gr <- obj$gr(obj$par)
  stopifnot(abs(val-f$nll)<.01,identical(names(obj$par),names(f$fit$par)))
  raw <- list(opt=list(par=obj$par,convergence=f$fit$convergence,message="Verified prior A"),
    parameters=obj$env$parList(obj$par),value=val,max_gradient=max(abs(gr)),
    sdr=list(pdHess=f$pdHess,cov.fixed=f$cov_fixed),pdHess=f$pdHess,
    quadrature_n=f$quadrature_n,skew=TRUE,zi=isTRUE(f$zi),k=f$k,spec=spec)
  wrapped <- joint_wrap(raw,stub,d,X)
  precision <- joint_precision(raw,d,X,max(61L,2L*f$quadrature_n-1L))
  direct <- joint_integrals(wrapped)
  dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
  write_result_table(precision,output_dir,"01_cached_precision_recheck.csv")
  write_result_table(direct,output_dir,"01_cached_direct_integrals.csv")
  writeLines(path,file.path(output_dir,"01_cached_source.txt"))
  valid <- isTRUE(candidate$primary_valid)&&basic_fit_valid(raw)&&isTRUE(precision$passed)&&all(direct$passed)
  list(primary=wrapped,primary_valid=valid,normal=NULL,normal_valid=FALSE,
    reproducible=isTRUE(candidate$reproducible),quadrature_ok=isTRUE(precision$passed),
    direct_integrals_ok=all(direct$passed),reasons=if(valid)character() else candidate$reasons,
    raw_fit=raw,spec=spec,source=path)
}

joint_pi0 <- function(candidate,nozi,d,X,output_dir) {
  f <- candidate$primary; b <- nozi$primary
  if(is.null(f)||is.null(b)||!isTRUE(nozi$primary_valid)||!isTRUE(f$zi))return(candidate)
  rv <- joint_row_values(f)
  if(max(rv$pi)>=1e-6)return(candidate)
  # Recheck the finite direction into the ZI model at constant small pi.
  p <- b$parameters; sp <- b$joint_spec; sp$zi <- TRUE; sp$zi_reads <- FALSE
  p$gamma_zi <- qlogis(1e-5)
  obj <- RTMB::MakeADFun(joint_objective(d,X,max(f$quadrature_n,b$quadrature_n),sp),p,
    random=joint_random(sp),map=joint_map(p,sp),silent=TRUE)
  val <- obj$fn(obj$par); score <- obj$gr(obj$par)[names(obj$par)=="gamma_zi"]/(1e-5*(1-1e-5))
  fq <- joint_total_objective(f); bq <- joint_total_objective(b)
  ok <- is.finite(score)&&score>=0&&abs(fq-bq)<.001&&fq>=bq-.001
  write_result_table(data.frame(max_pi=max(rv$pi),pi0_nll=b$nll,ZI_nll=f$nll,small_pi_nll=val,
    small_pi_score=unname(score),ZI_objective=fq,pi0_objective=bq,passed=ok),output_dir,"01_pi0_boundary_check.csv")
  if(ok){
    answer <- nozi; answer$primary$model_id <- f$model_id
    answer$boundary_note <- "额外零值概率经边界核查为0；按同一计数模型计参数和AIC"
    answer$spec$id <- candidate$spec$id; answer$primary$joint_spec$id <- candidate$spec$id
    return(answer)
  }
  candidate
}

joint_summary_row <- function(mc,model,n_source,d=NULL) {
  r <- new_comparison_row(mc,model,n_source)
  r$zi_reads_slope <- r$sd_reads_slope <- NA_real_
  r$pi_min <- r$pi_at_center <- r$pi_max <- NA_real_
  r$SD_min <- r$SD_at_center <- r$SD_max <- NA_real_
  r$variance_boundary <- ""; r$fit_source <- "new_fit"
  r$delta_AIC_vs_A <- r$delta_AIC_vs_best <- NA_real_
  r$diagnostic_flag_count <- NA_integer_
  if(!is.null(d)){
    r$n_used <- nrow(d); r$n_excluded <- n_source-nrow(d)
    r$n_repeat_groups <- nlevels(d$Repeat_group_ID); r$n_TS <- nlevels(d$TS)
    r$n_zero_MC <- sum(d$MC_reads==0); r$n_all_reads_MC <- sum(d$Other_reads==0)
    r$observed_mean_proportion <- mean(d$MC_reads/d$Ruegeria_reads_total)
  }
  r
}
joint_finish_row <- function(r,candidate) {
  r <- fill_fit_summary(r,candidate)
  f <- candidate$primary
  if(is.null(f))return(r)
  rv <- joint_row_values(f)
  if(isTRUE(f$joint_spec$sd_mean)||isTRUE(f$joint_spec$alpha_mean)) {
    center_X <- matrix(colMeans(f$X),nrow=1L,dimnames=list(NULL,colnames(f$X)))
    center <- joint_row_values(f,0,DO_c=0,TS=NULL,X=center_X)
  } else center <- joint_row_values(f,0)
  r$extra_zero_pi <- mean(rv$pi)
  r$pi_min <- min(rv$pi); r$pi_max <- max(rv$pi); r$pi_at_center <- center$pi
  r$SD_min <- min(rv$sd); r$SD_max <- max(rv$sd); r$SD_at_center <- center$sd
  r$zi_reads_slope <- if(isTRUE(f$joint_spec$zi_reads))f$parameters$zi_reads else 0
  r$sd_reads_slope <- if(isTRUE(f$joint_spec$sd_reads))f$parameters$sd_reads else 0
  r$variance_boundary <- paste(c(if(isTRUE(f$joint_spec$zero_repeat))"Repeat=0",
    if(isTRUE(f$joint_spec$zero_ts))"TS=0"),collapse="; ")
  if(!is.null(candidate$boundary_note))r$notes <- paste(Filter(nzchar,c(r$notes,candidate$boundary_note)),collapse="; ")
  r$zi_mode <- if(isTRUE(f$zi))if(isTRUE(f$joint_spec$zi_reads))"reads_dependent_pi" else "constant_pi" else "exact_pi0"
  r
}
joint_html <- function(summary,metrics,run_dir) {
  make_table <- function(v) {
    heads <- paste0("<tr>",paste0("<th>",html_escape(names(v)),"</th>",collapse=""),"</tr>")
    rows <- vapply(seq_len(nrow(v)),function(i)paste0("<tr>",paste0("<td>",
      vapply(v,function(a)html_escape(if(is.na(a[i]))"—" else as.character(a[i])),character(1)),
      "</td>",collapse=""),"</tr>"),character(1))
    paste0('<div class="scroll"><table><thead>',heads,'</thead><tbody>',paste(rows,collapse=""),'</tbody></table></div>')
  }
  top <- data.frame(MC=summary$MC_ID,最低AIC候选=summary$best_AIC_model,
    完成诊断且无主要提示=summary$adequate_models,建议=summary$next_step,
    最低AIC候选状态=summary$best_status,残差提示=summary$best_flags,check.names=FALSE)
  v <- data.frame(MC=metrics$MC_ID,候选=metrics$model,数值有效=metrics$numerical_valid,
    AIC=fmt_number(ifelse(metrics$numerical_valid,metrics$AIC,NA),2),
    ΔAIC较当前=fmt_number(metrics$delta_AIC_vs_A,2),ΔAIC较最低=fmt_number(metrics$delta_AIC_vs_best,2),
    实际零值=metrics$n_zero_MC,预测零值=fmt_number(metrics$expected_zero_mean,1),
    零值p=fmt_p(metrics$Zero_count_p),reads组零值BH=fmt_p(metrics$Zero_Reads_min_BH_p),
    零值reads斜率=fmt_number(metrics$zi_reads_slope),波动reads斜率=fmt_number(metrics$sd_reads_slope),
    pi最低=fmt_number(metrics$pi_min),pi中心=fmt_number(metrics$pi_at_center),pi最高=fmt_number(metrics$pi_max),
    SD最低=fmt_number(metrics$SD_min),SD中心=fmt_number(metrics$SD_at_center),SD最高=fmt_number(metrics$SD_max),
    DOβ=fmt_number(metrics$DO_beta),DO_BH=fmt_p(metrics$DO_p_BH),
    宿主BH=fmt_p(metrics$Host_p_BH),reads均值BH=fmt_p(metrics$Reads_p_BH),
    Uniformity=fmt_p(metrics$Uniformity_p),Dispersion=fmt_p(metrics$Dispersion_p),
    预测分位数=fmt_p(metrics$Predicted_quantiles_p),reads中心趋势BH=fmt_p(metrics$Reads_location_trend_BH_p),
    reads偏离趋势BH=fmt_p(metrics$Reads_deviation_trend_BH_p),
    方差边界=metrics$variance_boundary,状态=metrics$status,提示=metrics$diagnostic_flags,说明=metrics$notes,check.names=FALSE)
  html <- c('<!doctype html><html lang="zh-CN"><meta charset="utf-8"><title>MC8 reads零值和波动比较</title>',
    '<style>body{font:15px/1.6 Arial,"Microsoft YaHei",sans-serif;margin:24px;color:#193145}table{border-collapse:separate;border-spacing:0}th,td{padding:8px 12px;border:1px solid #dbe3e9;white-space:nowrap;text-align:right}th{background:#193c55;color:white;position:sticky;top:0;z-index:2}td:first-child{background:#edf3f7;font-weight:bold;position:sticky;left:0}tr:nth-child(even){background:#f6f8fa}.scroll{overflow:auto;max-height:70vh;margin-bottom:25px}td:nth-last-child(-n+3){white-space:normal;min-width:210px;text-align:left}</style>',
    '<h1>全部8个MC：reads对零值概率和样本波动的影响</h1>',
    '<p>A：共同pi＋共同SD；B：pi随reads改变；C：SD随reads改变；D：两项均随reads改变。平均丰度公式、df=3、样本、分母、偏正态家族一致。每个MC各自估计参数。</p>',
    '<p>最低AIC与残差是否合适分别列出。无效候选的推断留空。无零值时采用精确pi=0；经似然核实的方差边界明确标注。DO/Host/reads的BH分别跨8个计划MC，每个候选分别校正。</p>',
    '<h2>每个MC的重要结果</h2>',make_table(top),'<h2>四个候选的完整比较</h2>',make_table(v),
    '<p><a href="00_all_MC_joint_summary.csv">每个MC汇总</a> · <a href="01_all_model_metrics.csv">全部32行完整数值</a> · <a href="00_key_results.txt">中文摘要</a> · <a href="METHODS.txt">模型与检查说明</a></p>',
    '<p>pi是计数成分之外的额外零值概率；SD是中心化偏正态样本偏移的标准差。两种reads斜率均使用中心化log10 reads。总体预期比例=(1-pi_i)×E_count(p_i)。本轮不进行逐TS剔除。</p></html>')
  writeLines(enc2utf8(html),file.path(run_dir,"00_all_MC_joint_summary.html"))
}
joint_outputs <- function(table,MC_ids,run_dir,completed) {
  for(model in names(joint_specs)){
    idx <- which(table$model==model)
    for(term in c("DO","Host","Reads")){
      pp <- table[[paste0(term,"_p_raw")]][idx]; pp[!table$numerical_valid[idx]] <- NA_real_
      table[[paste0(term,"_p_BH")]][idx] <- adjust_valid_p(pp,n=length(MC_ids))
    }
  }
  summaries <- lapply(MC_ids,function(mc){
    idx <- which(table$MC_ID==mc); v <- table[idx,,drop=FALSE]
    valid <- which(v$numerical_valid & is.finite(v$AIC))
    a <- which(v$model=="A")
    if(length(valid)){
      best <- valid[which.min(v$AIC[valid])]
      table$delta_AIC_vs_best[idx[valid]] <<- v$AIC[valid]-v$AIC[best]
      if(length(a)&&v$numerical_valid[a])table$delta_AIC_vs_A[idx[valid]] <<- v$AIC[valid]-v$AIC[a]
      good <- which(v$numerical_valid & v$status=="数值有效；未检出主要残差偏差")
      data.frame(MC_ID=mc,best_AIC_model=v$model[best],best_AIC=v$AIC[best],
        valid_models=paste(v$model[valid],collapse=", "),adequate_models=paste(v$model[good],collapse=", "),
        best_status=v$status[best],best_flags=v$diagnostic_flags[best],
        next_step=if(v$status[best]=="数值有效；未检出主要残差偏差")"最低AIC候选可进入批次敏感性复核" else
          if(length(good))"最低AIC仍有偏差；与诊断较好候选权衡" else "仍需修复残差结构",stringsAsFactors=FALSE)
    } else data.frame(MC_ID=mc,best_AIC_model="",best_AIC=NA_real_,valid_models="",adequate_models="",
      best_status="没有数值有效候选",best_flags="",next_step="先修复数值拟合",stringsAsFactors=FALSE)
  })
  summary <- do.call(rbind,summaries)
  write_excel_csv(table,file.path(run_dir,"01_all_model_metrics.csv"))
  write_excel_csv(summary,file.path(run_dir,"00_all_MC_joint_summary.csv"))
  joint_html(summary,table,run_dir)
  txt <- capture.output({
    cat("全部MC：reads零值概率与样本波动的四候选比较\n")
    cat("A共同pi/共同SD；B reads->pi；C reads->SD；D两项均改变。\n")
    cat("完成MC：",completed,"/",length(MC_ids),"；数值有效候选：",sum(table$numerical_valid),"/",nrow(table),"。\n",sep="")
    cat("\n1. 每个MC的重要结论\n"); print(summary,row.names=FALSE)
    cat("\n2. 模型比较与零值预测\n")
    print(table[,c("MC_ID","model","numerical_valid","AIC","delta_AIC_vs_A","delta_AIC_vs_best","n_zero_MC",
      "expected_zero_mean","Zero_count_p","Zero_Reads_min_BH_p","variance_boundary")],row.names=FALSE)
    cat("\n3. reads零值与波动斜率、实际范围\n")
    print(table[,c("MC_ID","model","zi_reads_slope","sd_reads_slope","pi_min","pi_at_center","pi_max","SD_min","SD_max")],row.names=FALSE)
    cat("\n4. DO结果与残差状态\n")
    print(table[,c("MC_ID","model","DO_beta","DO_CI_low","DO_CI_high","DO_p_BH","Host_p_BH","Reads_p_BH","status")],row.names=FALSE)
    cat("\n5. 需要复核的项目\n")
    print(table[,c("MC_ID","model","diagnostic_flags","notes")],row.names=FALSE)
    cat("\n数值有效后才计算推断；AIC最低与诊断状态分别呈现。模拟2000次不重新拟合2000次。\n")
    cat("各MC保留零值；共用原始样本与分母。逐TS敏感性分析未在本轮执行。\n")
    cat("优先查看：00_all_MC_joint_summary.html；01_all_model_metrics.csv；00_key_results.txt。\n")
    cat("结果目录：",normalizePath(run_dir,winslash="/",mustWork=FALSE),"\n",sep="")
  })
  writeLines(enc2utf8(txt),file.path(run_dir,"00_key_results.txt"))
  invisible(list(table=table,summary=summary,text=txt))
}

run_all_MC_joint_reads <- function() {
  previous <- options(contrasts=c("contr.treatment","contr.poly"))
  on.exit(options(previous),add=TRUE)
  packages <- c("readxl","glmmTMB","RTMB","DHARMa","qgam","ggplot2")
  missing <- packages[!vapply(packages,requireNamespace,logical(1),quietly=TRUE)]
  if(length(missing))stop("缺少包：",paste(missing,collapse=", "))
  raw <- readxl::read_excel(input_file,sheet=input_sheet)
  MC_ids <- natural_MC_order(unique(as.character(raw$MC_ID[!is.na(raw$MC_ID)])))
  if(!is.null(MC_selection))MC_ids <- natural_MC_order(unique(MC_selection))
  shared <- build_shared_design(raw,MC_ids,host_reference)
  resuming <- nzchar(joint_resume_dir)
  run_dir <- if(resuming)joint_resume_dir else file.path(output_root,paste0("GLMM_joint_reads_MC8_",format(Sys.time(),"%Y%m%d_%H%M%S")))
  if(!resuming && dir.exists(run_dir))run_dir <- paste0(run_dir,"_",Sys.getpid())
  dir.create(run_dir,recursive=TRUE,showWarnings=FALSE)
  settings <- list(input_file=input_file,input_md5=unname(tools::md5sum(input_file)),MC_ids=MC_ids,
    n_sim=n_sim,seed=simulation_seed,mean_df=3L,specs=joint_specs,baseline_cache=joint_cache_dir,
    maximum_fit_nodes=maximum_fit_nodes,batch_sensitivity=FALSE,
    script_md5=unname(tools::md5sum(getOption("GLMM_JOINT_worker_source",file.path(project_dir,"scripts","GLMM_reads_diagnostics_full.R")))))
  if(resuming){
    old <- readRDS(file.path(run_dir,"01_run_settings.rds"))
    for(nm in c("input_md5","MC_ids","n_sim","seed","mean_df","specs","script_md5"))
      if(!identical(settings[[nm]],old[[nm]]))stop("恢复运行设置不一致：",nm)
  }
  saveRDS(settings,file.path(run_dir,"01_run_settings.rds"))
  saveRDS(shared,file.path(run_dir,"01_shared_model_design.rds"))
  write_excel_csv(shared$samples,file.path(run_dir,"01_shared_sample_design.csv"))
  writeLines(run_dir,file.path(tempdir(),"GLMM_joint_current_run.txt"))
  if(nzchar(Sys.getenv("GLMM_JOINT_POINTER")))writeLines(run_dir,Sys.getenv("GLMM_JOINT_POINTER"))
  writeLines("running",file.path(run_dir,"00_run_status.txt"))
  writeLines(capture.output(utils::sessionInfo()),file.path(run_dir,"sessionInfo.txt"))
  pool <- parallel::makePSOCKcluster(max(1L,min(as.integer(parallel_start_workers),parallel::detectCores())))
  on.exit(parallel::stopCluster(pool),add=TRUE)
  src <- getOption("GLMM_JOINT_worker_source",file.path(project_dir,"scripts","GLMM_reads_diagnostics_full.R"))
  parallel::clusterCall(pool,function(path){
    options(GLMM_ZI_definitions_only=TRUE,GLMM_JOINT_definitions_only=TRUE)
    Sys.setenv(OMP_NUM_THREADS=1,OPENBLAS_NUM_THREADS=1,MKL_NUM_THREADS=1)
    source(path,encoding="UTF-8"); requireNamespace("RTMB",quietly=TRUE); TRUE
  },src)
  rows <- list(); all_zeros <- list(); all_overall <- list(); all_coef <- list(); model_files <- list()
  for(i in seq_along(MC_ids)){
    mc <- MC_ids[i]; mc_dir <- file.path(run_dir,mc); dir.create(mc_dir,recursive=TRUE,showWarnings=FALSE)
    cat("\n========== ",i,"/",length(MC_ids)," ",mc," ==========\n",sep="");flush.console()
    prepared <- prepare_mc_shared(raw,mc,shared); d <- prepared$data_test; X <- prepared$X; stub <- make_model_stub(d)
    write_result_table(d,mc_dir,"01_model_input_rows.csv")
    write_result_table(prepared$sample_audit,mc_dir,"01_sample_audit.csv")
    candidates <- list(); mc_rows <- list()
    for(model in names(joint_specs)){
      spec <- joint_specs[[model]]; output_dir <- file.path(mc_dir,model)
      dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
      checkpoint <- file.path(output_dir,"fitted_candidate.rds")
      cat("  候选 ",model,"\n",sep="");flush.console()
      r <- joint_summary_row(mc,model,sum(as.character(raw$MC_ID)==mc,na.rm=TRUE),d)
      captured <- tryCatch({
        separation <- if(isTRUE(spec$zi_reads))joint_zero_reads_separation(d) else ""
        if(nzchar(separation)){
          candidate <- list(primary=NULL,primary_valid=FALSE,normal=NULL,normal_valid=FALSE,
            reproducible=FALSE,quadrature_ok=FALSE,direct_integrals_ok=FALSE,
            reasons=paste("reads零值完全分离",separation,sep=": "),spec=spec)
          writeLines(enc2utf8(candidate$reasons),file.path(output_dir,"01_nonidentifiable_zero_slope.txt"))
          cat("    ",candidate$reasons,"\n",sep="")
        } else if(file.exists(checkpoint)){
          candidate <- readRDS(checkpoint); cat("    恢复已保存候选。\n")
        } else if(model=="A"){
          candidate <- joint_verify_cached(mc,d,X,stub,output_dir)
          if(is.null(candidate))stop("当前A缓存缺失；请先核对缓存目录")
        } else if(!any(d$MC_reads==0)&&model%in%c("B","D")){
          # With no zeros, varying pi cannot increase the likelihood; exact pi=0.
          candidate <- candidates[[if(model=="B")"A" else "C"]]
          candidate$spec$id <- model; candidate$primary$joint_spec$id <- model
          candidate$primary$model_id <- paste0("JOINT_",model,"_PI0")
          candidate$boundary_note <- "无零值：精确pi=0，与相应计数模型相同，不另计零值参数"
        } else {
          seed_candidates <- Filter(function(c)!is.null(c$primary),candidates)
          seed_values <- vapply(seed_candidates,function(c)c$primary$nll,numeric(1))
          seed <- seed_candidates[[which.min(seed_values)]]$primary
          if(!any(d$MC_reads==0)){spec$zi <- FALSE;spec$zi_reads <- FALSE}
          # Start full-variance fits with finite SDs even when a prior candidate reached zero.
          p <- seed$parameters
          for(nm in c("log_sd_repeat","log_sd_ts"))if(!is.finite(p[[nm]]))p[[nm]] <- -6
          candidate <- joint_fit_pass(p,d,X,stub,spec,output_dir,pool,max(31L,seed$quadrature_n))
          candidate <- joint_variance_boundaries(candidate,d,X,stub,output_dir,pool)
          if(!is.null(candidate$primary)&&isTRUE(candidate$primary$zi)&&max(joint_row_values(candidate$primary)$pi)<1e-6){
            nozi <- if(!isTRUE(spec$sd_reads)&&!isTRUE(candidates$A$primary$zi))candidates$A else NULL
            if(is.null(nozi)){
              sp <- candidate$spec;sp$zi <- FALSE;sp$zi_reads <- FALSE
              nozi <- joint_fit_pass(candidate$primary$parameters,d,X,stub,sp,file.path(output_dir,"pi0_count_boundary"),pool)
              nozi <- joint_variance_boundaries(nozi,d,X,stub,file.path(output_dir,"pi0_count_boundary"),pool)
            }
            candidate <- joint_pi0(candidate,nozi,d,X,output_dir)
          }
        }
        # Check the same small-variance rule for A and expanded candidates.
        if(model=="A" && is.null(candidate$boundary_note))
          candidate <- joint_variance_boundaries(candidate,d,X,stub,output_dir,pool)
        saveRDS(candidate,checkpoint)
        candidates[[model]] <- candidate
        r <- joint_finish_row(r,candidate)
        if(model=="A")r$fit_source <- if(is.null(candidate$boundary_note))"verified_previous_A" else "verified_A_variance_boundary"
        if(r$numerical_valid){
          diagnostic_dir <- file.path(output_dir,"diagnostics");dir.create(diagnostic_dir,showWarnings=FALSE)
          diagnostic_file <- file.path(diagnostic_dir,"model_analysis_objects.rds")
          if(file.exists(diagnostic_file))diag <- readRDS(diagnostic_file) else {
            cat("    数值检查通过；2000次模拟及残差诊断。\n");flush.console()
            diag <- diagnose_skew_comparison_model(paste(mc,model),candidate$primary,d,n_sim,simulation_seed+i*100L,diagnostic_dir)
            diag <- zero_count_diagnostics(diag,d,diagnostic_dir)
            saveRDS(diag,diagnostic_file)
          }
          r <- fill_zero_diagnostic_row(r,diag)
          r$diagnostic_flag_count <- if(nzchar(r$diagnostic_flags))length(strsplit(r$diagnostic_flags,"; ",fixed=TRUE)[[1]]) else 0L
          all_zeros[[paste(mc,model)]] <- cbind(MC_ID=mc,model=model,diag$zero_group_table)
          all_overall[[paste(mc,model)]] <- cbind(MC_ID=mc,model=model,diag$overall_table)
          all_coef[[paste(mc,model)]] <- cbind(model=model,fixed_coefficient_table(candidate$primary,mc))
          rv <- joint_row_values(candidate$primary)
          write_result_table(data.frame(New_ID=d$New_ID,log_reads_c=d$log_reads_c,pi=rv$pi,
            Observation_SD=rv$sd,expected_proportion=diag$fitted_rows$expected_MC_proportion),output_dir,"03_row_specific_pi_and_SD.csv")
        }
        if(!is.null(candidate$primary))writeLines(enc2utf8(capture.output(base::summary(candidate$primary))),file.path(output_dir,"02_model_summary.txt"))
        TRUE
      },error=function(e){
        r$status <<- if(r$numerical_valid)"数值有效；诊断失败" else "处理失败"
        r$notes <<- paste(Filter(nzchar,c(r$notes,substr(conditionMessage(e),1,1000))),collapse="; ")
        writeLines(enc2utf8(conditionMessage(e)),file.path(output_dir,"00_error.txt"))
        cat("    已记录错误：",substr(conditionMessage(e),1,1000),"\n",sep="");FALSE
      })
      mc_rows[[model]] <- r
      cat("    ",model,": ",r$status,"；AIC=",round(r$AIC,3),"\n",sep="");flush.console()
      saveRDS(list(MC_ID=mc,candidates=candidates,metrics=do.call(rbind,mc_rows)),file.path(mc_dir,"fitted_comparison.rds"))
    }
    curves <- lapply(names(candidates),function(m){
      cc <- candidates[[m]]
      if(!isTRUE(cc$primary_valid))return(NULL)
      x <- seq(min(d$log_reads_c),max(d$log_reads_c),length.out=150)
      rv <- joint_row_values(cc$primary,x)
      data.frame(MC_ID=mc,model=m,log10_reads=x+shared$log_reads_center,pi=rv$pi,SD=rv$sd)
    })
    curves <- bind_flexible_tables(curves)
    if(nrow(curves)){
      write_result_table(curves,mc_dir,"05_zero_probability_and_SD_curves.csv")
      plot_data <- rbind(data.frame(curves[,c("model","log10_reads")],component="Extra-zero probability",value=curves$pi),
        data.frame(curves[,c("model","log10_reads")],component="Observation SD (logit)",value=curves$SD))
      fig <- ggplot2::ggplot(plot_data,ggplot2::aes(log10_reads,value,color=model))+
        ggplot2::geom_line(linewidth=.85)+ggplot2::facet_wrap(~component,scales="free_y",ncol=1)+
        ggplot2::theme_bw(base_size=12)+ggplot2::labs(title=mc,x="Total Ruegeria reads (log10)",y=NULL,color="Candidate")
      ggplot2::ggsave(file.path(mc_dir,"05_zero_probability_and_SD_curves.png"),fig,width=7,height=7,dpi=240,bg="white")
      ggplot2::ggsave(file.path(mc_dir,"05_zero_probability_and_SD_curves.pdf"),fig,width=7,height=7)
    }
    rows[[mc]] <- do.call(rbind,mc_rows);model_files[[mc]] <- file.path(mc_dir,"fitted_comparison.rds")
    table <- do.call(rbind,rows)
    for(pending in setdiff(MC_ids,names(rows)))for(model in names(joint_specs))
      table <- rbind(table,joint_summary_row(pending,model,sum(as.character(raw$MC_ID)==pending,na.rm=TRUE)))
    final <- joint_outputs(table,MC_ids,run_dir,length(rows))
    write_excel_csv(bind_flexible_tables(all_zeros),file.path(run_dir,"02_all_zero_group_diagnostics.csv"))
    write_excel_csv(bind_flexible_tables(all_overall),file.path(run_dir,"03_all_overall_diagnostics.csv"))
    write_excel_csv(bind_flexible_tables(all_coef),file.path(run_dir,"04_all_fixed_coefficients.csv"))
    saveRDS(list(metrics=final$table,summary=final$summary,model_files=model_files,completed=length(rows),settings=settings),
      file.path(run_dir,"GLMM_joint_reads_results.rds"))
    gc(FALSE)
  }
  writeLines(enc2utf8(c(
    "Fixed candidate set for all MCs: A constant pi/SD; B reads-dependent pi; C reads-dependent SD; D both.",
    "Count mean: DO_c + Host + ns(centered log10 Ruegeria reads,df=3) + Gaussian Repeat/TS + centered skew-normal Observation.",
    "logit(pi_i)=gamma0+gamma_reads*x_i; log(SD_i)=log_sd_observation+sd_reads*x_i; alpha common within MC.",
    "P(Y=0)=pi_i+(1-pi_i)*P_count(Y=0). The binomial denominator already accounts for ordinary count sampling zeros.",
    "Observation effects integrated by adaptive Gaussian quadrature; Repeat/TS integrated by Laplace. Six starts for new candidates.",
    "Require at least two valid starts reproducing best nll within .001, covariance positive, small gradient, higher-node and independent integral agreement.",
    "Near-zero Repeat/TS variance: exact zero-variance boundary is checked and used only if numerical checks pass and likelihood loss <=.001.",
    "No observed zeros -> exact pi=0; equivalent candidates use the same effective parameter count. Near-zero fitted pi also gets likelihood-direction verification.",
    "Strict separation of observed zeros versus positives by reads: the varying-pi logistic component has no finite ML solution; B/D are retained as invalid rows with empty inference/AIC.",
    "Candidate A reuses previous results only after exact record/response/denominator/design matching and fresh likelihood/integration checks. Diagnostics recomputed.",
    "All valid candidates simulate 2000 data sets conditional on fitted Repeat/TS, redrawing Observation, binomial and extra-zero indicators.",
    "Zero count, grouped zero counts, uniformity, dispersion, outliers, proportions, tails and residual trends checked. Diagnostically flagged candidates remain visible.",
    "DO, Host and mean-reads BH across all 8 planned MCs separately for each of A/B/C/D. Invalid inference blank; AIC ranked only for valid candidates within MC.",
    "AIC minimum does not automatically establish residual adequacy. No ordinary chi-square ZI boundary LRT, no leave-one-TS analysis this round.",
    "Inputs read-only; original total Ruegeria denominator and every zero retained. Missing DO exclusion audited identically.",
    "https://kaskr.r-universe.dev/RTMB/doc/manual.html",
    "https://glmmtmb.github.io/glmmTMB/articles/troubleshooting.html")),file.path(run_dir,"METHODS.txt"))
  stopifnot(identical(unname(tools::md5sum(input_file)),settings$input_md5))
  writeLines(c("complete",run_dir),file.path(run_dir,"00_run_status.txt"))
  cat("\n",paste(final$text,collapse="\n"),"\n",sep="")
  invisible(list(output_dir=run_dir,settings=settings,table=final$table,summary=final$summary))
}

# The active entry point is defined below, after the penalty-specific reporting.

# ===== Active experiment: D with three fixed zero-slope ridge penalties =====
# The same three strengths are fitted to every MC. The middle strength is a
# reference for sensitivity summaries, not a data-selected optimum.
zero_penalty_specs <- function() {
  s <- zero_penalty_prior_sds
  if (!is.numeric(s) || !length(s) || any(!is.finite(s)) || any(s <= 0) || anyDuplicated(s))
    stop("zero_penalty_prior_sds须为互不重复的有限正数。")
  if (length(zero_penalty_reference_sd) != 1L ||
      !zero_penalty_reference_sd %in% s) stop("中间候选标准差必须包含在候选列表中。")
  label <- names(s)
  if (is.null(label)) label <- paste0("strength",seq_along(s))
  specs <- lapply(seq_along(s),function(i) {
    lambda <- 1/(2*s[i]^2)
    id <- paste0("D_s",gsub(".","p",format(s[i],trim=TRUE,scientific=FALSE),fixed=TRUE),
      "_lambda",gsub(".","p",format(lambda,trim=TRUE,scientific=FALSE),fixed=TRUE))
    list(id=id,base_model="D",zi=TRUE,zi_reads=TRUE,sd_reads=TRUE,
      lambda_zi_reads=unname(lambda),prior_sd_zi_reads=unname(s[i]),strength_label=label[i],
      is_reference=isTRUE(unname(s[i])==zero_penalty_reference_sd))
  })
  names(specs) <- vapply(specs,`[[`,character(1),"id")
  specs
}

zero_penalty_summary_row <- function(mc,spec,n_source,d=NULL) {
  r <- joint_summary_row(mc,spec$id,n_source,d)
  r$base_model <- "D"; r$lambda_zi_reads <- spec$lambda_zi_reads
  r$prior_sd_zi_reads <- spec$prior_sd_zi_reads; r$strength_label <- spec$strength_label
  r$is_reference_strength <- spec$is_reference
  r$penalty_active <- FALSE; r$penalty_value <- r$objective_value <- r$raw_data_nll <- NA_real_
  r$zero_intercept <- NA_real_; r$zero_intercept_mode <- "not_fitted"
  r$zero_slope_separation <- ""; r$requested_penalty_target <- "zi_reads_slope_only"
  r$inference_method <- "point_estimates_only; formal intervals/p values not computed"
  r$CV_status <- "not_run"; r$lambda_selection <- "sensitivity_comparison; no automatic optimum"
  r
}

zero_penalty_finish_row <- function(r,candidate) {
  r <- joint_finish_row(r,candidate)
  f <- candidate$primary
  if (is.null(f)) return(r)
  r$raw_data_nll <- f$nll; r$logLik <- -f$nll
  r$penalty_value <- f$penalty_value; r$objective_value <- joint_total_objective(f)
  r$penalty_active <- isTRUE(f$zi) && isTRUE(f$joint_spec$zi_reads)
  r$zero_intercept <- if(isTRUE(f$zi))f$parameters$gamma_zi else -Inf
  r$zero_intercept_mode <- if(isTRUE(f$zi))"estimated_unpenalized_intercept" else "exact_pi0"
  # The inverse penalized Hessian is retained in the RDS for numerical checks;
  # it is not silently reported as validated frequentist Wald inference.
  blank <- c("DO_SE","DO_CI_low","DO_CI_high","DO_p_raw","DO_p_BH",
    "DO_OR_CI_low","DO_OR_CI_high","Host_Wald_chisq","Host_Wald_df","Host_p_raw","Host_p_BH",
    "Reads_Wald_chisq","Reads_Wald_df","Reads_p_raw","Reads_p_BH","alpha_SE",
    "AIC","BIC","delta_AIC_vs_A","delta_AIC_vs_best")
  for(nm in blank)r[[nm]] <- NA_real_
  if(isTRUE(r$numerical_valid))r$DO_OR <- exp(r$DO_beta)
  r
}

zero_penalty_seed <- function(mc,d,X,stub,output_dir) {
  # Old D is used only for initial parameters, never as a penalized result.
  path <- file.path(zero_penalty_seed_dir,mc,"D","fitted_candidate.rds")
  if(file.exists(path)) {
    candidate <- readRDS(path); f <- candidate$primary
    if(!is.null(f) && isTRUE(candidate$primary_valid)) {
      cols <- c("New_ID","Sample_ID","MC_ID","MC_reads","Other_reads","Ruegeria_reads_total",
        "DO_c","Host_type","Repeat_group_ID","TS","log_reads_c","mean_reads_s1","mean_reads_s2","mean_reads_s3")
      same <- nrow(f$data)==nrow(d) && all(vapply(cols,function(nm) {
        if(is.numeric(d[[nm]]))isTRUE(all.equal(as.numeric(d[[nm]]),as.numeric(f$data[[nm]]),tolerance=1e-12))
        else identical(as.character(d[[nm]]),as.character(f$data[[nm]]))
      },logical(1))) && identical(colnames(f$X),colnames(X)) &&
        isTRUE(all.equal(unname(f$X),unname(X),tolerance=1e-12))
      if(!same)stop("旧D起点与当前样本、响应、分母或设计不一致：",mc)
      writeLines(c("Verified data/design; initialization only",path),file.path(output_dir,"01_seed_source.txt"))
      return(list(parameters=f$parameters,nodes=f$quadrature_n,source=path))
    }
  }
  a <- joint_verify_cached(mc,d,X,stub,file.path(output_dir,"A_seed_verification"))
  if(!is.null(a) && isTRUE(a$primary_valid)) {
    writeLines(c("Verified A; initialization only",a$source),file.path(output_dir,"01_seed_source.txt"))
    return(list(parameters=a$primary$parameters,nodes=a$primary$quadrature_n,source=a$source))
  }
  p <- make_initial_parameters(d,X,stub)
  writeLines("Fresh normal starting values only; the fitted model remains skew-normal D.",
    file.path(output_dir,"01_seed_source.txt"))
  list(parameters=p,nodes=31L,source="fresh_initial_parameters")
}

zero_penalty_relabel_pi0 <- function(candidate,spec,note) {
  # Preserve any exact-zero Repeat/TS variance flags from the checked fit.
  actual <- utils::modifyList(candidate$primary$joint_spec,spec)
  actual$zi <- FALSE; actual$zi_reads <- FALSE
  candidate$spec <- actual; candidate$requested_spec <- spec
  candidate$primary$joint_spec <- actual
  candidate$primary$model_id <- paste0("JOINT_",spec$id,"_PI0")
  candidate$primary$penalty_requested <- TRUE
  candidate$primary$penalty_value <- 0
  candidate$primary$objective_value <- candidate$primary$nll
  if(!is.null(candidate$raw_fit)) {
    candidate$raw_fit$spec <- actual
    candidate$raw_fit$penalty_value <- 0
    candidate$raw_fit$data_nll <- candidate$raw_fit$value
  }
  candidate$boundary_note <- note
  candidate
}

zero_penalty_html <- function(summary,metrics,run_dir) {
  tab <- function(v) {
    heads <- paste0("<tr>",paste0("<th>",html_escape(names(v)),"</th>",collapse=""),"</tr>")
    body <- vapply(seq_len(nrow(v)),function(i)paste0("<tr>",paste0("<td>",
      vapply(v,function(a)html_escape(if(is.na(a[i]))"—" else as.character(a[i])),character(1)),
      "</td>",collapse=""),"</tr>"),character(1))
    paste0('<div class="scroll"><table><thead>',heads,'</thead><tbody>',paste(body,collapse=""),'</tbody></table></div>')
  }
  top <- data.frame(MC=summary$MC_ID,数值有效=summary$n_valid,已检查强度=summary$n_planned,
    无主要残差提示候选=summary$adequate_models,中间候选状态=summary$reference_status,
    有效候选DO范围=summary$DO_range_valid_candidates,复核事项=summary$next_step,check.names=FALSE)
  v <- data.frame(MC=metrics$MC_ID,候选=metrics$model,s=metrics$prior_sd_zi_reads,
    lambda=metrics$lambda_zi_reads,数值有效=metrics$numerical_valid,
    原始logLik=fmt_number(metrics$logLik,3),惩罚值=fmt_number(metrics$penalty_value,3),
    总目标=fmt_number(metrics$objective_value,3),零值reads斜率=fmt_number(metrics$zi_reads_slope),
    零值截距=fmt_number(metrics$zero_intercept),实际零值=metrics$n_zero_MC,
    模拟零值均值=fmt_number(metrics$expected_zero_mean,2),零值总数p=fmt_p(metrics$Zero_count_p),
    零值TS分配BH=fmt_p(metrics$Zero_TS_min_BH_p),零值reads分配BH=fmt_p(metrics$Zero_Reads_min_BH_p),
    DO点估计=fmt_number(metrics$DO_beta),波动reads斜率=fmt_number(metrics$sd_reads_slope),
    pi最低=fmt_number(metrics$pi_min),pi中心=fmt_number(metrics$pi_at_center),pi最高=fmt_number(metrics$pi_max),
    Uniformity=fmt_p(metrics$Uniformity_p),Dispersion=fmt_p(metrics$Dispersion_p),
    分位数=fmt_p(metrics$Predicted_quantiles_p),状态=metrics$status,残差提示=metrics$diagnostic_flags,
    说明=metrics$notes,check.names=FALSE)
  html <- c('<!doctype html><html lang="zh-CN"><meta charset="utf-8"><title>MC8 D模型零值斜率惩罚敏感性</title>',
    '<style>body{font:15px/1.6 Arial,"Microsoft YaHei",sans-serif;margin:24px;color:#193145}table{border-collapse:collapse}th,td{padding:8px 12px;border:1px solid #dbe3e9;white-space:nowrap;text-align:right}th{background:#193c55;color:white;position:sticky;top:0}tr:nth-child(even){background:#f6f8fa}.scroll{overflow:auto;max-height:70vh;margin:16px 0}td:nth-last-child(-n+3){white-space:normal;min-width:200px;text-align:left}</style>',
    '<h1>全部MC：D模型零值reads斜率的平方惩罚</h1>',
    paste0('<p>Q = 原始数据负对数似然 + lambda × gamma_reads²；lambda = 1/(2s²)。中间候选s=',
      zero_penalty_reference_sd,'、lambda=',1/(2*zero_penalty_reference_sd^2),'。每个MC使用同一套强度，各自估计参数。</p>'),
    '<p>保留D的reads均值样条、reads相关样本SD及中心化偏正态分布；仅惩罚额外零值reads斜率。本轮零值方程不含TS随机截距。</p>',
    '<p>本轮是强度敏感性比较，未运行分组交叉验证，未自动选择最优lambda。不同lambda的总目标不用于相互排名。普通AIC/BIC、科学效应的Wald区间/p值/BH留空；残差诊断p值仍用于检查模型。</p>',
    '<h2>每个MC的敏感性概览</h2>',tab(top),'<h2>全部惩罚候选</h2>',tab(v),
    '<p>无零值MC采用精确pi=0，并在三档设置中共用同一个计数模型；此时实际惩罚为0。分离记录保留为审计信息，正惩罚候选继续拟合并接受相同数值检查。</p>',
    '<p><a href="00_all_MC_penalty_summary.csv">每MC汇总</a> · <a href="01_all_model_metrics.csv">完整数值</a> · <a href="00_key_results.txt">中文摘要</a> · <a href="METHODS.txt">方法说明</a></p></html>')
  writeLines(enc2utf8(html),file.path(run_dir,"00_all_MC_penalty_summary.html"))
}

zero_penalty_outputs <- function(table,MC_ids,specs,run_dir,completed) {
  # A model which is numerically invalid never contributes scientific inference.
  infer <- c("DO_beta","DO_OR")
  table[!table$numerical_valid,infer] <- NA_real_
  summaries <- lapply(MC_ids,function(mc) {
    v <- table[table$MC_ID==mc,,drop=FALSE]
    valid <- which(v$numerical_valid)
    good <- which(v$numerical_valid & v$status=="数值有效；未检出主要残差偏差")
    ref <- which(v$is_reference_strength)
    b <- v$DO_beta[valid]; b <- b[is.finite(b)]
    data.frame(MC_ID=mc,n_planned=length(specs),n_valid=length(valid),
      valid_models=paste(v$model[valid],collapse=", "),adequate_models=paste(v$model[good],collapse=", "),
      reference_model=if(length(ref))v$model[ref[1]] else "",
      reference_status=if(length(ref))v$status[ref[1]] else "未设置",
      DO_range_valid_candidates=if(length(b))paste(format(range(b),digits=5),collapse=" to ") else "",
      DO_same_sign_valid_candidates=if(length(b)>1L)length(unique(sign(b)))==1L else NA,
      next_step=if(length(valid)<length(specs))"先复核失败候选；已完成结果继续保留" else
        if(!length(good))"继续定位零值、分位数或尾部偏差" else
          "结合强度敏感性，下一步规划分组预测验证与推断",stringsAsFactors=FALSE)
  })
  summary <- do.call(rbind,summaries)
  common <- names(specs)[vapply(names(specs),function(id) {
    v <- table[table$model==id,,drop=FALSE]
    nrow(v)==length(MC_ids) && all(v$numerical_valid) &&
      all(v$status=="数值有效；未检出主要残差偏差")
  },logical(1))]
  write_excel_csv(table,file.path(run_dir,"01_all_model_metrics.csv"))
  write_excel_csv(summary,file.path(run_dir,"00_all_MC_penalty_summary.csv"))
  zero_penalty_html(summary,table,run_dir)
  txt <- capture.output({
    cat("全部MC：D模型零值reads斜率的三档平方惩罚\n")
    cat("完成MC：",completed,"/",length(MC_ids),"；数值有效：",sum(table$numerical_valid),"/",nrow(table),"。\n",sep="")
    cat("本轮强度：\n")
    print(unique(table[,c("model","prior_sd_zi_reads","lambda_zi_reads","is_reference_strength")]),row.names=FALSE)
    cat("1. 每MC敏感性概览\n");print(summary,row.names=FALSE)
    cat("\n2. 原始似然、惩罚与总目标（不同lambda的Q不直接排名）\n")
    print(table[,c("MC_ID","model","lambda_zi_reads","numerical_valid","logLik","penalty_value","objective_value")],row.names=FALSE)
    cat("\n3. 参数、零值预测及批次分配\n")
    print(table[,c("MC_ID","model","zi_reads_slope","sd_reads_slope","pi_min","pi_at_center","pi_max",
      "n_zero_MC","expected_zero_mean","Zero_TS_min_BH_p","DO_beta","status")],row.names=FALSE)
    cat("\n4. 残差及边界复核事项\n");print(table[,c("MC_ID","model","diagnostic_flags","zero_slope_separation","notes")],row.names=FALSE)
    cat("\n在全部MC均数值有效且未检出主要残差偏差的强度：",if(length(common))paste(common,collapse=", ") else "目前没有", "\n",sep="")
    cat("本轮未做交叉验证，未选最优lambda，未给科学效应计算普通Wald区间/p值/BH。\n")
    cat("残差模拟条件于拟合的Repeat/TS，2000次模拟不重新拟合2000次。\n")
    cat("结果目录：",normalizePath(run_dir,winslash="/",mustWork=FALSE),"\n",sep="")
  })
  writeLines(enc2utf8(txt),file.path(run_dir,"00_key_results.txt"))
  invisible(list(table=table,summary=summary,text=txt,common_adequate_strengths=common))
}

# A single pool runs MC/lambda tasks. The no-zero MC uses one task for all
# strengths, since both its count fit and its seeded diagnostics are identical.
zero_penalty_worker_count <- function(n_tasks) {
  override <- Sys.getenv("GLMM_PARALLEL_WORKERS", unset="")
  requested <- if(nzchar(override))suppressWarnings(as.numeric(override)) else parallel_model_workers
  if(length(requested)!=1L || !is.finite(requested) || requested<1 || requested!=floor(requested))
    stop("parallel_model_workers / GLMM_PARALLEL_WORKERS 须为正整数。")
  cores <- parallel::detectCores(logical=TRUE)
  if(is.na(cores))cores <- 1L
  as.integer(max(1L,min(requested,max(1L,cores-2L),n_tasks)))
}

zero_penalty_clone_diagnostics <- function(diag,fit,source_dir,output_dir) {
  dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
  files <- list.files(source_dir,full.names=TRUE)
  files <- files[!dir.exists(files)]
  if(length(files) && !all(file.copy(files,output_dir,overwrite=TRUE)))
    stop("无法复制共用计数模型的诊断文件。")
  diag$fit_model <- fit
  diag$res_check$fittedModel <- fit
  diag$diagnostic_reuse <- list(source=source_dir,reason="Exact pi=0; identical count parameters, data and seed")
  writeLines(enc2utf8(c("共用已核查的精确 pi=0 计数模型，参数、数据、模拟种子相同。",
    "模拟、qgam和零值诊断共用一次计算；复制图中的模型标题标识原始共同计数拟合。",source_dir)),
    file.path(output_dir,"00_diagnostics_reused_from.txt"))
  writeLines(enc2utf8(capture.output(base::summary(fit))),file.path(output_dir,"01_model_summary.txt"))
  saveRDS(diag,file.path(output_dir,"model_analysis_objects.rds"))
  diag
}

zero_penalty_run_task <- function(job) {
  context <- .GLMM_penalty_task_contexts[[job$MC_ID]]
  specs <- .GLMM_penalty_task_specs
  mc <- job$MC_ID; d <- context$d; X <- context$X; stub <- context$stub
  mc_dir <- context$mc_dir; separation <- context$separation
  seed <- context$seed
  mc_rows <- all_zeros <- all_overall <- all_coef <- candidate_files <- list()
  exact_pi0 <- exact_diag <- NULL; exact_diag_dir <- ""
  for(model in job$models) {
    spec <- specs[[model]]; output_dir <- file.path(mc_dir,model)
    dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
    checkpoint <- file.path(output_dir,"fitted_candidate.rds")
    r <- zero_penalty_summary_row(mc,spec,context$n_source,d)
    r$zero_slope_separation <- separation
    started <- proc.time()[["elapsed"]]
    cat("\n",mc," / ",model,"；s=",spec$prior_sd_zi_reads,"；lambda=",spec$lambda_zi_reads,"\n",sep="")
    tryCatch({
      if(file.exists(checkpoint)) {
        candidate <- readRDS(checkpoint)
        if(!identical(candidate$requested_spec,spec))stop("检查点惩罚设置不一致：",model)
        cat("恢复本轮已保存候选。\n")
      } else if(!any(d$MC_reads==0) && !is.null(exact_pi0)) {
        candidate <- zero_penalty_relabel_pi0(exact_pi0,spec,
          "无观测零值：共用本轮已核查的D计数模型，精确pi=0，实际惩罚=0")
      } else {
        if(!is.null(seed$error))stop(seed$error)
        actual <- spec
        if(!any(d$MC_reads==0)) {actual$zi <- FALSE; actual$zi_reads <- FALSE}
        p <- seed$parameters
        for(nm in c("log_sd_repeat","log_sd_ts"))if(!is.finite(p[[nm]]))p[[nm]] <- -6
        # NULL is intentional: outer tasks share one pool, never nested pools.
        candidate <- joint_fit_pass(p,d,X,stub,actual,output_dir,pool=NULL,nodes=max(31L,seed$nodes))
        candidate <- joint_variance_boundaries(candidate,d,X,stub,output_dir,pool=NULL)
        f <- candidate$primary
        if(!is.null(f) && isTRUE(f$zi) && max(joint_row_values(f)$pi)<1e-6) {
          nozi_spec <- actual; nozi_spec$zi <- FALSE; nozi_spec$zi_reads <- FALSE
          nozi <- joint_fit_pass(f$parameters,d,X,stub,nozi_spec,file.path(output_dir,"pi0_count_boundary"),
            pool=NULL,nodes=f$quadrature_n)
          nozi <- joint_variance_boundaries(nozi,d,X,stub,file.path(output_dir,"pi0_count_boundary"),pool=NULL)
          candidate <- joint_pi0(candidate,nozi,d,X,output_dir)
        }
        candidate$requested_spec <- spec
        if(!is.null(candidate$primary) && !isTRUE(candidate$primary$zi))
          candidate <- zero_penalty_relabel_pi0(candidate,spec,
            if(!any(d$MC_reads==0))"无观测零值：精确pi=0，实际惩罚=0" else "经同一惩罚目标边界核查取精确pi=0")
      }
      saveRDS(candidate,checkpoint); candidate_files[[model]] <- checkpoint
      if(!any(d$MC_reads==0) && isTRUE(candidate$primary_valid))exact_pi0 <- candidate
      r <- zero_penalty_finish_row(r,candidate)
      r$fit_source <- if(!any(d$MC_reads==0))"checked_exact_pi0_count_D" else "new_penalized_D"
      if(r$numerical_valid) {
        diagnostic_dir <- file.path(output_dir,"diagnostics");dir.create(diagnostic_dir,showWarnings=FALSE)
        diagnostic_file <- file.path(diagnostic_dir,"model_analysis_objects.rds")
        if(file.exists(diagnostic_file))diag <- readRDS(diagnostic_file) else if(
          !any(d$MC_reads==0) && !is.null(exact_diag) && isTRUE(exact_diag$diagnostics_complete)) {
          diag <- zero_penalty_clone_diagnostics(exact_diag,candidate$primary,exact_diag_dir,diagnostic_dir)
          cat("共用同一精确pi=0模型的已完成诊断。\n")
        } else {
          diagnostic_started <- proc.time()[["elapsed"]]
          cat("数值检查通过；",n_sim,"次条件模拟及残差诊断。\n",sep="");flush.console()
          diag <- diagnose_skew_comparison_model(paste(mc,model),candidate$primary,d,n_sim,
            simulation_seed+context$index*100L,diagnostic_dir)
          diag <- zero_count_diagnostics(diag,d,diagnostic_dir)
          saveRDS(diag,diagnostic_file)
          cat("残差及零值诊断完成；",round(proc.time()[["elapsed"]]-diagnostic_started,1),"秒\n",sep="");flush.console()
        }
        if(!any(d$MC_reads==0)) {exact_diag <- diag; exact_diag_dir <- diagnostic_dir}
        r <- fill_zero_diagnostic_row(r,diag)
        r$diagnostic_flag_count <- if(nzchar(r$diagnostic_flags))length(strsplit(r$diagnostic_flags,"; ",fixed=TRUE)[[1]]) else 0L
        key <- paste(mc,model)
        all_zeros[[key]] <- cbind(MC_ID=mc,model=model,lambda=spec$lambda_zi_reads,diag$zero_group_table)
        all_overall[[key]] <- cbind(MC_ID=mc,model=model,lambda=spec$lambda_zi_reads,diag$overall_table)
        all_coef[[key]] <- cbind(model=model,lambda=spec$lambda_zi_reads,
          inference_method="point_estimates_only",fixed_coefficient_table(candidate$primary,mc))
        rv <- joint_row_values(candidate$primary)
        write_result_table(data.frame(New_ID=d$New_ID,log_reads_c=d$log_reads_c,pi=rv$pi,
          Observation_SD=rv$sd,expected_proportion=diag$fitted_rows$expected_MC_proportion),
          output_dir,"03_row_specific_pi_and_SD.csv")
      }
      if(!is.null(candidate$primary))writeLines(enc2utf8(capture.output(base::summary(candidate$primary))),
        file.path(output_dir,"02_model_summary.txt"))
    },error=function(e) {
      r$status <<- if(r$numerical_valid)"数值有效；诊断失败" else "处理失败"
      r$notes <<- paste(Filter(nzchar,c(r$notes,substr(conditionMessage(e),1,1000))),collapse="; ")
      writeLines(enc2utf8(conditionMessage(e)),file.path(output_dir,"00_error.txt"))
      cat("已记录错误：",conditionMessage(e),"\n",sep="")
    })
    mc_rows[[model]] <- r
    write_excel_csv(r,file.path(output_dir,"00_candidate_metrics.csv"))
    cat("完成：",r$status,"；耗时",round(proc.time()[["elapsed"]]-started,1),"秒\n",sep="");flush.console()
    gc(FALSE)
  }
  list(MC_ID=mc,models=job$models,metrics=do.call(rbind,mc_rows),candidate_files=candidate_files,
    zero_rows=all_zeros,overall_rows=all_overall,coefficient_rows=all_coef)
}

zero_penalty_task_worker <- function(job) {
  started <- Sys.time()
  log_connection <- file(job$log_file,open="at",encoding="UTF-8")
  output_depth <- sink.number(); sink(log_connection)
  on.exit({while(sink.number()>output_depth)sink();close(log_connection)},add=TRUE)
  writeLines(c("running",format(started),paste("PID",Sys.getpid())),job$status_file)
  result <- tryCatch(zero_penalty_run_task(job),error=function(e) {
    context <- .GLMM_penalty_task_contexts[[job$MC_ID]]
    rows <- lapply(job$models,function(m) {
      r <- zero_penalty_summary_row(job$MC_ID,.GLMM_penalty_task_specs[[m]],context$n_source,context$d)
      r$status <- "任务失败";r$notes <- conditionMessage(e);r
    })
    list(MC_ID=job$MC_ID,models=job$models,metrics=do.call(rbind,rows),candidate_files=list(),
      zero_rows=list(),overall_rows=list(),coefficient_rows=list(),error=conditionMessage(e))
  })
  result$elapsed_seconds <- as.numeric(difftime(Sys.time(),started,units="secs"))
  result$worker_pid <- Sys.getpid()
  # The master reads only results already returned by the worker. A temporary
  # sibling plus rename also protects recovery after interruption while saving.
  temporary <- paste0(job$result_file,".",Sys.getpid(),".tmp")
  saveRDS(result,temporary)
  if(!file.rename(temporary,job$result_file))stop("无法提交任务检查点：",job$result_file)
  writeLines(c("complete",format(Sys.time()),paste("elapsed_seconds",result$elapsed_seconds)),job$status_file)
  result
}

zero_penalty_dispatch <- function(pool,jobs,on_result) {
  n <- length(jobs);if(!n)return(invisible(NULL))
  # This small queue follows base R parallel's dynamicClusterApply protocol.
  # Receive one finished task, refill that worker immediately, then checkpoint
  # on the master. Statistical/plotting objects never share AD tapes.
  submit <- function(node,j)parallel:::sendCall(pool[[node]],zero_penalty_task_worker,list(jobs[[j]]),tag=j)
  active <- min(n,length(pool));next_job <- active+1L
  for(j in seq_len(active))submit(j,j)
  for(done in seq_len(n)) {
    received <- parallel:::recvOneResult(pool)
    if(inherits(received$value,"try-error"))stop(as.character(received$value))
    if(next_job<=n) {submit(received$node,next_job);next_job <- next_job+1L}
    on_result(received$tag,received$value)
  }
  invisible(NULL)
}

zero_penalty_collect_outputs <- function(results,contexts,MC_ids,specs,shared,run_dir,settings) {
  available <- Filter(Negate(is.null),results)
  metrics <- bind_flexible_tables(lapply(available,`[[`,"metrics"))
  zero_rows <- unlist(lapply(available,`[[`,"zero_rows"),recursive=FALSE)
  overall_rows <- unlist(lapply(available,`[[`,"overall_rows"),recursive=FALSE)
  coefficient_rows <- unlist(lapply(available,`[[`,"coefficient_rows"),recursive=FALSE)
  rows <- model_files <- list();completed_mc <- 0L
  for(mc in MC_ids) {
    context <- contexts[[mc]]
    current <- if(nrow(metrics))metrics[metrics$MC_ID==mc,,drop=FALSE] else data.frame()
    for(m in names(specs))if(!nrow(current) || !m%in%current$model)
      current <- bind_flexible_tables(list(current,zero_penalty_summary_row(mc,specs[[m]],context$n_source,context$d)))
    current <- current[match(names(specs),current$model),,drop=FALSE]
    rows[[mc]] <- current
    done_models <- if(nrow(metrics))metrics$model[metrics$MC_ID==mc] else character()
    if(all(names(specs)%in%done_models))completed_mc <- completed_mc+1L
    candidate_files <- unlist(lapply(available[vapply(available,function(a)identical(a$MC_ID,mc),logical(1))],
      `[[`,"candidate_files"),recursive=FALSE,use.names=TRUE)
    candidates <- lapply(candidate_files,readRDS)
    if(length(candidates)) {
      model_files[[mc]] <- file.path(context$mc_dir,"fitted_comparison.rds")
      saveRDS(list(MC_ID=mc,candidates=candidates,metrics=current),model_files[[mc]])
    }
    marker <- file.path(context$mc_dir,"05_curve_plots_complete.txt")
    if(all(names(specs)%in%done_models) && !file.exists(marker)) {
      curves <- bind_flexible_tables(lapply(names(candidates),function(m) {
        cc <- candidates[[m]];if(!isTRUE(cc$primary_valid))return(NULL)
        x <- seq(min(context$d$log_reads_c),max(context$d$log_reads_c),length.out=150)
        rv <- joint_row_values(cc$primary,x)
        data.frame(MC_ID=mc,model=m,log10_reads=x+shared$log_reads_center,pi=rv$pi,SD=rv$sd)
      }))
      if(nrow(curves)) {
        write_result_table(curves,context$mc_dir,"05_zero_probability_and_SD_curves.csv")
        plot_data <- rbind(data.frame(curves[,c("model","log10_reads")],component="Extra-zero probability",value=curves$pi),
          data.frame(curves[,c("model","log10_reads")],component="Observation SD (logit)",value=curves$SD))
        fig <- ggplot2::ggplot(plot_data,ggplot2::aes(log10_reads,value,color=model))+
          ggplot2::geom_line(linewidth=.85)+ggplot2::facet_wrap(~component,scales="free_y",ncol=1)+
          ggplot2::theme_bw(base_size=12)+ggplot2::labs(title=mc,x="Total Ruegeria reads (log10)",y=NULL,color="Penalty candidate")
        ggplot2::ggsave(file.path(context$mc_dir,"05_zero_probability_and_SD_curves.png"),fig,width=7,height=7,dpi=240,bg="white")
        ggplot2::ggsave(file.path(context$mc_dir,"05_zero_probability_and_SD_curves.pdf"),fig,width=7,height=7)
      }
      writeLines("complete",marker)
    }
  }
  final <- zero_penalty_outputs(do.call(rbind,rows),MC_ids,specs,run_dir,completed_mc)
  write_excel_csv(bind_flexible_tables(zero_rows),file.path(run_dir,"02_all_zero_group_diagnostics.csv"))
  write_excel_csv(bind_flexible_tables(overall_rows),file.path(run_dir,"03_all_overall_diagnostics.csv"))
  write_excel_csv(bind_flexible_tables(coefficient_rows),file.path(run_dir,"04_all_fixed_coefficients.csv"))
  saveRDS(list(metrics=final$table,summary=final$summary,model_files=model_files,completed=completed_mc,
    completed_models=nrow(metrics),settings=settings),file.path(run_dir,"GLMM_zero_penalty_results.rds"))
  final
}

run_all_MC_zero_penalty <- function() {
  previous <- options(contrasts=c("contr.treatment","contr.poly"))
  on.exit(options(previous),add=TRUE)
  packages <- c("readxl","glmmTMB","RTMB","DHARMa","qgam","ggplot2")
  missing <- packages[!vapply(packages,requireNamespace,logical(1),quietly=TRUE)]
  if(length(missing))stop("缺少包：",paste(missing,collapse=", "))
  specs <- zero_penalty_specs()
  raw <- readxl::read_excel(input_file,sheet=input_sheet)
  MC_ids <- natural_MC_order(unique(as.character(raw$MC_ID[!is.na(raw$MC_ID)])))
  if(!is.null(MC_selection))MC_ids <- natural_MC_order(unique(MC_selection))
  shared <- build_shared_design(raw,MC_ids,host_reference)
  src <- getOption("GLMM_JOINT_worker_source",file.path(project_dir,"scripts","GLMM_reads_diagnostics_full.R"))
  if(!file.exists(src))stop("工作进程脚本不存在：",src)
  resuming <- nzchar(zero_penalty_resume_dir)
  run_dir <- if(resuming)zero_penalty_resume_dir else file.path(output_root,
    paste0("GLMM_D_zero_penalty_MC8_",format(Sys.time(),"%Y%m%d_%H%M%S")))
  if(!resuming && dir.exists(run_dir))run_dir <- paste0(run_dir,"_",Sys.getpid())
  worker_count <- zero_penalty_worker_count(length(MC_ids)*length(specs))
  settings <- list(experiment="D_zero_slope_ridge_v1",input_file=input_file,
    input_md5=unname(tools::md5sum(input_file)),MC_ids=MC_ids,n_sim=n_sim,seed=simulation_seed,
    mean_df=3L,specs=specs,reference_sd=zero_penalty_reference_sd,seed_dir=zero_penalty_seed_dir,
    maximum_fit_nodes=maximum_fit_nodes,parallel_model_workers=worker_count,
    scheduling="flat_MC_lambda_pool_v2; six serial starts per task; fixed covariance only",
    batch_sensitivity=FALSE,zero_TS_random_effect=FALSE,CV=FALSE,
    script_md5=unname(tools::md5sum(src)))
  if(resuming) {
    old <- readRDS(file.path(run_dir,"01_run_settings.rds"))
    # Changing only the process count via the environment does not change a fit.
    for(nm in setdiff(names(settings),"parallel_model_workers"))
      if(!identical(settings[[nm]],old[[nm]]))stop("恢复运行设置不一致：",nm)
  }
  dir.create(run_dir,recursive=TRUE,showWarnings=FALSE)
  saveRDS(settings,file.path(run_dir,"01_run_settings.rds"))
  saveRDS(shared,file.path(run_dir,"01_shared_model_design.rds"))
  write_excel_csv(shared$samples,file.path(run_dir,"01_shared_sample_design.csv"))
  if(!resuming) {
    if(!file.copy(src,file.path(run_dir,"01_analysis_script_used.R")))stop("无法保存运行脚本快照")
    if(!file.copy(input_file,file.path(run_dir,"01_input_used_MC8.xlsx")))stop("无法保存输入快照")
  }
  manifest <- data.frame(role=c("script","input"),source=c(src,input_file),
    snapshot=c("01_analysis_script_used.R","01_input_used_MC8.xlsx"),
    md5=c(settings$script_md5,settings$input_md5),stringsAsFactors=FALSE)
  stopifnot(identical(unname(tools::md5sum(file.path(run_dir,manifest$snapshot))),manifest$md5))
  write_excel_csv(manifest,file.path(run_dir,"01_input_script_manifest.csv"))
  grid <- do.call(rbind,lapply(specs,function(s)data.frame(model=s$id,label=s$strength_label,
    prior_sd=s$prior_sd_zi_reads,lambda=s$lambda_zi_reads,reference=s$is_reference)))
  write_excel_csv(grid,file.path(run_dir,"01_penalty_strengths.csv"))
  writeLines("running",file.path(run_dir,"00_run_status.txt"))
  completed_run <- FALSE;pool <- NULL;contexts <- list();jobs <- list();results <- list();final <- NULL
  on.exit({
    if(!is.null(pool))try(parallel::stopCluster(pool),silent=TRUE)
    if(!completed_run) {
      # Recover only atomically committed tasks, after stopping the workers.
      if(length(results)<length(jobs))length(results) <- length(jobs)
      for(j in seq_along(jobs))if(is.null(results[[j]]) && file.exists(jobs[[j]]$result_file))
        results[j] <- list(tryCatch(readRDS(jobs[[j]]$result_file),error=function(e)NULL))
      if(length(contexts)==length(MC_ids) && length(results))try(
        zero_penalty_collect_outputs(results,contexts,MC_ids,specs,shared,run_dir,settings),silent=TRUE)
      writeLines("interrupted_or_failed; inspect partial outputs and per-task logs",file.path(run_dir,"00_run_status.txt"))
    }
  },add=TRUE)
  writeLines(run_dir,file.path(tempdir(),"GLMM_zero_penalty_current_run.txt"))
  if(nzchar(Sys.getenv("GLMM_ZERO_PENALTY_POINTER")))writeLines(run_dir,Sys.getenv("GLMM_ZERO_PENALTY_POINTER"))
  writeLines(capture.output(utils::sessionInfo()),file.path(run_dir,"sessionInfo.txt"))
  writeLines(enc2utf8(c(
    "D model; Q(theta) = total negative marginal data log likelihood + lambda * gamma_reads^2.",
    "lambda=1/(2*s^2), not divided by sample size; x=centered log10 total Ruegeria reads, not SD-standardized.",
    "Only gamma_reads is penalized. Original count mean, samples, denominator, zeros, spline df=3 and skew-normal family retained.",
    "D retains reads-dependent Observation SD. No zero-equation TS random effect, CV, automatic optimum or ordinary AIC/Wald/BH inference.",
    "Data-only nll = total penalized objective minus the fixed coefficient penalty; all three quantities saved separately.",
    "Same six starts, at least two valid starts within .001, original optimizer limits, gradient/fixed-Hessian/covariance and integration checks.",
    "Fixed covariance uses the same optimHess/chol/solve calculation as RTMB::sdreport; unused ADreport/random-effect SEs omitted.",
    "Near-zero pi and Repeat/TS variances checked against the SAME penalized objective. Old fits only initialize after data/design verification.",
    "Single flat MC/lambda process pool; no nested pool; one BLAS/OpenMP thread per process; base-R load-balanced socket queue.",
    "No observed zeros -> exact pi=0 and penalty=0. All strengths share one checked count fit and completed identical seeded diagnostics.",
    "Valid candidates retain 2000 (or configured n_sim) conditional simulations and all residual/zero-group checks.",
    "Each completed task updates global partial summaries. Each task has a log/status/result checkpoint; each start reports elapsed seconds.",
    "Inputs read-only. Source/input snapshots and settings stored. Resume requires matching script/input/settings.",
    "No sidebar, browser or interactive plot window opened.")),file.path(run_dir,"METHODS.txt"))
  for(i in seq_along(MC_ids)) {
    mc <- MC_ids[i];mc_dir <- file.path(run_dir,mc);dir.create(mc_dir,recursive=TRUE,showWarnings=FALSE)
    cat("准备 ",i,"/",length(MC_ids),"：",mc," 的共享数据、设计和起点。\n",sep="");flush.console()
    prepared <- prepare_mc_shared(raw,mc,shared);d <- prepared$data_test;X <- prepared$X;stub <- make_model_stub(d)
    write_result_table(d,mc_dir,"01_model_input_rows.csv")
    write_result_table(prepared$sample_audit,mc_dir,"01_sample_audit.csv")
    separation <- joint_zero_reads_separation(d)
    writeLines(enc2utf8(if(nzchar(separation))paste("Unpenalized separation audit:",separation,
      "Positive-penalty candidates continue to fitting.") else "No strict zero/positive reads separation detected."),
      file.path(mc_dir,"01_zero_reads_separation_audit.txt"))
    seed <- tryCatch(zero_penalty_seed(mc,d,X,stub,mc_dir),error=function(e)list(error=conditionMessage(e)))
    contexts[[mc]] <- list(index=i,d=d,X=X,stub=stub,mc_dir=mc_dir,separation=separation,seed=seed,
      n_source=sum(as.character(raw$MC_ID)==mc,na.rm=TRUE))
    groups <- if(any(d$MC_reads==0))lapply(names(specs),function(m)m) else list(names(specs))
    for(models in groups) {
      job_dir <- file.path(mc_dir,models[1]);dir.create(job_dir,recursive=TRUE,showWarnings=FALSE)
      jobs[[length(jobs)+1L]] <- list(MC_ID=mc,models=models,
        log_file=file.path(job_dir,"00_compute_log.txt"),status_file=file.path(job_dir,"00_task_status.txt"),
        result_file=file.path(job_dir,"00_task_result.rds"))
    }
  }
  results <- vector("list",length(jobs))
  for(j in seq_along(jobs))if(file.exists(jobs[[j]]$result_file))results[j] <- list(readRDS(jobs[[j]]$result_file))
  pending <- which(vapply(results,is.null,logical(1)))
  progress <- function() {
    data.frame(task=seq_along(jobs),MC_ID=vapply(jobs,`[[`,character(1),"MC_ID"),
      models=vapply(jobs,function(j)paste(j$models,collapse="; "),character(1)),
      state=vapply(seq_along(jobs),function(j)if(!is.null(results[[j]]))"complete" else
        if(file.exists(jobs[[j]]$status_file))readLines(jobs[[j]]$status_file,n=1L,warn=FALSE) else "queued",character(1)),
      elapsed_seconds=vapply(results,function(r)if(is.null(r))NA_real_ else r$elapsed_seconds,numeric(1)))
  }
  update_outputs <- function() {
    write_excel_csv(progress(),file.path(run_dir,"00_parallel_task_progress.csv"))
    zero_penalty_collect_outputs(results,contexts,MC_ids,specs,shared,run_dir,settings)
  }
  final <- update_outputs()
  cat("\n启动 ",min(worker_count,length(pending))," 个计算进程；",length(jobs)," 个独立 MC/lambda 任务（无零值 MC 的三档共用一个任务）。\n",
    "任务内保留六起点，不建立子进程池。详细进度见各候选的 00_compute_log.txt；任务完成后控制台及全局报告立即更新。\n",
    "输出：",run_dir,"\n",sep="");flush.console()
  if(length(pending)) {
    # Child libraries must see thread limits before R or its BLAS is loaded.
    thread_names <- c("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","BLIS_NUM_THREADS","VECLIB_MAXIMUM_THREADS")
    thread_old <- Sys.getenv(thread_names,unset=NA_character_)
    restore_threads <- function() {
      missing <- is.na(thread_old)
      if(any(missing))Sys.unsetenv(thread_names[missing])
      if(any(!missing))do.call(Sys.setenv,as.list(stats::setNames(thread_old[!missing],thread_names[!missing])))
    }
    do.call(Sys.setenv,as.list(stats::setNames(rep("1",length(thread_names)),thread_names)))
    pool <- tryCatch(parallel::makePSOCKcluster(min(worker_count,length(pending)),
      outfile=file.path(run_dir,"00_worker_messages.log")),finally=restore_threads())
    parallel::clusterCall(pool,function(path) {
      options(GLMM_ZI_definitions_only=TRUE,GLMM_JOINT_definitions_only=TRUE,contrasts=c("contr.treatment","contr.poly"))
      source(path,encoding="UTF-8");requireNamespace("RTMB",quietly=TRUE);TRUE
    },src)
    export <- new.env(parent=emptyenv())
    export$.GLMM_penalty_task_contexts <- contexts;export$.GLMM_penalty_task_specs <- specs
    parallel::clusterExport(pool,c(".GLMM_penalty_task_contexts",".GLMM_penalty_task_specs"),envir=export)
    zero_penalty_dispatch(pool,jobs[pending],function(index,result) {
      j <- pending[index];results[j] <<- list(result)
      cat("\n已完成任务 ",sum(!vapply(results,is.null,logical(1))),"/",length(jobs),"：",
        result$MC_ID," / ",paste(result$models,collapse=", "),"；",round(result$elapsed_seconds/60,1)," 分钟\n",sep="")
      print(result$metrics[,c("MC_ID","model","status","logLik","penalty_value","objective_value"),drop=FALSE],row.names=FALSE)
      flush.console()
      final <<- update_outputs()
    })
    parallel::stopCluster(pool);pool <- NULL
  }
  stopifnot(identical(unname(tools::md5sum(input_file)),settings$input_md5),
    identical(unname(tools::md5sum(src)),settings$script_md5),nrow(final$table)==length(MC_ids)*length(specs))
  writeLines(run_dir,file.path(output_root,"00_latest_GLMM_D_zero_penalty_MC8.txt"))
  completed_run <- TRUE
  writeLines(c("complete",run_dir),file.path(run_dir,"00_run_status.txt"))
  cat("\n",paste(final$text,collapse="\n"),"\n",sep="")
  invisible(list(output_dir=run_dir,settings=settings,table=final$table,summary=final$summary,
    common_adequate_strengths=final$common_adequate_strengths))
}

# ===== Active experiment: common structural screening suite, version 1 =====
# These overrides preserve the original D likelihood exactly when no new feature
# is active. Legacy functions remain available above for audit/reproduction.
suite_version <- "MC_abundance_structure_v2"
suite_reference_cache <- file.path(output_root, "GLMM_D_zero_penalty_MC8_20261012_160110")
suite_structure_reference <- file.path(output_root, "GLMM_structure_suite_MC8_20261012_185348")
suite_resume_dir <- Sys.getenv("GLMM_STRUCTURE_RESUME_DIR", unset = "")
suite_lambda_reads <- 0.125
suite_zero_DO_prior_sd <- 2       # logit extra-zero slope per +1 mg/L DO.
suite_zero_TS_prior_sd <- 1.5     # FIXED Gaussian SD, not an estimated variance.
suite_tail_weight <- 0.05        # Prespecified rare broad component.
suite_tail_SD_ratio <- 3         # Its SD / core SD; total SD remains sigma(x).
suite_run_top1_sensitivity <- FALSE # 上轮已做；本轮优先检验共同结构。
suite_save_plots <- TRUE
suite_include_previous_candidates <- FALSE # TRUE可同时重访上轮七个候选。
suite_worker_memory_gb <- 1.0 # 每个RTMB进程预留；启动时另外保留1.5 GB。

suite_legacy_objective <- joint_objective
suite_legacy_parameters <- joint_parameters
suite_legacy_random <- joint_random
suite_legacy_row_values <- joint_row_values
suite_legacy_integrals <- joint_integrals
suite_legacy_simulate <- simulate_centered_olre
suite_legacy_interior <- joint_interior
suite_legacy_starts <- joint_starts

suite_tail_components <- function(spec) {
  if (!isTRUE(spec$tail_mix)) return(list(weight = 1, scale = 1))
  w <- spec$tail_weight; k <- spec$tail_ratio
  if (!is.finite(w) || w <= 0 || w >= 1 || !is.finite(k) || k <= 1)
    stop("Invalid prespecified tail mixture.")
  list(weight = c(1-w, w), scale = c(1, k) / sqrt(1-w+w*k^2))
}

joint_penalty_value <- function(parameters, spec) {
  ans <- 0
  if (isTRUE(spec$zi) && isTRUE(spec$zi_reads) && !is.null(spec$lambda_zi_reads))
    ans <- ans + spec$lambda_zi_reads * parameters$zi_reads^2
  if (isTRUE(spec$zi) && isTRUE(spec$zi_DO))
    ans <- ans + parameters$zi_DO^2 / (2*spec$zero_DO_prior_sd^2)
  ans
}

joint_parameters <- function(p, d, spec) {
  original <- p
  p <- suite_legacy_parameters(p, d, spec)
  if (isTRUE(spec$sd_quadratic)) p$sd_reads2 <- if (is.null(original$sd_reads2)) 0 else original$sd_reads2
  if (isTRUE(spec$sd_mean)) p$sd_mean <- original$sd_mean %||% 0
  if (isTRUE(spec$alpha_mean)) p$alpha_mean <- original$alpha_mean %||% 0
  if (isTRUE(spec$zi) && isTRUE(spec$zi_DO))
    p$zi_DO <- if (is.null(original$zi_DO)) 0 else original$zi_DO
  if (isTRUE(spec$zi) && isTRUE(spec$zi_TS)) {
    p$b_zi_ts <- original$b_zi_ts
    if (is.null(p$b_zi_ts) || length(p$b_zi_ts) != nlevels(d$TS)) p$b_zi_ts <- rep(0, nlevels(d$TS))
  }
  p
}
joint_random <- function(spec) c(suite_legacy_random(spec),
  if (isTRUE(spec$zi) && isTRUE(spec$zi_TS)) "b_zi_ts")

joint_objective <- function(d, X, nodes, spec) {
  new_feature <- isTRUE(spec$sd_quadratic) || isTRUE(spec$zi_DO) ||
    isTRUE(spec$zi_TS) || isTRUE(spec$tail_mix) ||
    isTRUE(spec$sd_mean) || isTRUE(spec$alpha_mean)
  if (!new_feature) return(suite_legacy_objective(d, X, nodes, spec))
  y <- as.numeric(d$MC_reads); N <- as.numeric(d$Ruegeria_reads_total)
  x <- as.numeric(d$log_reads_c); DO <- as.numeric(d$DO_c)
  ir <- as.integer(d$Repeat_group_ID); it <- as.integer(d$TS)
  observed_logit <- qlogis((y+.5)/(N+1)); gh <- normal_quadrature(nodes)
  component <- suite_tail_components(spec)
  # Xc is fixed data; m = Xc beta is recomputed on every objective evaluation.
  # It excludes count Repeat/TS offsets and the Observation perturbation.
  Xc <- sweep(X, 2L, colMeans(X), "-")
  function(p) {
    eta <- as.vector(X %*% p$beta) + p$b_repeat[ir] + p$b_ts[it]
    m <- as.vector(Xc %*% p$beta)
    nll <- 0
    if (!isTRUE(spec$zero_repeat)) nll <- nll-sum(RTMB::dnorm(p$b_repeat,0,exp(p$log_sd_repeat),log=TRUE))
    if (!isTRUE(spec$zero_ts)) nll <- nll-sum(RTMB::dnorm(p$b_ts,0,exp(p$log_sd_ts),log=TRUE))
    if (isTRUE(spec$zi) && isTRUE(spec$zi_TS))
      nll <- nll-sum(RTMB::dnorm(p$b_zi_ts,0,spec$zero_TS_prior_sd,log=TRUE))
    ls <- p$log_sd_observation + if (isTRUE(spec$sd_reads)) p$sd_reads*x else 0
    if (isTRUE(spec$sd_quadratic)) ls <- ls+p$sd_reads2*x^2
    if (isTRUE(spec$sd_mean)) ls <- ls+p$sd_mean*m
    alpha <- p$alpha + if(isTRUE(spec$alpha_mean))p$alpha_mean*m else 0
    delta <- alpha/sqrt(1+alpha^2)
    for (k in seq_along(component$weight)) {
      omega <- exp(ls)*component$scale[k]/sqrt(1-2*delta^2/pi)
      xi <- -omega*delta*sqrt(2/pi); eta0 <- eta+xi
      zmode <- 10*tanh(((observed_logit-eta0)/omega)/10)
      for (iteration in 1:30) {
        probability <- RTMB::plogis(eta0+omega*zmode); az <- alpha*zmode
        mills <- exp(RTMB::dnorm(az,log=TRUE)-RTMB::pnorm(az,log.p=TRUE))
        curvature <- 1+omega^2*N*probability*(1-probability)+alpha^2*mills*(az+mills)
        step <- (omega*(y-N*probability)-zmode+alpha*mills)/curvature
        zmode <- zmode+3*tanh(step/3)
      }
      probability <- RTMB::plogis(eta0+omega*zmode); az <- alpha*zmode
      mills <- exp(RTMB::dnorm(az,log=TRUE)-RTMB::pnorm(az,log.p=TRUE))
      zsd <- 1/sqrt(1+omega^2*N*probability*(1-probability)+alpha^2*mills*(az+mills))
      for (j in seq_along(gh$nodes)) {
        z <- zmode+zsd*gh$nodes[j]
        term <- log(gh$weights[j])+RTMB::dbinom_robust(y,N,eta0+omega*z,log=TRUE)+
          RTMB::dnorm(z,0,1,log=TRUE)+log(2)+RTMB::pnorm(alpha*z,log.p=TRUE)-
          RTMB::dnorm(z,zmode,zsd,log=TRUE)
        if (j==1L) one_ll <- term else one_ll <- RTMB::logspace_add(one_ll,term)
      }
      one_ll <- one_ll+log(component$weight[k])
      if (k==1L) ll <- one_ll else ll <- RTMB::logspace_add(ll,one_ll)
    }
    if (isTRUE(spec$zi)) {
      g <- p$gamma_zi+if(isTRUE(spec$zi_reads))p$zi_reads*x else rep(0,length(x))
      if (isTRUE(spec$zi_DO)) g <- g+p$zi_DO*DO
      if (isTRUE(spec$zi_TS)) g <- g+p$b_zi_ts[it]
      lp <- -RTMB::logspace_add(0,-g); lq <- -RTMB::logspace_add(0,g)
      positive <- which(y>0); zero <- which(y==0)
      ll[positive] <- lq[positive]+ll[positive]
      ll[zero] <- RTMB::logspace_add(lp[zero],lq[zero]+ll[zero])
    }
    nll-sum(ll)+joint_penalty_value(p,spec)
  }
}

joint_row_values <- function(fit, x=fit$data$log_reads_c, DO_c=NULL, TS=NULL, X=NULL) {
  p <- fit$parameters; spec <- fit$joint_spec
  if (is.null(spec)) return(suite_legacy_row_values(fit,x))
  in_sample <- is.null(X) && identical(as.numeric(x),as.numeric(fit$data$log_reads_c))
  if (is.null(DO_c)) DO_c <- if(in_sample)fit$data$DO_c else rep(0,length(x))
  if (length(DO_c)==1L) DO_c <- rep(DO_c,length(x))
  if (length(DO_c)!=length(x)) stop("Prediction DO length does not match reads.")
  ls <- rep(p$log_sd_observation,length(x))
  if (isTRUE(spec$sd_reads)) ls <- ls+p$sd_reads*x
  if (isTRUE(spec$sd_quadratic)) ls <- ls+p$sd_reads2*x^2
  if(is.null(X)) {
    if(!in_sample && (isTRUE(spec$sd_mean)||isTRUE(spec$alpha_mean)))
      stop("Mean-linked prediction requires the full fixed-effect design X.")
    if(in_sample)X <- fit$X
  }
  if(!is.null(X)) {
    if(nrow(X)!=length(x)||!identical(colnames(X),colnames(fit$X)))
      stop("Prediction fixed-effect design does not match fitted columns.")
    m <- as.vector(sweep(X,2L,colMeans(fit$X),"-")%*%p$beta)
    fixed_eta <- as.vector(X%*%p$beta)
  } else {m <- rep(0,length(x));fixed_eta <- rep(NA_real_,length(x))}
  if(isTRUE(spec$sd_mean))ls <- ls+p$sd_mean*m
  alpha <- rep(p$alpha,length(x))
  if(isTRUE(spec$alpha_mean))alpha <- alpha+p$alpha_mean*m
  g <- rep(if(isTRUE(spec$zi))p$gamma_zi else -Inf,length(x))
  if (isTRUE(spec$zi) && isTRUE(spec$zi_reads)) g <- g+p$zi_reads*x
  if (isTRUE(spec$zi) && isTRUE(spec$zi_DO)) g <- g+p$zi_DO*DO_c
  if (isTRUE(spec$zi) && isTRUE(spec$zi_TS)) {
    if (is.null(TS) && in_sample) TS <- as.character(fit$data$TS)
    if (!is.null(TS)) {
      if (length(TS)==1L) TS <- rep(TS,length(x))
      index <- match(as.character(TS),levels(fit$data$TS))
      if (length(index)!=length(x) || anyNA(index)) stop("Unknown or unmatched prediction TS.")
      g <- g+p$b_zi_ts[index]
    } # New curve rows with TS=NULL explicitly use a zero TS offset.
  }
  list(sd=exp(ls),pi=plogis(g),alpha=alpha,m=m,fixed_eta=fixed_eta)
}

joint_integrals <- function(fit, rows=NULL, nodes=fit$quadrature_n) {
  if(is.null(rows) && isTRUE(fit$joint_spec$sd_mean)) {
    dd <- fit$data;v <- joint_row_values(fit);nn <- dd$Ruegeria_reads_total;yy <- dd$MC_reads
    rows <- unique(c(head(order(nn),3),head(order(-nn),3),head(order(yy/nn),4),
      head(order(-yy/nn),5),head(which(yy==0),4),head(order(v$sd),2),head(order(-v$sd),2),
      head(order(v$alpha),2),head(order(-v$alpha),2)))
  }
  if(!isTRUE(fit$joint_spec$alpha_mean) && !isTRUE(fit$joint_spec$tail_mix))
    return(suite_legacy_integrals(fit,rows,nodes))
  d <- fit$data;y <- d$MC_reads;N <- d$Ruegeria_reads_total
  rv <- joint_row_values(fit);parts <- suite_tail_components(fit$joint_spec)
  if(is.null(rows))rows <- unique(c(head(order(N),3),head(order(-N),3),
    head(order(y/N),4),head(order(-y/N),5),head(which(y==0),4),
    head(order(rv$alpha),2),head(order(-rv$alpha),2)))
  gh <- normal_quadrature(nodes)
  logsum <- function(v){a <- max(v);a+log(sum(exp(v-a)))}
  do.call(rbind,lapply(rows,function(i) {
    alpha <- rv$alpha[i];delta <- alpha/sqrt(1+alpha^2)
    counts <- lapply(seq_along(parts$weight),function(k) {
      omega <- rv$sd[i]*parts$scale[k]/sqrt(1-2*delta^2/pi)
      eta <- fit$eta_base[i]-omega*delta*sqrt(2/pi)
      z <- 10*tanh(((qlogis((y[i]+.5)/(N[i]+1))-eta)/omega)/10)
      for(j in 1:30) {
        pr <- plogis(eta+omega*z);az <- alpha*z
        mills <- exp(dnorm(az,log=TRUE)-pnorm(az,log.p=TRUE))
        step <- (omega*(y[i]-N[i]*pr)-z+alpha*mills)/
          (1+omega^2*N[i]*pr*(1-pr)+alpha^2*mills*(az+mills))
        z <- z+3*tanh(step/3)
      }
      pr <- plogis(eta+omega*z);az <- alpha*z
      mills <- exp(dnorm(az,log=TRUE)-pnorm(az,log.p=TRUE))
      zs <- 1/sqrt(1+omega^2*N[i]*pr*(1-pr)+alpha^2*mills*(az+mills))
      kernel <- function(a)dbinom(y[i],N[i],plogis(eta+omega*a),log=TRUE)+
        dnorm(a,log=TRUE)+log(2)+pnorm(alpha*a,log.p=TRUE)
      zz <- z+zs*gh$nodes
      quadrature <- logsum(log(gh$weights)+kernel(zz)-dnorm(zz,z,zs,log=TRUE))
      anchor <- kernel(z)
      direct <- integrate(function(a)exp(kernel(z+zs*a)-anchor)*zs,
        -Inf,Inf,rel.tol=1e-9,subdivisions=400L)
      c(quadrature=quadrature,direct=anchor+log(direct$value))
    })
    count_gh <- logsum(log(parts$weight)+vapply(counts,`[[`,numeric(1),"quadrature"))
    count_direct <- logsum(log(parts$weight)+vapply(counts,`[[`,numeric(1),"direct"))
    mix <- function(v) {
      if(!isTRUE(fit$zi))return(v)
      if(y[i]>0)return(log1p(-rv$pi[i])+v)
      logsum(c(log(rv$pi[i]),log1p(-rv$pi[i])+v))
    }
    a <- mix(count_gh);b <- mix(count_direct)
    data.frame(New_ID=d$New_ID[i],N=N[i],y=y[i],Observation_SD=rv$sd[i],
      alpha=alpha,m_i=rv$m[i],pi=rv$pi[i],log_likelihood_quadrature=a,
      log_likelihood_direct=b,absolute_difference=abs(a-b),
      count_absolute_difference=abs(count_gh-count_direct),
      passed=is.finite(b)&&abs(a-b)<1e-4&&abs(count_gh-count_direct)<1e-4)
  }))
}

simulate_centered_olre <- function(fit,n_sim,seed) {
  if (!isTRUE(fit$joint_spec$tail_mix) && !isTRUE(fit$joint_spec$alpha_mean))
    return(suite_legacy_simulate(fit,n_sim,seed))
  d <- fit$data; N <- d$Ruegeria_reads_total; rv <- joint_row_values(fit)
  parts <- suite_tail_components(fit$joint_spec)
  alpha <- rv$alpha; delta <- alpha/sqrt(1+alpha^2)
  set.seed(seed)
  sims <- replicate(n_sim,{
    component <- if(length(parts$weight)==1L)rep(1L,nrow(d)) else 1L+rbinom(nrow(d),1,parts$weight[2])
    omega <- rv$sd*parts$scale[component]/sqrt(1-2*delta^2/pi)
    u <- omega*(delta*(abs(rnorm(nrow(d)))-sqrt(2/pi))+sqrt(1-delta^2)*rnorm(nrow(d)))
    counts <- rbinom(nrow(d),N,plogis(fit$eta_base+u))
    if (isTRUE(fit$zi)) counts[rbinom(nrow(d),1,rv$pi)==1L] <- 0L
    counts
  })
  expected <- rowMeans(sims)
  res <- DHARMa::createDHARMa(simulatedResponse=sims,observedResponse=d$MC_reads,
    fittedPredictedResponse=expected,integerResponse=TRUE,seed=seed)
  res$fittedModel <- fit
  list(res_check=res,total_reads=N,expected_prop=expected/N,observation_sd=rv$sd,
    simulation_method="Conditional fitted count/zero TS and Repeat; redraw centered skew-normal with row-specific SD/alpha (and configured scale mixture), binomial count and extra-zero indicator")
}

joint_interior <- function(result,d) {
  if (!suite_legacy_interior(result,d)) return(FALSE)
  p <- result$parameters; sp <- result$spec
  ls <- p$log_sd_observation+if(isTRUE(sp$sd_reads))p$sd_reads*d$log_reads_c else 0
  if(isTRUE(sp$sd_quadratic))ls <- ls+p$sd_reads2*d$log_reads_c^2
  alpha <- rep(p$alpha,nrow(d))
  if(isTRUE(sp$sd_mean)||isTRUE(sp$alpha_mean)) {
    X <- result$abundance_X
    if(is.null(X))return(FALSE)
    m <- as.vector(sweep(X,2L,colMeans(X),"-")%*%p$beta)
    if(isTRUE(sp$sd_mean))ls <- ls+p$sd_mean*m
    if(isTRUE(sp$alpha_mean))alpha <- alpha+p$alpha_mean*m
  }
  all(is.finite(ls)) && min(ls)>-12 && max(ls)<6 &&
    all(is.finite(alpha)) && max(abs(alpha))<19.99 &&
    (!isTRUE(sp$sd_mean)||abs(p$sd_mean)<2.99) &&
    (!isTRUE(sp$alpha_mean)||abs(p$alpha_mean)<4.99) &&
    (!isTRUE(sp$sd_quadratic)||abs(p$sd_reads2)<3.99) &&
    (!isTRUE(sp$zi_DO)||abs(p$zi_DO)<9.99)
}
joint_starts <- function(seed,d,spec) {
  starts <- suite_legacy_starts(seed,d,spec)
  starts[[1]] <- joint_parameters(seed,d,spec) # Retain the verified warm optimum.
  for(j in seq_along(starts)) {
    if(isTRUE(spec$sd_quadratic))starts[[j]]$sd_reads2 <- c(
      if(is.null(seed$sd_reads2))0 else seed$sd_reads2,0,-.15,.15,-.5,.5)[j]
    if(isTRUE(spec$zi_DO))starts[[j]]$zi_DO <- c(
      if(is.null(seed$zi_DO))0 else seed$zi_DO,.1,-.5,.5,-1,1)[j]
    if(isTRUE(spec$sd_mean))starts[[j]]$sd_mean <- (seed$sd_mean %||% 0)+c(0,-.1,.1,-.3,.3,.5)[j]
    if(isTRUE(spec$alpha_mean))starts[[j]]$alpha_mean <- (seed$alpha_mean %||% 0)+c(0,-.1,.1,-.5,.5,1)[j]
    if(j%in%c(2L,3L) && (isTRUE(spec$sd_mean)||isTRUE(spec$alpha_mean))) {
      # Two nearby starts retain all verified parent parameters; three further
      # starts still explore different skew/read-width configurations.
      starts[[j]] <- joint_parameters(seed,d,spec)
      if(isTRUE(spec$sd_mean))starts[[j]]$sd_mean <- (seed$sd_mean %||% 0)+c(-.1,.1)[j-1L]
      if(isTRUE(spec$alpha_mean))starts[[j]]$alpha_mean <- (seed$alpha_mean %||% 0)+c(-.1,.1)[j-1L]
    }
    if(isTRUE(spec$sd_mean))starts[[j]]$sd_mean <- max(-2.5,min(2.5,starts[[j]]$sd_mean))
    if(isTRUE(spec$alpha_mean))starts[[j]]$alpha_mean <- max(-4,min(4,starts[[j]]$alpha_mean))
    if(j%in%c(4L,5L) && (isTRUE(spec$sd_mean)||isTRUE(spec$alpha_mean)))
      for(group in c("repeat","ts"))if(!isTRUE(spec[[paste0("zero_",group)]]) &&
        starts[[j]][[paste0("log_sd_",group)]]<log(.01))
        starts[[j]][[paste0("log_sd_",group)]] <- log(c(.1,.3)[j-3L])
  }
  starts
}

suite_specs <- function() {
  base <- list(zi=TRUE,zi_reads=TRUE,sd_reads=TRUE,sd_quadratic=FALSE,
    zi_DO=FALSE,zi_TS=FALSE,tail_mix=FALSE,sd_mean=FALSE,alpha_mean=FALSE,lambda_zi_reads=suite_lambda_reads,
    zero_DO_prior_sd=suite_zero_DO_prior_sd,zero_TS_prior_sd=suite_zero_TS_prior_sd,
    tail_weight=suite_tail_weight,tail_ratio=suite_tail_SD_ratio)
  changes <- list(D0=list(),Z_DO=list(zi_DO=TRUE),Z_TS=list(zi_TS=TRUE),
    SD_Q=list(sd_quadratic=TRUE),Z_DO_SD_Q=list(zi_DO=TRUE,sd_quadratic=TRUE),
    Z_TS_SD_Q=list(zi_TS=TRUE,sd_quadratic=TRUE),TAIL=list(tail_mix=TRUE),
    SD_M=list(zi_TS=TRUE,sd_quadratic=TRUE,sd_mean=TRUE),
    ALPHA_M=list(zi_TS=TRUE,sd_quadratic=TRUE,alpha_mean=TRUE),
    SD_ALPHA_M=list(zi_TS=TRUE,sd_quadratic=TRUE,sd_mean=TRUE,alpha_mean=TRUE))
  if(!isTRUE(suite_include_previous_candidates))
    changes <- changes[c("D0","Z_TS_SD_Q","SD_M","ALPHA_M","SD_ALPHA_M")]
  lapply_names <- lapply(names(changes),function(id)utils::modifyList(base,c(list(id=id),changes[[id]])))
  names(lapply_names) <- names(changes); lapply_names
}
suite_effective_spec <- function(spec,d) {
  if (!any(d$MC_reads==0)) {spec$zi <- spec$zi_reads <- spec$zi_DO <- spec$zi_TS <- FALSE}
  spec
}
suite_statistical_spec <- function(spec) {
  spec$id <- NULL; spec$zero_repeat <- spec$zero_ts <- NULL
  if(!isTRUE(spec$zi)) {spec$lambda_zi_reads <- 0;spec$zero_DO_prior_sd <- spec$zero_TS_prior_sd <- NULL}
  if(!isTRUE(spec$zi_DO))spec$zero_DO_prior_sd <- NULL
  if(!isTRUE(spec$zi_TS))spec$zero_TS_prior_sd <- NULL
  if(!isTRUE(spec$tail_mix))spec$tail_weight <- spec$tail_ratio <- NULL
  # Canonical false flags permit a checked v1 bridge for unchanged models.
  if(!isTRUE(spec$sd_mean))spec$sd_mean <- NULL
  if(!isTRUE(spec$alpha_mean))spec$alpha_mean <- NULL
  spec
}
suite_same_data <- function(fit,d,X) {
  cols <- c("New_ID","Sample_ID","MC_ID","MC_reads","Other_reads","Ruegeria_reads_total",
    "DO_c","Host_type","Repeat_group_ID","TS","log_reads_c",paste0("mean_reads_s",1:3))
  !is.null(fit) && nrow(fit$data)==nrow(d) && all(cols%in%names(fit$data)) &&
    all(vapply(cols,function(nm)if(is.numeric(d[[nm]]))
      isTRUE(all.equal(as.numeric(d[[nm]]),as.numeric(fit$data[[nm]]),tolerance=1e-12)) else
      identical(as.character(d[[nm]]),as.character(fit$data[[nm]])),logical(1))) &&
    identical(colnames(fit$X),colnames(X)) && isTRUE(all.equal(unname(fit$X),unname(X),tolerance=1e-12)) &&
    identical(levels(fit$data$TS),levels(d$TS)) && identical(levels(fit$data$Repeat_group_ID),levels(d$Repeat_group_ID))
}
suite_atomic_save <- function(x,path) {
  tmp <- paste0(path,".",Sys.getpid(),".tmp");saveRDS(x,tmp)
  if(!file.rename(tmp,path)) {unlink(tmp);stop("Cannot commit checkpoint: ",path)}
}
suite_key <- function(context,spec) digest::digest(list(version=suite_version,
  data=context$d,X=context$X,spec=suite_statistical_spec(spec),engine=.GLMM_suite_settings$engine),algo="sha256")

suite_cache_candidate <- function(context,spec,key,output_dir) {
  checks <- list(); old <- NULL; source <- ""
  dirs <- .GLMM_suite_settings$cache_dirs
  for(dir in dirs) {
    path <- file.path(dir,context$MC_ID,spec$id,"fitted_candidate.rds")
    if(!file.exists(path))next
    cc <- tryCatch(readRDS(path),error=function(e)NULL)
    ok <- !is.null(cc) && identical(cc$suite_cache_key,key) && isTRUE(cc$primary_valid) && suite_same_data(cc$primary,context$d,context$X)
    checks[[length(checks)+1L]] <- data.frame(path=path,eligible=ok,reason=if(ok)"exact suite key/data" else "suite key/data/validity mismatch")
    if(ok){old <- cc;source <- path;break}
  }
  # Verified v1 controls can cross the engine/version change. New m-linked
  # models never reuse a parent result; those files provide initial values only.
  if(is.null(old) && !isTRUE(spec$sd_mean) && !isTRUE(spec$alpha_mean)) {
    bridge_dirs <- unique(c(suite_structure_reference,dirs))
    for(dir in bridge_dirs) {
      aliases <- unique(c(spec$id,if(!any(context$d$MC_reads==0) &&
        isTRUE(spec$sd_quadratic))"SD_Q"))
      for(alias in aliases) {
        path <- file.path(dir,context$MC_ID,alias,"fitted_candidate.rds")
        if(!file.exists(path))next
        cc <- tryCatch(readRDS(path),error=function(e)NULL)
        ss <- tryCatch(readRDS(file.path(dir,"01_run_settings.rds")),error=function(e)NULL)
        f <- cc$primary
        requested_old <- cc$requested_spec %||% cc$spec
        snapshot <- file.path(dir,"01_analysis_script_used.R")
        session_file <- file.path(dir,"sessionInfo.txt")
        tokens <- if(file.exists(session_file))unlist(strsplit(paste(
          readLines(session_file,warn=FALSE),collapse=" "),"[[:space:]]+")) else character()
        versions_ok <- all(vapply(c("RTMB","TMB","Matrix"),function(pkg)
          paste0(pkg,"_",as.character(utils::packageVersion(pkg)))%in%tokens,logical(1)))
        ok <- !is.null(f) && !is.null(ss) && isTRUE(cc$primary_valid) &&
          isTRUE(cc$reproducible) && isTRUE(cc$quadrature_ok) && isTRUE(cc$direct_integrals_ok) &&
          identical(ss$experiment,"MC_common_structure_v1") && versions_ok &&
          identical(ss$input_md5,.GLMM_suite_settings$input_md5) &&
          file.exists(snapshot) && identical(unname(tools::md5sum(snapshot)),ss$script_md5) &&
          identical(digest::digest(file=snapshot,algo="sha256"),
            "d5dd3cf39f585a78f966a1cb779ce262ebd8ae6c0113d3e351f668932c24ed96") &&
          identical(unname(tools::md5sum(file.path(dir,"01_input_used_MC8.xlsx"))),ss$input_md5) &&
          suite_same_data(f,context$d,context$X) &&
          identical(suite_statistical_spec(suite_effective_spec(requested_old,context$d)),
            suite_statistical_spec(spec))
        reason <- "v1 source/data/spec/package/validity mismatch"
        if(ok) {
          actual <- f$joint_spec # Preserve checked pi=0 and group-variance boundaries.
          actual$id <- spec$id;actual$sd_mean <- actual$alpha_mean <- FALSE
          pp <- joint_parameters(f$parameters,context$d,actual)
          obj <- RTMB::MakeADFun(joint_objective(context$d,context$X,f$quadrature_n,actual),
            pp,random=joint_random(actual),map=joint_map(pp,actual),silent=TRUE)
          ok <- identical(names(obj$par),names(f$fit$par))
          if(ok) {
            val <- obj$fn(f$fit$par);gr <- obj$gr(f$fit$par)
            ok <- is.finite(val)&&abs(val-joint_total_objective(f))<.01&&
              all(is.finite(gr))&&max(abs(gr))<.02
            write_result_table(data.frame(old_objective=joint_total_objective(f),
              current_objective=val,difference=val-joint_total_objective(f),
              max_gradient=max(abs(gr)),passed=ok),output_dir,"01_v1_objective_identity.csv")
          }
          rm(obj);gc(FALSE)
          if(ok) {cc$primary$joint_spec <- actual;cc$spec <- actual;old <- cc;source <- path}
          reason <- if(ok)"verified v1 control; objective/gradient checked; no refit" else "v1 objective/gradient mismatch"
        }
        checks[[length(checks)+1L]] <- data.frame(path=path,eligible=ok,reason=reason)
        if(!is.null(old))break
      }
      if(!is.null(old))break
    }
  }
  # Only the exact legacy D reference has a deliberate cross-version bridge.
  if(is.null(old) && identical(spec$id,"D0") && suite_lambda_reads==.125) {
    dir <- suite_reference_cache; path <- file.path(dir,context$MC_ID,"D_s2_lambda0p125","fitted_candidate.rds")
    if(file.exists(path)) {
      cc <- tryCatch(readRDS(path),error=function(e)NULL)
      settings <- tryCatch(readRDS(file.path(dir,"01_run_settings.rds")),error=function(e)NULL)
      f <- cc$primary
      ok <- !is.null(f) && isTRUE(cc$primary_valid) && isTRUE(cc$reproducible) &&
        isTRUE(cc$quadrature_ok) && isTRUE(cc$direct_integrals_ok) &&
        suite_same_data(f,context$d,context$X) && !is.null(settings) &&
        identical(settings$input_md5,.GLMM_suite_settings$input_md5) &&
        identical(settings$experiment,"D_zero_slope_ridge_v1") &&
        identical(unname(tools::md5sum(file.path(dir,"01_analysis_script_used.R"))),settings$script_md5) &&
        identical(unname(tools::md5sum(file.path(dir,"01_input_used_MC8.xlsx"))),settings$input_md5) &&
        isTRUE(f$joint_spec$sd_reads) &&
        (isTRUE(f$zi)&&isTRUE(all.equal(f$joint_spec$lambda_zi_reads,suite_lambda_reads)) ||
          !any(context$d$MC_reads==0)&&!isTRUE(f$zi))
      reason <- "legacy provenance/data/spec mismatch"
      if(ok) {
        # Evaluate once, without refitting: the old optimum must have the same
        # objective and a small gradient under today's likelihood implementation.
        actual <- utils::modifyList(f$joint_spec,spec)
        pp <- joint_parameters(f$parameters,context$d,actual)
        obj <- RTMB::MakeADFun(joint_objective(context$d,context$X,f$quadrature_n,actual),pp,
          random=joint_random(actual),map=joint_map(pp,actual),silent=TRUE)
        ok <- identical(names(obj$par),names(f$fit$par))
        if(ok) {
          val <- obj$fn(f$fit$par); gr <- obj$gr(f$fit$par)
          ok <- is.finite(val)&&abs(val-joint_total_objective(f))<.01&&all(is.finite(gr))&&max(abs(gr))<.02
          write_result_table(data.frame(old_objective=joint_total_objective(f),current_objective=val,
            difference=val-joint_total_objective(f),max_gradient=max(abs(gr)),passed=ok),output_dir,"01_legacy_objective_identity.csv")
        }
        reason <- if(ok)"legacy exact D; inherited six-start/integration checks; objective/gradient rechecked without refit" else "legacy objective/gradient mismatch"
        if(ok){cc$primary$joint_spec <- actual;cc$spec <- actual;old <- cc;source <- path}
        rm(obj);gc(FALSE)
      }
      checks[[length(checks)+1L]] <- data.frame(path=path,eligible=ok,reason=reason)
    }
  }
  write_result_table(bind_flexible_tables(checks),output_dir,"01_cache_audit.csv")
  if(!is.null(old)) {old$suite_cache_key <- key;old$reuse_source <- source}
  old
}

suite_clone_diagnostics <- function(diag,fit,source_dir,output_dir) {
  dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
  files <- list.files(source_dir,full.names=TRUE);files <- files[!dir.exists(files)]
  same_dir <- identical(normalizePath(source_dir,winslash="/"),normalizePath(output_dir,winslash="/"))
  if(!same_dir&&length(files)&&!all(file.copy(files,output_dir,overwrite=TRUE)))stop("Cannot copy compatible diagnostics.")
  diag$fit_model <- fit;diag$res_check$fittedModel <- fit
  writeLines(c("Reused the identical fitted distribution, data rows, simulation count and seed.",
    "Existing plot titles retain their original candidate identifier.",source_dir),file.path(output_dir,"00_diagnostics_reused_from.txt"))
  diag
}

suite_diagnostics <- function(candidate,context,output_dir) {
  f <- candidate$primary;path <- file.path(output_dir,"diagnostics");dir.create(path,showWarnings=FALSE)
  old <- NULL
  if(nzchar(candidate$reuse_source %||% "")) {
    source_dir <- file.path(dirname(candidate$reuse_source),"diagnostics")
    file <- file.path(source_dir,"model_analysis_objects.rds")
    if(file.exists(file)) {
      a <- tryCatch(readRDS(file),error=function(e)NULL)
      ss <- tryCatch(readRDS(file.path(dirname(dirname(dirname(candidate$reuse_source))),"01_run_settings.rds")),error=function(e)NULL)
      seed_ok <- !is.null(ss) && identical(as.integer(ss$seed),as.integer(simulation_seed)) &&
        identical(as.character(ss$MC_ids),as.character(.GLMM_suite_settings$MC_ids)) &&
        !is.null(a) && identical(as.character(a$res_check$randomState$kind),as.character(RNGkind()))
      versions_ok <- FALSE
      session_file <- file.path(dirname(dirname(dirname(candidate$reuse_source))),"sessionInfo.txt")
      if(file.exists(session_file)) {
        tokens <- unlist(strsplit(paste(readLines(session_file,warn=FALSE),collapse=" "),"[[:space:]]+"))
        versions_ok <- all(vapply(c("DHARMa","qgam"),function(p)
          paste0(p,"_",as.character(utils::packageVersion(p)))%in%tokens,logical(1)))
      }
      ok <- !is.null(a) && isTRUE(a$diagnostics_complete) && seed_ok &&
        versions_ok &&
        ncol(a$res_check$simulatedResponse)==n_sim &&
        identical(as.numeric(a$res_check$observedResponse),as.numeric(context$d$MC_reads)) &&
        suite_same_data(a$fit_model,context$d,context$X) &&
        isTRUE(all.equal(a$fit_model$parameters,f$parameters,tolerance=1e-12))
      if(ok)old <- suite_clone_diagnostics(a,f,source_dir,path)
    }
  }
  if(is.null(old)) {
    old <- diagnose_skew_comparison_model(paste(context$MC_ID,f$joint_spec$id),f,context$d,n_sim,
      simulation_seed+context$index*100L,path)
    old <- zero_count_diagnostics(old,context$d,path)
  }
  suite_atomic_save(old,file.path(path,"model_analysis_objects.rds"));old
}
`%||%` <- function(x,y) if(is.null(x))y else x

suite_calibration <- function(diag,d) {
  sim <- diag$res_check$simulatedResponse; N <- d$Ruegeria_reads_total
  z <- diag$zero_group_table; ts <- z[z$grouping=="TS",,drop=FALSE]
  zero_rmse <- sqrt(mean(((ts$observed-ts$simulated_mean)/ts$n)^2))
  bands <- cut(rank(N,ties.method="average"),
    unique(stats::quantile(rank(N,ties.method="average"),seq(0,1,.25))),include.lowest=TRUE)
  dev <- 2*abs(diag$res_check$scaledResiduals-.5)
  reads_error <- sqrt(mean((tapply(dev,bands,mean)-.5)^2))
  residual <- as.numeric(diag$res_check$scaledResiduals)
  ordered <- sort(residual);n <- length(ordered)
  ks_distance <- max(seq_len(n)/n-ordered,ordered-(seq_len(n)-1)/n)
  ranks <- rank(diag$fitted_rows$expected_MC_proportion,ties.method="average")
  breaks <- unique(quantile(ranks,seq(0,1,.25)))
  predicted_band <- if(length(breaks)>1L)cut(ranks,breaks,include.lowest=TRUE) else factor(rep("all",n))
  predicted_errors <- unlist(lapply(split(residual,predicted_band),function(r)
    as.numeric(quantile(r,c(.25,.5,.75)))-c(.25,.5,.75)))
  dispersion_error <- if(is.finite(diag$dispersion_ratio)&&diag$dispersion_ratio>0)abs(log(diag$dispersion_ratio)) else NA_real_
  stat <- function(y)c(Q95=unname(quantile(y/N,.95)),Q99=unname(quantile(y/N,.99)),
    Maximum=max(y/N),Above25=mean(y/N>.25),Above50=mean(y/N>.5))
  observed <- stat(d$MC_reads); simulated <- apply(sim,2,stat)
  lo <- apply(simulated,1,quantile,.025);hi <- apply(simulated,1,quantile,.975)
  excess <- pmax(observed-hi,lo-observed,0)/pmax(hi-lo,.005)
  tail <- data.frame(statistic=names(observed),observed=observed,
    simulated_mean=rowMeans(simulated),simulated_95_low=lo,simulated_95_high=hi,
    standardized_outside_interval=excess,row.names=NULL)
  # Positive-count conditional ranks reuse existing simulations. Removing only
  # the observed zeros from unconditional DHARMa ranks would change the target.
  positive <- d$MC_reads>0;positive_n <- rowSums(sim>0)
  positive_pred <- rowSums(sim)/pmax(positive_n,1)/N
  less <- rowSums(sweep(sim,1,d$MC_reads,"<") & sim>0)
  equal <- rowSums(sweep(sim,1,d$MC_reads,"=="))
  midpoint_rank <- (less+.5*(equal+1))/(positive_n+1)
  positive_ranks <- rank(positive_pred[positive],ties.method="average")
  breaks <- unique(quantile(positive_ranks,seq(0,1,.25)))
  band <- if(length(breaks)>1L)cut(positive_ranks,breaks,include.lowest=TRUE) else
    factor(rep("All",sum(positive)))
  indices <- which(positive)
  positive_bands <- do.call(rbind,lapply(split(indices,band),function(j) {
    q <- quantile(midpoint_rank[j],c(.25,.5,.75))
    data.frame(n=length(j),observed_mean_proportion=mean(d$MC_reads[j]/N[j]),
      predicted_positive_mean=mean(positive_pred[j]),residual_Q25=q[1],
      residual_Q50=q[2],residual_Q75=q[3],upper_5pct_fraction=mean(midpoint_rank[j]>.95),
      lower_5pct_fraction=mean(midpoint_rank[j]<.05),min_positive_sim=min(positive_n[j]),row.names=NULL)
  }))
  positive_bands$prediction_band <- seq_len(nrow(positive_bands))
  pq <- as.matrix(positive_bands[,c("residual_Q25","residual_Q50","residual_Q75")])
  p_error <- sweep(pq,2,c(.25,.5,.75),"-")
  list(zero_TS_RMSE=zero_rmse,reads_deviation_RMSE=reads_error,
    dispersion_log_error=dispersion_error,residual_KS_D=ks_distance,
    predicted_quantile_RMSE=sqrt(mean(predicted_errors^2)),
    tail_interval_excess=sqrt(mean(excess^2)),tail=tail,
    positive_quantile_RMSE=sqrt(mean(p_error^2)),
    positive_center_RMSE=sqrt(mean((positive_bands$residual_Q50-.5)^2)),
    positive_width_RMSE=sqrt(mean((positive_bands$residual_Q75-positive_bands$residual_Q25-.5)^2)),
    positive_upper_tail_RMSE=sqrt(mean((positive_bands$upper_5pct_fraction-.05)^2)),
    positive_min_sim=min(positive_n[positive]),positive_bands=positive_bands)
}

suite_metric_row <- function(context,spec) {
  r <- zero_penalty_summary_row(context$MC_ID,list(id=spec$id,lambda_zi_reads=spec$lambda_zi_reads,
    prior_sd_zi_reads=sqrt(1/(2*spec$lambda_zi_reads)),strength_label="fixed_reference",is_reference=spec$id=="D0"),
    context$n_source,context$d)
  r$base_model <- "common_structure_suite";r$fit_source <- "pending";r$reuse_source <- ""
  r$zi_DO_slope <- r$sd_reads_quadratic <- NA_real_
  r$zero_TS_fixed_prior_SD <- if(isTRUE(spec$zi_TS))spec$zero_TS_prior_sd else NA_real_
  r$tail_mixture <- isTRUE(spec$tail_mix)
  r$zero_TS_RMSE <- r$reads_deviation_RMSE <- r$tail_interval_excess <- NA_real_
  r$dispersion_log_error <- r$residual_KS_D <- r$predicted_quantile_RMSE <- NA_real_
  r$diagnostic_flag_count <- NA_integer_;r$sensitivity_only <- FALSE;r$excluded_New_ID <- ""
  r$diagnostic_n_sim <- n_sim;r$CV_status <- "not_run; structural_screening_only"
  for(nm in c("sd_mean_slope","alpha_mean_slope","alpha_row_min","alpha_row_max","m_min","m_max",
    "positive_quantile_RMSE","positive_center_RMSE","positive_width_RMSE",
    "positive_upper_tail_RMSE","positive_min_sim"))r[[nm]] <- NA_real_
  r$lambda_selection <- "reference fixed; previous sensitivity retained separately"
  r$requested_penalty_target <- "zero_reads; zero_DO if present; Gaussian zero_TS prior integrated"
  r
}

suite_run_one <- function(job) {
  context <- .GLMM_suite_contexts[[job$MC_ID]]
  specs <- .GLMM_suite_specs
  requested <- specs[[job$models[1]]]
  actual <- suite_effective_spec(requested,context$d)
  if(isTRUE(job$sensitivity)) {
    drop <- which.max(context$d$MC_reads/context$d$Ruegeria_reads_total)
    excluded <- as.character(context$d$New_ID[drop])
    write_result_table(context$d[drop,,drop=FALSE],job$output_dir,"01_sensitivity_excluded_row.csv")
    context$d <- context$d[-drop,,drop=FALSE];context$X <- context$X[-drop,,drop=FALSE]
    context$stub <- make_model_stub(context$d)
    requested <- specs$D0;requested$id <- "SENS_TOP1"
    actual <- suite_effective_spec(requested,context$d)
    if(qr(context$X)$rank<ncol(context$X)) {
      r <- suite_metric_row(context,requested);r$sensitivity_only <- TRUE;r$excluded_New_ID <- excluded
      r$status <- "敏感性未拟合：删除最高比例行后固定设计秩不足";return(list(metrics=r,candidate_files=list()))
    }
  }
  id <- requested$id; out <- job$output_dir
  key <- suite_key(context,actual)
  candidate <- tryCatch(suite_cache_candidate(context,actual,key,out),error=function(e) {
    writeLines(enc2utf8(conditionMessage(e)),file.path(out,"01_cache_rejected_error.txt"));NULL
  })
  cached <- !is.null(candidate)
  seed_source <- "fresh_start"
  fit_nodes <- 31L
  if(is.null(candidate)) {
    base_file <- file.path(context$mc_dir,"D0","fitted_candidate.rds")
    if(file.exists(base_file)) {
      baseline <- readRDS(base_file);seed <- baseline$primary$parameters
      seed_source <- if(isTRUE(baseline$primary_valid))"verified_D_start" else "previous_unverified_start"
      if(isTRUE(baseline$primary_valid))fit_nodes <- max(fit_nodes,baseline$primary$quadrature_n)
    } else seed <- context$initial
    if(actual$id %in% c("Z_DO_SD_Q","Z_TS_SD_Q")) {
      parent <- file.path(context$mc_dir,"SD_Q","fitted_candidate.rds")
      if(file.exists(parent)) {
        cc <- readRDS(parent)
        if(isTRUE(cc$primary_valid)) {
          seed <- cc$primary$parameters;seed_source <- "verified_SD_Q_start"
          fit_nodes <- max(fit_nodes,cc$primary$quadrature_n)
        }
      }
      zero_parent <- file.path(context$mc_dir,if(isTRUE(actual$zi_DO))"Z_DO" else "Z_TS","fitted_candidate.rds")
      if(file.exists(zero_parent)) {
        cc <- readRDS(zero_parent)
        if(isTRUE(cc$primary_valid))for(nm in c("zi_DO","b_zi_ts"))
          if(!is.null(cc$primary$parameters[[nm]]))seed[[nm]] <- cc$primary$parameters[[nm]]
      }
    }
    if(isTRUE(actual$sd_mean)||isTRUE(actual$alpha_mean)) {
      # Always prefer a checked full-data parent with the same zero/width base.
      parent <- file.path(context$mc_dir,"Z_TS_SD_Q","fitted_candidate.rds")
      if(!file.exists(parent) && !any(context$d$MC_reads==0))
        parent <- file.path(context$mc_dir,"SD_Q","fitted_candidate.rds")
      if(file.exists(parent)) {
        cc <- readRDS(parent)
        if(isTRUE(cc$primary_valid)) {
          seed <- cc$primary$parameters;seed_source <- "verified_Z_TS_SD_Q_start"
          fit_nodes <- max(fit_nodes,cc$primary$quadrature_n)
        }
      }
      if(identical(actual$id,"SD_ALPHA_M")) {
        for(parent_id in c("SD_M","ALPHA_M")) {
          parent <- file.path(context$mc_dir,parent_id,"fitted_candidate.rds")
          if(file.exists(parent)) {
            cc <- readRDS(parent)
            if(isTRUE(cc$primary_valid)) {
              if(parent_id=="SD_M")seed <- cc$primary$parameters
              else seed$alpha_mean <- cc$primary$parameters$alpha_mean
              fit_nodes <- max(fit_nodes,cc$primary$quadrature_n)
              seed_source <- paste(seed_source,parent_id,sep="+")
            }
          }
        }
      }
      # The added slopes start at zero when absent, nesting the verified parent.
      seed$sd_mean <- seed$sd_mean %||% 0
      seed$alpha_mean <- seed$alpha_mean %||% 0
    }
    if(is.null(seed))stop("No usable starting parameters.")
    for(nm in c("log_sd_repeat","log_sd_ts"))if(!is.finite(seed[[nm]]))seed[[nm]] <- -6
    candidate <- joint_fit_pass(seed,context$d,context$X,context$stub,actual,out,pool=NULL,nodes=fit_nodes)
    candidate <- joint_variance_boundaries(candidate,context$d,context$X,context$stub,out,pool=NULL)
    if(!is.null(candidate$primary) && isTRUE(candidate$primary$zi) &&
        max(joint_row_values(candidate$primary)$pi)<1e-6) {
      nozi <- actual;nozi$zi <- nozi$zi_reads <- nozi$zi_DO <- nozi$zi_TS <- FALSE
      boundary <- joint_fit_pass(candidate$primary$parameters,context$d,context$X,context$stub,
        nozi,file.path(out,"pi0_boundary"),pool=NULL,nodes=candidate$primary$quadrature_n)
      boundary <- joint_variance_boundaries(boundary,context$d,context$X,context$stub,file.path(out,"pi0_boundary"),pool=NULL)
      candidate <- joint_pi0(candidate,boundary,context$d,context$X,out)
    }
    candidate$reuse_source <- ""
  }
  candidate$requested_spec <- requested;candidate$suite_cache_key <- key
  candidate_file <- file.path(out,"fitted_candidate.rds")
  suite_atomic_save(candidate,candidate_file)
  r <- suite_metric_row(context,requested)
  r <- zero_penalty_finish_row(r,candidate)
  r$fit_source <- if(cached)"verified_cache" else paste0("new_fit_from_",seed_source)
  r$reuse_source <- candidate$reuse_source %||% ""
  if(!is.null(candidate$primary)) {
    p <- candidate$primary$parameters
    r$zi_DO_slope <- p$zi_DO %||% 0;r$sd_reads_quadratic <- p$sd_reads2 %||% 0
    rv <- joint_row_values(candidate$primary)
    r$sd_mean_slope <- p$sd_mean %||% 0;r$alpha_mean_slope <- p$alpha_mean %||% 0
    r$alpha_row_min <- min(rv$alpha);r$alpha_row_max <- max(rv$alpha)
    r$m_min <- min(rv$m);r$m_max <- max(rv$m)
    writeLines(enc2utf8(capture.output(base::summary(candidate$primary))),file.path(out,"02_model_summary.txt"))
    write_result_table(data.frame(TS=levels(context$d$TS),count_TS_offset=p$b_ts,
      extra_zero_TS_offset=p$b_zi_ts %||% rep(0,nlevels(context$d$TS))),out,"04_TS_offsets.csv")
  }
  if(isTRUE(job$sensitivity)) {
    r$sensitivity_only <- TRUE;r$excluded_New_ID <- excluded;r$diagnostic_n_sim <- 0L
    r$diagnostic_status <- "not_run; point-estimate influence refit only"
    r$status <- if(isTRUE(r$numerical_valid))"敏感性拟合数值有效；仅比较DO点估计" else r$status
  } else if(isTRUE(candidate$primary_valid)) tryCatch({
    diag <- suite_diagnostics(candidate,context,out)
    r <- fill_zero_diagnostic_row(r,diag)
    r$diagnostic_flag_count <- if(nzchar(r$diagnostic_flags))length(strsplit(r$diagnostic_flags,"; ",fixed=TRUE)[[1]]) else 0L
    cal <- suite_calibration(diag,context$d)
    for(nm in c("zero_TS_RMSE","reads_deviation_RMSE","tail_interval_excess",
      "dispersion_log_error","residual_KS_D","predicted_quantile_RMSE","positive_quantile_RMSE",
      "positive_center_RMSE","positive_width_RMSE","positive_upper_tail_RMSE","positive_min_sim"))r[[nm]] <- cal[[nm]]
    write_result_table(cal$tail,out,"04_tail_calibration_effect_sizes.csv")
    write_result_table(cal$positive_bands,out,"06_positive_prediction_band_calibration.csv")
    rows <- diag$fitted_rows
    quantiles <- t(apply(sweep(diag$res_check$simulatedResponse,1,context$d$Ruegeria_reads_total,"/"),1,quantile,c(.025,.975)))
    rows$row_prediction_95_low <- quantiles[,1];rows$row_prediction_95_high <- quantiles[,2]
    extreme <- rows[order(rows$MC_reads/rows$Ruegeria_reads_total,decreasing=TRUE),,drop=FALSE]
    write_result_table(head(extreme,10),out,"05_extreme_row_audit.csv")
    rv <- joint_row_values(candidate$primary)
    write_result_table(data.frame(New_ID=context$d$New_ID,log_reads_c=context$d$log_reads_c,
      DO_c=context$d$DO_c,TS=context$d$TS,fixed_eta=rv$fixed_eta,m_i=rv$m,
      extra_zero_pi=rv$pi,Observation_SD=rv$sd,alpha_i=rv$alpha),out,"03_row_specific_components.csv")
    rm(diag);gc(FALSE)
  },error=function(e) {
    r$status <<- "数值有效；诊断或诊断汇总失败"
    r$diagnostic_status <<- "failed"
    r$notes <<- paste(r$notes,conditionMessage(e),sep="; ")
    writeLines(enc2utf8(conditionMessage(e)),file.path(out,"00_diagnostic_error.txt"))
  })
  # Deduplicated no-zero candidates get the same fit AND diagnostic values.
  rows <- list();files <- list()
  for(model in job$models) {
    rr <- r;rr$model <- model
    if(length(job$models)>1L && model!=id)rr$fit_source <- paste0("exact_no_zero_alias_of_",id)
    rows[[model]] <- rr;files[[model]] <- candidate_file
  }
  metrics <- do.call(rbind,rows)
  write_excel_csv(metrics,file.path(out,"00_candidate_metrics.csv"))
  list(MC_ID=job$MC_ID,models=job$models,metrics=metrics,candidate_files=files)
}

suite_task_worker <- function(job) {
  started <- Sys.time();con <- file(job$log_file,open="at",encoding="UTF-8")
  depth <- sink.number();sink(con)
  on.exit({while(sink.number()>depth)sink();close(con)},add=TRUE)
  writeLines(c("running",format(started),paste("PID",Sys.getpid())),job$status_file)
  cat(job$MC_ID," / ",paste(job$models,collapse=", "),"\n",sep="");flush.console()
  ans <- tryCatch(suite_run_one(job),error=function(e) {
    context <- .GLMM_suite_contexts[[job$MC_ID]]
    rr <- lapply(job$models,function(m) {
      spec <- if(m=="SENS_TOP1")utils::modifyList(.GLMM_suite_specs$D0,list(id=m)) else .GLMM_suite_specs[[m]]
      r <- suite_metric_row(context,spec);r$status <- "任务失败";r$notes <- conditionMessage(e)
      r$sensitivity_only <- isTRUE(job$sensitivity);r
    })
    writeLines(enc2utf8(conditionMessage(e)),file.path(job$output_dir,"00_error.txt"))
    list(MC_ID=job$MC_ID,models=job$models,metrics=do.call(rbind,rr),candidate_files=list())
  })
  ans$elapsed_seconds <- as.numeric(difftime(Sys.time(),started,units="secs"))
  ans$worker_pid <- Sys.getpid();ans$run_signature <- .GLMM_suite_settings$run_signature
  suite_atomic_save(ans,job$result_file)
  writeLines(c("complete",format(Sys.time())),job$status_file);ans
}

suite_dispatch <- function(pool,jobs,on_result,send_call=parallel:::sendCall,
  receive_result=parallel:::recvOneResult) {
  if(!length(jobs))return(invisible(NULL))
  submit <- function(node,j)send_call(pool[[node]],suite_task_worker,list(jobs[[j]]),tag=j)
  # Per-MC dependency graph: a slow model no longer blocks the next stage for
  # every other MC. Dependencies absent from this pending list are completed.
  keys <- unlist(lapply(jobs,function(j)paste(j$MC_ID,j$models,sep="::")))
  owners <- rep(seq_along(jobs),lengths(lapply(jobs,`[[`,"models")))
  dependencies <- lapply(jobs,function(j) {
    dep <- match(paste(j$MC_ID,j$depends %||% character(),sep="::"),keys)
    unique(owners[dep[!is.na(dep)]])
  })
  launched <- completed <- rep(FALSE,length(jobs));idle <- seq_along(pool)
  while(!all(completed)) {
    ready <- which(!launched & vapply(dependencies,function(d)all(completed[d]),logical(1)))
    if(length(ready))ready <- ready[order(vapply(jobs[ready],`[[`,integer(1),"phase"))]
    while(length(idle)&&length(ready)) {
      node <- idle[1];j <- ready[1];idle <- idle[-1];ready <- ready[-1]
      submit(node,j);launched[j] <- TRUE
    }
    if(!any(launched & !completed))stop("No runnable tasks: dependency graph has a cycle.")
    received <- receive_result(pool)
    if(inherits(received$value,"try-error"))stop(as.character(received$value))
    completed[received$tag] <- TRUE;idle <- c(idle,received$node)
    on_result(received$tag,received$value)
  }
}

suite_worker_count <- function(n_tasks) {
  requested <- zero_penalty_worker_count(n_tasks)
  memory <- c(total=NA_real_,available=NA_real_)
  if(.Platform$OS.type=="windows") {
    command <- paste0("$m=Get-CimInstance Win32_OperatingSystem -ErrorAction Stop; ",
      "Write-Output ([string]($m.TotalVisibleMemorySize/1048576)+','+",
      "[string]($m.FreePhysicalMemory/1048576))")
    info <- tryCatch(suppressWarnings(system2("powershell.exe",
      c("-NoProfile","-NonInteractive","-Command",shQuote(command)),
      stdout=TRUE,stderr=FALSE,timeout=8)),error=function(e)character())
    if(length(info)) {
      values <- suppressWarnings(as.numeric(strsplit(tail(info,1),",",fixed=TRUE)[[1]]))
      if(length(values)==2L && all(is.finite(values)) && all(values>0))memory[] <- values
    }
  }
  cap <- requested
  if(Sys.getenv("GLMM_MEMORY_GUARD","1")!="0") {
    cap <- if(all(is.finite(memory)))max(1L,floor(min(memory["total"]-4,
      memory["available"]-1.5)/suite_worker_memory_gb)) else min(8L,requested)
  }
  count <- as.integer(max(1L,min(requested,cap)))
  attr(count,"audit") <- data.frame(requested_after_core_task_cap=requested,
    total_memory_GB=memory["total"],available_memory_GB=memory["available"],
    reserved_per_worker_GB=suite_worker_memory_gb,selected_workers=as.integer(count),
    threads_per_worker=1L,memory_guard=Sys.getenv("GLMM_MEMORY_GUARD","1")!="0",
    row.names=NULL)
  count
}

suite_problem_audit <- function(contexts,run_dir) {
  rows <- lapply(contexts,function(ctx) {
    d <- ctx$d
    do.call(rbind,lapply(levels(d$TS),function(ts) {
      a <- d[d$TS==ts,,drop=FALSE]
      data.frame(MC_ID=ctx$MC_ID,TS=ts,n=nrow(a),n_zero=sum(a$MC_reads==0),
        zero_fraction=mean(a$MC_reads==0),DO_c=mean(a$DO_c),
        distinct_DO_values=length(unique(a$DO_c)),median_reads=median(a$Ruegeria_reads_total))
    }))
  })
  write_excel_csv(do.call(rbind,rows),file.path(run_dir,"06_TS_DO_zero_design_audit.csv"))
  d <- contexts[[1]]$d;design <- model.matrix(~DO_c+TS,d)
  writeLines(enc2utf8(c(paste("Independent observed TS levels:",nlevels(d$TS)),
    paste("Rank / columns for DO + fixed TS:",qr(design)$rank,"/",ncol(design)),
    "DO is shared within TS in this dataset. Rows within TS do not add independent DO levels.",
    "Zero-DO and zero-TS are competing descriptions; this run first refines the common structure and defers grouped validation.",
    "Zero-TS SD is fixed at the configured value; seven observed TS offsets are shrunk and Laplace-integrated.")),
    file.path(run_dir,"06_design_interpretation.txt"))
}

suite_relative_gain <- function(base,current,floor) {
  if(!is.finite(base)||!is.finite(current))return(NA_real_)
  (base-current)/max(base,floor)
}

suite_report <- function(results,contexts,specs,run_dir,complete=FALSE) {
  tab <- bind_flexible_tables(lapply(Filter(Negate(is.null),results),`[[`,"metrics"))
  if(!nrow(tab))return(invisible(NULL))
  measures <- c(zero="zero_TS_RMSE",reads="reads_deviation_RMSE",tail="tail_interval_excess",
    dispersion="dispersion_log_error",quantile="predicted_quantile_RMSE",
    positive="positive_quantile_RMSE",positive_center="positive_center_RMSE",
    positive_width="positive_width_RMSE",positive_tail="positive_upper_tail_RMSE")
  floors <- c(zero=.005,reads=.02,tail=.05,dispersion=.1,quantile=.02,
    positive=.02,positive_center=.02,positive_width=.02,positive_tail=.02)
  for(reference in c("D0","Z_TS_SD_Q")) {
    for(nm in names(measures))tab[[paste0(nm,"_gain_vs_",reference)]] <- NA_real_
    tab[[paste0("DO_change_vs_",reference)]] <- NA_real_
    for(mc in names(contexts)) {
      ii <- which(tab$MC_ID==mc);b <- tab[ii[tab$model[ii]==reference],,drop=FALSE]
      if(nrow(b)!=1L||!isTRUE(b$numerical_valid))next
      for(j in ii) {
        if(isTRUE(tab$numerical_valid[j])&&!isTRUE(tab$sensitivity_only[j]))
          for(nm in names(measures))tab[[paste0(nm,"_gain_vs_",reference)]][j] <-
            suite_relative_gain(b[[measures[nm]]],tab[[measures[nm]]][j],floors[nm])
        tab[[paste0("DO_change_vs_",reference)]][j] <- tab$DO_beta[j]-b$DO_beta
      }
    }
  }
  directions <- lapply(names(contexts),function(mc) {
    a <- tab[tab$MC_ID==mc & !tab$sensitivity_only,,drop=FALSE]
    b <- a[a$model=="Z_TS_SD_Q",,drop=FALSE]
    if(nrow(b)!=1L)return(data.frame(MC_ID=mc,baseline_flags="pending",next_step="等待共同参考模型。"))
    newer <- a[a$model%in%c("SD_M","ALPHA_M","SD_ALPHA_M"),,drop=FALSE]
    improved <- function(column)newer$model[newer$numerical_valid &
      is.finite(newer[[column]]) & newer[[column]]>=.20]
    p <- improved("positive_gain_vs_Z_TS_SD_Q")
    c <- improved("positive_center_gain_vs_Z_TS_SD_Q")
    w <- improved("positive_width_gain_vs_Z_TS_SD_Q")
    t <- improved("tail_gain_vs_Z_TS_SD_Q")
    failed <- newer$model[!newer$numerical_valid]
    remaining <- newer[newer$numerical_valid & newer$diagnostic_status=="ok",,drop=FALSE]
    notes <- c(if(!isTRUE(b$numerical_valid))"共同参考数值无效，先核查起点/积分。",
      if(length(p))paste0("非零条件分位数误差改善至少20%：",paste(p,collapse=", "),"。"),
      if(length(c))paste0("非零中心校准改善：",paste(c,collapse=", "),"。"),
      if(length(w))paste0("非零散布校准改善：",paste(w,collapse=", "),"。"),
      if(length(t))paste0("整体高比例尾部校准改善：",paste(t,collapse=", "),"。"),
      if(length(failed))paste0("数值检查未通过：",paste(failed,collapse=", "),"；先修复拟合。"),
      if(nrow(remaining))paste(vapply(seq_len(nrow(remaining)),function(i)paste0(
        remaining$model[i],"：",if(nzchar(remaining$diagnostic_flags[i]))remaining$diagnostic_flags[i] else
          "未检出主要残差偏差"),character(1)),collapse="；"))
    if(!length(notes))notes <- "等待新候选，或尚无达到20%筛查阈值的改善；检查中心、散布及上尾分别如何变化。"
    data.frame(MC_ID=mc,baseline_flags=b$diagnostic_flags,next_step=paste(notes,collapse=" "))
  })
  direction_table <- bind_flexible_tables(directions)
  main <- tab[!tab$sensitivity_only,,drop=FALSE]
  finite_mean <- function(v)if(any(is.finite(v)))mean(v[is.finite(v)]) else NA_real_
  overview <- do.call(rbind,lapply(names(specs),function(model) {
    a <- main[main$model==model,,drop=FALSE]
    eligible <- a$numerical_valid & a$diagnostic_status=="ok"
    data.frame(model=model,completed_MC=nrow(a),numeric_valid_MC=sum(a$numerical_valid),
      complete_diagnostics_MC=sum(eligible),planned_MC=length(contexts),
      no_major_flags_MC=sum(eligible & !nzchar(a$diagnostic_flags)),
      all_MC_numerically_valid=nrow(a)==length(contexts)&&all(a$numerical_valid),
      all_MC_no_major_flags=nrow(a)==length(contexts)&&all(eligible & !nzchar(a$diagnostic_flags)),
      paired_positive_MC=sum(is.finite(a$positive_gain_vs_Z_TS_SD_Q)),
      positive_gain_equal_MC=finite_mean(a$positive_gain_vs_Z_TS_SD_Q),
      center_gain_equal_MC=finite_mean(a$positive_center_gain_vs_Z_TS_SD_Q),
      width_gain_equal_MC=finite_mean(a$positive_width_gain_vs_Z_TS_SD_Q),
      zero_gain_equal_MC=finite_mean(a$zero_gain_vs_Z_TS_SD_Q),
      tail_gain_equal_MC=finite_mean(a$tail_gain_vs_Z_TS_SD_Q),
      maximum_absolute_DO_change=if(any(is.finite(a$DO_change_vs_Z_TS_SD_Q)))
        max(abs(a$DO_change_vs_Z_TS_SD_Q),na.rm=TRUE) else NA_real_)
  }))
  write_excel_csv(tab,file.path(run_dir,"01_all_candidate_metrics.csv"))
  write_excel_csv(overview,file.path(run_dir,"02_common_structure_comparison.csv"))
  write_excel_csv(direction_table,file.path(run_dir,"00_next_step_directions.csv"))
  if(any(tab$sensitivity_only)) {
    influence <- tab[tab$sensitivity_only,c("MC_ID","numerical_valid","excluded_New_ID","DO_beta","DO_change_vs_D0")]
    write_excel_csv(influence,file.path(run_dir,"07_top1_influence_summary.csv"))
  }
  write_excel_csv(main[,c("MC_ID","model","numerical_valid","diagnostic_status","diagnostic_flags",
    "positive_quantile_RMSE","positive_center_RMSE","positive_width_RMSE","positive_upper_tail_RMSE",
    "zero_TS_RMSE","tail_interval_excess","DO_beta","sd_mean_slope","alpha_mean_slope",
    "alpha_row_min","alpha_row_max","m_min","m_max")],file.path(run_dir,"11_abundance_structure_comparison.csv"))
  text <- c(if(complete)"丰度相关共同结构筛查完成。" else "运行中：已完成任务的部分结果。",
    "D0保留历史参照；新增结构的改善统一相对于Z_TS_SD_Q计算。",
    "m_i=(X_i-colMeans(X_training)) beta；不含Repeat/TS/Observation偏移，随beta在完整似然中联合更新。",
    "SD_M：log sigma=a0+a1*x+a2*x^2+a3*m；ALPHA_M：alpha_i=alpha0+c1*m；SD_ALPHA_M同时加入两项。",
    "每行扰动均中心化、标准化至均值0和SD sigma_i；Observation通过积分处理，不进行mi/ui外层交替插件拟合。",
    "保留六起点、Hessian/梯度、高节点及独立积分核查；数值失败不参与有效候选比较。",
    "非零校准使用条件Y>0的预测分布和midpoint ranks，按正计数预测均值分成四组；这是描述性误差，无新增p值。",
    "若positive_min_sim很小，非零条件误差的Monte Carlo精度不足，需先补充模拟复核。",
    "增益=(参考误差-当前误差)/max(参考误差,预设下限)，不是所有增益都等同普通百分比；平均只包含有效配对MC。",
    "20%仅为误差筛查阈值；同时查看零值、中心、散布和尾部有无恶化，不按诊断p>0.05自动选模型。",
    "DO_beta表示计数logit中心的直接系数；m相关宽度/偏斜亦随DO变化，完整比例效应需结合输出的DO预测曲线。",
    "SD_at_center表示x=0且m=0的参考宽度；每个实际样本的m_i、sigma_i、alpha_i见03_row_specific_components.csv。",
    "本轮先完成共同结构筛查；未做交叉验证、最优lambda选择或普通Wald科学效应区间/p值/BH。",
    capture.output(print(overview,row.names=FALSE)),unlist(lapply(seq_len(nrow(direction_table)),function(i)
      c("",paste0(direction_table$MC_ID[i],"；参考提示：",direction_table$baseline_flags[i]),direction_table$next_step[i]))))
  writeLines(enc2utf8(text),file.path(run_dir,"00_next_steps.txt"))
  table_html <- function(a) {
    if(!nrow(a))return("<p>Pending</p>")
    for(nm in names(a))if(is.numeric(a[[nm]]))a[[nm]] <- signif(a[[nm]],5)
    paste0("<table><tr>",paste0("<th>",html_escape(names(a)),"</th>",collapse=""),"</tr>",
      paste(apply(a,1,function(r)paste0("<tr>",paste0("<td>",html_escape(r),"</td>",collapse=""),"</tr>")),collapse=""),"</table>")
  }
  html <- paste0("<!doctype html><html><head><meta charset='utf-8'><title>Abundance-linked MC structures</title>",
    "<style>body{font-family:Arial,'Microsoft YaHei',sans-serif;margin:24px}table{border-collapse:collapse;margin:20px 0}th,td{border:1px solid #ddd;padding:8px;text-align:left;vertical-align:top}th{background:#eee}pre{white-space:pre-wrap}.scroll{overflow-x:auto}</style></head><body>",
    "<h1>丰度相关共同结构筛查</h1><p>中心、宽度及偏斜联合拟合；复用兼容参考结果；非零诊断单独汇总。</p>",
    table_html(direction_table),"<div class='scroll'>",table_html(overview),"</div><pre>",
    html_escape(paste(text,collapse="\n")),"</pre></body></html>")
  writeLines(enc2utf8(html),file.path(run_dir,"00_next_steps.html"))
  invisible(list(metrics=tab,overview=overview,directions=direction_table))
}

suite_population_prediction <- function(fit,X,x,DO_c) {
  p <- fit$parameters;sp <- fit$joint_spec;parts <- suite_tail_components(sp)
  gh <- normal_quadrature(61L);gr <- normal_quadrature(21L)
  rv <- joint_row_values(fit,x,DO_c=DO_c,TS=NULL,X=X)
  alpha <- rv$alpha;delta <- alpha/sqrt(1+alpha^2)
  weights <- sweep(2*pnorm(outer(alpha,gh$nodes)),2L,gh$weights,"*")
  weights <- weights/rowSums(weights)
  rr_sd <- sqrt(exp(2*p$log_sd_repeat)+exp(2*p$log_sd_ts))
  eta <- as.vector(X%*%p$beta);count <- numeric(length(x))
  for(k in seq_along(parts$weight))for(j in seq_along(gr$nodes)) {
    omega <- rv$sd*parts$scale[k]/sqrt(1-2*delta^2/pi)
    xi <- -omega*delta*sqrt(2/pi)
    prob <- outer(eta+xi+rr_sd*gr$nodes[j],rep(1,length(gh$nodes)))+outer(omega,gh$nodes)
    count <- count+parts$weight[k]*gr$weights[j]*rowSums(plogis(prob)*weights)
  }
  if(!isTRUE(sp$zi))zero <- rep(0,length(x)) else {
    g <- p$gamma_zi+if(isTRUE(sp$zi_reads))p$zi_reads*x else rep(0,length(x))
    if(isTRUE(sp$zi_DO))g <- g+p$zi_DO*DO_c
    if(isTRUE(sp$zi_TS))zero <- as.vector(plogis(outer(g,rep(1,length(gh$nodes)))+
      outer(rep(sp$zero_TS_prior_sd,length(g)),gh$nodes))%*%gh$weights) else zero <- plogis(g)
  }
  list(expected_proportion=(1-zero)*count,count_proportion=count,extra_zero_pi=zero,
    SD=rv$sd,alpha=rv$alpha,m=rv$m)
}

suite_curves <- function(results,contexts,shared,run_dir) {
  available <- Filter(Negate(is.null),results)
  # Use the per-result MC/model identity rather than concatenated list names.
  tables <- list();turns <- list()
  for(ctx in contexts) {
    selected <- available[vapply(available,function(a)identical(a$MC_ID,ctx$MC_ID),logical(1))]
    map <- unlist(lapply(selected,`[[`,"candidate_files"),recursive=FALSE,use.names=TRUE)
    map <- map[names(map)%in%names(.GLMM_suite_specs)]
    for(model in names(map)) {
      candidate <- readRDS(map[[model]])
      if(!isTRUE(candidate$primary_valid))next
      fit <- candidate$primary
      x <- seq(min(ctx$d$log_reads_c),max(ctx$d$log_reads_c),length.out=100L)
      nd <- ctx$d[rep(1L,length(x)),,drop=FALSE]
      nd$Host_type <- factor(rep(host_reference,length(x)),levels=levels(ctx$d$Host_type))
      nd$DO_c <- 0;nd$log_reads_c <- x
      nd <- set_reads_basis_columns(nd,x,shared$mean_basis)
      X <- model.matrix(~DO_c+Host_type+mean_reads_s1+mean_reads_s2+mean_reads_s3,nd)
      value <- suite_population_prediction(fit,X,x,nd$DO_c)
      basis <- predict(shared$mean_basis,newx=x);b0 <- as.numeric(predict(shared$mean_basis,newx=0))
      index <- match(paste0("mean_reads_s",1:3),colnames(fit$X))
      curve <- data.frame(MC_ID=ctx$MC_ID,model=model,log10_reads=x+shared$log_reads_center,
        Ruegeria_reads_total=10^(x+shared$log_reads_center),
        mean_reads_contribution_logit=as.vector(sweep(basis,2,b0)%*%fit$parameters$beta[index]),
        Observation_SD_logit=value$SD,Observation_alpha=value$alpha,m_i=value$m,
        extra_zero_pi_DOc0=value$extra_zero_pi,
        expected_proportion_DOc0=value$expected_proportion)
      tables[[paste(ctx$MC_ID,model)]] <- curve
      a2 <- fit$parameters$sd_reads2 %||% 0;a1 <- fit$parameters$sd_reads %||% 0
      mean_linked <- isTRUE(fit$joint_spec$sd_mean) && abs(fit$parameters$sd_mean)>1e-10
      vertex <- if(!mean_linked && abs(a2)>1e-8)-a1/(2*a2) else NA_real_
      turns[[paste(ctx$MC_ID,model)]] <- data.frame(MC_ID=ctx$MC_ID,model=model,a1=a1,a2=a2,
        vertex_centered_log10_reads=vertex,
        vertex_reads=if(is.finite(vertex)&&abs(vertex+shared$log_reads_center)<300)10^(vertex+shared$log_reads_center) else NA_real_,
        vertex_inside_observed_range=is.finite(vertex)&&vertex>=min(x)&&vertex<=max(x),
        vertex_type=if(mean_linked)"mean_linked_not_pure_quadratic" else
          if(!is.finite(vertex))"none" else if(a2>0)"minimum_SD" else "maximum_SD",
        grid_minimum_SD_reads=10^(x[which.min(value$SD)]+shared$log_reads_center),
        grid_maximum_SD_reads=10^(x[which.max(value$SD)]+shared$log_reads_center),
        grid_extrema_are_reference_host_DOc0=TRUE)
      do <- seq(min(ctx$d$DO_c),max(ctx$d$DO_c),length.out=80L)
      xm <- rep(median(ctx$d$log_reads_c),length(do))
      nd <- ctx$d[rep(1L,length(do)),,drop=FALSE]
      nd$Host_type <- factor(rep(host_reference,length(do)),levels=levels(ctx$d$Host_type))
      nd$DO_c <- do;nd <- set_reads_basis_columns(nd,xm,shared$mean_basis)
      X <- model.matrix(~DO_c+Host_type+mean_reads_s1+mean_reads_s2+mean_reads_s3,nd)
      pred <- suite_population_prediction(fit,X,xm,do)
      do_curve <- data.frame(MC_ID=ctx$MC_ID,model=model,DO_c=do,
        expected_proportion=pred$expected_proportion,count_component_proportion=pred$count_proportion,
        extra_zero_pi=pred$extra_zero_pi,Observation_SD=pred$SD,alpha_i=pred$alpha,m_i=pred$m,host=host_reference,
        reference_Ruegeria_reads=10^(xm+shared$log_reads_center))
      write_result_table(do_curve,ctx$mc_dir,paste0("09_DO_curve_",model,".csv"))
    }
    cur <- bind_flexible_tables(tables[startsWith(names(tables),paste0(ctx$MC_ID," "))])
    if(nrow(cur)) {
      write_result_table(cur,ctx$mc_dir,"09_reads_component_curves.csv")
      if(isTRUE(suite_save_plots)) {
        components <- c(mean_reads_contribution_logit="Mean reads contribution (logit)",
          Observation_SD_logit="Observation SD (logit)",Observation_alpha="Observation shape alpha",
          m_i="Centered fixed-effect predictor m",extra_zero_pi_DOc0="Extra-zero probability (DO_c = 0)",
          expected_proportion_DOc0="Expected MC proportion (DO_c = 0)")
        long <- do.call(rbind,lapply(names(components),function(nm)data.frame(
          model=cur$model,log10_reads=cur$log10_reads,component=unname(components[nm]),value=cur[[nm]])))
        fig <- ggplot2::ggplot(long,ggplot2::aes(log10_reads,value,color=model))+
          ggplot2::geom_line(linewidth=.7)+ggplot2::facet_wrap(~component,scales="free_y",ncol=2)+
          ggplot2::theme_bw(base_size=11)+ggplot2::labs(title=ctx$MC_ID,x="Total Ruegeria reads (log10)",y=NULL,
            subtitle="Fitted point curves; Acropora reference; Gaussian group effects integrated for proportion curves")
        ggplot2::ggsave(file.path(ctx$mc_dir,"09_reads_component_comparison.png"),fig,width=10,height=9,dpi=220,bg="white")
      }
    }
  }
  write_excel_csv(bind_flexible_tables(turns),file.path(run_dir,"08_SD_turning_points.csv"))
  write_excel_csv(bind_flexible_tables(tables),file.path(run_dir,"09_all_reads_component_curves.csv"))
}

run_all_MC_structure_suite <- function() {
  previous <- options(contrasts=c("contr.treatment","contr.poly"));on.exit(options(previous),add=TRUE)
  packages <- c("readxl","glmmTMB","RTMB","DHARMa","qgam","ggplot2","digest")
  missing <- packages[!vapply(packages,requireNamespace,logical(1),quietly=TRUE)]
  if(length(missing))stop("缺少包：",paste(missing,collapse=", "))
  stopifnot(is.finite(suite_lambda_reads),suite_lambda_reads>0,suite_zero_DO_prior_sd>0,
    suite_zero_TS_prior_sd>0,n_sim>=250L,is.finite(suite_worker_memory_gb),suite_worker_memory_gb>0)
  specs <- suite_specs();raw <- readxl::read_excel(input_file,sheet=input_sheet)
  MC_ids <- natural_MC_order(unique(as.character(raw$MC_ID[!is.na(raw$MC_ID)])))
  if(!is.null(MC_selection))MC_ids <- natural_MC_order(unique(MC_selection))
  shared <- build_shared_design(raw,MC_ids,host_reference)
  src <- getOption("GLMM_JOINT_worker_source",file.path(project_dir,"scripts","GLMM_reads_diagnostics_full.R"))
  if(!file.exists(src))stop("Worker source does not exist: ",src)
  resuming <- nzchar(suite_resume_dir)
  run_dir <- if(resuming)suite_resume_dir else file.path(output_root,paste0("GLMM_structure_suite_MC8_",format(Sys.time(),"%Y%m%d_%H%M%S")))
  if(!resuming&&dir.exists(run_dir))run_dir <- paste0(run_dir,"_",Sys.getpid())
  functions <- c("joint_objective","suite_legacy_objective","joint_parameters","joint_random",
    "joint_map","joint_fit","joint_penalty_value","joint_row_values","joint_integrals",
    "simulate_centered_olre","suite_tail_components","joint_interior","joint_starts")
  engine <- digest::digest(list(version=suite_version,script_md5=unname(tools::md5sum(src)),
    functions=lapply(functions,function(nm)body(get(nm))),
    packages=vapply(c("RTMB","TMB","Matrix"),function(p)as.character(utils::packageVersion(p)),character(1))),algo="sha256")
  settings <- list(experiment=suite_version,input_file=input_file,input_md5=unname(tools::md5sum(input_file)),
    MC_ids=MC_ids,n_sim=n_sim,seed=simulation_seed,rng_kind=RNGkind(),specs=specs,mean_df=3L,
    maximum_fit_nodes=maximum_fit_nodes,script_md5=unname(tools::md5sum(src)),engine=engine,
    top1_sensitivity=suite_run_top1_sensitivity,reference_cache=suite_reference_cache,
    structure_reference=suite_structure_reference,
    save_plots=suite_save_plots,CV=FALSE)
  settings$run_signature <- digest::digest(settings,algo="sha256")
  if(resuming) {
    old <- readRDS(file.path(run_dir,"01_run_settings.rds"))
    if(!identical(old$run_signature,settings$run_signature))stop("Resume input/script/statistical settings differ.")
  }
  previous_dirs <- list.dirs(output_root,full.names=TRUE,recursive=FALSE)
  previous_dirs <- sort(previous_dirs[startsWith(basename(previous_dirs),"GLMM_structure_suite_MC8_")],decreasing=TRUE)
  settings$cache_dirs <- unique(c(if(resuming)run_dir,suite_structure_reference,
    head(setdiff(previous_dirs,run_dir),10L)))
  dir.create(run_dir,recursive=TRUE,showWarnings=FALSE)
  saveRDS(settings,file.path(run_dir,"01_run_settings.rds"));saveRDS(shared,file.path(run_dir,"01_shared_model_design.rds"))
  if(!resuming) {
    if(!file.copy(src,file.path(run_dir,"01_analysis_script_used.R")) ||
      !file.copy(input_file,file.path(run_dir,"01_input_used_MC8.xlsx")))stop("Cannot save source/input snapshots.")
  }
  stopifnot(identical(unname(tools::md5sum(file.path(run_dir,"01_analysis_script_used.R"))),settings$script_md5),
    identical(unname(tools::md5sum(file.path(run_dir,"01_input_used_MC8.xlsx"))),settings$input_md5))
  write_excel_csv(data.frame(model=names(specs),zi_DO=vapply(specs,`[[`,logical(1),"zi_DO"),
    zi_TS=vapply(specs,`[[`,logical(1),"zi_TS"),sd_quadratic=vapply(specs,`[[`,logical(1),"sd_quadratic"),
    tail_mix=vapply(specs,`[[`,logical(1),"tail_mix"),sd_mean=vapply(specs,`[[`,logical(1),"sd_mean"),
    alpha_mean=vapply(specs,`[[`,logical(1),"alpha_mean")),file.path(run_dir,"01_candidate_design.csv"))
  writeLines(capture.output(utils::sessionInfo()),file.path(run_dir,"sessionInfo.txt"))
  writeLines("running",file.path(run_dir,"00_run_status.txt"))
  writeLines(run_dir,file.path(output_root,"00_current_GLMM_structure_suite_MC8.txt"))
  methods <- c(
    "Common structural screening; default D0/Z_TS_SD_Q/SD_M/ALPHA_M/SD_ALPHA_M share mean df=3/data/denominator.",
    "D0 and Z_TS_SD_Q reuse compatible 20261012_185348 fits and diagnostics after source/data/spec/package/objective/gradient checks.",
    "New cache: exact data/design/effective statistical spec/engine key. Legacy bridge: provenance, objective and gradient identity; no refit.",
    "Mean: DO_c + Host + ns(centered log10 total Ruegeria reads,df=3) + count Repeat/TS Gaussian intercepts.",
    "Zero: reads slope ridge; Z_DO adds a N(0, configured SD^2) fixed DO slope penalty.",
    "Z_TS adds independent Gaussian zero-TS offsets with FIXED configured SD=1.5 by default, Laplace-integrated. This SD is not estimated.",
    "SD_Q: log(sigma)=a0+a1*x+a2*x^2. a2 is not ridge-penalized; its optimizer bounds and observed-range log SD are checked.",
    "m_i=(X_i-colMeans(X_training)) beta is computed inside the objective with automatic derivatives, excluding count Repeat/TS and Observation offsets.",
    "SD_M adds sd_mean*m_i to log sigma; ALPHA_M adds alpha_mean*m_i to alpha; SD_ALPHA_M adds both. Beta and new slopes are jointly optimized.",
    "New slopes initialize at zero on a checked Z_TS_SD_Q parent; the combined candidate also uses checked single-feature fits as starts.",
    "Each row is centered and standardized: E[u_i]=0, SD[u_i]=sigma_i. New slopes are unpenalized with optimizer bounds sd_mean +/-3, alpha_mean +/-5.",
    "Fit validity also requires all observed row alpha within +/-19.99 and log SD within (-12,6); boundary hits are reported, not silently clipped.",
    "TAIL: 95%/5% centered skew-normal scale mixture with broad/core SD ratio=3; normalization retains total SD sigma(x). These constants are prespecified.",
    "TAIL is a directional broad-component sensitivity candidate, not an estimated outlier class or a Student-t model.",
    "Same six starts, optimizer/Hessian/gradient, higher-node and independent integral validation as the prior analysis.",
    "Conditional 2000 simulations by default: retain fitted count Repeat/TS and fitted zero-TS offsets; redraw Observation/mixture/count/zero indicators.",
    "No-zero MC: exact pi=0; identical zero-only variants are deduplicated, including their diagnostics.",
    "Top1 sensitivity is off by default because the previous run already did it; all main-analysis observations remain included.",
    "Calibrations: equal-TS zero-rate RMSE; reads-quartile deviation RMSE (Uniform target 0.5); tail interval excess; abs(log dispersion ratio); residual KS distance; predicted-quartile residual quantile RMSE.",
    "A 20% error reduction is a prespecified screening aid, not a significance threshold or an automatic final-model selector.",
    "Positive-count conditional midpoint ranks use saved simulations with Y>0, grouped by conditional positive mean; separate center/width/upper-tail errors need no new simulation or qgam fit.",
    "Original logLik, fixed penalties, objective, diagnostics and DO changes are separate; no ordinary AIC/Wald scientific inference or Q ranking.",
    "DO is shared within seven TS. This run concentrates on common structural adequacy and does not do CV.",
    "Population proportion curves integrate both count group distributions and Observation; zero-TS prior also integrated. No confidence intervals.",
    "Target 12-worker flat process pool, limited by CPU/tasks/available memory; one BLAS/OpenMP thread per process; no nested pools.",
    "Per-MC dependency scheduling releases tasks immediately when their own parents finish; no global phase barrier. Six starts and integral checks are retained.",
    "Memory budget uses 1 GB per worker plus available-memory reserve; 00_parallel_memory_budget.csv records requested and actual workers.",
    "Resume: set GLMM_STRUCTURE_RESUME_DIR to an interrupted suite directory; exact statistical/script/input settings required.")
  writeLines(methods,file.path(run_dir,"METHODS.txt"))
  # Preserve earlier lambda evidence without claiming that it applies unchanged to new structures.
  old_settings <- tryCatch(readRDS(file.path(suite_reference_cache,"01_run_settings.rds")),error=function(e)NULL)
  if(!is.null(old_settings)&&identical(old_settings$input_md5,settings$input_md5)) {
    file.copy(file.path(suite_reference_cache,"01_all_model_metrics.csv"),file.path(run_dir,"10_previous_D_lambda_sensitivity.csv"),overwrite=TRUE)
  }
  contexts <- jobs <- results <- list();pool <- NULL;finished <- FALSE
  on.exit({
    if(!is.null(pool))try(parallel::stopCluster(pool),silent=TRUE)
    if(!finished) {
      if(length(results)<length(jobs))length(results) <- length(jobs)
      for(j in seq_along(jobs))if(is.null(results[[j]])&&file.exists(jobs[[j]]$result_file))
        results[j] <- list(tryCatch(readRDS(jobs[[j]]$result_file),error=function(e)NULL))
      if(length(contexts))try(suite_report(results,contexts,specs,run_dir,FALSE),silent=TRUE)
      writeLines("interrupted_or_failed; inspect partial reports and task logs",file.path(run_dir,"00_run_status.txt"))
    }
  },add=TRUE)
  add_job <- function(ctx,models,phase,sensitivity=FALSE) {
    out <- file.path(ctx$mc_dir,models[1]);dir.create(out,recursive=TRUE,showWarnings=FALSE)
    parent_map <- list(Z_TS_SD_Q="D0",Z_DO="D0",Z_TS="D0",SD_Q="D0",TAIL="D0",
      Z_DO_SD_Q=c("SD_Q","Z_DO"),SD_M="Z_TS_SD_Q",ALPHA_M="Z_TS_SD_Q",
      SD_ALPHA_M=c("SD_M","ALPHA_M"))
    depends <- if(sensitivity)"D0" else parent_map[[models[1]]] %||% character()
    if(identical(models[1],"Z_TS_SD_Q")&&isTRUE(suite_include_previous_candidates))
      depends <- c("SD_Q","Z_TS")
    depends <- intersect(depends,names(specs))
    list(MC_ID=ctx$MC_ID,models=models,phase=phase,depends=depends,sensitivity=sensitivity,output_dir=out,
      log_file=file.path(out,"00_compute_log.txt"),status_file=file.path(out,"00_task_status.txt"),
      result_file=file.path(out,"00_task_result.rds"))
  }
  for(i in seq_along(MC_ids)) {
    mc <- MC_ids[i];prepared <- prepare_mc_shared(raw,mc,shared);d <- prepared$data_test;X <- prepared$X
    mc_dir <- file.path(run_dir,mc);dir.create(mc_dir,showWarnings=FALSE)
    write_result_table(d,mc_dir,"01_model_input_rows.csv");write_result_table(prepared$sample_audit,mc_dir,"01_sample_audit.csv")
    initial <- list(beta=rep(0,ncol(X)),log_sd_repeat=log(.25),log_sd_ts=log(.25),
      log_sd_observation=log(.75),alpha=0,b_repeat=rep(0,nlevels(d$Repeat_group_ID)),b_ts=rep(0,nlevels(d$TS)))
    initial$beta[1] <- qlogis((sum(d$MC_reads)+.5)/(sum(d$Ruegeria_reads_total)+1))
    ctx <- list(MC_ID=mc,index=i,d=d,X=X,stub=make_model_stub(d),mc_dir=mc_dir,initial=initial,
      n_source=sum(as.character(raw$MC_ID)==mc,na.rm=TRUE))
    contexts[[mc]] <- ctx
    keys <- vapply(specs,function(s)digest::digest(suite_statistical_spec(suite_effective_spec(s,d))),character(1))
    groups <- split(names(specs),factor(keys,levels=unique(keys)))
    for(models in groups)jobs[[length(jobs)+1L]] <- add_job(ctx,models,
      if("D0"%in%models)1L else if(models[1]%in%c("SD_M","ALPHA_M"))3L else
        if(models[1]=="SD_ALPHA_M")4L else 2L)
    if(isTRUE(suite_run_top1_sensitivity))jobs[[length(jobs)+1L]] <- add_job(ctx,"SENS_TOP1",2L,TRUE)
  }
  suite_problem_audit(contexts,run_dir)
  results <- vector("list",length(jobs))
  for(j in seq_along(jobs))if(file.exists(jobs[[j]]$result_file)) {
    rr <- readRDS(jobs[[j]]$result_file)
    if(!identical(rr$run_signature,settings$run_signature))stop("Task checkpoint settings mismatch.")
    results[j] <- list(rr)
  }
  assign(".GLMM_suite_contexts",contexts,envir=.GlobalEnv)
  assign(".GLMM_suite_specs",specs,envir=.GlobalEnv)
  assign(".GLMM_suite_settings",settings,envir=.GlobalEnv)
  update <- function() {
    write_excel_csv(data.frame(MC_ID=vapply(jobs,`[[`,character(1),"MC_ID"),
      models=vapply(jobs,function(j)paste(j$models,collapse="; "),character(1)),phase=vapply(jobs,`[[`,integer(1),"phase"),
      depends_on=vapply(jobs,function(j)paste(j$depends,collapse="; "),character(1)),
      state=vapply(seq_along(jobs),function(j)if(!is.null(results[[j]]))"complete" else if(file.exists(jobs[[j]]$status_file))
        readLines(jobs[[j]]$status_file,n=1,warn=FALSE) else "queued",character(1))),file.path(run_dir,"00_parallel_task_progress.csv"))
    suite_report(results,contexts,specs,run_dir,FALSE)
  }
  update()
  pending <- which(vapply(results,is.null,logical(1)))
  gc(FALSE)
  worker_count <- suite_worker_count(max(1L,length(pending)))
  write_excel_csv(attr(worker_count,"audit"),file.path(run_dir,"00_parallel_memory_budget.csv"))
  cat("结构筛查：",length(MC_ids)," 个MC，",length(specs),"个共同候选；",length(jobs)," 个去重后任务。\n",
    "目标",parallel_model_workers,"进程；按CPU/内存/任务数实际分配",as.integer(worker_count),"进程。\n",
    "逐MC释放依赖任务；复用已核查参考，联合拟合m相关宽度/偏斜。输出：",run_dir,"\n",sep="");flush.console()
  if(length(pending)) {
    thread_names <- c("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","BLIS_NUM_THREADS","VECLIB_MAXIMUM_THREADS")
    thread_old <- Sys.getenv(thread_names,unset=NA_character_)
    do.call(Sys.setenv,as.list(stats::setNames(rep("1",length(thread_names)),thread_names)))
    pool <- tryCatch(parallel::makePSOCKcluster(worker_count,outfile=file.path(run_dir,"00_worker_messages.log")),finally={
      if(any(is.na(thread_old)))Sys.unsetenv(thread_names[is.na(thread_old)])
      if(any(!is.na(thread_old)))do.call(Sys.setenv,as.list(stats::setNames(thread_old[!is.na(thread_old)],thread_names[!is.na(thread_old)])))
    })
    parallel::clusterCall(pool,function(path) {
      options(GLMM_ZI_definitions_only=TRUE,GLMM_JOINT_definitions_only=TRUE,GLMM_STRUCTURE_definitions_only=TRUE,
        contrasts=c("contr.treatment","contr.poly"))
      source(path,encoding="UTF-8");requireNamespace("RTMB",quietly=TRUE);TRUE
    },src)
    parallel::clusterExport(pool,c(".GLMM_suite_contexts",".GLMM_suite_specs",".GLMM_suite_settings"),envir=.GlobalEnv)
    suite_dispatch(pool,jobs[pending],function(index,result) {
        j <- pending[index];results[j] <<- list(result)
        cat("完成 ",sum(!vapply(results,is.null,logical(1))),"/",length(jobs),"：",result$MC_ID," / ",
          paste(result$models,collapse=", "),"；",round(result$elapsed_seconds/60,1),"分钟\n",sep="");flush.console()
        update()
    })
    parallel::stopCluster(pool);pool <- NULL
  }
  final <- suite_report(results,contexts,specs,run_dir,TRUE)
  suite_curves(results,contexts,shared,run_dir)
  audits <- lapply(jobs,function(job) {
    path <- file.path(job$output_dir,"01_cache_audit.csv")
    if(!file.exists(path))return(NULL)
    a <- tryCatch(read.csv(path,fileEncoding="UTF-8-BOM",stringsAsFactors=FALSE),error=function(e)NULL)
    if(is.null(a)||!nrow(a))return(NULL)
    cbind(MC_ID=job$MC_ID,model=job$models[1],a)
  })
  write_excel_csv(bind_flexible_tables(audits),file.path(run_dir,"03_cache_reuse_audit.csv"))
  saveRDS(list(settings=settings,metrics=final$metrics,overview=final$overview,directions=final$directions,
    task_results=lapply(jobs,`[[`,"result_file")),file.path(run_dir,"GLMM_structure_suite_results.rds"))
  stopifnot(identical(unname(tools::md5sum(input_file)),settings$input_md5),
    identical(unname(tools::md5sum(src)),settings$script_md5),
    nrow(final$metrics[!final$metrics$sensitivity_only,,drop=FALSE])==length(MC_ids)*length(specs))
  writeLines(c("complete",run_dir),file.path(run_dir,"00_run_status.txt"))
  writeLines(run_dir,file.path(output_root,"00_latest_GLMM_structure_suite_MC8.txt"));finished <- TRUE
  cat("\n结构筛查完成。先查看 00_next_steps.html 或 00_next_steps.txt。\n结果目录：",run_dir,"\n",sep="")
  invisible(list(output_dir=run_dir,settings=settings,metrics=final$metrics,summary=final$overview,directions=final$directions))
}

# Saved-result review. These functions never call a fitter or simulation routine.
review_csv <- function(path) read.csv(path,fileEncoding="UTF-8-BOM",stringsAsFactors=FALSE,check.names=FALSE)
review_spec <- function(x) {
  x <- suite_statistical_spec(x)
  x[sort(names(x))]
}
review_source_index <- function(path) {
  settings <- readRDS(file.path(path,"01_run_settings.rds"))
  files <- file.path(path,c("01_analysis_script_used.R","01_input_used_MC8.xlsx"))
  status <- readLines(file.path(path,"00_run_status.txt"),n=1,warn=FALSE)
  if(!identical(status,"complete") || !all(file.exists(files)) ||
     !identical(unname(tools::md5sum(files[1])),settings$script_md5) ||
     !identical(unname(tools::md5sum(files[2])),settings$input_md5) ||
     !identical(unname(tools::md5sum(input_file)),settings$input_md5))
    stop("结果来源未完成、快照不一致或当前输入已变更：",path)
  shared <- readRDS(file.path(path,"01_shared_model_design.rds"))
  list(path=path,settings=settings,shared=shared,
    metrics=review_csv(file.path(path,"01_all_candidate_metrics.csv")))
}
review_worker_count <- function(n_tasks) {
  override <- Sys.getenv("GLMM_PARALLEL_WORKERS",unset="")
  request <- if(nzchar(override))suppressWarnings(as.numeric(override)) else review_parallel_workers
  if(length(request)!=1L||!is.finite(request)||request<1||request!=floor(request))stop("复核并行数应为正整数。")
  cores <- parallel::detectCores(logical=TRUE);if(is.na(cores))cores <- 1L
  request <- max(1L,min(request,max(1L,cores-2L),n_tasks))
  # Reuse the existing memory allocator with a review-specific reservation.
  allocator <- suite_worker_count
  scope <- new.env(parent=environment(allocator))
  scope$suite_worker_memory_gb <- review_worker_memory_gb
  scope$zero_penalty_worker_count <- function(n_tasks)as.integer(min(request,n_tasks))
  environment(allocator) <- scope
  allocator(n_tasks)
}
review_predict <- function(fit,shared,hosts,DO,x) {
  n <- max(length(hosts),length(DO),length(x))
  hosts <- rep(hosts,length.out=n);DO <- rep(DO,length.out=n);x <- rep(x,length.out=n)
  nd <- fit$data[rep(1L,n),,drop=FALSE]
  nd$Host_type <- factor(hosts,levels=levels(fit$data$Host_type))
  if(anyNA(nd$Host_type))stop("预测宿主未出现在拟合设计中。")
  nd$DO_c <- DO;nd$log_reads_c <- x
  nd <- set_reads_basis_columns(nd,x,shared$mean_basis)
  X <- model.matrix(~DO_c+Host_type+mean_reads_s1+mean_reads_s2+mean_reads_s3,nd)
  if(!identical(colnames(X),colnames(fit$X)))stop("预测设计列与拟合不一致。")
  out <- suite_population_prediction(fit,X,x,DO)
  if(any(!is.finite(out$expected_proportion))||any(out$expected_proportion<0|out$expected_proportion>1))
    stop("保存参数的预测超出有效范围。")
  out
}
review_mean_bands <- function(observed,predicted,positive_only=FALSE) {
  keep <- is.finite(observed)&is.finite(predicted)
  if(positive_only)keep <- keep&observed>0
  observed <- observed[keep];predicted <- predicted[keep]
  if(!length(observed))return(data.frame())
  ranks <- rank(predicted,ties.method="average")
  breaks <- unique(quantile(ranks,seq(0,1,.25)))
  bands <- if(length(breaks)>1L)cut(ranks,breaks,include.lowest=TRUE) else factor(rep("All",length(ranks)))
  out <- lapply(split(seq_along(ranks),bands,drop=TRUE),function(j)data.frame(
    n=length(j),observed_mean=mean(observed[j]),predicted_mean=mean(predicted[j])))
  ans <- do.call(rbind,out);rownames(ans) <- NULL;ans$prediction_band <- seq_len(nrow(ans));ans
}
review_one_MC <- function(task) {
  previous <- options(contrasts=c("contr.treatment","contr.poly"));on.exit(options(previous))
  mc <- task$MC_ID;shared <- task$sources[[1]]$shared
  samples <- shared$samples;basis <- as.matrix(samples[,paste0("mean_reads_s",1:3)])
  audits <- metrics <- curves <- contrasts <- hosts <- bands <- positive <- zeros <- tails <- groups <- list()
  manifest <- list();reference_data <- reference_X <- NULL
  track <- function(path) {
    manifest[[length(manifest)+1L]] <<- data.frame(path=normalizePath(path,winslash="/"),md5=unname(tools::md5sum(path)))
  }
  for(model in task$models) {
    chosen <- NULL
    for(source in task$sources) {
      file <- file.path(source$path,mc,model,"fitted_candidate.rds")
      if(!file.exists(file))next
      track(file);cc <- tryCatch(readRDS(file),error=function(e)NULL)
      row <- source$metrics[source$metrics$MC_ID==mc&source$metrics$model==model,,drop=FALSE]
      reason <- "saved numerical checks failed or fit missing"
      f <- cc$primary
      ok <- !is.null(f)&&nrow(row)==1L&&isTRUE(row$numerical_valid)&&
        isTRUE(cc$primary_valid)&&isTRUE(cc$reproducible)&&isTRUE(cc$quadrature_ok)&&isTRUE(cc$direct_integrals_ok)
      if(ok) {
        matched <- match(f$data$New_ID,samples$New_ID)
        sd <- samples[matched,,drop=FALSE]
        want <- suite_effective_spec(suite_specs()[[model]],f$data)
        have <- suite_effective_spec(cc$requested_spec %||% cc$spec,f$data)
        same <- !anyNA(matched)&&!anyDuplicated(f$data$New_ID)&&nrow(f$data)==nrow(samples)&&
          all(vapply(c("DO_c","Host_type","Repeat_group_ID","TS","Ruegeria_reads_total","log_reads_c"),
            function(nm)isTRUE(all.equal(as.character(f$data[[nm]]),as.character(sd[[nm]]))),logical(1)))&&
          isTRUE(all.equal(unname(as.matrix(f$data[,paste0("mean_reads_s",1:3)])),unname(basis[matched,,drop=FALSE]),tolerance=1e-12))&&
          identical(review_spec(have),review_spec(want))
        if(!is.null(reference_data))same <- same&&suite_same_data(f,reference_data,reference_X)
        ok <- isTRUE(same);reason <- if(ok)"saved checks passed; same data/design/spec; no refit" else "data/design/spec mismatch"
      }
      audits[[length(audits)+1L]] <- data.frame(MC_ID=mc,model=model,path=file,eligible=ok,reason=reason)
      if(ok) {chosen <- list(fit=f,row=row,source=source$path,file=file);break}
    }
    if(is.null(chosen)) {
      metrics[[model]] <- data.frame(MC_ID=mc,model=model,numerical_valid=FALSE,
        diagnostic_status="not_available",diagnostic_flags="",source_run="",review_note="待核查；未自动拟合")
      next
    }
    f <- chosen$fit;d <- f$data
    if(is.null(reference_data)){reference_data <- d;reference_X <- f$X}
    row <- chosen$row
    keep <- intersect(c("MC_ID","model","n_used","DO_beta","numerical_valid","diagnostic_status","diagnostic_flags",
      "Uniformity_p","Dispersion_p","Dispersion_ratio","Predicted_quantiles_p","Tail_min_BH_p","Zero_TS_min_BH_p"),names(row))
    row <- row[,keep,drop=FALSE];row$source_run <- chosen$source
    row$reused_previous_run <- !identical(chosen$source,task$sources[[1]]$path)
    row$review_note <- "沿用保存的数值检查；本次未重新优化参数"
    diag_file <- file.path(dirname(chosen$file),"diagnostics","model_analysis_objects.rds")
    if(file.exists(diag_file)) {
      track(diag_file);diag <- readRDS(diag_file);sim <- diag$res_check$simulatedResponse
      valid_diag <- isTRUE(diag$diagnostics_complete)&&is.matrix(sim)&&nrow(sim)==nrow(d)&&ncol(sim)>=250L&&
        identical(as.character(diag$fitted_rows$New_ID),as.character(d$New_ID))&&
        isTRUE(all.equal(diag$fit_model$parameters,f$parameters,tolerance=1e-10,check.attributes=FALSE))&&
        identical(as.numeric(diag$res_check$observedResponse),as.numeric(d$MC_reads))&&
        all(is.finite(sim))&&all(sim>=0)&&all(sweep(sim,1,d$Ruegeria_reads_total,"<="))&&
        all(sim==floor(sim))&&all(rowSums(sim>0)[d$MC_reads>0]>0)
      if(valid_diag) {
        cal <- suite_calibration(diag,d)
        for(nm in c("zero_TS_RMSE","positive_quantile_RMSE","positive_center_RMSE","positive_width_RMSE",
          "positive_upper_tail_RMSE","positive_min_sim","tail_interval_excess"))row[[nm]] <- cal[[nm]]
        row$diagnostic_n_sim <- ncol(sim)
        row$positive_MC_precision_flag <- cal$positive_min_sim<review_min_positive_sim
        all_bands <- review_mean_bands(d$MC_reads/d$Ruegeria_reads_total,rowMeans(sim)/d$Ruegeria_reads_total)
        bands[[model]] <- cbind(MC_ID=mc,model=model,all_bands)
        positive[[model]] <- cbind(MC_ID=mc,model=model,cal$positive_bands)
        z <- diag$zero_group_table;z <- z[z$grouping=="TS",,drop=FALSE]
        z$observed_rate <- z$observed/z$n;z$predicted_rate <- z$simulated_mean/z$n
        zeros[[model]] <- cbind(MC_ID=mc,model=model,z)
        tails[[model]] <- cbind(MC_ID=mc,model=model,cal$tail)
        for(variable in c("Host_type","TS","reads_band")) {
          g <- lapply(split(seq_len(nrow(d)),as.character(d[[variable]])),function(j)data.frame(
            group=as.character(d[[variable]][j[1]]),n=length(j),
            observed_mean=mean(d$MC_reads[j]/d$Ruegeria_reads_total[j]),
            predicted_mean=mean(rowMeans(sim[j,,drop=FALSE])/d$Ruegeria_reads_total[j]),
            mean_DHARMa_residual=mean(diag$res_check$scaledResiduals[j])))
          groups[[paste(model,variable)]] <- cbind(MC_ID=mc,model=model,grouping=variable,do.call(rbind,g))
        }
        rm(sim,diag)
      } else {row$diagnostic_status <- "saved_diagnostic_mismatch";row$review_note <- "拟合有效；保存诊断不匹配，未重做模拟"}
    } else {row$diagnostic_status <- "saved_diagnostics_missing";row$review_note <- "拟合有效；缺少保存诊断，未重做模拟"}
    metrics[[model]] <- row
    xref <- median(samples$log_reads_c)
    do <- seq(min(samples$DO_c),max(samples$DO_c),length.out=80L)
    pred <- review_predict(f,shared,host_reference,do,xref)
    offset <- if("DO_mg_L"%in%names(d))mean(d$DO_mg_L-d$DO_c) else 0
    curve <- data.frame(MC_ID=mc,model=model,DO_c=do,DO_display=do+offset,
      DO_axis=if("DO_mg_L"%in%names(d))"Dissolved oxygen (mg/L)" else "Centered dissolved oxygen",
      expected_proportion=pred$expected_proportion,host=host_reference,
      reference_reads=10^(xref+shared$log_reads_center),source_run=chosen$source)
    curves[[model]] <- curve
    contrasts[[model]] <- data.frame(MC_ID=mc,model=model,DO_low=curve$DO_display[1],DO_high=tail(curve$DO_display,1),
      predicted_low=curve$expected_proportion[1],predicted_high=tail(curve$expected_proportion,1),
      change_percentage_points=100*(tail(curve$expected_proportion,1)-curve$expected_proportion[1]),
      host=host_reference,reference_reads=curve$reference_reads[1],source_run=chosen$source)
    hs <- levels(d$Host_type);hp <- review_predict(f,shared,hs,0,xref)$expected_proportion
    support <- lapply(hs,function(h) {
      j <- which(as.character(d$Host_type)==h)
      data.frame(Host=h,n=length(j),n_TS=length(unique(d$TS[j])),
        reference_DO_in_range=min(d$DO_c[j])<=0&&max(d$DO_c[j])>=0,
        reference_reads_in_range=min(d$log_reads_c[j])<=xref&&max(d$log_reads_c[j])>=xref)
    })
    htab <- do.call(rbind,support);htab$expected_proportion <- hp
    htab$difference_vs_Acropora_pp <- 100*(hp-hp[match(host_reference,hs)])
    htab$reference_DO_c <- 0;htab$reference_reads <- 10^(xref+shared$log_reads_center)
    hosts[[model]] <- cbind(MC_ID=mc,model=model,htab,source_run=chosen$source)
    rm(f,chosen);gc(FALSE)
  }
  answer <- lapply(list(metrics=metrics,DO_curves=curves,DO_contrasts=contrasts,host_predictions=hosts,
    all_calibration=bands,positive_calibration=positive,zero_TS=zeros,tail_calibration=tails,
    group_calibration=groups,source_audit=audits,source_manifest=manifest),bind_flexible_tables)
  answer$MC_ID <- mc;answer
}
review_plots <- function(tables,out) {
  if(!isTRUE(review_save_plots))return(invisible(NULL))
  if(!requireNamespace("ggplot2",quietly=TRUE))stop("作图需要ggplot2。")
  colors <- c(D0="#E64B35",Z_TS_SD_Q="#4DBBD5",SD_M="#3C5488")
  theme <- ggplot2::theme_classic(base_family="Arial",base_size=11)+
    ggplot2::theme(strip.background=ggplot2::element_blank(),strip.text=ggplot2::element_text(face="bold"),
      legend.position="top",plot.caption=ggplot2::element_text(hjust=0,size=9),axis.line=ggplot2::element_line(linewidth=.4))
  save <- function(p,name,w=12,h=7) {
    ggplot2::ggsave(file.path(out,paste0(name,".png")),p,width=w,height=h,dpi=300,bg="white")
    ggplot2::ggsave(file.path(out,paste0(name,".pdf")),p,width=w,height=h,bg="white",device=grDevices::cairo_pdf)
  }
  for(nm in names(tables))if("MC_ID"%in%names(tables[[nm]]))
    tables[[nm]]$MC_ID <- factor(tables[[nm]]$MC_ID,levels=natural_MC_order(unique(tables$metrics$MC_ID)))
  a <- tables$DO_curves
  if(nrow(a)) {
    p <- ggplot2::ggplot(a,ggplot2::aes(DO_display,100*expected_proportion,color=model))+
      ggplot2::geom_line(linewidth=.7)+ggplot2::facet_wrap(~MC_ID,ncol=4,scales="free_y")+
      ggplot2::scale_color_manual(values=colors,drop=FALSE)+theme+
      ggplot2::labs(title="Dissolved oxygen associations across three structures",x=unique(a$DO_axis)[1],
        y="Expected MC proportion (%)",color="Structure",caption=paste0(
          "Reference host: ",host_reference,"; fixed total Ruegeria reads. Curves integrate group and observation distributions.\n",
          "Point estimates only. Unavailable fits are omitted; exact sources are recorded in the CSV tables."))
    save(p,"02_DO_structure_comparison")
  }
  h <- tables$host_predictions
  if(nrow(h)) {
    levels_h <- unique(h$Host);h$host_position <- match(h$Host,levels_h)+
      c(-.22,0,.22)[match(h$model,c("D0","Z_TS_SD_Q","SD_M"))]
    h$support <- ifelse(h$n<10|!h$reference_DO_in_range|!h$reference_reads_in_range,"Limited support","Reference within observed range")
    p <- ggplot2::ggplot(h,ggplot2::aes(100*expected_proportion,host_position,color=model,shape=support))+
      ggplot2::geom_point(size=1.8)+ggplot2::facet_wrap(~MC_ID,ncol=4,scales="free_x")+
      ggplot2::scale_y_continuous(breaks=seq_along(levels_h),labels=levels_h)+
      ggplot2::scale_color_manual(values=colors,drop=FALSE)+ggplot2::scale_shape_manual(values=c("Limited support"=1,"Reference within observed range"=16))+
      theme+ggplot2::labs(title="Host predictions across three structures",x="Expected MC proportion (%)",y=NULL,
        color="Structure",shape=NULL,caption="DO is centered at zero; total Ruegeria reads are fixed. Open symbols flag n < 10 or a reference outside the host's observed range.\nPoint estimates only; sample counts and TS coverage are recorded in 04_host_predictions.csv.")
    save(p,"03_host_structure_comparison",15,11)
  }
  for(kind in c("all_calibration","positive_calibration")) {
    b <- tables[[kind]];if(!nrow(b))next
    if(kind=="positive_calibration") {b$observed_mean <- b$observed_mean_proportion;b$predicted_mean <- b$predicted_positive_mean}
    p <- ggplot2::ggplot(b,ggplot2::aes(100*predicted_mean,100*observed_mean,color=model))+
      ggplot2::geom_abline(slope=1,intercept=0,linetype=2,color="grey45")+
      ggplot2::geom_line(ggplot2::aes(group=model),linewidth=.5)+ggplot2::geom_point(size=2)+
      ggplot2::facet_wrap(~MC_ID,ncol=4,scales="free")+ggplot2::scale_color_manual(values=colors,drop=FALSE)+theme+
      ggplot2::labs(title=if(kind=="all_calibration")"Calibration: all observations" else "Calibration: positive counts",
        x="Predicted mean MC proportion (%)",y="Observed mean MC proportion (%)",color="Structure",
        caption="Points summarize up to four prediction bands within each structure; band membership can differ between structures.\nPositive-count predictions are conditional on Y > 0. All summaries reuse saved simulations; dashed line indicates equality.")
    save(p,if(kind=="all_calibration")"04_all_observation_calibration" else "05_positive_count_calibration")
  }
  z <- tables$zero_TS
  if(nrow(z)) {
    observed <- unique(z[,c("MC_ID","group","observed_rate")])
    p <- ggplot2::ggplot(z,ggplot2::aes(group,100*predicted_rate,color=model,group=model))+
      ggplot2::geom_line(linewidth=.5)+ggplot2::geom_point(size=1.5)+
      ggplot2::geom_point(data=observed,ggplot2::aes(group,100*observed_rate),inherit.aes=FALSE,shape=4,color="black",size=2.3)+
      ggplot2::facet_wrap(~MC_ID,ncol=4,scales="free_y")+ggplot2::scale_color_manual(values=colors,drop=FALSE)+theme+
      ggplot2::labs(title="Zero-count calibration by sampling batch",x="Sampling batch",y="Zero-count frequency (%)",
        color="Structure",caption="Black crosses: observed frequency. Colored points: mean frequency in saved conditional simulations.")
    save(p,"06_zero_TS_calibration")
  }
}
review_report <- function(tables,out) {
  metrics <- tables$metrics;pending <- list()
  for(mc in natural_MC_order(unique(metrics$MC_ID)))for(model in review_models) {
    r <- metrics[metrics$MC_ID==mc&metrics$model==model,,drop=FALSE]
    if(!nrow(r)||!isTRUE(r$numerical_valid)) {
      action <- if(mc=="MC6"&&model=="SD_M")"若此敏感性对照影响判断，再定点核查pi=0边界及零值TS数值稳定性" else "先核对保存结果与来源；确有必要再定点补跑"
      pending[[length(pending)+1L]] <- data.frame(MC_ID=mc,model=model,status="unavailable",action=action,automatic_refit=FALSE)
    } else if(isTRUE(r$reused_previous_run)) {
      pending[[length(pending)+1L]] <- data.frame(MC_ID=mc,model=model,status="valid_previous_fit_reused",
        action="已通过保存数值检查和数据/设计/结构核对；本次复用，未重新计算梯度或拟合",automatic_refit=FALSE)
    }
    if(nrow(r)&&isTRUE(r$numerical_valid)&&!identical(r$diagnostic_status,"ok"))
      pending[[length(pending)+1L]] <- data.frame(MC_ID=mc,model=model,status=r$diagnostic_status,
        action="先核对保存诊断；本次未补充模拟",automatic_refit=FALSE)
  }
  pending <- bind_flexible_tables(pending)
  if(!nrow(pending))pending <- data.frame(MC_ID=character(),model=character(),status=character(),action=character(),automatic_refit=logical())
  write_excel_csv(pending,file.path(out,"00_pending_checks.csv"))
  lines <- c("三个结构的已有结果精简复核完成。",
    "D0=简单参照；Z_TS_SD_Q=主要工作候选；SD_M=丰度相关宽度敏感性对照。",
    "本次只读历史结果，未重新拟合、未生成新随机模拟、未做CV、未选lambda。",
    "数值有效性沿用保存的多起点/梯度/Hessian/积分检查，并核对快照、数据、设计和结构。",
    "无效拟合不产生科学预测；兼容的上一轮有效拟合标注source_run和reused_previous_run。",
    "DO和宿主图是积分后的模型预测点估计；无正式效应区间/p值/BH。",
    "DO对比使用相同宿主、reads及观测DO范围；宿主图使用DO_c=0和相同reads。",
    "宿主n<10或参考条件超出该宿主范围以空心符号标记；该阈值仅用于展示支持度。",
    "校准复用保存模拟；正计数预测与分位数条件于Y>0。各模型分组成员可以不同。",
    "04_host_predictions.csv包含宿主样本量、TS覆盖和参考条件支持度。",
    "00_pending_checks.csv仅列出缺口和历史有效拟合复用情况，不自动触发补跑。",
    "先结合偏差大小和效应稳定性判断，不以全部诊断p>0.05作为停止条件。")
  text <- c(lines,"",capture.output(print(metrics[,intersect(c("MC_ID","model","numerical_valid","diagnostic_flags","reused_previous_run","positive_quantile_RMSE"),names(metrics)),drop=FALSE],row.names=FALSE)),
    "",capture.output(print(tables$DO_contrasts[,intersect(c("MC_ID","model","change_percentage_points"),names(tables$DO_contrasts)),drop=FALSE],row.names=FALSE)),
    "",capture.output(print(pending,row.names=FALSE)))
  writeLines(enc2utf8(text),file.path(out,"00_next_steps.txt"))
  esc <- function(x) {x <- gsub("&","&amp;",as.character(x),fixed=TRUE);x <- gsub("<","&lt;",x,fixed=TRUE);gsub(">","&gt;",x,fixed=TRUE)}
  figures <- list.files(out,pattern="^[0-9].*\\.png$",full.names=FALSE)
  html <- paste0("<!doctype html><html><head><meta charset='utf-8'><title>Three-structure review</title>",
    "<style>body{font-family:Arial,'Microsoft YaHei',sans-serif;margin:24px}pre{white-space:pre-wrap}img{max-width:100%;height:auto}</style></head><body>",
    "<h1>三个结构的已有结果精简复核</h1><pre>",esc(paste(text,collapse="\n")),"</pre>",
    paste0("<figure><img src='",figures,"'><figcaption>",figures,"</figcaption></figure>",collapse=""),"</body></html>")
  writeLines(enc2utf8(html),file.path(out,"00_next_steps.html"));pending
}
run_saved_MC_review <- function() {
  stopifnot(identical(review_models,c("D0","Z_TS_SD_Q","SD_M")),
    is.finite(review_worker_memory_gb),review_worker_memory_gb>0)
  previous <- options(contrasts=c("contr.treatment","contr.poly"));on.exit(options(previous),add=TRUE)
  sources <- lapply(review_source_dirs,review_source_index)
  shared <- sources[[1]]$shared
  for(s in sources)if(!identical(s$settings$input_md5,sources[[1]]$settings$input_md5)||
    !isTRUE(all.equal(s$shared,shared,tolerance=1e-12)))stop("历史来源的共同设计不一致。")
  ids <- natural_MC_order(shared$MC_ids)
  if(!is.null(MC_selection)){if(any(!MC_selection%in%ids))stop("MC_selection超出保存结果。");ids <- natural_MC_order(unique(MC_selection))}
  tasks <- lapply(ids,function(mc)list(MC_ID=mc,sources=sources,models=review_models))
  out <- file.path(output_root,paste0("GLMM_saved_review_MC8_",format(Sys.time(),"%Y%m%d_%H%M%S")))
  if(dir.exists(out))out <- paste0(out,"_",Sys.getpid())
  dir.create(out,recursive=TRUE,showWarnings=FALSE)
  complete <- FALSE;pool <- NULL
  on.exit({if(!is.null(pool))try(parallel::stopCluster(pool),silent=TRUE)
    if(!complete)writeLines("incomplete",file.path(out,"00_run_status.txt"))},add=TRUE)
  writeLines("running",file.path(out,"00_run_status.txt"))
  worker_count <- review_worker_count(length(tasks))
  write_excel_csv(attr(worker_count,"audit"),file.path(out,"00_parallel_memory_budget.csv"))
  cat("读取已有结果：",length(ids),"个MC；三个结构；",as.integer(worker_count),"个轻量复核进程。\n输出：",out,"\n",sep="");flush.console()
  if(as.integer(worker_count)>1L) {
    src <- getOption("GLMM_JOINT_worker_source",file.path(project_dir,"scripts","GLMM_reads_diagnostics_full.R"))
    thread_names <- c("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","BLIS_NUM_THREADS")
    old <- Sys.getenv(thread_names,unset=NA_character_)
    do.call(Sys.setenv,as.list(setNames(rep("1",length(thread_names)),thread_names)))
    pool <- tryCatch(parallel::makePSOCKcluster(as.integer(worker_count),outfile=file.path(out,"00_worker_messages.log")),finally={
      if(anyNA(old))Sys.unsetenv(thread_names[is.na(old)])
      if(any(!is.na(old)))do.call(Sys.setenv,as.list(setNames(old[!is.na(old)],thread_names[!is.na(old)])))})
    parallel::clusterCall(pool,function(path) {
      options(GLMM_ZI_definitions_only=TRUE,GLMM_JOINT_definitions_only=TRUE,GLMM_STRUCTURE_definitions_only=TRUE,GLMM_REVIEW_definitions_only=TRUE)
      source(path,encoding="UTF-8");TRUE
    },src)
    parallel::clusterExport(pool,c("review_min_positive_sim","host_reference"),envir=environment(run_saved_MC_review))
    result <- parallel::parLapply(pool,tasks,review_one_MC)
    parallel::stopCluster(pool);pool <- NULL
  } else result <- lapply(tasks,review_one_MC)
  names_out <- setdiff(names(result[[1]]),"MC_ID")
  tables <- setNames(lapply(names_out,function(nm)bind_flexible_tables(lapply(result,`[[`,nm))),names_out)
  files <- c(metrics="01_three_structure_metrics.csv",DO_curves="02_DO_curves.csv",DO_contrasts="03_DO_range_contrasts.csv",
    host_predictions="04_host_predictions.csv",all_calibration="05_all_observation_calibration.csv",
    positive_calibration="06_positive_count_calibration.csv",zero_TS="07_zero_TS_calibration.csv",tail_calibration="08_tail_calibration.csv",
    group_calibration="09_group_calibration.csv",source_audit="10_saved_fit_audit.csv",source_manifest="11_source_manifest.csv")
  for(nm in names(files))write_excel_csv(tables[[nm]],file.path(out,files[[nm]]))
  manifest <- tables$source_manifest
  if(nrow(manifest)&&!identical(unname(tools::md5sum(manifest$path)),as.character(manifest$md5)))stop("复核期间来源文件发生变化。")
  stopifnot(nrow(tables$metrics)==length(ids)*length(review_models))
  review_plots(tables,out);pending <- review_report(tables,out)
  saveRDS(list(source_dirs=review_source_dirs,models=review_models,metrics=tables$metrics,pending=pending,
    refitted=FALSE,new_simulations=FALSE,CV=FALSE),file.path(out,"GLMM_saved_review_results.rds"))
  src <- getOption("GLMM_JOINT_worker_source",file.path(project_dir,"scripts","GLMM_reads_diagnostics_full.R"))
  if(file.exists(src))file.copy(src,file.path(out,"01_review_script_used.R"))
  writeLines(capture.output(sessionInfo()),file.path(out,"sessionInfo.txt"))
  writeLines(c("complete",out),file.path(out,"00_run_status.txt"))
  writeLines(out,file.path(output_root,"00_latest_GLMM_saved_review_MC8.txt"));complete <- TRUE
  cat("\n精简复核完成。先查看00_next_steps.html；待核查项见00_pending_checks.csv。\n结果目录：",out,"\n",sep="")
  invisible(list(output_dir=out,metrics=tables$metrics,pending=pending))
}

if(!isTRUE(getOption("GLMM_REVIEW_definitions_only",FALSE)) &&
   !isTRUE(getOption("GLMM_STRUCTURE_definitions_only",FALSE)) &&
   !isTRUE(getOption("GLMM_JOINT_definitions_only",FALSE)) &&
   !isTRUE(getOption("GLMM_ZI_definitions_only",FALSE))) {
  all_MC_saved_review_results <- run_saved_MC_review()
  all_MC_saved_review_metrics <- all_MC_saved_review_results$metrics
}
