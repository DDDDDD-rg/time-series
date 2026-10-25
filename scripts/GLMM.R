# ============================================================
# Z_TS_SD_Q：8 个 MC 的统一结构完整重拟合（由用户运行）。
# logit(p) = DO_c + Host + ns(centered log10 total Ruegeria reads, df=3)
#            + Repeat 随机截距 + TS 随机截距 + 中心化偏正态观测扰动。
# log sigma(x) = a0 + a1*x + a2*x^2；alpha 为每个 MC 的常数。
# logit(pi) = gamma0 + gamma_reads*x + z_TS；z_TS ~ N(0,1.5^2)。
# 零值 reads 斜率惩罚：lambda*gamma_reads^2，lambda=0.125。
# 8 个 MC 使用相同结构、样本规则和样条设计，各自估计参数。
# 无观测零值时 pi=0；近零概率及随机效应方差边界重新检查。
# 每个 MC 重新优化、重新模拟2000次；历史结果仅可提供初始参数。
# 保留6起点、Hessian/梯度、多节点/独立积分、零值/正值/尾部诊断。
# 默认最多8个独立进程，内部单线程；实际数量按CPU和可用内存限制。
# beta：固定惩罚/先验的近似Wald检验及BH；不做CV、lambda选择或重拟合bootstrap。
# 顶部参数可修改。主入口位于文件末尾；运行结果写入独立新目录。
# 此脚本内嵌已有验证的计算函数，不需要source其他分析脚本。
# ============================================================
project_dir <- "C:/Users/xzy/Desktop/CUHK_Study/GNBF6010_Research_Projects/Research/Msc main project-time series"
input_file <- file.path(project_dir,"data","GLMM_input_MC8_20261011.xlsx")
input_sheet <- "Input_long"
MC_selection <- c("MC0","MC6","MC10","MC15","MC21","MC25","MC507","MC520")
host_reference <- "Acropora sp."
output_root <- file.path(project_dir,"results")
n_sim <- 2000L
simulation_seed <- 123L
quadrature_nodes <- 31L
quadrature_check_nodes <- 61L
maximum_fit_nodes <- 121L
alpha_starts <- c(-6,-2,2,6)
parallel_model_workers <- 8L # GLMM_PARALLEL_WORKERS可覆盖，MC级任务最多8个。
suite_worker_memory_gb <- 1.0 # 另保留1.5GB可用内存；审计记录实际进程数。
suite_save_plots <- TRUE
suite_lambda_reads <- 0.125
suite_zero_TS_prior_sd <- 1.5
suite_version <- "Z_TS_SD_Q_Wald_essential_v2"
use_previous_initial_values <- TRUE # 只用初值，不复用最终拟合/诊断。
warm_start_dirs <- file.path(output_root,c(
  "GLMM_structure_suite_MC8_20261012_235153",
  "GLMM_structure_suite_MC8_20261012_185348"))
# 中断后把环境变量设为同一运行目录；只跳过已经完整完成的MC。
main_resume_dir <- Sys.getenv("GLMM_Z_TS_SD_Q_RESUME_DIR",unset="")

# ===== Embedded fitting and diagnostic engine =====
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

valid_p <- function(p) is.finite(p) & p >= 0 & p <= 1

adjust_valid_p <- function(p, method = "BH", n = length(p)) {
  answer <- rep(NA_real_, length(p))
  keep <- valid_p(p)
  if (any(keep)) answer[keep] <- stats::p.adjust(p[keep], method = method, n = n)
  answer
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
    cat("Legacy summary omits penalized Wald tests. See 00_beta_significance.csv or 02_beta_Wald_approximate.csv for explicitly labelled approximate inference.\n")
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

formula.mc_olre_fit <- function(x, ...) x$formula

write_excel_csv <- function(data, path) {
  con <- file(path, open = "wt", encoding = "UTF-8")
  on.exit(close(con))
  cat("\ufeff", file = con)
  utils::write.table(data, con, sep = ",", row.names = FALSE,
    col.names = TRUE, na = "", qmethod = "double")
  invisible(path)
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

html_escape <- function(x) {
  x <- gsub("&", "&amp;", as.character(x), fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub('"', "&quot;", x, fixed = TRUE)
  x
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

zero_penalty_worker_count <- function(n_tasks) {
  override <- Sys.getenv("GLMM_PARALLEL_WORKERS", unset="")
  requested <- if(nzchar(override))suppressWarnings(as.numeric(override)) else parallel_model_workers
  if(length(requested)!=1L || !is.finite(requested) || requested<1 || requested!=floor(requested))
    stop("parallel_model_workers / GLMM_PARALLEL_WORKERS 须为正整数。")
  cores <- parallel::detectCores(logical=TRUE)
  if(is.na(cores))cores <- 1L
  as.integer(max(1L,min(requested,max(1L,cores-2L),n_tasks)))
}

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
  r$base_model <- "Z_TS_SD_Q";r$fit_source <- "pending";r$reuse_source <- ""
  r$zi_DO_slope <- r$sd_reads_quadratic <- NA_real_
  r$zero_TS_fixed_prior_SD <- if(isTRUE(spec$zi_TS))spec$zero_TS_prior_sd else NA_real_
  r$tail_mixture <- isTRUE(spec$tail_mix)
  r$zero_TS_RMSE <- r$reads_deviation_RMSE <- r$tail_interval_excess <- NA_real_
  r$dispersion_log_error <- r$residual_KS_D <- r$predicted_quantile_RMSE <- NA_real_
  r$diagnostic_flag_count <- NA_integer_;r$sensitivity_only <- FALSE;r$excluded_New_ID <- ""
  r$diagnostic_n_sim <- n_sim;r$CV_status <- "not_run; single_structure_refit"
  for(nm in c("sd_mean_slope","alpha_mean_slope","alpha_row_min","alpha_row_max","m_min","m_max",
    "positive_quantile_RMSE","positive_center_RMSE","positive_width_RMSE",
    "positive_upper_tail_RMSE","positive_min_sim"))r[[nm]] <- NA_real_
  r$lambda_selection <- "fixed at 0.125; not selected in this run"
  r$requested_penalty_target <- "zero_reads ridge; fixed-SD Gaussian zero_TS prior integrated"
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
  seed_source <- context$seed_source %||% "fresh_start"
  fit_nodes <- max(quadrature_nodes, context$initial_nodes %||% quadrature_nodes)
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

# Approximate inference conditional on the specified zero-component penalty/prior.
# The covariance includes nuisance-parameter uncertainty through the full inverse
# penalized Hessian. No fitting or Hessian recomputation is performed here.
main_beta_inference <- function(candidate, mc) {
  f <- candidate$primary; terms <- colnames(f$X); b <- as.numeric(f$parameters$beta)
  ans <- data.frame(MC_ID=mc,term=terms,estimate=b,SE=NA_real_,z=NA_real_,
    CI95_low=NA_real_,CI95_high=NA_real_,p_approx=NA_real_,
    inference_status="numerical_or_covariance_check_failed",
    method="Approximate penalized-Hessian Wald; fixed lambda/prior",stringsAsFactors=FALSE)
  ix <- which(names(f$fit$par)=="beta"); V <- f$cov_fixed
  if(!isTRUE(candidate$primary_valid)||!isTRUE(f$pdHess)||length(ix)!=length(b)||
     !is.matrix(V)||nrow(V)!=length(f$fit$par)||ncol(V)!=nrow(V))return(ans)
  V <- V[ix,ix,drop=FALSE]
  if(any(!is.finite(V))||!isSymmetric(V,tol=1e-7)||
     inherits(try(chol(V),silent=TRUE),"try-error"))return(ans)
  ans$SE <- sqrt(diag(V));ans$z <- b/ans$SE
  ans$CI95_low <- b-qnorm(.975)*ans$SE;ans$CI95_high <- b+qnorm(.975)*ans$SE
  ans$p_approx <- 2*pnorm(-abs(ans$z));ans$inference_status <- "ok_approximate"
  ans
}

main_beta_family <- function(term) {
  ifelse(term=="DO_c","DO_across_MC8",
    ifelse(grepl("^Host_type",term),"Host_contrasts_across_MC8",
      ifelse(grepl("^mean_reads_s",term),"Reads_basis_across_MC8","Intercept_across_MC8")))
}

main_write_inference <- function(tables,contexts,tab,run_dir,complete) {
  # Include missing/failed MC terms in the pre-specified BH families.
  full <- do.call(rbind,lapply(names(contexts),function(mc)
    data.frame(MC_ID=mc,term=colnames(contexts[[mc]]$X),stringsAsFactors=FALSE)))
  got <- bind_flexible_tables(tables)
  if(nrow(got))full <- merge(full,got,by=c("MC_ID","term"),all.x=TRUE,sort=FALSE)
  else {full$estimate <- full$SE <- full$z <- full$CI95_low <- full$CI95_high <- full$p_approx <- NA_real_
    full$inference_status <- NA_character_;full$method <- "Approximate penalized-Hessian Wald; fixed lambda/prior"}
  full$inference_status[is.na(full$inference_status)] <- "not_available"
  full$BH_family <- main_beta_family(full$term);full$q_BH <- NA_real_
  for(fam in unique(full$BH_family)) {
    j <- which(full$BH_family==fam);ok <- j[is.finite(full$p_approx[j])]
    if(length(ok))full$q_BH[ok] <- p.adjust(full$p_approx[ok],"BH",n=length(j))
  }
  full$BH_significance <- ifelse(is.na(full$q_BH),"",ifelse(full$q_BH<.001,"***",
    ifelse(full$q_BH<.01,"**",ifelse(full$q_BH<.05,"*","ns"))))
  full$report_stage <- if(complete)"final_available_results" else "interim_BH_may_change"
  full$diagnostic_flags <- tab$diagnostic_flags[match(full$MC_ID,tab$MC_ID)]
  tab$inference_method <- "Approximate penalized-Hessian Wald; see 00_beta_significance.csv"
  write_excel_csv(full,file.path(run_dir,"00_beta_significance.csv"))
  do <- full[full$term=="DO_c",,drop=FALSE]
  write_excel_csv(do,file.path(run_dir,"00_DO_significance.csv"))
  for(nm in c("estimate","SE","p_approx","q_BH","CI95_low","CI95_high")) {
    target <- switch(nm,estimate="DO_beta",SE="DO_SE",p_approx="DO_p_approx",q_BH="DO_q_BH",
      CI95_low="DO_CI95_low",CI95_high="DO_CI95_high")
    tab[[target]] <- do[[nm]][match(tab$MC_ID,do$MC_ID)]
  }
  tab
}

# One shared simulation bank; omit qgam fitting and bootstrap outlier tests.
diagnose_skew_comparison_model <- function(id,fit,d,n_sim,seed,output_dir) {
  started <- Sys.time();dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
  simulation <- simulate_centered_olre(fit,n_sim,seed);res <- simulation$res_check
  tests <- list(Uniformity=DHARMa::testUniformity(res,plot=FALSE),
    Dispersion=DHARMa::testDispersion(res,plot=FALSE))
  overall <- do.call(rbind,lapply(names(tests),function(nm)data.frame(
    check=nm,p_value=tests[[nm]]$p.value,status="ok")))
  rows <- cbind(d,expected_MC_proportion=simulation$expected_prop,
    DHARMa_residual=res$scaledResiduals,absolute_deviation=2*abs(res$scaledResiduals-.5))
  write_result_table(rows,output_dir,"09_fitted_rows_and_residuals.csv")
  writeLines(c("Skipped by design: qgam quantile regressions, repeated group KS/variance tests,",
    "bootstrap outlier test and duplicate global abundance/tail tests.",
    "Retained: uniformity, dispersion, grouped zeros, positive and upper-tail calibration.",
    "Calibration tables use the same simulation bank; no simulation refits."),
    file.path(output_dir,"00_diagnostic_scope.txt"))
  if(isTRUE(suite_save_plots)) {
    grDevices::png(file.path(output_dir,"10_residual_overview.png"),width=1800,height=650,res=160)
    tryCatch({par(mfrow=c(1,3),mar=c(4,4,2,1))
      plot(ppoints(nrow(d)),sort(res$scaledResiduals),xlab="Uniform quantile",ylab="Residual quantile");abline(0,1,col=2)
      plot(simulation$expected_prop,res$scaledResiduals,xlab="Expected MC proportion",ylab="Scaled residual");abline(h=.5,col=2)
      plot(d$log_reads_c,res$scaledResiduals,xlab="Centered log10 total reads",ylab="Scaled residual");abline(h=.5,col=2)
    },finally=grDevices::dev.off())
  }
  write_result_table(data.frame(stage="simulation_and_basic_checks",seconds=as.numeric(
    difftime(Sys.time(),started,units="secs"))),output_dir,"00_diagnostic_timing.csv")
  list(fit_model=fit,res_check=res,overall_table=overall,
    tail_table=data.frame(p_BH_within_six=numeric()),
    reads_diagnostics=list(main_table=data.frame(diagnostic=c("Residual_location","Absolute_deviation"),p_BH=NA_real_)),
    fitted_rows=rows,dispersion_ratio=unname(tests$Dispersion$statistic),
    diagnostics_complete=all(is.finite(overall$p_value)),
    diagnostic_status=if(all(is.finite(overall$p_value)))"ok" else "partial",
    simulation_method=simulation$simulation_method,plots=list())
}
# ===== Dedicated Z_TS_SD_Q entry and reporting =====
suite_specs <- function() list(Z_TS_SD_Q=list(id="Z_TS_SD_Q",zi=TRUE,
  zi_reads=TRUE,sd_reads=TRUE,sd_quadratic=TRUE,zi_DO=FALSE,zi_TS=TRUE,
  tail_mix=FALSE,sd_mean=FALSE,alpha_mean=FALSE,lambda_zi_reads=suite_lambda_reads,
  zero_DO_prior_sd=2,zero_TS_prior_sd=suite_zero_TS_prior_sd,
  tail_weight=.05,tail_ratio=3))

# A new run always refits. Completed checkpoints are handled only by the runner.
suite_cache_candidate <- function(context,spec,key,output_dir) NULL

# Always simulate from this run's new fitted model; no legacy diagnostic reuse.
suite_diagnostics <- function(candidate,context,output_dir) {
  f <- candidate$primary;path <- file.path(output_dir,"diagnostics")
  dir.create(path,recursive=TRUE,showWarnings=FALSE)
  write_excel_csv(main_beta_inference(candidate,context$MC_ID),
    file.path(output_dir,"02_beta_Wald_approximate.csv"))
  ans <- diagnose_skew_comparison_model(paste(context$MC_ID,"Z_TS_SD_Q"),
    f,context$d,n_sim,simulation_seed+context$index*100L,path)
  ans <- zero_count_diagnostics(ans,context$d,path)
  suite_atomic_save(ans,file.path(path,"model_analysis_objects.rds"));ans
}

main_initial_values <- function(ctx,spec) {
  ctx$seed_source <- "fresh_start";ctx$seed_path <- "";ctx$seed_md5 <- ""
  ctx$initial_nodes <- quadrature_nodes
  if(!isTRUE(use_previous_initial_values))return(ctx)
  for(dir in warm_start_dirs) {
    path <- file.path(dir,ctx$MC_ID,"Z_TS_SD_Q","fitted_candidate.rds")
    if(!file.exists(path))next
    old <- tryCatch(readRDS(path),error=function(e)NULL)
    if(is.null(old)||!isTRUE(old$primary_valid)||!suite_same_data(old$primary,ctx$d,ctx$X))next
    requested <- old$requested_spec %||% old$spec
    if(!identical(suite_statistical_spec(suite_effective_spec(requested,ctx$d)),
      suite_statistical_spec(suite_effective_spec(spec,ctx$d))))next
    ctx$initial <- old$primary$parameters
    ctx$initial_nodes <- max(quadrature_nodes,old$primary$quadrature_n)
    ctx$seed_source <- "saved_parameters_only_refitted"
    ctx$seed_path <- path;ctx$seed_md5 <- unname(tools::md5sum(path));break
  }
  ctx
}

main_report <- function(results,contexts,jobs,run_dir,complete=FALSE) {
  rows <- lapply(seq_along(jobs),function(j) {
    if(!is.null(results[[j]]))return(results[[j]]$metrics)
    r <- suite_metric_row(contexts[[jobs[[j]]$MC_ID]],suite_specs()[[1]])
    r$status <- "未完成";r
  })
  tab <- bind_flexible_tables(rows)
  tab$positive_MC_precision_flag <- ifelse(is.finite(tab$positive_min_sim),tab$positive_min_sim<50,NA)
  write_excel_csv(tab,file.path(run_dir,"00_all_MC_summary.csv"))
  fields <- c("MC_ID","model","status","numerical_valid","convergence","pdHess","max_gradient",
    "multistart_reproducible","quadrature_ok","direct_integrals_ok","diagnostic_status",
    "Uniformity_p","Dispersion_p","Dispersion_ratio","Outliers_p","Predicted_quantiles_p",
    "Reads_quantiles_p","Mean_proportion_p","Q90_proportion_p","Zero_count_p",
    "Tail_min_BH_p","Zero_TS_min_BH_p","diagnostic_flags","zero_TS_RMSE",
    "positive_quantile_RMSE","positive_center_RMSE","positive_width_RMSE",
    "positive_upper_tail_RMSE","positive_min_sim","positive_MC_precision_flag","tail_interval_excess")
  checks <- tab[,intersect(fields,names(tab)),drop=FALSE]
  write_excel_csv(checks,file.path(run_dir,"00_model_checks.csv"))
  states <- vapply(seq_along(jobs),function(j) {
    if(!is.null(results[[j]]))return("completed")
    if(file.exists(jobs[[j]]$status_file))return(readLines(jobs[[j]]$status_file,n=1,warn=FALSE))
    "queued"
  },character(1))
  write_excel_csv(data.frame(MC_ID=names(contexts),state=states),file.path(run_dir,"00_parallel_task_progress.csv"))
  coefficients <- inference <- list()
  for(result in Filter(Negate(is.null),results)) {
    file <- result$candidate_files[["Z_TS_SD_Q"]]
    if(is.null(file)||!file.exists(file))next
    candidate <- readRDS(file)
    if(!is.null(candidate$primary))inference[[result$MC_ID]] <- main_beta_inference(candidate,result$MC_ID)
    if(!isTRUE(candidate$primary_valid))next
    f <- candidate$primary;p <- f$parameters
    values <- c(stats::setNames(as.numeric(p$beta),colnames(f$X)),
      SD_Repeat=exp(p$log_sd_repeat),SD_TS=exp(p$log_sd_ts),
      a0=p$log_sd_observation,a1=p$sd_reads,a2=p$sd_reads2,alpha=p$alpha,
      gamma0=if(isTRUE(f$zi))p$gamma_zi else NA_real_,
      gamma_reads=if(isTRUE(f$zi))p$zi_reads else NA_real_)
    values <- c(values,stats::setNames(p$b_ts,paste0("count_",levels(f$data$TS))),
      stats::setNames(if(isTRUE(f$zi))p$b_zi_ts else rep(NA_real_,nlevels(f$data$TS)),paste0("zero_",levels(f$data$TS))))
    coefficients[[result$MC_ID]] <- data.frame(MC_ID=result$MC_ID,term=names(values),estimate=unname(values),
      extra_zero_active=isTRUE(f$zi),numerical_valid=TRUE,stringsAsFactors=FALSE)
  }
  cf <- bind_flexible_tables(coefficients)
  if(nrow(cf)) {
    write_excel_csv(cf,file.path(run_dir,"00_model_coefficients_long.csv"))
    wide <- data.frame(term=unique(cf$term),stringsAsFactors=FALSE)
    for(mc in names(contexts)) {
      one <- cf[cf$MC_ID==mc,,drop=FALSE];wide[[mc]] <- one$estimate[match(wide$term,one$term)]
    }
    write_excel_csv(wide,file.path(run_dir,"00_model_coefficients_by_MC.csv"))
  }
  tab <- main_write_inference(inference,contexts,tab,run_dir,complete)
  write_excel_csv(tab,file.path(run_dir,"00_all_MC_summary.csv"))
  txt <- c(if(complete)"Z_TS_SD_Q：8个MC运行流程结束。" else "Z_TS_SD_Q：运行中，以下是部分结果。",
    "每个MC重新拟合同一结构；历史有效参数仅作起点。",
    "均值：DO_c + Host + ns(log10(reads)中心化,df=3) + Repeat/TS随机截距。",
    "log sigma=a0+a1*x+a2*x^2；观测层面中心化偏正态，alpha为常数。",
    "额外零值：gamma0+gamma_reads*x+z_TS；lambda=0.125；zero-TS先验SD=1.5。",
    "模拟条件于本次拟合的Repeat/TS及zero-TS；2000次模拟不重新拟合2000次。",
    "beta输出SE、95%区间、近似Wald p及BH q；固定lambda/先验，用惩罚Hessian逆矩阵。DO仅7个TS，采用大样本近似。",
    "BH：DO跨8MC；Host全部对照系数跨8MC；reads样条全部基系数跨8MC；截距单独一组。未完成MC仍计入计划检验数。",
    "星号依据BH q（<0.05/0.01/0.001）；样条单个基系数不表示整体reads效应。CV、lambda选择、重拟合bootstrap未运行。",
    "numerical_valid只描述数值检查；残差偏差及校准误差另列。",
    "positive_min_sim<50标记正值条件模拟精度不足；不以全部残差p>0.05作为合格标准。",
    "主要输出：00_DO_significance.csv、00_beta_significance.csv、00_model_coefficients_by_MC.csv、00_model_checks.csv。",
    capture.output(print(tab[,intersect(c("MC_ID","status","DO_beta","DO_p_approx","DO_q_BH","diagnostic_flags"),names(tab)),drop=FALSE],row.names=FALSE)))
  writeLines(enc2utf8(txt),file.path(run_dir,"00_key_results.txt"))
  writeLines(enc2utf8(paste0("<!doctype html><meta charset='utf-8'><title>Z_TS_SD_Q</title>",
    "<style>body{font-family:Arial,'Microsoft YaHei';margin:24px}pre{white-space:pre-wrap}</style><pre>",
    html_escape(paste(txt,collapse="\n")),"</pre>")),file.path(run_dir,"00_key_results.html"))
  tab
}

run_Z_TS_SD_Q_MC8 <- function() {
  old_options <- options(contrasts=c("contr.treatment","contr.poly"))
  on.exit(options(old_options),add=TRUE)
  packages <- c("readxl","glmmTMB","RTMB","TMB","Matrix","DHARMa","ggplot2","digest")
  missing <- packages[!vapply(packages,requireNamespace,logical(1),quietly=TRUE)]
  if(length(missing))stop("缺少包，请先安装：",paste(missing,collapse=", "))
  expected <- c("MC0","MC6","MC10","MC15","MC21","MC25","MC507","MC520")
  stopifnot(identical(MC_selection,expected),n_sim>=250L,n_sim==as.integer(n_sim),
    suite_lambda_reads>0,suite_zero_TS_prior_sd>0,quadrature_nodes>=15L,
    maximum_fit_nodes>=quadrature_nodes,is.finite(suite_worker_memory_gb),suite_worker_memory_gb>0)
  raw <- readxl::read_excel(input_file,sheet=input_sheet)
  if(!"MC_ID"%in%names(raw)||!all(expected%in%as.character(raw$MC_ID)))stop("输入未包含完整的8个MC。")
  shared <- build_shared_design(raw,expected,host_reference)
  src <- getOption("GLMM_JOINT_worker_source",file.path(project_dir,"scripts","GLMM.R"))
  if(!file.exists(src))stop("无法找到本脚本：",src)
  specs <- suite_specs()
  versions <- vapply(packages,function(p)as.character(utils::packageVersion(p)),character(1))
  settings <- list(experiment=suite_version,input_file=input_file,input_md5=unname(tools::md5sum(input_file)),
    script_md5=unname(tools::md5sum(src)),MC_ids=expected,specs=specs,n_sim=n_sim,seed=simulation_seed,
    rng_kind=RNGkind(),R_version=R.version.string,mean_df=3L,host_reference=host_reference,versions=versions,
    quadrature_nodes=quadrature_nodes,maximum_fit_nodes=maximum_fit_nodes,
    initial_values_only=use_previous_initial_values,warm_start_dirs=warm_start_dirs,
    save_plots=suite_save_plots,CV=FALSE,beta_inference="approximate_penalized_Hessian_Wald",diagnostics="essential_v2")
  settings$engine <- digest::digest(settings,algo="sha256")
  settings$run_signature <- settings$engine
  settings$cache_dirs <- character()
  resuming <- nzchar(main_resume_dir)
  out <- if(resuming)main_resume_dir else file.path(output_root,paste0("GLMM_Z_TS_SD_Q_MC8_",format(Sys.time(),"%Y%m%d_%H%M%S")))
  if(!resuming&&dir.exists(out))out <- paste0(out,"_",Sys.getpid())
  if(resuming) {
    old <- readRDS(file.path(out,"01_run_settings.rds"))
    if(!identical(old$run_signature,settings$run_signature))stop("恢复失败：脚本、数据、包版本或统计设置已改变。")
  }
  dir.create(out,recursive=TRUE,showWarnings=FALSE)
  if(!resuming) {
    saveRDS(settings,file.path(out,"01_run_settings.rds"));saveRDS(shared,file.path(out,"01_shared_model_design.rds"))
    if(!file.copy(src,file.path(out,"01_analysis_script_used.R"))||
       !file.copy(input_file,file.path(out,"01_input_used_MC8.xlsx")))stop("无法保存输入/脚本快照。")
  }
  writeLines(capture.output(sessionInfo()),file.path(out,"sessionInfo.txt"))
  contexts <- jobs <- list();results <- vector("list",length(expected));pool <- NULL;finished <- FALSE
  on.exit({if(!is.null(pool))try(parallel::stopCluster(pool),silent=TRUE)
    if(!finished)writeLines("interrupted_or_failed; see task checkpoints",file.path(out,"00_run_status.txt"))},add=TRUE)
  for(i in seq_along(expected)) {
    mc <- expected[i];prepared <- prepare_mc_shared(raw,mc,shared);d <- prepared$data_test;X <- prepared$X
    mc_dir <- file.path(out,mc);dir.create(mc_dir,showWarnings=FALSE)
    checkpoint <- file.path(mc_dir,"01_fit_context.rds")
    if(resuming&&file.exists(checkpoint))ctx <- readRDS(checkpoint) else {
      initial <- list(beta=rep(0,ncol(X)),log_sd_repeat=log(.25),log_sd_ts=log(.25),
        log_sd_observation=log(.75),alpha=0,b_repeat=rep(0,nlevels(d$Repeat_group_ID)),b_ts=rep(0,nlevels(d$TS)))
      initial$beta[1] <- qlogis((sum(d$MC_reads)+.5)/(sum(d$Ruegeria_reads_total)+1))
      ctx <- list(MC_ID=mc,index=i,d=d,X=X,stub=make_model_stub(d),mc_dir=mc_dir,
        initial=initial,n_source=sum(as.character(raw$MC_ID)==mc,na.rm=TRUE))
      ctx <- main_initial_values(ctx,specs[[1]])
      saveRDS(ctx,checkpoint)
    }
    contexts[[mc]] <- ctx
    write_result_table(d,mc_dir,"01_model_input_rows.csv")
    write_result_table(prepared$sample_audit,mc_dir,"01_sample_audit.csv")
    task_dir <- file.path(mc_dir,"Z_TS_SD_Q");dir.create(task_dir,showWarnings=FALSE)
    jobs[[i]] <- list(MC_ID=mc,models="Z_TS_SD_Q",phase=1L,depends=character(),sensitivity=FALSE,
      output_dir=task_dir,log_file=file.path(task_dir,"00_compute_log.txt"),
      status_file=file.path(task_dir,"00_task_status.txt"),result_file=file.path(task_dir,"00_task_result.rds"))
    if(resuming&&file.exists(jobs[[i]]$result_file)) {
      rr <- readRDS(jobs[[i]]$result_file)
      if(!identical(rr$run_signature,settings$run_signature))stop("任务检查点设置不一致。")
      # Failed/incomplete diagnostics retry; valid fits with residual flags are complete.
      if(all(rr$metrics$numerical_valid%in%TRUE)&&all(rr$metrics$diagnostic_status%in%"ok"))results[i] <- list(rr)
    }
  }
  write_excel_csv(do.call(rbind,lapply(contexts,function(ctx)data.frame(MC_ID=ctx$MC_ID,
    seed_source=ctx$seed_source,path=ctx$seed_path,md5=ctx$seed_md5,final_fit_reused=FALSE))),file.path(out,"01_initial_values_audit.csv"))
  suite_problem_audit(contexts,out)
  assign(".GLMM_suite_contexts",contexts,envir=.GlobalEnv)
  assign(".GLMM_suite_specs",specs,envir=.GlobalEnv)
  assign(".GLMM_suite_settings",settings,envir=.GlobalEnv)
  main_report(results,contexts,jobs,out)
  writeLines("running",file.path(out,"00_run_status.txt"))
  writeLines(out,file.path(output_root,"00_current_GLMM_Z_TS_SD_Q_MC8.txt"))
  pending <- which(vapply(results,is.null,logical(1)))
  workers <- suite_worker_count(max(1L,length(pending)))
  write_excel_csv(attr(workers,"audit"),file.path(out,"00_parallel_memory_budget.csv"))
  cat("Z_TS_SD_Q：8个MC；待运行",length(pending),"个；并行进程",as.integer(workers),"，每进程1线程。\n结果目录：",out,"\n",sep="");flush.console()
  if(length(pending)) {
    vars <- c("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","BLIS_NUM_THREADS","VECLIB_MAXIMUM_THREADS")
    old_threads <- Sys.getenv(vars,unset=NA_character_)
    do.call(Sys.setenv,as.list(stats::setNames(rep("1",length(vars)),vars)))
    pool <- tryCatch(parallel::makePSOCKcluster(as.integer(workers),outfile=file.path(out,"00_worker_messages.log")),finally={
      if(anyNA(old_threads))Sys.unsetenv(vars[is.na(old_threads)])
      if(any(!is.na(old_threads)))do.call(Sys.setenv,as.list(stats::setNames(old_threads[!is.na(old_threads)],vars[!is.na(old_threads)])))
    })
    parallel::clusterCall(pool,function(path,rng) {
      options(GLMM_MAIN_definitions_only=TRUE,contrasts=c("contr.treatment","contr.poly"))
      do.call(RNGkind,as.list(rng))
      source(path,encoding="UTF-8");requireNamespace("RTMB",quietly=TRUE);TRUE
    },src,settings$rng_kind)
    parallel::clusterExport(pool,c(".GLMM_suite_contexts",".GLMM_suite_specs",".GLMM_suite_settings"),envir=.GlobalEnv)
    suite_dispatch(pool,jobs[pending],function(index,result) {
      j <- pending[index];results[j] <<- list(result)
      main_report(results,contexts,jobs,out)
      cat("完成 ",sum(!vapply(results,is.null,logical(1))),"/8：",result$MC_ID,
        "；",round(result$elapsed_seconds/60,1),"分钟\n",sep="");flush.console()
    })
    parallel::stopCluster(pool);pool <- NULL
  }
  metrics <- main_report(results,contexts,jobs,out,TRUE)
  curve_ok <- tryCatch({suite_curves(results,contexts,shared,out);TRUE},error=function(e) {
    writeLines(enc2utf8(conditionMessage(e)),file.path(out,"00_curve_error.txt"));FALSE
  })
  saveRDS(list(settings=settings,metrics=metrics,task_results=lapply(jobs,`[[`,"result_file")),file.path(out,"GLMM_Z_TS_SD_Q_results.rds"))
  stopifnot(nrow(metrics)==8L,identical(unname(tools::md5sum(input_file)),settings$input_md5),
    identical(unname(tools::md5sum(src)),settings$script_md5))
  status <- if(all(metrics$numerical_valid%in%TRUE)&&all(metrics$diagnostic_status%in%"ok")&&curve_ok)"complete" else "complete_with_issues"
  writeLines(c(status,out),file.path(out,"00_run_status.txt"))
  writeLines(out,file.path(output_root,"00_latest_GLMM_Z_TS_SD_Q_MC8.txt"));finished <- TRUE
  cat("\n运行流程结束：",status,"。先查看00_key_results.html和00_model_checks.csv。\n结果目录：",out,"\n",sep="")
  invisible(list(output_dir=out,metrics=metrics))
}

if(!isTRUE(getOption("GLMM_MAIN_definitions_only",FALSE))) {
  all_MC_Z_TS_SD_Q_results <- run_Z_TS_SD_Q_MC8()
}
