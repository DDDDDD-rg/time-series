# Definitions only. Requires glmm_MC19_pilot_helpers.R definitions in the same environment.
# LOTO = leave one sampling time out; influence only, never a selection or P-value rule.
# Primary estimates and inference always retain all 98 coral samples including TS4/TS6.

mc19_loto_version <- function() "MC19_coral_nitrate_LOTO_v1"

mc19_loto_latest <- function(cfg) file.path(cfg$results_root, "00_latest_GLMM_MC19_LOTO.txt")

mc19_loto_scripts <- function() c("GLMM_leave_one_TS_out.R", "plot_GLMM_MC19_LOTO.R",
  "glmm_MC19_loto_helpers.R", "glmm_MC19_pilot_helpers.R")

mc19_loto_prepare <- function(cfg, source) {
  entries <- lapply(source$manifest$model_files, readRDS)
  center <- source$manifest$center; scale <- source$manifest$scale
  if (length(center) != 1L || length(scale) != 1L ||
      !is.finite(center) || !is.finite(scale) || scale <= 0) stop("原运行标准化参数无效。")
  reference <- entries[[paste(cfg$mc_ids[1], cfg$model_ids[1], sep = "__")]]$data
  shared <- c("New_ID", "TS", "E_raw", "E_z", "Ruegeria_reads_total")
  required <- c(shared, "MC_ID", "MC_reads", "Other_reads")
  if (!all(required %in% names(reference))) stop("原模型缺少必要样本字段。")
  ts <- levels(reference$TS)
  if (nrow(reference) != 98L || !identical(ts, paste0("TS", 1:7)) ||
      anyDuplicated(reference$New_ID)) stop("本入口预期原98样本、TS1–TS7。")
  if (!identical(cfg$mc_ids, c("MC10", "MC6", "MC95", "MC376")) ||
      !identical(cfg$model_ids, c("BB_GLM", "BB_TS")) ||
      !identical(cfg$environment, "Nitrate_mgN_L")) stop("本入口仅支持原4MC/硝酸盐/两种结构。")
  data_by_mc <- list()
  for (mc in cfg$mc_ids) for (model in cfg$model_ids) {
    key <- paste(mc, model, sep = "__"); entry <- entries[[key]]; d <- entry$data
    baseline <- source$results[source$results$MC_ID == mc & source$results$model == model, , drop = FALSE]
    if (!identical(entry$version, cfg$version) || !identical(entry$MC_ID, mc) ||
        !identical(entry$model, model) || !isTRUE(entry$result$numerical_ok) ||
        !isTRUE(baseline$numerical_ok) || is.null(entry$fit) ||
        !all(required %in% names(d)) || nrow(d) != 98L || !identical(levels(d$TS), ts) ||
        anyNA(d[, required]) || anyDuplicated(d$New_ID)) stop("原完整模型不适合直接比较：", key)
    for (field in shared) {
      same <- if (is.numeric(reference[[field]])) mc19_same(reference[[field]], d[[field]]) else
        identical(as.character(reference[[field]]), as.character(d[[field]]))
      if (!same) stop("完整模型间样本/环境/分母不一致：", key, " / ", field)
    }
    if (!mc19_same(d$E_z, (d$E_raw - center) / scale) ||
        !mc19_same(d$MC_reads + d$Other_reads, d$Ruegeria_reads_total) ||
        !all(d$MC_ID == mc) ||
        !mc19_same(entry$result$beta_per_TS_SD, baseline$beta_per_TS_SD)) {
      stop("原计数/标准化/结果表与模型不一致：", key)
    }
    if (is.null(data_by_mc[[mc]])) data_by_mc[[mc]] <- d else {
      if (!mc19_same(data_by_mc[[mc]]$MC_reads, d$MC_reads) ||
          !mc19_same(data_by_mc[[mc]]$Other_reads, d$Other_reads)) stop("两种结构计数不一致：", mc)
    }
  }
  if (any(vapply(split(reference$E_raw, reference$TS), function(x)
      length(unique(x)) != 1L, logical(1)))) stop("同一个TS的硝酸盐数值不唯一。")
  folds <- do.call(rbind, lapply(ts, function(left) {
    keep <- as.character(reference$TS) != left
    d <- reference[keep, , drop = FALSE]
    data.frame(left_out_TS = left, n_left_out = sum(!keep), n_train = sum(keep),
      n_train_TS = length(unique(d$TS)), train_E_low = min(d$E_raw), train_E_high = max(d$E_raw))
  }))
  common_range <- c(max(folds$train_E_low), min(folds$train_E_high))
  if (any(!is.finite(common_range)) || common_range[1] >= common_range[2]) {
    stop("7轮没有共同覆盖的硝酸盐区间，不能按同一范围比较预测变化。")
  }
  folds$common_E_low <- common_range[1]; folds$common_E_high <- common_range[2]
  if (nrow(folds) != 7L || any(folds$n_train_TS != 6L) ||
      length(cfg$mc_ids) * length(cfg$model_ids) * nrow(folds) != 56L) stop("留一任务数不等于56。")
  list(entries = entries, data_by_mc = data_by_mc, folds = folds,
    common_range = common_range, center = center, scale = scale, ts = ts)
}

mc19_loto_add_comparison <- function(entry, full_result, left_out, prepared) {
  row <- if (left_out == "None") full_result else entry$result
  # 原有变化量按每个训练集自身范围计算，改名避免与共同范围混淆。
  names(row)[names(row) == "fixed_part_change_pp"] <- "training_range_fixed_change_pp"
  names(row)[names(row) == "marginal_change_pp"] <- "training_range_marginal_change_pp"
  row$fit_scope <- if (left_out == "None") "Full" else "Leave_one_TS"
  row$left_out_TS <- left_out
  row$n_left_out <- if (left_out == "None") 0L else
    prepared$folds$n_left_out[match(left_out, prepared$folds$left_out_TS)]
  row$full_beta <- full_result$beta_per_TS_SD
  row$full_SE <- full_result$SE
  row$beta_change <- row$beta_change_in_full_SE <- row$SE_ratio <- NA_real_
  row$sign_reversed <- NA
  row$common_E_low <- prepared$common_range[1]
  row$common_E_high <- prepared$common_range[2]
  row$common_fixed_change_pp <- row$common_marginal_change_pp <- NA_real_
  row$common_marginal_p_low <- row$common_marginal_p_high <- NA_real_
  row$prediction_error <- ""
  row$comparison_usable <- isTRUE(row$numerical_ok) &&
    all(is.finite(c(row$beta_per_TS_SD, row$SE, row$full_beta, row$full_SE))) &&
    row$SE > 0 && row$full_SE > 0
  if (row$comparison_usable) {
    row$beta_change <- row$beta_per_TS_SD - row$full_beta
    row$beta_change_in_full_SE <- row$beta_change / row$full_SE
    row$SE_ratio <- row$SE / row$full_SE
    row$sign_reversed <- row$beta_per_TS_SD * row$full_beta < 0
    predictions <- mc19_capture(function() {
      beta <- glmmTMB::fixef(entry$fit)$cond
      z <- (prepared$common_range - prepared$center) / prepared$scale
      eta <- unname(beta["(Intercept)"]) + unname(beta["E_z"]) * z
      marginal <- mc19_marginal_probability(eta, row$TS_SD_logit)
      if (any(!is.finite(marginal))) stop("共同范围预测非有限。")
      c(fixed_change = 100 * diff(stats::plogis(eta)),
        marginal_change = 100 * diff(marginal), p_low = marginal[1], p_high = marginal[2])
    })
    if (!nzchar(predictions$error)) {
      row$common_fixed_change_pp <- unname(predictions$value["fixed_change"])
      row$common_marginal_change_pp <- unname(predictions$value["marginal_change"])
      row$common_marginal_p_low <- unname(predictions$value["p_low"])
      row$common_marginal_p_high <- unname(predictions$value["p_high"])
    } else row$prediction_error <- predictions$error
    if (length(predictions$warnings)) row$warnings <- paste(
      c(row$warnings, predictions$warnings)[nzchar(c(row$warnings, predictions$warnings))], collapse = " | ")
  }
  row$P_value <- NA_real_; row$q_value <- NA_real_
  row$inference_status <- "LOTO_sensitivity_not_calibrated"
  row$source_pilot_diagnostic_status <- full_result$diagnostic_status
  row$diagnostic_status <- if (left_out == "None") "source_pilot_only" else "not_simulated_in_LOTO"
  row
}

mc19_loto_summarise <- function(estimates, cfg) {
  min_or_na <- function(x) if (any(is.finite(x))) min(x[is.finite(x)]) else NA_real_
  max_or_na <- function(x) if (any(is.finite(x))) max(x[is.finite(x)]) else NA_real_
  rows <- list()
  for (mc in cfg$mc_ids) for (model in cfg$model_ids) {
    part <- estimates[estimates$MC_ID == mc & estimates$model == model, , drop = FALSE]
    full <- part[part$fit_scope == "Full", , drop = FALSE]
    loo <- part[part$fit_scope == "Leave_one_TS", , drop = FALSE]
    good <- loo$comparison_usable %in% TRUE
    usable <- loo[good, , drop = FALSE]
    flips <- usable$left_out_TS[usable$sign_reversed %in% TRUE]
    review <- loo$gradient_review_flag %in% TRUE | loo$large_coefficient_flag %in% TRUE |
      loo$TS_variance_boundary %in% TRUE
    largest <- if (nrow(usable)) {
      change <- abs(usable$beta_change)
      paste(usable$left_out_TS[change >= max(change) - 1e-12], collapse = ";")
    } else ""
    rows[[paste(mc, model, sep = "__")]] <- data.frame(MC_ID = mc, model = model,
      environment = cfg$environment, full_beta = full$beta_per_TS_SD, full_SE = full$SE,
      full_Wald95_low_pilot = full$Wald95_low_pilot,
      full_Wald95_high_pilot = full$Wald95_high_pilot,
      n_planned_folds = 7L, n_finished_folds = nrow(loo), n_usable_folds = sum(good),
      unusable_TS = paste(loo$left_out_TS[!good], collapse = ";"),
      numerical_review_TS = paste(loo$left_out_TS[review], collapse = ";"),
      prediction_unavailable_TS = paste(loo$left_out_TS[
        !is.finite(loo$common_marginal_change_pp)], collapse = ";"),
      sign_reversal_count = length(flips), sign_reversal_TS = paste(flips, collapse = ";"),
      beta_min = min_or_na(usable$beta_per_TS_SD), beta_max = max_or_na(usable$beta_per_TS_SD),
      largest_abs_beta_change = max_or_na(abs(usable$beta_change)),
      largest_abs_change_in_full_SE = max_or_na(abs(usable$beta_change_in_full_SE)),
      largest_shift_TS = largest,
      SE_ratio_min = min_or_na(usable$SE_ratio), SE_ratio_max = max_or_na(usable$SE_ratio),
      common_E_low = full$common_E_low, common_E_high = full$common_E_high,
      full_common_marginal_change_pp = full$common_marginal_change_pp,
      LOTO_common_marginal_change_pp_min = min_or_na(usable$common_marginal_change_pp),
      LOTO_common_marginal_change_pp_max = max_or_na(usable$common_marginal_change_pp),
      TS_SD_min = min_or_na(usable$TS_SD_logit), TS_SD_max = max_or_na(usable$TS_SD_logit),
      P_value = NA_real_, q_value = NA_real_, inference_status = "LOTO_sensitivity_not_calibrated",
      evaluation = "Review changes and uncertainty; no automatic stability cutoff")
  }
  result <- do.call(rbind, rows); rownames(result) <- NULL
  result
}

mc19_loto_run <- function(cfg, source_results_dir) {
  mc19_require("glmmTMB")
  source <- mc19_read_run(cfg, source_results_dir)
  # 只读原8个模型及冻结样本，不重新读取或改写原始input。
  source_paths <- c(file.path(source$support, "run_manifest.rds"),
    file.path(source$run_dir, "00_result.csv"), unname(source$manifest$model_files))
  source_hashes <- mc19_hash(source_paths)
  prepared <- mc19_loto_prepare(cfg, source)
  scripts <- file.path(cfg$project_dir, "scripts", "GLMM", mc19_loto_scripts())
  script_hashes <- mc19_hash(scripts)
  tag <- paste0("GLMM_MC19_nitrate_LOTO_", format(Sys.time(), "%Y%m%d_%H%M%S"))
  run_dir <- file.path(cfg$results_root, tag)
  if (dir.exists(run_dir)) stop("LOTO运行目录已存在，请稍后重试。")
  support <- file.path(run_dir, "supporting_files")
  dir.create(file.path(support, "models"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(support, "scripts_used"), showWarnings = FALSE)
  run_dir <- normalizePath(run_dir, winslash = "/", mustWork = TRUE)
  on.exit(mc19_write_text(capture.output(utils::sessionInfo()),
    file.path(support, "sessionInfo.txt")), add = TRUE)
  for (script in scripts) if (!file.copy(script,
      file.path(support, "scripts_used", basename(script)), overwrite = FALSE)) stop("脚本快照保存失败。")
  mc19_write_csv(source_hashes, file.path(support, "source_fingerprints.csv"))
  mc19_write_csv(script_hashes, file.path(support, "script_fingerprints.csv"))
  mc19_write_csv(source$results, file.path(support, "source_full_results.csv"))
  mc19_write_csv(data.frame(model_key = names(source$manifest$model_files),
    source_model_file = unname(source$manifest$model_files)), file.path(support, "full_model_references.csv"))
  mc19_write_csv(prepared$folds, file.path(support, "fold_plan.csv"))
  mc19_write_csv(do.call(rbind, prepared$data_by_mc), file.path(support, "frozen_full_input.csv"))
  mc19_write_csv(data.frame(environment = cfg$environment, center = prepared$center,
    scale = prepared$scale, common_E_low = prepared$common_range[1],
    common_E_high = prepared$common_range[2], rule = "Reuse original seven-TS center/SD in all folds"),
    file.path(support, "frozen_scaling_and_range.csv"))
  manifest <- list(version = mc19_loto_version(), pilot_version = cfg$version, cfg = cfg,
    source_run_dir = source$run_dir, source_fingerprints = source_hashes,
    source_model_files = source$manifest$model_files,
    center = prepared$center, scale = prepared$scale, common_range = prepared$common_range,
    folds = prepared$folds, model_files = character(), analysis_complete = FALSE,
    expected_new_fits = 56L, reused_full_fits = 8L, new_full_fits = 0L,
    source_averaging = FALSE, simulations = FALSE, final_inference = FALSE)
  saveRDS(manifest, file.path(support, "loto_manifest.rds"))
  mc19_write_text("running_56_sensitivity_fits", file.path(support, "run_status.txt"))
  mc19_write_text(c("逐TS留一：原8个完整模型只读复用，新增4MC×2结构×7轮=56次拟合。",
    paste("原完整结果：", source$run_dir),
    "原始input、完整模型、诊断、结果和PDF不改写；不重拟合完整模型，不模拟。",
    "每轮从同一份98样本开始，临时留出一个TS；BB_TS仅使用剩余6个TS。",
    "E_z固定使用原7个TS的中心与SD；共同硝酸盐范围取7个训练集覆盖区间的交集。",
    "00_result.csv为8行汇总；supporting_files/all_estimates.csv为8个完整估计+56个留一估计。",
    "记录β变化、SE比、方向翻转、拟合和边界状态；不自动判稳健/淘汰，不按比例变化设门槛。",
    "beta_change_in_full_SE仅为完整模型SE的参考尺度，不是正式检验统计量。",
    "Wald95仍为先导近似区间，P/q留空；本分析不是TS层级P值校准或独立验证。",
    "边际占比变化为共同范围的点预测差；图上不为此差值构造未经校准的区间。"),
    file.path(run_dir, "00_readme.txt"))
  rows <- list()
  for (mc in cfg$mc_ids) for (model in cfg$model_ids) {
    key <- paste(mc, model, sep = "__")
    full_result <- source$results[source$results$MC_ID == mc & source$results$model == model, , drop = FALSE]
    rows[[paste0(key, "__Full")]] <- mc19_loto_add_comparison(
      prepared$entries[[key]], full_result, "None", prepared)
  }
  estimates <- do.call(rbind, rows); rownames(estimates) <- NULL
  mc19_write_csv(estimates, file.path(support, "all_estimates.csv"))
  new_fits <- 0L
  for (mc in cfg$mc_ids) for (model in cfg$model_ids) for (left in prepared$ts) {
    key <- paste(mc, model, sep = "__")
    full_result <- source$results[source$results$MC_ID == mc & source$results$model == model, , drop = FALSE]
    full_data <- prepared$data_by_mc[[mc]]
    train <- droplevels(full_data[as.character(full_data$TS) != left, , drop = FALSE])
    fold <- prepared$folds[prepared$folds$left_out_TS == left, , drop = FALSE]
    if (nrow(train) != fold$n_train || nlevels(train$TS) != 6L ||
        left %in% as.character(train$TS) ||
        !mc19_same(train$E_z, (train$E_raw - prepared$center) / prepared$scale)) {
      stop("留一数据或冻结标准化检查失败：", key, " / ", left)
    }
    cat("留一拟合 ", new_fits + 1L, "/56：", key, "，留出", left,
      "，训练样本", nrow(train), "\n", sep = ""); flush.console()
    entry <- mc19_fit_one(mc, model, train, cfg) # 同一Beta-binomial及两种均值结构。
    entry$left_out_TS <- left
    entry$loto_version <- mc19_loto_version()
    row <- mc19_loto_add_comparison(entry, full_result, left, prepared)
    file_key <- paste(key, paste0("minus_", left), sep = "__")
    model_path <- file.path(support, "models", paste0(file_key, ".rds"))
    saveRDS(entry, model_path)
    manifest$model_files[file_key] <- model_path
    rows[[file_key]] <- row
    new_fits <- new_fits + 1L
    estimates <- do.call(rbind, rows); rownames(estimates) <- NULL
    mc19_write_csv(estimates, file.path(support, "all_estimates.csv"))
    manifest$n_attempted_fits <- new_fits
    saveRDS(manifest, file.path(support, "loto_manifest.rds"))
  }
  if (new_fits != 56L || nrow(estimates) != 64L ||
      !identical(mc19_hash(source_paths)$md5, source_hashes$md5) ||
      !identical(mc19_hash(scripts)$md5, script_hashes$md5)) {
    stop("任务数、原结果或当前脚本发生变化，本批不发布为latest。")
  }
  summary <- mc19_loto_summarise(estimates, cfg)
  mc19_write_csv(summary, file.path(run_dir, "00_result.csv"))
  manifest$analysis_complete <- TRUE
  manifest$n_attempted_fits <- new_fits
  manifest$finished_at <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
  saveRDS(manifest, file.path(support, "loto_manifest.rds"))
  any_review <- any(!estimates$comparison_usable) || any(estimates$gradient_review_flag) ||
    any(estimates$large_coefficient_flag) || any(estimates$TS_variance_boundary %in% TRUE) ||
    any(nzchar(estimates$prediction_error))
  mc19_write_text(if (any_review) "LOTO_complete_with_fit_review_flags" else
    "LOTO_complete_stability_review_pending", file.path(support, "run_status.txt"))
  mc19_write_text(run_dir, mc19_loto_latest(cfg))
  cat("56次留一拟合结束。结果：", run_dir, "\n")
  invisible(run_dir)
}

mc19_loto_read_run <- function(cfg, run_dir = "") {
  if (!nzchar(run_dir)) {
    latest <- mc19_loto_latest(cfg)
    if (!file.exists(latest)) stop("尚无LOTO结果，请先由用户运行GLMM_leave_one_TS_out.R。")
    run_dir <- trimws(readLines(latest, warn = FALSE, encoding = "UTF-8")[1])
  }
  if (is.na(run_dir) || !nzchar(run_dir)) stop("LOTO运行目录为空。")
  run_dir <- normalizePath(run_dir, winslash = "/", mustWork = TRUE)
  support <- file.path(run_dir, "supporting_files")
  manifest <- readRDS(file.path(support, "loto_manifest.rds"))
  estimates <- utils::read.csv(file.path(support, "all_estimates.csv"), stringsAsFactors = FALSE)
  summary <- utils::read.csv(file.path(run_dir, "00_result.csv"), stringsAsFactors = FALSE)
  if (!identical(manifest$version, mc19_loto_version()) || !isTRUE(manifest$analysis_complete) ||
      manifest$n_attempted_fits != 56L || length(manifest$model_files) != 56L ||
      nrow(estimates) != 64L || nrow(summary) != 8L ||
      sum(estimates$fit_scope == "Full") != 8L || sum(estimates$fit_scope == "Leave_one_TS") != 56L ||
      anyDuplicated(paste(estimates$MC_ID, estimates$model, estimates$left_out_TS))) {
    stop("不是完整的当前LOTO结果。")
  }
  for (mc in cfg$mc_ids) for (model in cfg$model_ids) {
    part <- estimates[estimates$MC_ID == mc & estimates$model == model, , drop = FALSE]
    if (nrow(part) != 8L || !setequal(part$left_out_TS, c("None", paste0("TS", 1:7)))) {
      stop("LOTO结果缺少某轮次：", mc, " / ", model)
    }
  }
  list(run_dir = run_dir, support = support, manifest = manifest, estimates = estimates, summary = summary)
}

mc19_loto_axis_range <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(c(-1, 1))
  limits <- range(c(0, x))
  width <- diff(limits)
  if (width < 1e-8) return(limits + c(-.1, .1))
  limits + c(-1, 1) * width * .08
}

mc19_loto_panel <- function(part, metric, limits) {
  order <- match(c("None", paste0("TS", 1:7)), part$left_out_TS)
  part <- part[order, , drop = FALSE]
  y <- rev(seq_len(nrow(part)))
  valid <- part$comparison_usable %in% TRUE
  values <- if (metric == "beta") part$beta_per_TS_SD else part$common_marginal_change_pp
  valid <- valid & is.finite(values)
  if (metric == "beta") valid <- valid & is.finite(part$Wald95_low_pilot) &
    is.finite(part$Wald95_high_pilot)
  labels <- ifelse(part$left_out_TS == "None", paste0("Full (", part$n_samples, ")"),
    paste0("- ", part$left_out_TS, " (", part$n_samples, ")"))
  full_value <- values[part$left_out_TS == "None"]
  graphics::plot(NA, xlim = limits, ylim = c(.4, 8.6), yaxt = "n", ylab = "",
    xlab = if (metric == "beta") "Nitrate beta per original seven-TS SD" else
      "Marginal predicted change (percentage points)",
    main = paste(part$MC_ID[1], part$model[1], if (metric == "beta") "beta" else "common range"))
  graphics::axis(2, at = y, labels = labels, las = 1, cex.axis = .8)
  graphics::abline(v = 0, lty = 3, col = "grey55")
  if (length(full_value) == 1L && is.finite(full_value)) graphics::abline(v = full_value, lty = 2)
  if (metric == "beta" && any(valid)) graphics::segments(part$Wald95_low_pilot[valid],
    y[valid], part$Wald95_high_pilot[valid], y[valid], col = "grey35")
  if (any(valid)) graphics::points(values[valid], y[valid],
    pch = ifelse(part$left_out_TS[valid] == "None", 18, 16),
    cex = ifelse(part$left_out_TS[valid] == "None", 1.4, .9))
  if (any(!valid)) graphics::text(rep(limits[1], sum(!valid)), y[!valid], "Unavailable", adj = 0, cex = .7)
  review <- valid & (part$gradient_review_flag %in% TRUE | part$large_coefficient_flag %in% TRUE |
    part$TS_variance_boundary %in% TRUE)
  if (any(review)) graphics::points(values[review], y[review], pch = 1, cex = 1.7)
  # Interval and point-estimate definitions are on the guide page; keep axis titles unobscured.
}

mc19_loto_plots <- function(cfg, run_dir = "") {
  # 只读取CSV与本批清单，不读取fit对象，不调用拟合或模拟。
  run <- mc19_loto_read_run(cfg, run_dir)
  figure_dir <- file.path(cfg$figures_root, basename(run$run_dir))
  support <- file.path(figure_dir, "supporting_files")
  dir.create(support, recursive = TRUE, showWarnings = FALSE)
  mc19_pdf(file.path(figure_dir, "00_LOTO_stability.pdf"), function() {
    mc19_note_page("Leave-one-TS-out: nitrate pilot", c(
      "Source: eight saved full-data models; no full-data refits.",
      "New fits: four MCs x two structures x seven omitted TS = 56.",
      "Full-data samples remain the main analysis; omitted TS are sensitivity subsets.",
      "Every fold uses the original seven-TS nitrate center and SD.",
      paste0("Common nitrate range: ", format(run$manifest$common_range[1], digits = 4), " to ",
        format(run$manifest$common_range[2], digits = 4), " mg N/L."),
      "Top panels: beta with approximate Wald 95% intervals. Bottom: marginal changes.",
      "Compare slopes within the same structure; GLMM slopes are conditional.",
      "A surrounding open circle marks a numerical/variance-boundary review flag.",
      "Sign reversals near zero need interpretation together with size and uncertainty.",
      "No automatic stability cutoff; P/q calibration is still pending."))
    for (mc in cfg$mc_ids) {
      part <- run$estimates[run$estimates$MC_ID == mc, , drop = FALSE]
      good <- part$comparison_usable %in% TRUE
      beta_range <- mc19_loto_axis_range(c(part$Wald95_low_pilot[good], part$Wald95_high_pilot[good]))
      change_range <- mc19_loto_axis_range(part$common_marginal_change_pp[good])
      graphics::par(mfrow = c(2, 2), mar = c(5, 6.5, 3, 1), oma = c(1, 0, 0, 0))
      for (model in cfg$model_ids) mc19_loto_panel(part[part$model == model, , drop = FALSE], "beta", beta_range)
      for (model in cfg$model_ids) mc19_loto_panel(part[part$model == model, , drop = FALSE], "change", change_range)
      graphics::mtext("Diamond: full model. Dots: omit one TS. Numbers in labels are training sample counts.",
        outer = TRUE, side = 1, line = .1, cex = .75)
    }
  })
  mc19_write_csv(run$estimates, file.path(support, "plot_data.csv"))
  file.copy(file.path(run$run_dir, "00_result.csv"), file.path(support, "result_snapshot.csv"), overwrite = TRUE)
  for (script in c("plot_GLMM_MC19_LOTO.R", "glmm_MC19_loto_helpers.R")) file.copy(
    file.path(cfg$project_dir, "scripts", "GLMM", script), file.path(support, script), overwrite = TRUE)
  mc19_write_text(c(paste("LOTO results:", run$run_dir),
    paste("Original full models:", run$manifest$source_run_dir),
    "Plotting reads saved tables only; no model fitting or simulation."), file.path(support, "README.txt"))
  mc19_write_text(figure_dir, file.path(cfg$figures_root, "00_latest_GLMM_MC19_LOTO.txt"))
  mc19_write_text("plots_complete", file.path(run$support, "plot_status.txt"))
  cat("逐TS留一合并图目录：", figure_dir, "\n")
  invisible(figure_dir)
}
