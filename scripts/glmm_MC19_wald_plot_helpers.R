# Saved-result plots only. No model fitting, response simulation or P/q recalculation.

mc19_wald_note <- function(title, lines) {
  graphics::par(oma = c(0, 0, 0, 0))
  mc19_note_page(title, unlist(lapply(lines, strwrap, width = 105), use.names = FALSE))
}

mc19_wald_forest <- function(r, title) {
  good <- r$numerical_ok %in% TRUE & is.finite(r$CI95_low) & is.finite(r$CI95_high)
  graphics::par(mfrow = c(1, 1), mar = c(5.5, 8, 4, 2), oma = c(0, 0, 0, 0))
  if (!any(good)) return(mc19_wald_note(title, "No numerically usable effect intervals; all model rows retained in the result table."))
  limits <- range(c(0, r$CI95_low[good], r$CI95_high[good]))
  if (diff(limits) == 0) limits <- limits + c(-1, 1)
  y <- rev(seq_len(nrow(r)))
  graphics::plot(NA, xlim = limits, ylim = c(.3, nrow(r) + .7), yaxt = "n",
    xlab = "Nutrient coefficient per SD across 7 TS (logit scale)", ylab = "", main = title)
  graphics::abline(v = 0, lty = 3, col = "grey50")
  graphics::axis(2, at = y, labels = r$MC_ID, las = 1, cex.axis = .85)
  graphics::segments(r$CI95_low[good], y[good], r$CI95_high[good], y[good])
  graphics::points(r$beta_per_TS_SD[good], y[good],
    pch = ifelse(is.finite(r$beta_q_value[good]) & r$beta_q_value[good] < .05, 19, 1))
  if (any(!good)) graphics::text(rep(limits[1], sum(!good)), y[!good], "Numerical review", adj = 0, cex = .7)
  graphics::mtext("Wald 95% intervals; filled points: BH q < .05 across all 76 TS-intercept tests.", side = 1, line = 3.7, cex = .8)
}

mc19_wald_matrix <- function(results, cfg, model) {
  graphics::par(mfrow = c(1, 1), mar = c(6, 6, 5, 1), oma = c(0, 0, 0, 0))
  graphics::plot(NA, xlim = c(.5, 4.5), ylim = c(.5, 19.5), xaxt = "n", yaxt = "n",
    xlab = "Nutrient", ylab = "", main = paste(model, "- full-data beta and BH q"))
  graphics::axis(1, at = 1:4, labels = unname(cfg$short_names[cfg$environments]))
  graphics::axis(2, at = 19:1, labels = cfg$mc_ids, las = 1, cex.axis = .8)
  for (i in seq_along(cfg$mc_ids)) for (j in seq_along(cfg$environments)) {
    r <- results[results$MC_ID == cfg$mc_ids[i] & results$environment == cfg$environments[j] &
      results$model == model, , drop = FALSE]
    y <- 20 - i
    graphics::rect(j - .5, y - .5, j + .5, y + .5, border = "grey70")
    label <- if (isTRUE(r$numerical_ok) && is.finite(r$beta_P_value))
      sprintf("%+.2f; q=%s%s", r$beta_per_TS_SD, formatC(r$beta_q_value, format = "g", digits = 3),
        if (is.finite(r$beta_q_value) && r$beta_q_value < .05) " *" else "") else "Review / NA"
    graphics::text(j, y, label, cex = .7)
  }
  graphics::mtext("Each cell: beta per nutrient SD; q from one family of 76 TS-intercept Wald tests. * q < .05.",
    side = 1, line = 4.3, cex = .8)
}

mc19_wald_curves <- function(mc, environment, entries, grids, cfg) {
  graphics::par(mfrow = c(1, 1), mar = c(5, 4.5, 4, 1), oma = c(2, 0, 0, 0))
  for (model in cfg$model_ids) {
    key <- mc19_wald_key(mc, environment, model)
    entry <- entries[[key]]; d <- entry$data; grid <- grids[[key]]
    upper <- max(c(d$MC_proportion, if (!is.null(grid)) grid$Wald95_high_pilot), na.rm = TRUE)
    pchs <- c(1, 2, 0, 5, 6, 16, 17)
    graphics::plot(d$E_raw, 100 * d$MC_proportion, pch = pchs[as.integer(d$TS)],
      xlab = unname(cfg$axis_labels[environment]), ylab = "MC proportion (%)",
      ylim = c(0, max(1, 105 * upper)), main = paste(mc, model))
    if (!is.null(grid)) {
      graphics::polygon(c(grid$E_raw, rev(grid$E_raw)),
        100 * c(grid$Wald95_low_pilot, rev(grid$Wald95_high_pilot)),
        col = "grey85", border = NA)
      graphics::lines(grid$E_raw, 100 * grid$fixed_part, lwd = 2)
      graphics::lines(grid$E_raw, 100 * grid$marginal, lwd = 2, lty = 2)
      graphics::points(d$E_raw, 100 * d$MC_proportion, pch = pchs[as.integer(d$TS)])
      conditional <- as.numeric(stats::predict(entry$fit, type = "response", re.form = NULL))
      first <- !duplicated(d$TS)
      graphics::points(d$E_raw[first], 100 * conditional[first], pch = 4, cex = 1.4, lwd = 2)
    } else graphics::mtext("Prediction unavailable; inspect numerical status", side = 3, cex = .65)
    graphics::legend("topleft", legend = levels(d$TS), pch = pchs, bty = "n", cex = .65, ncol = 2)
    graphics::legend("topright", c("Fixed part (b_TS=0)", "Marginal over TS", "Wald 95%: fixed part"),
      lty = c(1, 2, NA), pch = c(NA, NA, 15), col = c("black", "black", "grey85"), bty = "n", cex = .6)
  }
  graphics::mtext("All 98 original samples retained. Crosses: conditional TS fits. Bands: fixed-part uncertainty, not prediction intervals.",
    outer = TRUE, side = 1, line = .3, cex = .72)
}

mc19_wald_diagnostic_page <- function(entry, diagnosis, cfg) {
  d <- entry$data
  graphics::par(mfrow = c(2, 2), mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0))
  DHARMa::plotQQunif(diagnosis$dharma, testUniformity = FALSE, testOutliers = FALSE, testDispersion = FALSE)
  graphics::mtext(paste(entry$MC_ID, cfg$short_names[entry$environment], entry$model), side = 3, outer = TRUE, line = .4)
  graphics::plot(d$E_raw, diagnosis$dharma$scaledResiduals,
    xlab = unname(cfg$axis_labels[entry$environment]), ylab = "Simulation-scaled residual", pch = 16)
  graphics::abline(h = .5, lty = 3)
  graphics::boxplot(split(diagnosis$dharma$scaledResiduals, d$TS), col = "white", xlab = "TS", ylab = "Simulation-scaled residual")
  graphics::abline(h = .5, lty = 3)
  graphics::plot(d$Ruegeria_reads_total, diagnosis$dharma$scaledResiduals, log = "x",
    xlab = "Total Ruegeria reads (log scale)", ylab = "Simulation-scaled residual", pch = 16)
  graphics::abline(h = .5, lty = 3)
}

mc19_wald_plots <- function(cfg, run_dir = "") {
  mc19_require("glmmTMB")
  run <- mc19_wald_read_run(cfg, run_dir)
  entries <- lapply(run$manifest$model_files, readRDS)
  figure_dir <- file.path(cfg$figures_root, basename(run$run_dir))
  support <- file.path(figure_dir, "supporting_files")
  dir.create(file.path(support, "diagnostics"), recursive = TRUE, showWarnings = FALSE)
  grids <- list(); notices <- list()
  for (key in names(entries)) {
    entry <- entries[[key]]
    scaling <- run$manifest$scaling[run$manifest$scaling$environment == entry$environment, ]
    prediction <- mc19_capture(function() mc19_predict_grid(entry, scaling))
    if (!nzchar(prediction$error) && !is.null(prediction$value)) {
      prediction$value$environment <- entry$environment
      grids[key] <- list(prediction$value)
    }
    if (nzchar(prediction$error) || length(prediction$warnings)) notices[[key]] <- data.frame(
      key = key, error = prediction$error, warnings = paste(prediction$warnings, collapse = " | "))
  }
  mc19_pdf(file.path(figure_dir, "00_summary.pdf"), function() {
    mc19_wald_note("MC19 x four nutrients: full-data Wald analysis", c(
      "All 98 original coral samples, TS1-TS7 including TS4 and TS6. No sample averaging.",
      "19 MCs x 4 nutrients = 76 TS-random-intercept models; one BH family of 76 beta tests.",
      "Every formula: cbind(MC_reads, Other_reads) ~ E_z + (1 | TS).",
      "Each fit uses one nutrient: beta-binomial/logit; constant dispersion; no extra zero-inflation.",
      "Beta P comes directly from glmmTMB summary (Wald z test); no bootstrap or LOTO is run.",
      "Wald inference is asymptotic; nutrient variation is represented by seven sampling times.",
      "NO3 high values are concentrated in TS4 and TS6; both are retained in all main fits.",
      "Below detection: NO2 four of seven TS; DIP three of seven TS; input uses L/2 with flags.",
      "Model/diagnostic issues remain visible. No model is selected by significance, AIC or LOTO flips."))
    for (model in cfg$model_ids) mc19_wald_matrix(run$results, cfg, model)
  }, width = 12, height = 10)
  for (j in seq_along(cfg$environments)) {
    environment <- cfg$environments[j]; short <- unname(cfg$short_names[environment])
    r <- run$results[run$results$environment == environment, , drop = FALSE]
    mc19_pdf(file.path(figure_dir, sprintf("%02d_%s_effects_and_curves.pdf", j, short)), function() {
      mc19_wald_note(paste(short, "- effects and curves"), c(
        "Full-data beta/SE and Wald 95% intervals; q adjusted across all 76 TS-intercept tests.",
        "Solid curves: fixed part at b_TS=0; dashed: integrate over the fitted TS-intercept distribution.",
        "Shading describes fixed-part mean uncertainty, not individual-sample prediction intervals.",
        "TS-random-intercept structure only. All biological samples, including extremes, remain visible."))
      for (model in cfg$model_ids) mc19_wald_forest(r[r$model == model, , drop = FALSE], paste(short, model))
      for (mc in cfg$mc_ids) mc19_wald_curves(mc, environment, entries, grids, cfg)
    }, width = 12, height = 9)
  }
  overview_path <- file.path(run$support, "diagnostics", "diagnostic_overview.csv")
  overview <- if (file.exists(overview_path)) utils::read.csv(overview_path, stringsAsFactors = FALSE) else NULL
  if (!is.null(overview)) mc19_require("DHARMa")
  for (environment in cfg$environments) {
    short <- unname(cfg$short_names[environment])
    mc19_pdf(file.path(support, "diagnostics", paste0(short, "_diagnostics.pdf")), function() {
      mc19_wald_note(paste(short, "- diagnostic supplement"), c(
        "Saved unconditional response simulations only; TS random intercepts are redrawn.",
        "Diagnostic P values check distributional fit, not nutrient-beta significance.",
        "Central/zero/TS patterns and tails are reviewed together; no requirement that every test P exceeds .05."))
      for (mc in cfg$mc_ids) for (model in cfg$model_ids) {
        key <- mc19_wald_key(mc, environment, model); entry <- entries[[key]]
        row <- if (!is.null(overview)) overview[mc19_wald_key(overview$MC_ID, overview$environment, overview$model) == key, ] else NULL
        path <- file.path(run$support, "diagnostics", paste0(key, ".rds"))
        diagnosis <- if (file.exists(path)) readRDS(path) else NULL
        usable <- !is.null(row) && nrow(row) == 1L &&
          row$status %in% c("review_flags", "completed_no_automatic_flags") && !is.null(diagnosis) &&
          identical(diagnosis$n_sim, as.integer(row$n_sim)) && identical(diagnosis$seed, as.integer(row$seed)) &&
          identical(diagnosis$model_md5, unname(tools::md5sum(run$manifest$model_files[[key]])))
        if (usable && isTRUE(entry$result$numerical_ok)) {
          mc19_wald_diagnostic_page(entry, diagnosis, cfg)
          mc19_plot_simulation_summary(entry, diagnosis)
        } else mc19_wald_note(key, "Saved diagnostics unavailable or numerical review needed; row retained.")
      }
    }, width = 12, height = 9)
  }
  if (length(grids)) mc19_write_csv(do.call(rbind, grids), file.path(support, "prediction_curves.csv"))
  if (length(notices)) mc19_write_csv(do.call(rbind, notices), file.path(support, "plot_notices.csv"))
  file.copy(file.path(run$run_dir, "00_result.csv"), file.path(support, "result_snapshot.csv"), overwrite = TRUE)
  for (script in cfg$scripts) file.copy(file.path(cfg$project_dir, "scripts", script),
    file.path(support, script), overwrite = TRUE)
  mc19_write_text(c(paste("Source:", run$run_dir),
    "Five main PDFs: overview plus one combined effects/curves PDF per nutrient.",
    "Diagnostic PDFs and all coordinates/scripts are under supporting_files.",
    "Replot only reads saved models, P/q, and diagnostics. No fitting/simulation/P recalculation."),
    file.path(support, "README.txt"))
  mc19_write_text(figure_dir, file.path(cfg$figures_root, basename(cfg$latest_file)))
  cat("Wald合并图目录：", figure_dir, "\n", sep = "")
  invisible(figure_dir)
}
