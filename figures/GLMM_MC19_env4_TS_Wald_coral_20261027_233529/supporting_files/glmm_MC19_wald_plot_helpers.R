# Saved-result plots only. No model fitting, response simulation or P/q recalculation.

mc19_wald_note <- function(title, lines) {
  graphics::par(oma = c(0, 0, 0, 0))
  mc19_note_page(title, unlist(lapply(lines, strwrap, width = 105), use.names = FALSE))
}

# Main figures use one flat canvas: no page title, all labels in English.
mc19_main_heatmap <- function(results, cfg) {
  graphics::par(mar = rep(0, 4), family = "sans", xpd = NA)
  graphics::plot.new()
  graphics::plot.window(c(0, 1), c(0, 1), xaxs = "i", yaxs = "i")
  beta <- matrix(NA_real_, length(cfg$mc_ids), length(cfg$environments))
  q <- beta
  for (i in seq_along(cfg$mc_ids)) for (j in seq_along(cfg$environments)) {
    r <- results[results$MC_ID == cfg$mc_ids[i] & results$environment == cfg$environments[j], ]
    if (isTRUE(r$numerical_ok)) beta[i, j] <- r$beta_per_TS_SD
    q[i, j] <- r$beta_q_value
  }
  limit <- max(abs(beta), na.rm = TRUE)
  if (!is.finite(limit) || limit == 0) limit <- 1
  palette <- grDevices::colorRampPalette(c("#083DCC", "#FFFFFF", "#D00000"))(201)
  colour <- function(v) if (!is.finite(v)) "#D9D9D9" else
    palette[1L + round(200 * (v + limit) / (2 * limit))]
  left <- .20; right <- .71; bottom <- .09; top <- .96
  w <- (right - left) / 4; h <- (top - bottom) / 19
  for (i in seq_len(19)) for (j in seq_len(4)) {
    x <- left + (j - 1) * w; y <- top - i * h
    graphics::rect(x, y, x + w, y + h, col = colour(beta[i, j]), border = "white", lwd = .8)
    mark <- if (!is.finite(beta[i, j])) "NA" else
      if (is.finite(q[i, j]) && q[i, j] < .05) "*" else ""
    graphics::text(x + w / 2, y + h / 2, mark, cex = 1.15,
      col = if (is.finite(beta[i, j]) && abs(beta[i, j]) > .65 * limit) "white" else "black")
  }
  graphics::text(left - .02, top - (seq_len(19) - .5) * h,
    sub("MC2s1", "MC2S1", cfg$mc_ids, fixed = TRUE), adj = 1, cex = .88)
  graphics::text(left + (seq_len(4) - .5) * w, bottom - .024,
    expression(plain(NO)[3]*"-N", plain(NO)[2]*"-N", plain(NH)[4]*"-N", plain(DIP)), cex = .95)
  breaks <- seq(-limit, limit, length.out = 202)
  for (k in 1:201) graphics::rect(.80, .28 + (k - 1) / 201 * .43,
    .835, .28 + k / 201 * .43, col = palette[k], border = NA)
  ticks <- c(-limit, 0, limit)
  graphics::text(.85, .28 + (ticks + limit) / (2 * limit) * .43,
    formatC(ticks, format = "f", digits = 2), adj = 0, cex = .8)
  graphics::text(.8175, .75, "Beta", cex = .92)
  graphics::text(.8175, .727, "per nutrient SD", cex = .65)
  graphics::text(.5, .025, "* BH q < 0.05 (76 tests)", cex = .8)
}

# One panel per significant MC x nutrient combination, ordered by the analysis plan.
mc19_significant_curve_rows <- function(results, cfg) {
  selected <- results[results$model == "BB_TS" & results$numerical_ok %in% TRUE &
    is.finite(results$beta_q_value) & results$beta_q_value < .05, , drop = FALSE]
  selected <- selected[order(match(selected$MC_ID, cfg$mc_ids),
    match(selected$environment, cfg$environments)), , drop = FALSE]
  selected$key <- mc19_wald_key(selected$MC_ID, selected$environment, selected$model)
  if (anyDuplicated(selected$key) || anyNA(match(selected$MC_ID, cfg$mc_ids)) ||
      anyNA(match(selected$environment, cfg$environments))) stop("Invalid significant curve keys.")
  # One MC per row; all its significant nutrients stay on the same page.
  mc_row <- match(selected$MC_ID, unique(selected$MC_ID))
  selected$page <- (mc_row - 1L) %/% 3L + 1L
  selected$row <- (mc_row - 1L) %% 3L + 1L
  selected$column <- as.integer(ave(seq_len(nrow(selected)), selected$MC_ID, FUN=seq_along))
  n_cols <- max(2L, as.integer(table(selected$MC_ID)))
  selected$panel <- (selected$row - 1L) * n_cols + selected$column
  rownames(selected) <- NULL
  selected
}

# Dimensions in inches: fixed axis/row/footer clearances, independent of page count.
mc19_curve_page_size <- function(selected) {
  n_cols <- max(2L, as.integer(table(selected$MC_ID)))
  rows <- if (nrow(selected)) max(selected$row) else 1L
  list(n_cols=n_cols, rows=rows, width=7.8*n_cols+.4,
    height=.8 + 2.8*rows + 1.0*(rows-1L) + 1.4)
}

mc19_curve_limits <- function(all_panels, entries, grids) {
  answer <- lapply(unique(all_panels$MC_ID), function(mc) {
    keys <- all_panels$key[all_panels$MC_ID == mc]
    full <- max(vapply(keys, function(k)
      max(100*c(entries[[k]]$data$MC_proportion, grids[[k]]$Wald95_high_pilot)), numeric(1)))
    focus <- max(vapply(keys, function(k) max(
      stats::quantile(100*entries[[k]]$data$MC_proportion, .85),
      100*grids[[k]]$Wald95_high_pilot), numeric(1)))
    zoomed <- is.finite(focus) && focus > 0 && full > 2.5*focus
    list(ymax=max(full,.001)*1.08, full=max(full,.001)*1.08,
      zoom=max(focus,.001)*1.12, zoomed=zoomed)
  })
  names(answer) <- unique(all_panels$MC_ID)
  answer
}

# A dedicated right-hand zoom area; even its white frame stays outside the main plot.
mc19_curve_panel_boxes <- function(row, column, width, height) {
  ytop <- 1-(.8+(row-1L)*3.8)/height
  list(main=c((.8+(column-1L)*7.8)/width,(5.05+(column-1L)*7.8)/width,
      ytop-2.8/height,ytop),
    zoom=c((5.68+(column-1L)*7.8)/width,(7.65+(column-1L)*7.8)/width,
      ytop-1.45/height,ytop-.10/height))
}

mc19_main_curves <- function(results, entries, grids, cfg, environment_table, panels,
    page_height=NULL, page_width=NULL) {
  all_panels <- mc19_significant_curve_rows(results,cfg)
  size <- mc19_curve_page_size(all_panels)
  width <- if (is.null(page_width)) size$width else page_width
  height <- if (is.null(page_height)) mc19_curve_page_size(panels)$height else page_height
  graphics::par(mar=rep(0,4),family="sans",xpd=NA)
  graphics::plot.new(); graphics::plot.window(c(0,1),c(0,1),xaxs="i",yaxs="i")
  if (!nrow(panels)) {
    graphics::text(.5,.5,"No associations with BH q < 0.05.",cex=.9)
  } else {
    ts_levels <- as.character(environment_table$TS)
    colours <- c("#247BA0","#39A59C","#6B9A45","#D89C31","#D26A56",
      "#8E6BB1","#BD6E99","#555555")[seq_along(ts_levels)]
    graphics::legend(.5,1-.08/height,legend=ts_levels,xjust=.5,yjust=1,
      horiz=TRUE,bty="n",pch=21,pt.bg=colours,col="white",cex=.87)
    limits <- mc19_curve_limits(all_panels,entries,grids)
    panel <- function(r,box,label,zoom=FALSE) {
      d <- entries[[r$key]]$data; g <- grids[[r$key]]
      span <- diff(range(d$E_raw)); xr <- range(d$E_raw)+c(-.04,.04)*span
      limit <- limits[[r$MC_ID]]
      ymax <- if (zoom) limit$zoom else limit$ymax
      px <- function(x) box[1]+(x-xr[1])/diff(xr)*diff(box[1:2])
      py <- function(y) box[3]+y/ymax*diff(box[3:4])
      if (zoom) graphics::rect(box[1]-.30/width,box[3]-.23/height,
        box[2]+.04/width,box[4]+.24/height,col="white",border="#B9C0C7",lwd=.7)
      graphics::clip(box[1],box[2],box[3],box[4])
      yt <- pretty(c(0,ymax),if(zoom) 2 else 4); yt <- yt[yt>=0 & yt<=ymax]
      graphics::segments(box[1],py(yt),box[2],py(yt),col="#ECEFF1",lwd=.7)
      graphics::polygon(px(c(g$E_raw,rev(g$E_raw))),
        py(100*c(g$Wald95_low_pilot,rev(g$Wald95_high_pilot))),border=NA,col="#D9DDE1")
      offsets <- ave(seq_len(nrow(d)),d$TS,FUN=function(v) seq(-.009,.009,length.out=length(v)))
      ts_colour <- colours[match(as.character(d$TS),ts_levels)]
      censored <- as.character(d$TS) %in% environment_table$TS[
        environment_table[[paste0(r$environment,"_censored")]] %in% TRUE]
      graphics::points(px(d$E_raw+offsets*span),py(100*d$MC_proportion),pch=21,
        bg=ifelse(censored,"white",grDevices::adjustcolor(ts_colour,alpha.f=.50)),
        col=ifelse(censored,ts_colour,"white"),lwd=.6,cex=if(zoom) .42 else .80)
      graphics::lines(px(g$E_raw),py(100*g$fixed_part),lwd=if(zoom) 1.5 else 2.5,col="#20272D")
      graphics::clip(0,1,0,1)
      graphics::segments(c(box[1],box[1]),c(box[4],box[3]),
        c(box[1],box[2]),c(box[3],box[3]),col="#8C969D",lwd=.8)
      xt <- pretty(xr,if(zoom) 2 else 4);xt <- xt[xt>=xr[1] & xt<=xr[2]]
      graphics::text(px(xt),box[3]-if(zoom) .10/height else .18/height,
        formatC(xt,format="f",digits=3),cex=if(zoom) .60 else .78)
      graphics::text(box[1]-if(zoom) .04/width else .10/width,py(yt),
        format(yt,trim=TRUE,digits=2),adj=1,cex=if(zoom) .60 else .80)
      if (zoom) {
        graphics::text(mean(box[1:2]),box[4]+.11/height,"Zoom",cex=.65)
      } else {
        graphics::text(box[1],box[4]+.25/height,
          paste(label,sub("MC2s1","MC2S1",r$MC_ID,fixed=TRUE)),adj=0,font=2,cex=1)
        graphics::text(box[2],box[4]+.25/height,
          sprintf("Beta = %.2f; q = %s",r$beta_per_TS_SD,formatC(r$beta_q_value,format="g",digits=2)),
          adj=1,cex=.70)
        graphics::text(box[1]-.58/width,mean(box[3:4]),"MC proportion (%)",srt=90,cex=.88)
        graphics::text(mean(box[1:2]),box[3]-.48/height,
          unname(cfg$axis_labels[r$environment]),cex=.88)
      }
    }
    for (k in seq_len(nrow(panels))) {
      r <- panels[k,,drop=FALSE]
      boxes <- mc19_curve_panel_boxes(r$row,r$column,width,height)
      panel(r,boxes$main,LETTERS[k])
      if (limits[[r$MC_ID]]$zoomed) panel(r,boxes$zoom,"",zoom=TRUE)
    }
  }
  if (nrow(panels)) {
    graphics::text(.5,.62/height,"Line: b_TS = 0; band: 95% CI; open points: below detection (L/2).",cex=.70)
    graphics::text(.5,.39/height,"Main panels: full range; right-hand panels: y-axis zoom (same x-axis).",cex=.70)
  }
  graphics::text(.5,.16/height,sprintf("%s; beta per nutrient SD; BH q (%d tests). n = %d; %d TS.",
    cfg$analysis_scope,cfg$q_family_size,results$n_samples[1],results$n_TS[1]),cex=.70)
}

mc19_export_significant_curves <- function(stem,results,entries,grids,cfg,environment_table,selected) {
  pages <- if(nrow(selected)) split(seq_len(nrow(selected)),selected$page) else list(integer())
  size <- mc19_curve_page_size(selected)
  grDevices::pdf(paste0(stem,".pdf"),width=size$width,height=size$height,
    family="sans",pointsize=16,useDingbats=FALSE)
  tryCatch({
    for(ix in pages) mc19_main_curves(results,entries,grids,cfg,environment_table,
      selected[ix,,drop=FALSE],page_height=size$height,page_width=size$width)
  },finally=grDevices::dev.off())
  for(i in seq_along(pages)) {
    ix <- pages[[i]]; height <- mc19_curve_page_size(selected[ix,,drop=FALSE])$height
    grDevices::png(sprintf("%s_page_%02d.png",stem,i),width=size$width,height=height,
      units="in",res=300,pointsize=16,type="cairo",bg="white")
    tryCatch(mc19_main_curves(results,entries,grids,cfg,environment_table,
      selected[ix,,drop=FALSE],page_height=height,page_width=size$width),finally=grDevices::dev.off())
  }
}

mc19_main_export <- function(stem, draw, width, height) {
  grDevices::pdf(paste0(stem, ".pdf"), width = width, height = height,
    family = "sans", pointsize = 16, useDingbats = FALSE)
  tryCatch(draw(), finally = grDevices::dev.off())
  grDevices::png(paste0(stem, ".png"), width = width, height = height,
    units = "in", res = 300, pointsize = 16, type = "cairo", bg = "white")
  tryCatch(draw(), finally = grDevices::dev.off())
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
  selected <- mc19_significant_curve_rows(run$results, cfg)
  curve_keys <- selected$key
  for (key in curve_keys) {
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
  mc19_write_csv(if (length(notices)) do.call(rbind, notices) else
    data.frame(key = character(), error = character(), warnings = character()),
    file.path(support, "plot_notices.csv"))
  if (!all(curve_keys %in% names(grids))) stop("Significant-model prediction unavailable; inspect plot_notices.csv.")
  mc19_write_csv(selected, file.path(support, "significant_curve_index.csv"))
  axis_limits <- mc19_curve_limits(selected, entries, grids)
  axis_table <- if (length(axis_limits)) do.call(rbind, lapply(names(axis_limits), function(mc) {
    x <- axis_limits[[mc]]
    data.frame(MC_ID=mc, main_ymax_percent=x$ymax, zoom_ymax_percent=x$zoom, zoom_enabled=x$zoomed)
  })) else data.frame(MC_ID=character(),main_ymax_percent=numeric(),
    zoom_ymax_percent=numeric(),zoom_enabled=logical())
  mc19_write_csv(axis_table, file.path(support, "curve_axis_ranges.csv"))
  environment_table <- utils::read.csv(file.path(run$support, "input_environment_TS.csv"), stringsAsFactors = FALSE)
  # Archive replaced curve figures and known legacy main PDFs; retain all originals locally.
  old <- file.path(figure_dir, c("00_summary.pdf", "01_NO3_effects_and_curves.pdf",
    "02_NO2_effects_and_curves.pdf", "03_NH4_effects_and_curves.pdf", "04_DIP_effects_and_curves.pdf",
    "02_MC98_MC80_MC15_nutrient_curves.pdf", "02_MC98_MC80_MC15_nutrient_curves.png"))
  old <- c(old, list.files(figure_dir,
    pattern = "^02_all_significant_nutrient_curves([.]pdf|_page_[0-9]+[.]png)$", full.names = TRUE))
  old <- old[file.exists(old)]
  if (length(old)) {
    archive <- file.path(support, "previous_main_figures", format(Sys.time(), "%Y%m%d_%H%M%S"))
    dir.create(archive, recursive = TRUE, showWarnings = FALSE)
    if (!all(file.rename(old, file.path(archive, basename(old))))) stop("Could not archive previous main PDFs.")
  }
  mc19_main_export(file.path(figure_dir, "01_MC19_nutrient_heatmap"),
    function() mc19_main_heatmap(run$results, cfg), width = 7, height = 8.5)
  mc19_export_significant_curves(file.path(figure_dir, "02_all_significant_nutrient_curves"),
    run$results, entries, grids, cfg, environment_table, selected)
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
  prediction_table <- if (length(grids)) do.call(rbind, grids) else
    data.frame(MC_ID = character(), environment = character(), E_raw = numeric(),
      fixed_part = numeric(), Wald95_low_pilot = numeric(), Wald95_high_pilot = numeric())
  mc19_write_csv(prediction_table, file.path(support, "prediction_curves.csv"))
  mc19_write_csv(if (length(notices)) do.call(rbind, notices) else
    data.frame(key = character(), error = character(), warnings = character()),
    file.path(support, "plot_notices.csv"))
  file.copy(file.path(run$run_dir, "00_result.csv"), file.path(support, "result_snapshot.csv"), overwrite = TRUE)
  for (script in cfg$scripts) file.copy(file.path(cfg$project_dir, "scripts", "GLMM", script),
    file.path(support, script), overwrite = TRUE)
  mc19_write_text(c(paste("Source:", run$run_dir),
    "Main: 01_MC19_nutrient_heatmap (PDF/PNG); 02_all_significant_nutrient_curves.pdf and page PNGs.",
    "Curves: every numerically valid MC x nutrient with saved BH q < 0.05; one MC/row, up to three MCs/page.",
    "Panel/page mapping: supporting_files/significant_curve_index.csv. No q-value recalculation.",
    "Saved diagnostics, predictions, scripts and previous main PDFs: supporting_files/.",
    "Replot only reads saved models, P/q, and diagnostics. No fitting/simulation/P recalculation."),
    file.path(support, "README.txt"))
  mc19_write_text(figure_dir, file.path(cfg$figures_root, basename(cfg$latest_file)))
  cat("Wald合并图目录：", figure_dir, "\n", sep = "")
  invisible(figure_dir)
}
