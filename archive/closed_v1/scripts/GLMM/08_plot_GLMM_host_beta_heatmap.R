# Historical MC8 only. Current MC19 pilot figures: plot_GLMM_MC19.R.
# Figure A: host coefficients from the saved Z_TS_SD_Q MC8 model.
# Saved-result plotting only: no refitting, simulation or BH recalculation.
# English labels; Arial; white background; no title/subtitle.
# Match 06_plot_pairwise_permanova.R: high-contrast viridis, separate beta/q/stars.
# Run this script yourself when ready. Output goes directly under figures.

# ---- 1. Settings ------------------------------------------------------------
project_dir <- "C:/Users/xzy/Desktop/CUHK_Study/GNBF6010_Research_Projects/Research/Msc main project-time series"
results_dir <- file.path(project_dir, "results", "GLMM_Z_TS_SD_Q_MC8_20261013_185720")
coefficient_file <- file.path(results_dir, "00_beta_significance.csv")
checks_file <- file.path(results_dir, "00_model_checks.csv")
style_file <- file.path(project_dir, "scripts", "PERMANOVA", "permanova_plot_style.R")
output_dir <- file.path(project_dir, "figures",
  paste0("GLMM_host_beta_heatmap_", format(Sys.time(), "%Y%m%d_%H%M%S")))

mc_order <- c("MC0", "MC6", "MC10", "MC15", "MC21", "MC25", "MC507", "MC520")
reference_host <- "Acropora sp."
# Keep every GLMM category, including categories absent from the PERMANOVA subset.
host_order <- c("Acropora sp.", "Montipora sp.", "Platygyra sp.", "CCA",
  "Diadema sp.", "Drupella sp.", "Haliclona sp.", "Holothuria sp.",
  "Perna sp.", "Sea anemone", "Shrimp", "Trochus sp.", "Cuttlefish",
  "Sargassum", "Seawater", "Sediment")

font_family <- "Arial"
dpi <- 300
save_pdf <- FALSE                 # TRUE adds a vector PDF (requires Cairo).
show_q_values <- TRUE             # beta in center; q on lower line; stars above.
show_caption <- FALSE             # Full caption stays in A_figure_caption.txt; more space for cells.
panel_tag <- NULL                 # No top-left panel letter.
beta_limits <- NULL               # Automatic symmetric range; or e.g. c(-2.5,2.5).
figure_width <- 24                # Keep the canvas unchanged in this font revision.
figure_height <- 15
# geom_text sizes are millimetres; theme sizes are points.
beta_text_size <- 9.0
q_text_size <- 7.0
star_text_size <- 10.0
axis_text_size <- 30

plot_GLMM_host_heatmap <- function() {
  for (pkg in c("ggplot2", "viridisLite", "scales")) {
    if (!requireNamespace(pkg, quietly = TRUE)) stop("Required package: ", pkg)
  }
  if (!file.exists(style_file)) stop("Missing shared display settings: ", style_file)
  # Source display definitions only; this file performs no statistical analysis.
  display <- new.env(parent = globalenv())
  sys.source(style_file, envir = display, keep.source = FALSE)
  fill_colours <- display$pv_continuous
  if (save_pdf && !capabilities("cairo")) {
    stop("Cairo is needed for Arial PDF export. Set save_pdf <- FALSE for PNG only.")
  }
  if (anyDuplicated(mc_order) || anyDuplicated(host_order)) stop("Duplicate axis labels.")
  if (length(beta_limits) && (length(beta_limits) != 2L ||
      any(!is.finite(beta_limits)) || beta_limits[1] >= 0 ||
      beta_limits[2] <= 0 || abs(sum(beta_limits)) > 1e-10)) {
    stop("beta_limits must be symmetric about zero, e.g. c(-2.5, 2.5).")
  }
  read_saved <- function(path) {
    if (!file.exists(path)) stop("Missing saved result: ", path)
    read.csv(path, stringsAsFactors = FALSE, check.names = FALSE,
      fileEncoding = "UTF-8-BOM", na.strings = c("NA", ""))
  }
  source_hashes <- tools::md5sum(c(coefficient_file, checks_file, style_file))
  coefficients <- read_saved(coefficient_file)
  checks <- read_saved(checks_file)
  required <- c("MC_ID", "term", "estimate", "SE", "CI95_low", "CI95_high",
    "p_approx", "q_BH", "BH_family", "inference_status", "report_stage")
  if (length(setdiff(required, names(coefficients)))) stop("Missing coefficient columns.")
  if (!all(c("MC_ID", "numerical_valid") %in% names(checks))) stop("Missing check columns.")
  if (anyDuplicated(checks$MC_ID) || !setequal(checks$MC_ID, mc_order)) {
    stop("The check table must contain the planned eight MCs exactly once.")
  }
  d <- coefficients[startsWith(coefficients$term, "Host_type"), , drop = FALSE]
  if (!nrow(d) || !setequal(unique(d$MC_ID), mc_order)) stop("Incomplete MC identifiers.")
  d$Host <- substring(d$term, nchar("Host_type") + 1L)
  if (anyNA(d$Host) || any(!nzchar(d$Host)) || reference_host %in% d$Host ||
      anyDuplicated(d[c("MC_ID", "Host")] )) stop("Invalid or duplicate host contrasts.")
  if (!all(d$BH_family == "Host_contrasts_across_MC8") ||
      !all(d$report_stage == "final_available_results")) {
    stop("Expected final saved host-contrast inference with the original BH family.")
  }
  numeric_cols <- c("estimate", "SE", "CI95_low", "CI95_high", "p_approx", "q_BH")
  if (!all(vapply(d[numeric_cols], is.numeric, logical(1)))) stop("Non-numeric inference columns.")
  for (field in c("p_approx", "q_BH")) {
    x <- d[[field]]
    if (any(!is.na(x) & (!is.finite(x) | x < 0 | x > 1))) stop("Invalid ", field)
  }
  hosts <- unique(d$Host)
  expected <- expand.grid(MC_ID = mc_order, Host = hosts, stringsAsFactors = FALSE)
  # Missing inference remains a row with NA fields; missing planned rows are errors.
  if (nrow(d) != nrow(expected)) stop("Saved host contrasts are not a complete planned grid.")
  groups <- c(host_order[host_order %in% c(reference_host, hosts)],
    sort(setdiff(c(reference_host, hosts), host_order)))
  if (!reference_host %in% groups) groups <- c(reference_host, groups)
  family_size <- nrow(expected)
  numeric_valid <- toupper(as.character(checks$numerical_valid)) == "TRUE"
  numeric_valid[is.na(numeric_valid)] <- FALSE
  valid_map <- stats::setNames(numeric_valid, checks$MC_ID)
  d$Available <- unname(valid_map[d$MC_ID]) &
    !is.na(d$inference_status) & d$inference_status == "ok_approximate" &
    is.finite(d$estimate) & is.finite(d$SE) & d$SE > 0 &
    is.finite(d$CI95_low) & is.finite(d$CI95_high) & is.finite(d$q_BH)
  if (any(d$Available & (d$CI95_low > d$estimate | d$CI95_high < d$estimate))) {
    stop("A saved confidence interval does not contain its estimate.")
  }
  d$Cell_status <- ifelse(d$Available, "estimated", "unavailable")
  # A reference coefficient is fixed by model coding; it is not a fitted test.
  ref <- d[rep(1L, length(mc_order)), , drop = FALSE]
  ref$MC_ID <- mc_order; ref$Host <- reference_host; ref$term <- "(reference)"
  ref[numeric_cols] <- NA_real_; ref$Available <- FALSE
  ref$Cell_status <- ifelse(unname(valid_map[mc_order]), "reference", "unavailable")
  ref$inference_status <- "reference_not_tested"
  d <- rbind(d, ref)
  valid <- d[d$Cell_status == "estimated", , drop = FALSE]
  unavailable <- d[d$Cell_status == "unavailable", , drop = FALSE]
  reference <- d[d$Cell_status == "reference", , drop = FALSE]
  if (!nrow(valid)) stop("There are no valid saved host coefficients to plot.")
  limits <- beta_limits
  if (is.null(limits)) {
    bound <- max(.5, ceiling(max(abs(valid$estimate)) * 2) / 2)
    limits <- c(-bound, bound)
  }
  if (any(valid$estimate < limits[1] | valid$estimate > limits[2])) {
    stop("beta_limits would clip estimated coefficients; enlarge the limits.")
  }
  stars <- function(q) ifelse(q < .001, "***", ifelse(q < .01, "**", ifelse(q < .05, "*", "")))
  format_q <- function(q) vapply(q, function(x) {
    if (x < .001) formatC(x, format = "e", digits = 1) else formatC(x, format = "f", digits = 3)
  }, character(1))
  valid$Beta_label <- sprintf("%.3f", valid$estimate)
  valid$Q_label <- paste0("q=", format_q(valid$q_BH))
  valid$Stars_label <- stars(valid$q_BH)
  valid$Text_colour <- display$pv_text_colour(valid$estimate, limits, fill_colours)
  failed_mc <- names(valid_map)[!valid_map]
  caption <- paste0(
    "Fill/center: host coefficient beta (logit scale), relative to Acropora sp.; lower line: saved BH q.\n",
    "Stars above beta: * q < 0.05; ** q < 0.01; *** q < 0.001. BH family: ", family_size,
    " planned host contrasts across eight MCs.\n",
    "Positive beta: higher count-component composition odds than Acropora; negative beta: lower odds.\n",
    "Adjusted for DO, reads spline and Repeat/TS effects; approximate penalized-Hessian Wald inference.\n",
    "Ref.: fixed reference, not tested. NA: numerical fit or coefficient inference unavailable.")
  if (length(failed_mc)) caption <- paste0(caption, " Numerical checks failed: ", paste(failed_mc, collapse = ", "), ".")
  suppressPackageStartupMessages(library(ggplot2))
  plot <- ggplot() +
    geom_tile(data = reference, aes(x = Host, y = MC_ID), fill = "#F1F1F1",
      colour = "white", linewidth = .55) +
    geom_tile(data = unavailable, aes(x = Host, y = MC_ID), fill = "#DADADA",
      colour = "white", linewidth = .55) +
    geom_tile(data = valid, aes(x = Host, y = MC_ID, fill = estimate),
      colour = "white", linewidth = .55) +
    geom_text(data = valid, aes(x = Host, y = MC_ID, label = Beta_label, colour = Text_colour),
      position = position_nudge(y = 0), family = font_family, size = beta_text_size, show.legend = FALSE) +
    geom_text(data = valid, aes(x = Host, y = MC_ID, label = Stars_label, colour = Text_colour),
      position = position_nudge(y = .32), family = font_family, size = star_text_size, fontface = "bold", show.legend = FALSE) +
    geom_text(data = reference, aes(x = Host, y = MC_ID), label = "Ref.",
      family = font_family, size = beta_text_size, colour = "#444444") +
    geom_text(data = unavailable, aes(x = Host, y = MC_ID), label = "NA",
      family = font_family, size = beta_text_size, colour = "#444444") +
    scale_colour_identity() +
    scale_fill_gradientn(colours = fill_colours, limits = limits,
      breaks = seq(limits[1], limits[2], length.out = 5L),
      labels = function(x) sprintf("%.2f", x), name = expression(Host~beta),
      guide = guide_colourbar(title.position = "top", barwidth = grid::unit(90, "mm"),
        barheight = grid::unit(7, "mm"), ticks.colour = "black", frame.colour = "black")) +
    scale_x_discrete(limits = groups, drop = FALSE, expand = expansion(add = .5)) +
    scale_y_discrete(limits = rev(mc_order), drop = FALSE, expand = expansion(add = .5)) +
    coord_fixed(clip = "off") +
    labs(title = NULL, subtitle = NULL, x = "Host / sample category", y = "MC",
      tag = panel_tag, caption = if (show_caption) caption else NULL) +
    theme_classic(base_size = 16, base_family = font_family) +
    theme(axis.line = element_line(colour = "black", linewidth = .45),
      axis.ticks = element_line(colour = "black", linewidth = .4),
      axis.ticks.length = grid::unit(2, "mm"), axis.text = element_text(colour = "black", size = axis_text_size),
      axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
      axis.title = element_text(colour = "black", size = 30),
      axis.title.x = element_text(margin = margin(t = 10)), axis.title.y = element_text(margin = margin(r = 10)),
      plot.title = element_blank(), plot.subtitle = element_blank(),
      plot.tag = element_blank(),
      plot.caption = element_text(size = 13, colour = "#333333", hjust = 0,
        lineheight = 1.15, margin = margin(t = 16)), plot.caption.position = "plot",
      legend.position = "top", legend.justification = "right",
      legend.title = element_text(size = 24), legend.text = element_text(size = 22),
      legend.margin = margin(t = 2, b = 8), panel.grid = element_blank(),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background = element_rect(fill = "white", colour = NA),
      plot.margin = margin(t = 12, r = 20, b = 18, l = 18))
  if (show_q_values) plot <- plot +
    geom_text(data = valid, aes(x = Host, y = MC_ID, label = Q_label, colour = Text_colour),
      position = position_nudge(y = -.32), family = font_family, size = q_text_size, show.legend = FALSE)
  if (!identical(unname(tools::md5sum(names(source_hashes))), unname(source_hashes))) {
    stop("Saved results or display settings changed while preparing this figure.")
  }
  if (dir.exists(output_dir)) stop("Output directory exists; choose a fresh output_dir.")
  if (!dir.create(output_dir, recursive = TRUE)) stop("Could not create output directory.")
  png_device <- if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_png else "png"
  ggsave(file.path(output_dir, "A_GLMM_host_beta_heatmap.png"), plot,
    device = png_device, width = figure_width, height = figure_height, units = "in", dpi = dpi, bg = "white")
  if (save_pdf) ggsave(file.path(output_dir, "A_GLMM_host_beta_heatmap.pdf"), plot,
    device = grDevices::cairo_pdf, width = figure_width, height = figure_height, units = "in", bg = "white")
  # Export labels and exact saved statistics so every plotted cell is traceable.
  audit <- d
  audit$Stars <- ifelse(audit$Available, stars(audit$q_BH), "")
  write.csv(audit, file.path(output_dir, "A_heatmap_cell_values.csv"), row.names = FALSE, fileEncoding = "UTF-8")
  write.csv(data.frame(Source = names(source_hashes), MD5 = unname(source_hashes)),
    file.path(output_dir, "figure_sources.csv"), row.names = FALSE, fileEncoding = "UTF-8")
  writeLines(caption, file.path(output_dir, "A_figure_caption.txt"), useBytes = TRUE)
  saveRDS(list(plot = plot, cells = audit, beta_limits = limits, host_order = groups, MC_order = mc_order),
    file.path(output_dir, "A_plot_objects.rds"))
  message("Saved: ", output_dir, "\nSignificant host contrasts: ", sum(valid$q_BH < .05),
    "; estimated: ", nrow(valid), "; planned: ", family_size)
  invisible(list(plot = plot, output_dir = output_dir, cells = audit))
}

# Source/Run from top to bottom to generate Figure A. Definition-only loading is
# available to inspect code without reading results or generating figures.
if (!isTRUE(getOption("GLMM_HOST_HEATMAP_definitions_only", FALSE))) {
  GLMM_host_beta_heatmap <- plot_GLMM_host_heatmap()
}
