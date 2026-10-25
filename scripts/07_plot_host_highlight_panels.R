# One host highlighted per panel, using SAVED NMDS coordinates only.
# Run from top to bottom in RStudio. Default layout: 3 columns; row count follows the number of types.
# No distance calculation, NMDS fitting, PERMANOVA, or BH recalculation.
(function() {
  # ---- Settings --------------------------------------------------------------
  project_dir <- "C:/Users/xzy/Desktop/CUHK_Study/GNBF6010_Research_Projects/Research/Msc main project-time series"
  latest_pointer <- file.path(project_dir,"results","00_latest_PERMANOVA_MC8.txt")
  if (!file.exists(latest_pointer) && is.null(getOption("mc.results_root"))) stop("Run run_permanova.R first.")
  results_root <- getOption("mc.results_root", if(file.exists(latest_pointer))
    readLines(latest_pointer,n=1L,warn=FALSE,encoding="UTF-8") else "")
  analysis_scope <- "all_hosts" # "coral" gives the three coral hosts in one row.
  n_columns <- 3L
  show_ellipses <- FALSE        # Default: highlight points only, with no ellipses.
  ellipse_fill_alpha <- 0.07
  ellipse_linewidth <- 0.4
  ellipse_line_alpha <- 0.65
  save_pdf <- TRUE
  dpi <- 300
  host_order <- c("Acropora sp.", "Montipora sp.", "Platygyra sp.", "CCA",
    "Diadema sp.", "Drupella sp.", "Haliclona sp.", "Holothuria sp.",
    "Perna sp.", "Sea anemone", "Shrimp", "Trochus sp.", "Seawater", "Sediment")
  output_dir <- file.path(project_dir, "figures", paste0("host_highlight_",
    analysis_scope, "_", format(Sys.time(), "%Y%m%d_%H%M%S")))

  stopifnot(analysis_scope %in% c("all_hosts", "coral"), n_columns >= 1L)
  required <- c("ggplot2", "viridisLite")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly=TRUE)]
  if (length(missing)) stop("Required packages: ", paste(missing, collapse=", "))
  suppressPackageStartupMessages(library(ggplot2))
  source(file.path(project_dir, "scripts", "permanova_plot_style.R"), local=TRUE)

  # ---- Read and match the saved coordinates to the completed analysis --------
  result_dir <- file.path(results_root, analysis_scope)
  coordinate_file <- file.path(result_dir, "plot_ordinations.rds")
  analysis_file <- file.path(result_dir, "host_time_analysis.rds")
  if (!file.exists(coordinate_file)) stop(
    "Saved coordinates are missing: ", coordinate_file,
    ". Run 04_plot_permanova.R for this scope first; this script never refits NMDS.")
  if (!file.exists(analysis_file)) stop("Missing completed analysis: ", analysis_file)
  ordinations <- readRDS(coordinate_file)
  analysis <- readRDS(analysis_file)
  data <- as.data.frame(ordinations$host$data)
  original <- as.data.frame(analysis$data)
  if (analysis_scope=="all_hosts" && !all(c("Seawater","Sediment") %in% original$Host_type))
    stop("These saved results exclude environmental samples; run the updated analysis first.")
  selected_mc <- c("MC6", "MC507", "MC21", "MC520", "MC10", "MC15", "MC0", "MC25")
  needed <- c("Repeat_group_ID", "Host_type", "TS", "NMDS1", "NMDS2", selected_mc)
  if (!all(needed %in% names(data)) || !nrow(data)) stop("Incomplete saved coordinates.")
  if (!all(setdiff(needed, c("NMDS1", "NMDS2")) %in% names(original)))
    stop("Incomplete saved analysis data.")
  if (anyNA(data[, needed]) || anyDuplicated(data$Repeat_group_ID) ||
      anyDuplicated(original$Repeat_group_ID) ||
      !setequal(as.character(data$Repeat_group_ID), as.character(original$Repeat_group_ID)))
    stop("Coordinate IDs and analysis IDs must match one-to-one.")
  original <- original[match(data$Repeat_group_ID, original$Repeat_group_ID), , drop=FALSE]
  for (column in c("Host_type", "TS", selected_mc)) {
    same <- if (column %in% selected_mc)
      isTRUE(all.equal(data[[column]], original[[column]], tolerance=1e-12)) else
      identical(as.character(data[[column]]), as.character(original[[column]]))
    if (!same) stop("Saved coordinates belong to different data: ", column)
  }
  result <- analysis$results[analysis$results$Factor == "Host_type", , drop=FALSE]
  if (nrow(result) != 1L || result$Status != "ok" || result$Scope != analysis_scope ||
      !identical(ordinations$MC, selected_mc) ||
      !setequal(as.character(ordinations$input_md5), as.character(result$Input_MD5)))
    stop("Saved result scope, MC set, or source fingerprint does not match.")
  stress <- ordinations$host$fit$stress
  if (length(stress) != 1L || !is.finite(stress) ||
      !all(is.finite(as.matrix(data[, c("NMDS1", "NMDS2")]))))
    stop("Non-finite saved NMDS coordinates or stress.")

  groups <- unique(as.character(data$Host_type))
  groups <- c(host_order[host_order %in% groups], sort(setdiff(groups, host_order)))
  palette <- pv_subset(pv_hosts, groups)
  if ("Trochus sp." %in% groups) palette["Trochus sp."] <- "#111111"
  counts <- table(factor(data$Host_type, levels=groups))
  strip_labels <- setNames(sprintf("%s (n = %d)", groups, as.integer(counts)), groups)
  foreground <- data
  foreground$Panel <- factor(as.character(data$Host_type), levels=groups)
  # Each panel has every OTHER host in grey, drawn below its highlighted host.
  background <- do.call(rbind, lapply(groups, function(host) {
    d <- data[as.character(data$Host_type) != host, c("NMDS1", "NMDS2"), drop=FALSE]
    d$Panel <- factor(rep(host, nrow(d)), levels=groups)
    d
  }))
  padded_range <- function(x) {
    limits <- range(x)
    padding <- diff(limits) * 0.05
    if (padding == 0) padding <- 0.05
    limits + c(-padding, padding)
  }
  base_font <- if (isTRUE(capabilities("cairo"))) "Arial" else "sans"
  caption <- sprintf("Stress = %.3f; n = %d; R\u00b2 = %.3f; BH p = %s",
    stress, nrow(data), result$R2_raw, formatC(result$P_BH_raw, format="g", digits=3))

  # ---- Build the panels: fixed coordinates, fixed scales, equal aspect --------
  # Compute ellipse paths once, then include their full extent in every panel.
  ellipse_paths <- data.frame(NMDS1=numeric(), NMDS2=numeric(),
    Host_type=character(), Panel=factor(levels=groups))
  skipped <- character()
  if (show_ellipses) {
    split_rows <- split(seq_len(nrow(foreground)), foreground$Panel)
    eligible <- vapply(split_rows, function(index) {
      xy <- as.matrix(foreground[index, c("NMDS1", "NMDS2"), drop=FALSE])
      nrow(xy) >= 4L && qr(scale(xy, center=TRUE, scale=FALSE))$rank == 2L
    }, logical(1))
    skipped <- names(split_rows)[!eligible]
    paths <- lapply(names(split_rows)[eligible], function(host) {
      d <- foreground[split_rows[[host]], , drop=FALSE]
      ellipse_plot <- ggplot(d, aes(NMDS1, NMDS2)) +
        stat_ellipse(type="norm", level=0.95, segments=180)
      path <- ggplot_build(ellipse_plot)$data[[1L]]
      if (!nrow(path) || !all(is.finite(as.matrix(path[, c("x", "y")]))))
        stop("Could not compute a finite ellipse for: ", host)
      data.frame(NMDS1=path$x, NMDS2=path$y, Host_type=host,
        Panel=factor(rep(host, nrow(path)), levels=groups))
    })
    if (length(paths)) ellipse_paths <- do.call(rbind, paths)
  }
  x_limits <- padded_range(c(data$NMDS1, ellipse_paths$NMDS1))
  y_limits <- padded_range(c(data$NMDS2, ellipse_paths$NMDS2))
  plot <- ggplot()
  if (nrow(ellipse_paths)) plot <- plot +
    geom_polygon(data=ellipse_paths, aes(NMDS1, NMDS2, fill=Host_type, group=Panel),
      colour=NA, alpha=ellipse_fill_alpha, show.legend=FALSE)
  plot <- plot + geom_point(data=background, aes(NMDS1, NMDS2),
    shape=16, colour="#D5D5D5", size=1.35, alpha=0.65)
  if (nrow(ellipse_paths)) plot <- plot +
    geom_path(data=ellipse_paths, aes(NMDS1, NMDS2, colour=Host_type, group=Panel),
      linewidth=ellipse_linewidth, alpha=ellipse_line_alpha, show.legend=FALSE)
  plot <- plot +
    geom_point(data=foreground, aes(NMDS1, NMDS2, fill=Host_type),
      shape=21, colour="white", size=3.0, stroke=0.25, alpha=1) +
    scale_fill_manual(values=palette, guide="none") +
    scale_colour_manual(values=palette, guide="none") +
    facet_wrap(~Panel, ncol=n_columns, scales="fixed",
      labeller=as_labeller(strip_labels), drop=FALSE) +
    coord_fixed(ratio=1, xlim=x_limits, ylim=y_limits, expand=FALSE) +
    labs(title=NULL, subtitle=NULL, x="NMDS1", y="NMDS2", caption=caption) +
    theme_classic(base_size=12, base_family=base_font) +
    theme(axis.text=element_text(colour="black", size=10),
      axis.line=element_line(colour="black", linewidth=0.4),
      axis.ticks=element_line(colour="black", linewidth=0.35),
      strip.background=element_blank(),
      strip.text=element_text(colour="black", size=12, margin=margin(b=7)),
      panel.spacing=grid::unit(7, "mm"), legend.position="none",
      plot.title=element_blank(), plot.subtitle=element_blank(),
      plot.caption=element_text(hjust=0, size=11, margin=margin(t=12)),
      plot.background=element_rect(fill="white", colour=NA),
      plot.margin=margin(12,14,12,12))

  # ---- Save when the user runs this script ------------------------------------
  if (dir.exists(output_dir)) stop("Output directory already exists: ", output_dir)
  if (!dir.create(output_dir, recursive=TRUE)) stop("Could not create output directory.")
  n_rows <- ceiling(length(groups)/n_columns)
  figure_width <- 12.5
  figure_height <- max(4.0, 3.4*n_rows+0.7)
  stem <- paste0("host_highlight_", analysis_scope, "_", n_columns, "x", n_rows)
  png_device <- if (requireNamespace("ragg", quietly=TRUE)) ragg::agg_png else "png"
  ggsave(file.path(output_dir, paste0(stem, ".png")), plot, device=png_device,
    width=figure_width, height=figure_height, units="in", dpi=dpi, bg="white")
  if (save_pdf) {
    device <- if (isTRUE(capabilities("cairo"))) grDevices::cairo_pdf else grDevices::pdf
    ggsave(file.path(output_dir, paste0(stem, ".pdf")), plot, device=device,
      width=figure_width, height=figure_height, units="in", bg="white")
  }
  write.csv(data[, c("Repeat_group_ID", "Host_type", "TS", "NMDS1", "NMDS2")],
    file.path(output_dir, "saved_coordinates_used.csv"), row.names=FALSE, fileEncoding="UTF-8")
  full_caption <- c(
    "Host distributions in the shared NMDS ordination of eight selected MCs.",
    "Each panel highlights one host; grey points show all other hosts. One point represents one Repeat_group mean.",
    "All panels reuse exactly the same saved NMDS coordinates, axis limits and equal aspect ratio. The full point range is retained.",
    "Coordinates were previously calculated from Bray-Curtis distances after square-root transformation of the MC proportions.",
    "Host PERMANOVA adjusts for TS with permutations within TS. NMDS coordinates are not host-adjusted.",
    "The figure-wide R-squared and BH-adjusted P refer to the overall host test, not individual panels. Host and TS omnibus tests form the two-test BH family.",
    sprintf("PERMDISP P (raw / Lingoes, unadjusted): %.4g / %.4g.",
      result$P_disp_raw, result$P_disp_lingoes),
    if (any(c(result$P_disp_raw, result$P_disp_lingoes) <= 0.05, na.rm=TRUE))
      "Within-host dispersion differences were detected; interpret the host composition test together with this diagnostic." else
      "The dispersion diagnostic is reported with the saved analysis.",
    if (show_ellipses) paste0("Only the highlighted host has a 95% data ellipse (normal approximation). ",
      "The matching-colour fill has alpha = ", ellipse_fill_alpha,
      "; the thin outline has linewidth = ", ellipse_linewidth,
      " and alpha = ", ellipse_line_alpha, ". These are not confidence regions for group centres. ",
      "Shared axis limits include all points and all computed ellipses.") else "Ellipses are omitted.",
    if (length(skipped)) paste("No ellipse (<4 points or collinear):", paste(skipped, collapse=", ")),
    paste("Saved ordination:", coordinate_file), paste("Saved analysis:", analysis_file),
    "This script performs saved-result plotting only; it does not refit or rerun statistical tests.")
  writeLines(enc2utf8(full_caption), file.path(output_dir, "figure_caption.txt"), useBytes=TRUE)
  write.csv(data.frame(Source=c(coordinate_file, analysis_file),
    MD5=unname(tools::md5sum(c(coordinate_file, analysis_file)))),
    file.path(output_dir, "figure_sources.csv"), row.names=FALSE)
  message("Host highlight panels saved in: ", output_dir)
})()
