# 当前MC19 coral方案停用本历史入口；原代码留作查阅，原文件另有备份。
stop("本轮没有两两检验结果，不绘制旧MC8两两热图。MC逐样本热图在04中输出。请使用run_permanova.R或04_plot_permanova.R。",call.=FALSE)

# Five pairwise PERMANOVA heatmaps, using SAVED results only.
# Run this file from top to bottom when ready to create the figures.
# No model fitting, distance calculation, permutations, NMDS, or BH recalculation.
# All figure text is English. Style: the supplied PM11 seasonal/time-series figures.
# No titles/subtitles. R-squared, q, stars, and daggers have separate cell positions.

# ---- 1. Settings --------------------------------------------------------------
project_dir <- "C:/Users/xzy/Desktop/CUHK_Study/GNBF6010_Research_Projects/Research/Msc main project-time series"
source(file.path(project_dir, "scripts", "PERMANOVA", "permanova_plot_style.R"), local=TRUE)
# Read the latest completed run, or an explicitly selected mc.results_root.
latest_pointer <- file.path(project_dir,"results","00_latest_PERMANOVA_MC8.txt")
  if (!file.exists(latest_pointer) && is.null(getOption("mc.results_root"))) stop("Run run_permanova.R first.")
  results_root <- getOption("mc.results_root", if(file.exists(latest_pointer))
    readLines(latest_pointer,n=1L,warn=FALSE,encoding="UTF-8") else "")
# Save directly beneath figures, with one timestamped folder per plotting run.
output_dir <- getOption("mc.pairwise_figure_dir",file.path(project_dir, "figures",
                        paste0("pairwise_heatmaps_", format(Sys.time(), "%Y%m%d_%H%M%S"))))

font_family <- "Arial"
dpi <- 300
save_pdf <- FALSE             # TRUE also saves a vector PDF for each panel.
show_dispersion_flags <- TRUE # Dagger: saved BH-adjusted PERMDISP P <= 0.05.
r2_limits <- c(0, 0.60)       # Shared scale; the current maximum is about 0.531.
fill_colours <- pv_continuous # Shared high-contrast, linear 0-0.60 colour scale.

# Original host names are retained; coral genera appear first.
host_order <- c("Acropora sp.", "Montipora sp.", "Platygyra sp.", "CCA",
                "Diadema sp.", "Drupella sp.", "Haliclona sp.",
                "Holothuria sp.", "Perna sp.", "Sea anemone", "Shrimp", "Trochus sp.", "Seawater", "Sediment")
ts_order <- paste0("TS", 1:8)

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Package ggplot2 is required. Install it before running this script.")
}
if (save_pdf && !capabilities("cairo")) {
  stop("Vector PDF export requires Cairo support; set save_pdf = FALSE for PNG only.")
}
suppressPackageStartupMessages(library(ggplot2))

specs <- list(
  list(tag = "A", scope = "all_hosts", factor = "Host_type", table = "host",
       axis = "Host",
       file = "A_all_hosts_host_pairwise", width = 14.5, height = 14.0),
  list(tag = "B", scope = "coral", factor = "Host_type", table = "host",
       axis = "Host",
       file = "B_coral_host_pairwise", width = 9.0, height = 9.0),
  list(tag = "C", scope = "all_hosts", factor = "TS", table = "TS",
       axis = "Sampling batch (TS)",
       file = "C_all_hosts_TS_pairwise", width = 11.2, height = 10.7),
  list(tag = "D", scope = "coral", factor = "TS", table = "TS",
       axis = "Sampling batch (TS)",
       file = "D_coral_TS_pairwise", width = 11.2, height = 10.7)
)

# ---- 2. Read and validate the saved comparison tables ---------------------
pair_key <- function(a, b) {
  vapply(seq_along(a), function(i) paste(sort(c(a[i], b[i])), collapse = " | "),
         character(1))
}

read_comparisons <- function(spec) {
  path <- file.path(results_root, spec$scope, "pairwise",
                    paste0("pairwise_", spec$table, "_results.csv"))
  if (!file.exists(path)) stop("Missing saved result: ", path)
  d <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE,
                fileEncoding = "UTF-8-BOM", na.strings = c("NA", ""))
  required <- c("Factor", "Scope", "Group1", "Group2", "N_MC", "MC_set",
                "Input_MD5", "Block_factor", "Analysis_unit", "BH_family_size",
                "R2_raw", "P_BH_raw", "Status", "P_disp_BH_raw", "P_disp_BH_lingoes")
  missing <- setdiff(required, names(d))
  if (length(missing)) stop("Missing columns in ", path, ": ", paste(missing, collapse = ", "))
  if (!nrow(d) || anyNA(d[, c("Group1", "Group2", "Status")])) {
    stop("Empty or incomplete comparison identifiers: ", path)
  }
  expected_block <- if (spec$factor != "TS") "TS" else "Host_type"
  expected_mc <- "MC6;MC507;MC21;MC520;MC10;MC15;MC0;MC25"
  if (!isTRUE(all(d$Factor == spec$factor)) || !isTRUE(all(d$Scope == spec$scope)) ||
      !isTRUE(all(d$Block_factor == expected_block)) || !isTRUE(all(d$N_MC == 8L)) ||
      !isTRUE(all(d$MC_set == expected_mc)) ||
      !isTRUE(all(d$Analysis_unit == "Repeat_group\u5747\u503c"))) {
    stop("Unexpected factor, scope, permutation block, MC set, or analysis unit: ", path)
  }
  groups <- unique(c(d$Group1, d$Group2))
  if(spec$scope=="all_hosts" && spec$factor=="Host_type" &&
      !all(c("Seawater","Sediment") %in% groups)) stop("Old results exclude environmental samples; run the updated analysis first.")
  planned <- choose(length(groups), 2L)
  if (any(d$Group1 == d$Group2) || anyDuplicated(pair_key(d$Group1, d$Group2)) ||
      nrow(d) != planned || !isTRUE(all(d$BH_family_size == planned))) {
    stop("The saved table is not one complete planned pairwise family: ", path)
  }
  d$Estimated <- d$Status == "ok"
  ok <- d$Estimated
  if (!any(ok) || any(!is.finite(d$R2_raw[ok])) ||
      any(!is.finite(d$P_BH_raw[ok])) ||
      any(d$R2_raw[ok] < 0 | d$R2_raw[ok] > 1) ||
      any(d$P_BH_raw[ok] < 0 | d$P_BH_raw[ok] > 1)) {
    stop("Invalid R-squared or BH-adjusted P among estimated comparisons: ", path)
  }
  preferred <- switch(spec$factor, Host_type=host_order, Sample_type=c("HOST","Seawater","Sediment"), TS=ts_order)
  groups <- c(preferred[preferred %in% groups], sort(setdiff(groups, preferred)))
  list(data = d, groups = groups, source = path, family_size = planned)
}

specs[[length(specs)+1L]] <- list(tag="E",scope="sample_type",factor="Sample_type",
  table="sample_type",axis="Sample type",file="E_HOST_seawater_sediment_pairwise",width=9,height=9)
saved <- lapply(specs, read_comparisons)
fingerprints <- unique(unlist(lapply(saved, function(x) x$data$Input_MD5)))
if (length(fingerprints) != 1L || is.na(fingerprints) || !nzchar(fingerprints)) {
  stop("The comparison tables do not identify the same source-data fingerprint.")
}
all_r2 <- unlist(lapply(saved, function(x) x$data$R2_raw[x$data$Estimated]))
r2_limits[2] <- max(r2_limits[2], ceiling(max(all_r2)*10)/10)
if (length(r2_limits) != 2L || any(!is.finite(r2_limits)) ||
    r2_limits[1] >= r2_limits[2] || min(all_r2) < r2_limits[1] ||
    max(all_r2) > r2_limits[2]) {
  stop("r2_limits must include every estimated R-squared across all panels.")
}

# P_BH_raw is ALREADY BH-adjusted. Here 'raw' identifies uncorrected Bray-Curtis
# distances, rather than Lingoes-corrected distances; it does not mean unadjusted P.
# Unavailable comparisons stay NA. They are never converted to nonsignificance.
q_stars <- function(q) {
  ifelse(q <= 0.001, "***", ifelse(q <= 0.01, "**", ifelse(q <= 0.05, "*", "")))
}
q_label <- function(q) {
  vapply(q, function(x) {
    if (x < 0.001) sub("\\.0e", "e", formatC(x, format = "e", digits = 1)) else
      formatC(x, format = "f", digits = 3)
  }, character(1))
}
wrap_caption <- function(x, width) {
  lines <- strsplit(x, "\n", fixed = TRUE)[[1]]
  paste(vapply(lines, function(z) paste(strwrap(z, width = width), collapse = "\n"),
               character(1)), collapse = "\n")
}

# ---- 3. Build the five heatmaps ------------------------------------------------
make_heatmap <- function(spec, input) {
  d <- input$data
  groups <- input$groups
  n <- length(groups)
  index1 <- match(d$Group1, groups)
  index2 <- match(d$Group2, groups)
  d$Column <- groups[pmin(index1, index2)]
  d$Row <- groups[pmax(index1, index2)]
  valid <- d[d$Estimated, , drop = FALSE]
  unavailable <- d[!d$Estimated, , drop = FALSE]
  diagonal <- data.frame(Column = groups, Row = groups)

  # Dispersion has its own saved BH family. The dagger is separate from stars.
  dispersion_flag <- (is.finite(valid$P_disp_BH_raw) & valid$P_disp_BH_raw <= 0.05) |
    (is.finite(valid$P_disp_BH_lingoes) & valid$P_disp_BH_lingoes <= 0.05)
  valid$R2_label <- sprintf("%.3f", valid$R2_raw)
  valid$Q_label <- paste0("q=", q_label(valid$P_BH_raw))
  valid$Stars_label <- q_stars(valid$P_BH_raw)
  valid$Dispersion_label <- ifelse(show_dispersion_flags & dispersion_flag, "\u2020", "")
  valid$Text_colour <- pv_text_colour(valid$R2_raw, r2_limits, fill_colours)
  # geom_text sizes are in millimetres. Keep each label short; symbols do not
  # lengthen the numeric line. The larger coral canvas also leaves wider cells.
  text_sizes <- if (n >= 10) c(r2 = 3.8, q = 2.9, symbol = 2.9) else if (n >= 6)
    c(r2 = 4.0, q = 3.0, symbol = 3.1) else c(r2 = 4.2, q = 3.2, symbol = 3.3)
  design_note <- if (spec$factor != "TS") {
    if(spec$factor=="Sample_type") "Sample-type comparisons adjust for TS, with profile permutations within shared TS." else
      "Host comparisons adjust for TS, with profile permutations within shared TS."
  } else {
    "TS comparisons adjust for host, with profile permutations within shared hosts."
  }
  caption <- paste0(
    "Fill/center: marginal PERMANOVA R\u00b2; lower line: BH-adjusted P (q).\n",
    "Stars above R\u00b2: * q \u2264 0.05; ** q \u2264 0.01; *** q \u2264 0.001. BH family: ",
    input$family_size, " planned pairs.\n",
    "Eight MC proportions: Repeat_group means, square-root transformation, Bray\u2013Curtis distances.\n",
    design_note)
  if (show_dispersion_flags) {
    caption <- paste0(caption, "\n\u2020 (upper right): BH-adjusted PERMDISP P \u2264 0.05 for raw or Lingoes-corrected distances.")
  }
  if (nrow(unavailable)) {
    caption <- paste0(caption, "\nNA: comparison unavailable; see the saved source table for the reason.")
  }
  caption <- wrap_caption(caption, width = if (n >= 10) 132L else if (n >= 6) 110L else 85L)

  ggplot() +
    geom_tile(data = diagonal, aes(x = Column, y = Row), fill = "#F1F1F1",
              colour = "white", linewidth = 0.55) +
    geom_tile(data = valid, aes(x = Column, y = Row, fill = R2_raw),
              colour = "white", linewidth = 0.55) +
    geom_tile(data = unavailable, aes(x = Column, y = Row), fill = "#DADADA",
              colour = "white", linewidth = 0.55) +
    geom_text(data = valid, aes(x = Column, y = Row, label = R2_label, colour = Text_colour),
              position = position_nudge(y = 0.02), family = font_family,
              size = unname(text_sizes["r2"]), show.legend = FALSE) +
    geom_text(data = valid, aes(x = Column, y = Row, label = Q_label, colour = Text_colour),
              position = position_nudge(y = -0.26), family = font_family,
              size = unname(text_sizes["q"]), show.legend = FALSE) +
    geom_text(data = valid, aes(x = Column, y = Row, label = Stars_label, colour = Text_colour),
              position = position_nudge(y = 0.31), family = font_family,
              size = unname(text_sizes["symbol"]), show.legend = FALSE) +
    geom_text(data = valid, aes(x = Column, y = Row, label = Dispersion_label, colour = Text_colour),
              position = position_nudge(x = 0.35, y = 0.31), family = font_family,
              size = unname(text_sizes["symbol"]), show.legend = FALSE) +
    geom_text(data = unavailable, aes(x = Column, y = Row), label = "NA",
              family = font_family, size = if (n >= 10) 3.7 else 4.1, colour = "#444444") +
    scale_colour_identity() +
    scale_fill_gradientn(colours = fill_colours, limits = r2_limits,
                          breaks = seq(r2_limits[1], r2_limits[2], length.out = 7L),
                          labels = function(x) sprintf("%.2f", x),
                          name = expression(PERMANOVA~R^2),
                          guide = guide_colourbar(title.position = "top",
                            barwidth = grid::unit(75, "mm"), barheight = grid::unit(6, "mm"),
                            ticks.colour = "black", frame.colour = "black")) +
    scale_x_discrete(limits = groups, drop = FALSE, expand = expansion(add = 0.5)) +
    scale_y_discrete(limits = rev(groups), drop = FALSE, expand = expansion(add = 0.5)) +
    coord_fixed(clip = "off") +
    labs(title = NULL, subtitle = NULL, x = spec$axis, y = spec$axis,
          tag = spec$tag, caption = caption) +
    theme_classic(base_size = 14, base_family = font_family) +
    theme(
      axis.line = element_line(colour = "black", linewidth = 0.45),
      axis.ticks = element_line(colour = "black", linewidth = 0.4),
      axis.ticks.length = grid::unit(2, "mm"),
      axis.text = element_text(colour = "black", size = 13),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      axis.title = element_text(colour = "black", size = 15),
      axis.title.x = element_text(margin = margin(t = 10)),
      axis.title.y = element_text(margin = margin(r = 10)),
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      plot.tag = element_text(size = 27, face = "bold", colour = "black"),
      plot.tag.position = "topleft",
      plot.caption = element_text(size = 10, colour = "#333333", hjust = 0,
                                   lineheight = 1.15, margin = margin(t = 16)),
      plot.caption.position = "plot",
      legend.position = "top", legend.justification = "right",
      legend.title = element_text(size = 12), legend.text = element_text(size = 10),
      legend.margin = margin(t = 2, b = 8),
      panel.grid = element_blank(),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background = element_rect(fill = "white", colour = NA),
      plot.margin = margin(t = 12, r = 20, b = 18, l = 18)
    )
}

plots <- Map(make_heatmap, specs, saved)
names(plots) <- vapply(specs, function(x) x$file, character(1))

# ---- 4. Save only when the user runs this script -------------------------------
if (dir.exists(output_dir)) stop("Output directory already exists; choose a new output_dir.")
if (!dir.create(output_dir, recursive = TRUE)) stop("Could not create output directory: ", output_dir)
png_device <- if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_png else "png"
manifest <- vector("list", length(specs))
for (i in seq_along(specs)) {
  s <- specs[[i]]
  png_path <- file.path(output_dir, paste0(s$file, ".png"))
  ggsave(filename = png_path, plot = plots[[i]], device = png_device,
          width = s$width, height = s$height, units = "in", dpi = dpi, bg = "white")
  if (save_pdf) {
    ggsave(filename = file.path(output_dir, paste0(s$file, ".pdf")), plot = plots[[i]],
            device = grDevices::cairo_pdf, width = s$width, height = s$height,
            units = "in", bg = "white")
  }
  manifest[[i]] <- data.frame(
    Panel = s$tag, Scope = s$scope, Factor = s$factor, Source_CSV = saved[[i]]$source,
    Source_MD5 = unname(tools::md5sum(saved[[i]]$source)), Input_MD5 = fingerprints,
    Planned_pairs = saved[[i]]$family_size, Estimated_pairs = sum(saved[[i]]$data$Estimated),
    BH_significant_pairs = sum(saved[[i]]$data$P_BH_raw[saved[[i]]$data$Estimated] <= 0.05),
    PNG = png_path, R2_min = r2_limits[1], R2_max = r2_limits[2],
    stringsAsFactors = FALSE
  )
}
write.csv(do.call(rbind, manifest), file.path(output_dir, "figure_sources.csv"),
           row.names = FALSE, fileEncoding = "UTF-8")
writeLines(c("Saved-result plotting only; no statistical analyses were rerun.",
              paste0("Results: ", results_root),
              "R-squared: R2_raw. BH-adjusted P: P_BH_raw (used directly).",
              "Dispersion flags: P_disp_BH_raw or P_disp_BH_lingoes <= 0.05.",
              "Each panel is a separate BH family, including planned unavailable pairs.",
              "No titles/subtitles. Cells: R-squared at center, q below, stars above, dagger at upper right.",
              paste0("Font: ", font_family, "; PNG resolution: ", dpi, " dpi.")),
            file.path(output_dir, "README.txt"), useBytes = TRUE)
message(length(specs)," heatmaps saved in: ", output_dir)
