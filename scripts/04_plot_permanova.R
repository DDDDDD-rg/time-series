# 最终展示图：宿主、时间、干湿季，以及BH校正后显著的环境指标。
# 先完成02、03和05，再从头运行本文件；不重跑PERMANOVA或置换检验。
# 椭圆描述组内分布，不表示该组显著，也不用于判断两组是否显著不同。

# ===== 1. 设置和读取结果 =====
project_dir <- "C:/Users/xzy/Desktop/CUHK_Study/GNBF6010_Research_Projects/Research/Msc main project-time series"
source(file.path(project_dir, "scripts", "permanova_plot_style.R"), local=TRUE)
analysis_scope <- getOption("mc.analysis_scope", "all_hosts")  # 与01/05范围一致
random_seed <- 20261004
ellipse_level <- 0.95
show_pairwise_letters <- FALSE # Pairwise results are shown in the separate heatmaps.
pairwise_alpha <- 0.05
make_interactive_3d <- TRUE   # FALSE keeps the existing 2D-only output.
nmds_3d_workers <- 4L        # Windows-compatible PSOCK workers for random starts.
nmds_3d_reuse_saved <- TRUE # Reuse compatible saved coordinates for display-only edits.
nmds_3d_trymax <- 200L
latest_result_pointer <- file.path(project_dir, "results", "00_latest_PERMANOVA_MC8.txt")
latest_figure_pointer <- file.path(project_dir, "figures", "00_latest_PERMANOVA_MC8.txt")
default_results_root <- if (file.exists(latest_result_pointer))
  readLines(latest_result_pointer, n = 1L, warn = FALSE, encoding = "UTF-8") else
  file.path(project_dir, "results", "PERMANOVA_MC8")
default_figures_root <- if (file.exists(latest_figure_pointer))
  readLines(latest_figure_pointer, n = 1L, warn = FALSE, encoding = "UTF-8") else
  file.path(project_dir, "figures", "PERMANOVA_MC8")
result_dir <- file.path(getOption("mc.results_root", default_results_root), analysis_scope)
figure_dir <- getOption("mc.figure_dir", file.path(default_figures_root,
                         paste0(analysis_scope, "_", format(Sys.time(), "%Y%m%d_%H%M%S"))))

packages <- c("ggplot2", "vegan")
if (make_interactive_3d) packages <- c(packages,"plotly","htmlwidgets","htmltools")
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("请先安装R包：", paste(missing, collapse = ", "))
suppressPackageStartupMessages(library(ggplot2))
html_3d_selfcontained <- FALSE
if (make_interactive_3d) {
  stopifnot(length(nmds_3d_workers)==1L, is.finite(nmds_3d_workers), nmds_3d_workers>=1L,
            length(nmds_3d_trymax)==1L, is.finite(nmds_3d_trymax), nmds_3d_trymax>=30L)
  html_3d_selfcontained <- requireNamespace("rmarkdown",quietly=TRUE) && rmarkdown::pandoc_available()
  if (!html_3d_selfcontained) message("HTML export: Pandoc unavailable; HTML and local _files folders will be saved together.")
}
stopifnot(analysis_scope %in% c("all_hosts", "coral"))
files <- c(
  file.path(result_dir, "host_time_analysis.rds"),
  file.path(result_dir, "repeat_groups", "season_environment_analysis.rds")
)
if (!all(file.exists(files))) stop("缺少结果RDS，请先完成对应范围的02和03。")
host_saved <- readRDS(files[1])
season_env_saved <- readRDS(files[2])
host_data <- as.data.frame(host_saved$data)
if (analysis_scope=="all_hosts" && !all(c("Seawater","Sediment") %in% host_data$Host_type))
  stop("Saved results exclude environmental samples. Run the updated run_permanova.R first.")
season_data <- as.data.frame(season_env_saved$season_input)
env_data <- as.data.frame(season_env_saved$environment_input)
host_result <- host_saved$results[host_saved$results$Factor == "Host_type", , drop = FALSE]
time_result <- host_saved$results[host_saved$results$Factor == "TS", , drop = FALSE]
season_result <- season_env_saved$season
env_results <- season_env_saved$environment
stopifnot(nrow(host_result) == 1, nrow(time_result) == 1, nrow(season_result) == 1,
          host_result$Status == "ok", time_result$Status == "ok", season_result$Status == "ok",
          all(host_saved$results$Scope == analysis_scope), all(season_result$Scope == analysis_scope),
          all(env_results$Scope == analysis_scope))
# 季节和环境图必须读取取消第二次平均后的结果，防止混用旧均值结果。
stopifnot("Repeat_group_ID" %in% names(season_data),
          "Repeat_group_ID" %in% names(env_data),
          anyDuplicated(season_data$Repeat_group_ID) == 0,
          anyDuplicated(env_data$Repeat_group_ID) == 0,
          all(season_result$Analysis_unit == "Repeat_group均值"),
          all(env_results$Analysis_unit == "Repeat_group均值"))
mc_cols <- grep("^MC", names(host_data), value = TRUE)
stopifnot(identical(mc_cols, grep("^MC", names(season_data), value = TRUE)),
          identical(mc_cols, grep("^MC", names(env_data), value = TRUE)), length(mc_cols) == 8L)
host_data$Host_type <- droplevels(factor(host_data$Host_type))
host_data$TS <- droplevels(factor(host_data$TS, levels = paste0("TS", 1:8)))

# 字母取自05的BH两两比较；先检查它与02使用同一批Repeat_group数据。
letter_displays <- list()
if (show_pairwise_letters) {
  source(file.path(project_dir, "scripts", "pairwise_letters.R"), encoding = "UTF-8", local = TRUE)
  pairwise_file <- file.path(result_dir, "pairwise", "pairwise_host_time_analysis.rds")
  if (!file.exists(pairwise_file)) stop("缺少两两结果，请先运行同范围的05；或将show_pairwise_letters设为FALSE。")
  pairwise_saved <- readRDS(pairwise_file)
  if (!identical(pairwise_saved$scope, analysis_scope)) stop("05与04的分析范围不同。")
  if (!all(c("Host_type", "TS") %in% names(pairwise_saved$families))) {
    stop("字母图需要宿主和TS两类结果，请在05选择两类比较后运行。")
  }
  comparison_data <- as.data.frame(pairwise_saved$data)
  if (anyDuplicated(comparison_data$Repeat_group_ID) ||
      !setequal(comparison_data$Repeat_group_ID, host_data$Repeat_group_ID)) {
    stop("02与05的Repeat_group不一致，请更新对应结果后绘图。")
  }
  comparison_data <- comparison_data[match(host_data$Repeat_group_ID,
                                            comparison_data$Repeat_group_ID), , drop = FALSE]
  if (!all(c(mc_cols, "Host_type", "TS", "n_records") %in% names(comparison_data))) {
    stop("05缺少与02匹配的组成或元数据列。")
  }
  for (column in c(mc_cols, "Host_type", "TS", "n_records")) {
    same <- if (column %in% mc_cols)
      isTRUE(all.equal(host_data[[column]], comparison_data[[column]], tolerance = 1e-12)) else
      identical(as.character(host_data[[column]]), as.character(comparison_data[[column]]))
    if (!same) stop("02与05的数据内容不一致：", column, "。请更新结果后绘图。")
  }
  for (effect in c("Host_type", "TS")) {
    family_results <- pairwise_saved$families[[effect]]$results
    if (!all(family_results$Factor == effect) ||
        !all(family_results$Scope == analysis_scope) ||
        !all(family_results$Analysis_unit == "Repeat_group均值")) {
      stop("05的因素、范围或分析单位与04不一致：", effect)
    }
    letter_displays[[effect]] <- derive_pairwise_letters(
      family_results, levels(host_data[[effect]]), pairwise_alpha
    )
  }
}

plots <- list()
ellipse_notes <- character()
overview_captions <- character()

# ===== 2. 共用工具：坐标、图注、椭圆和保存 =====
p_text <- function(p) ifelse(is.finite(p), formatC(p, format = "g", digits = 3), "NA")
stars <- function(p) ifelse(!is.finite(p), "NA", ifelse(p <= 0.001, "***",
                         ifelse(p <= 0.01, "**", ifelse(p <= 0.05, "*", "ns"))))
wrap_caption <- function(text) paste(vapply(strsplit(text, "\n", fixed = TRUE)[[1]],
  function(line) paste(strwrap(line, width = 110), collapse = "\n"), character(1)), collapse = "\n")
plot_theme <- pv_theme()

make_nmds <- function(data, reuse_saved=TRUE) {
  if (reuse_saved && exists("host_nmds",envir=environment(make_nmds),inherits=FALSE)) {
    candidate <- host_nmds
    index <- match(as.character(data$Repeat_group_ID),as.character(candidate$data$Repeat_group_ID))
    if (nrow(data)==nrow(candidate$data) && !anyNA(index) && !anyDuplicated(index) &&
        all(vapply(c("Host_type","TS",mc_cols),function(column)
          isTRUE(all.equal(as.character(data[[column]]),
            as.character(candidate$data[[column]][index]))),logical(1)))) {
      candidate$data <- cbind(data,candidate$data[index,c("NMDS1","NMDS2"),drop=FALSE])
      return(candidate)
    }
  }
  # Reuse matching 2D coordinates too when only display settings have changed.
  cached_file <- file.path(result_dir,"plot_ordinations.rds")
  if (reuse_saved && file.exists(cached_file)) {
    cached <- tryCatch(readRDS(cached_file),error=function(e) NULL)
    if (!is.null(cached) && identical(cached$MC,mc_cols) &&
        setequal(as.character(cached$input_md5),as.character(unique(host_result$Input_MD5)))) {
      for (key in c("host","season","environment")) {
        candidate <- cached[[key]]
        if (is.null(candidate) || is.null(candidate$data)) next
        d <- candidate$data
        if (!all(c("Repeat_group_ID","Host_type","TS",mc_cols,"NMDS1","NMDS2") %in% names(d)) ||
            !identical(as.character(d$Repeat_group_ID),as.character(data$Repeat_group_ID))) next
        same <- vapply(c("Host_type","TS",mc_cols),function(column) {
          if (column %in% mc_cols) isTRUE(all.equal(d[[column]],data[[column]],tolerance=1e-12)) else
            identical(as.character(d[[column]]),as.character(data[[column]]))
        },logical(1))
        if (!all(same) || !identical(as.integer(candidate$fit$ndim),2L) ||
            !is.finite(candidate$fit$stress) ||
            !all(is.finite(as.matrix(d[,c("NMDS1","NMDS2")])))) next
        candidate$data <- cbind(data,d[,c("NMDS1","NMDS2"),drop=FALSE])
        message("Reusing saved 2D NMDS coordinates: ",key)
        return(candidate)
      }
    }
  }
  matrix <- as.matrix(data[, mc_cols, drop = FALSE])
  stopifnot(!anyNA(matrix), all(is.finite(matrix)), all(matrix >= 0), all(rowSums(matrix) > 0))
  rownames(matrix) <- if ("Repeat_group_ID" %in% names(data)) as.character(data$Repeat_group_ID) else
    paste(data$TS, data$Host_type, sep = "|")
  stopifnot(anyDuplicated(rownames(matrix)) == 0)
  distance <- vegan::vegdist(sqrt(matrix), method = "bray")
  set.seed(random_seed)
  fit <- vegan::metaMDS(distance, k = 2, trymax = 100, trace = FALSE)
  xy <- as.data.frame(vegan::scores(fit, display = "sites", choices = 1:2))
  stopifnot(identical(rownames(xy), attr(distance, "Labels")), is.finite(fit$stress))
  list(fit = fit, data = cbind(data, xy))
}

# stat_ellipse需要足够的点和二维变异：不足4点或坐标共线的组仅保留散点。
ellipse_data <- function(data, group) {
  groups <- split(seq_len(nrow(data)), as.character(data[[group]]))
  eligible <- vapply(groups, function(index) {
    xy <- as.matrix(data[index, c("NMDS1", "NMDS2"), drop = FALSE])
    nrow(xy) >= 4 && qr(scale(xy, center = TRUE, scale = FALSE))$rank == 2
  }, logical(1))
  skipped <- names(groups)[!eligible]
  list(data = data[unlist(groups[eligible], use.names = FALSE), , drop = FALSE], skipped = skipped)
}

result_caption <- function(result, use_bh = TRUE) {
  p <- if (use_bh) result$P_BH_raw else result$P_raw
  ling_p <- if (use_bh) result$P_BH_lingoes else result$P_lingoes
  dispersion <- c(result$P_disp_raw, result$P_disp_lingoes)
  caution <- if (any(is.finite(dispersion) & dispersion <= 0.05)) " Dispersion differences detected: interpret with caution." else ""
  paste0("PERMANOVA (raw Bray): R2 = ", sprintf("%.3f", result$R2_raw),
    if (use_bh) "; BH p = " else "; p = ", p_text(p), " (", stars(p), "). Lingoes ",
    if (use_bh) "BH p = " else "p = ", p_text(ling_p), ".\n",
    "PERMDISP p (raw / Lingoes, unadjusted) = ", p_text(result$P_disp_raw), " / ",
    p_text(result$P_disp_lingoes), ".", caution)
}

group_plot <- function(ordination, group, title, legend, result, unit, colours = NULL,
                       use_bh = TRUE, letter_display = NULL,
                       show_ellipses = identical(group, "Season")) {
  data <- ordination$data
  data$Group <- droplevels(factor(data[[group]]))
  # Compact host/TS overviews: identical circles, colour-coded groups, thin ellipses.
  if (group %in% c("Host_type", "TS")) {
    palette <- pv_subset(if (group == "Host_type") pv_hosts else pv_ts, levels(data$Group))
    ellipse <- ellipse_data(data, "Group")
    if (length(ellipse$skipped)) ellipse_notes <<- c(ellipse_notes,
      paste(group, "No ellipse (<4 points or collinear):", paste(ellipse$skipped, collapse=", ")))
    overview_captions[group] <<- paste0(
      if (group == "Host_type") "Host overview. Colour represents host. " else
        "Sampling-batch overview. Colour represents sampling batch. ",
      "Each point is one Repeat_group mean; eight MC relative abundances were square-root transformed before Bray-Curtis NMDS. ",
      "All points and the full coordinate range are retained. Ellipses show ",
      round(100 * ellipse_level), "% data ellipses (normal approximation), not confidence regions for group centres. ",
      if (group == "Host_type") "The host PERMANOVA adjusts for TS, with permutations within TS. " else
        "The TS PERMANOVA adjusts for host, with permutations within host. ",
      "NMDS coordinates show unadjusted Bray-Curtis distances. ",
      result_caption(result, use_bh),
      " Pairwise comparisons are reported in the separate BH-adjusted heatmaps.",
      if (length(ellipse$skipped)) paste0(" No ellipse (<4 points or collinear): ",
        paste(ellipse$skipped, collapse=", "), ".") else "")
    return(ggplot(data, aes(NMDS1, NMDS2)) +
      stat_ellipse(data=ellipse$data, aes(colour=Group, group=Group),
        type="norm", level=ellipse_level, linewidth=0.4, alpha=0.65, show.legend=FALSE) +
      geom_point(aes(fill=Group), shape=21, colour="white", size=2.6, stroke=0.25, alpha=1) +
      scale_fill_manual(values=palette, name=legend) +
      scale_colour_manual(values=palette, guide="none") +
      guides(fill=guide_legend(override.aes=list(shape=21, colour="white", size=3.2, alpha=1))) +
      coord_equal() +
      labs(title=NULL, subtitle=NULL,
        caption=sprintf("Stress = %.3f; n = %d; R\u00b2 = %.3f; BH p = %s",
          ordination$fit$stress, nrow(data), result$R2_raw, p_text(result$P_BH_raw))) +
      plot_theme)
  }
  colours <- pv_subset(switch(group, Host_type=pv_hosts, TS=pv_ts, Season=pv_season), levels(data$Group))
  data$Shape_group <- if (group %in% c("Host_type", "TS"))
    droplevels(factor(data$Host_type)) else data$Group
  shapes <- pv_subset(if (group %in% c("Host_type", "TS")) pv_host_shapes else
    c(Dry=16, Wet=17), levels(data$Shape_group))
  legend_labels <- setNames(levels(data$Group), levels(data$Group))
  letter_note <- ""
  if (!is.null(letter_display)) {
    table <- letter_display$table
    labels <- table$Letters[match(levels(data$Group), table$Group)]
    labels[is.na(labels)] <- "n/a"
    legend_labels[] <- paste0(levels(data$Group), " [", labels, "]")
    letter_note <- paste0("\nLegend letters: pairwise PERMANOVA, BH across ",
      letter_display$family_size, " comparisons; alpha = ", letter_display$alpha,
      ". Shared letters = difference not detected; no shared letters = significant.")
    if (length(letter_display$excluded)) {
      letter_note <- paste0(letter_note, "\nn/a = incomplete comparisons; letters apply to labelled groups.")
    }
  }
  note <- ""
  ellipse_note <- if (show_ellipses) paste0(round(100 * ellipse_level),
    "% data ellipses (normal approximation); ellipses are not significance tests.") else
    if (group == "TS") "Colour = sampling batch; point shape = host. Ellipses omitted." else
      "Colour and point shape = host. Ellipses omitted."
  if (show_ellipses) {
    ellipse <- ellipse_data(data, "Group")
    note <- if (length(ellipse$skipped)) paste0("\nNo ellipse (<4 points or collinear): ",
      paste(ellipse$skipped, collapse = ", "), ".") else ""
    if (length(ellipse$skipped)) ellipse_notes <<- c(ellipse_notes, paste(title, note))
  }
  shape_labels <- if (group == "TS") names(shapes) else unname(legend_labels[names(shapes)])
  plot <- ggplot(data, aes(NMDS1, NMDS2, colour = Group, shape = Shape_group)) +
    geom_point(size = 3.2, alpha = 1, stroke = 0.85) + coord_equal() +
    scale_shape_manual(values=shapes, limits=names(shapes), labels=shape_labels) +
    labs(title = NULL, colour = legend, shape = if (group == "TS") "Host" else legend,
      subtitle = sprintf("%d selected MCs; stress = %.3f; %d %s", length(mc_cols), ordination$fit$stress, nrow(data), unit),
      caption = wrap_caption(paste0(result_caption(result, use_bh), "\n",
        ellipse_note,
        sprintf("\nStress = %.3f; n = %d %s.", ordination$fit$stress, nrow(data), unit),
        letter_note, note))) + plot_theme
  if (show_ellipses) plot <- plot +
    stat_ellipse(data=ellipse$data, aes(NMDS1, NMDS2, colour=Group, group=Group),
      inherit.aes=FALSE, type="norm", level=ellipse_level,
      linewidth=0.85, alpha=0.9, show.legend=FALSE)
  if (group == "TS") plot <- plot + guides(
    colour=guide_legend(order=1, override.aes=list(shape=16, size=3.5, alpha=1)),
    shape=guide_legend(order=2, override.aes=list(colour="black", size=3.5, alpha=1)))
  if (!is.null(colours)) {
    plot <- plot + scale_colour_manual(values = colours, limits = names(colours),
                                       labels = unname(legend_labels[names(colours)]))
  } else {
    plot <- plot + scale_colour_discrete(limits = levels(data$Group),
                                         labels = unname(legend_labels))
  }
  plot
}

# Definitions sourced by the updated 04 script; no fitting until called.
fit_interactive_nmds_3d <- function(data, mc, seed, workers=4L, trymax=200L) {
  community <- as.matrix(data[, mc, drop=FALSE])
  ids <- as.character(data$Repeat_group_ID)
  stopifnot(nrow(data)>7L, !anyNA(community), all(is.finite(community)),
    all(community>=0), all(rowSums(community)>0), !anyNA(ids), !anyDuplicated(ids))
  rownames(community) <- ids
  distance <- vegan::vegdist(sqrt(community), method="bray")
  cores <- parallel::detectCores(logical=FALSE)
  if (!is.finite(cores)) cores <- 1L
  workers <- max(1L, min(as.integer(workers), as.integer(cores)))
  cluster <- NULL
  if (workers>1L) {
    cluster <- parallel::makePSOCKcluster(workers)
    on.exit(parallel::stopCluster(cluster), add=TRUE)
    parallel::clusterEvalQ(cluster, library(vegan))
    parallel::clusterSetRNGStream(cluster, iseed=seed)
  }
  set.seed(seed)
  fit <- vegan::metaMDS(distance, k=3, try=30, trymax=trymax, trace=FALSE,
    parallel=if (is.null(cluster)) 1L else cluster)
  xyz <- as.data.frame(vegan::scores(fit, display="sites", choices=1:3))
  stopifnot(identical(rownames(xyz), ids),
    identical(names(xyz), c("NMDS1","NMDS2","NMDS3")),
    all(is.finite(as.matrix(xyz))), is.finite(fit$stress))
  if (is.null(fit$converged) || fit$converged<1L)
    message("3D NMDS: best solution was not repeated; inspect the saved diagnostics before interpretation.")
  list(fit=fit, data=cbind(data,xyz), workers=workers, distance=distance)
}

nmds_3d_hover <- function(data, mc) {
  escape <- function(x) as.character(htmltools::htmlEscape(as.character(x)))
  text <- paste0("Host: ",escape(data$Host_type),"<br>TS: ",escape(data$TS),
    "<br>Repeat group: ",escape(data$Repeat_group_ID))
  if ("DO_mg_L" %in% names(data)) text <- paste0(text,"<br>DO: ",
    ifelse(is.finite(data$DO_mg_L),sprintf("%.2f mg/L",data$DO_mg_L),"not available"))
  text <- paste0(text,"<br>NMDS1: ",sprintf("%.4f",data$NMDS1),
    "<br>NMDS2: ",sprintf("%.4f",data$NMDS2),"<br>NMDS3: ",sprintf("%.4f",data$NMDS3))
  for (name in mc) text <- paste0(text,"<br>",escape(name),": ",sprintf("%.2f%%",100*data[[name]]))
  text
}

# Axes are rendered as geometry at fixed data coordinates, not auto-selected box edges.
nmds_3d_fixed_axes <- function(widget, data) {
  coordinates <- as.matrix(data[,c("NMDS1","NMDS2","NMDS3")])
  bounds <- apply(coordinates,2,range)
  span <- pmax(bounds[2,]-bounds[1,],0.01)
  low <- bounds[1,]-0.04*span
  high <- bounds[2,]+0.04*span
  annotations <- list()
  for (dimension in 1:3) {
    start <- low; end <- low; end[dimension] <- high[dimension]
    widget <- plotly::add_trace(widget,x=c(start[1],end[1]),y=c(start[2],end[2]),
      z=c(start[3],end[3]),type="scatter3d",mode="lines",showlegend=FALSE,
      hoverinfo="skip",name=paste0("fixed_axis_",dimension),
      line=list(color="#333333",width=3),inherit=FALSE)
    ticks <- pretty(bounds[,dimension],n=4)
    ticks <- ticks[ticks>=low[dimension] & ticks<=high[dimension]]
    offset_axis <- if (dimension==1L) 2L else 1L
    for (value in ticks) {
      position <- low; position[dimension] <- value
      tip <- position; tip[offset_axis] <- tip[offset_axis]-0.012*span[offset_axis]
      widget <- plotly::add_trace(widget,x=c(position[1],tip[1]),y=c(position[2],tip[2]),
        z=c(position[3],tip[3]),type="scatter3d",mode="lines",showlegend=FALSE,
        hoverinfo="skip",name=paste0("fixed_tick_",dimension,"_",value),
        line=list(color="#333333",width=2),inherit=FALSE)
      label_position <- tip
      label_position[offset_axis] <- label_position[offset_axis]-0.025*span[offset_axis]
      annotations[[length(annotations)+1L]] <- list(x=label_position[1],y=label_position[2],
        z=label_position[3],text=format(value,trim=TRUE,digits=3),
        showarrow=FALSE,font=list(size=11,color="#333333"),captureevents=FALSE)
    }
    label_position <- low
    label_position[dimension] <- high[dimension]+0.05*span[dimension]
    annotations[[length(annotations)+1L]] <- list(x=label_position[1],y=label_position[2],
      z=label_position[3],text=paste0("NMDS",dimension),showarrow=FALSE,
      font=list(size=14,color="black"),captureevents=FALSE)
  }
  list(widget=widget, annotations=annotations)
}

make_interactive_nmds_3d <- function(ordination, group, palette, result, mc) {
  data <- ordination$data
  # Saturated categorical colours; mapping stays fixed for both analysis scopes.
  if (group=="Host_type") palette <- c("Acropora sp."="#0057B8",
    "Montipora sp."="#E64B00","Platygyra sp."="#008F39","CCA"="#D10087",
    "Diadema sp."="#B59A00","Drupella sp."="#6422C8","Haliclona sp."="#00A0B0",
    "Holothuria sp."="#C40023","Perna sp."="#00B9F2","Sea anemone"="#717600",
    "Shrimp"="#713B00","Trochus sp."="#111111",
    "Seawater"="#005F73","Sediment"="#9B5DE5")
  values <- as.character(data[[group]])
  groups <- if (is.factor(data[[group]])) levels(droplevels(data[[group]])) else unique(values)
  groups <- names(palette)[names(palette) %in% groups]
  stopifnot(setequal(groups,unique(values)))
  hover <- nmds_3d_hover(data,mc)
  widget <- plotly::plot_ly()
  axes <- nmds_3d_fixed_axes(widget,data)
  widget <- axes$widget
  for (name in groups) {
    index <- which(values==name)
    widget <- plotly::add_trace(widget,x=data$NMDS1[index],y=data$NMDS2[index],
      z=data$NMDS3[index],type="scatter3d",mode="markers",showlegend=TRUE,
      name=paste0(name," (n = ",length(index),")"),legendgroup=name,
      text=hover[index],hoverinfo="text",
      marker=list(symbol="circle",size=3,opacity=1,color=unname(palette[name]),
        line=list(color="white",width=0.25)),inherit=FALSE)
  }
  axis <- function(values) {
    limits <- range(values)
    pad <- max(diff(limits)*0.18,0.02)
    list(range=limits+c(-pad,pad),autorange=FALSE,visible=FALSE,
      showspikes=FALSE,showbackground=FALSE)
  }
  summary <- sprintf("Stress (3D) = %.3f; n = %d; R\u00b2 = %.3f; BH p = %s",
    ordination$fit$stress,nrow(data),result$R2_raw,
    formatC(result$P_BH_raw,format="g",digits=3))
  widget <- plotly::layout(widget,title=list(text=""),uirevision="NMDS_3D_fixed_axes",
    paper_bgcolor="white",font=list(family="Arial",size=14,color="black"),
    margin=list(l=0,r=20,t=15,b=100),
    scene=list(xaxis=axis(data$NMDS1),yaxis=axis(data$NMDS2),zaxis=axis(data$NMDS3),
      aspectmode="data",dragmode="turntable",uirevision="NMDS_3D_fixed_axes",
      annotations=axes$annotations,
      camera=list(eye=list(x=1.5,y=1.5,z=1.1),up=list(x=0,y=0,z=1),
        projection=list(type="orthographic"))),
    legend=list(title=list(text=if(group=="Host_type") "Host" else "Sampling batch"),
      x=1.01,y=1,itemclick="toggle",itemdoubleclick="toggleothers",itemsizing="constant"),
    annotations=list(
      list(text=summary,x=0,y=-0.08,xref="paper",yref="paper",xanchor="left",showarrow=FALSE),
      list(text="Hover legend: highlight without fading. Click: toggle multiple groups; Show all: reset.",
        x=0,y=-0.14,xref="paper",yref="paper",xanchor="left",showarrow=FALSE,font=list(size=12))))
  widget <- plotly::layout(widget,
    updatemenus=list(list(type="buttons",direction="right",x=1.01,y=1.10,
      xanchor="left",yanchor="top",showactive=FALSE,
      buttons=list(list(label="Show all",name="nmds_show_all",method="skip",args=list())))))
  widget <- htmlwidgets::onRender(widget,"function(el,x,options){
  var axes=[],samples=[],original={},selected=new Set(),hovered=null;
  var restoring=false,pending=false,dirty=false;
  el.data.forEach(function(t,i){
    if(/^fixed_(axis|tick)_/.test(t.name||'')) axes.push(i);
    if(t.type==='scatter3d' && t.mode==='markers'){
      samples.push(i);
      original[i]={color:t.marker.color,size:t.marker.size,opacity:t.marker.opacity,
        lineColor:t.marker.line.color,lineWidth:t.marker.line.width};
    }
  });
  // Coalesce rapid pointer/click changes; never restyle coordinates or the camera.
  function render(){
    if(pending){dirty=true;return;}
    pending=true;dirty=false;
    var update={'marker.color':[],'marker.size':[],'marker.opacity':[],
      'marker.line.color':[],'marker.line.width':[],visible:[]};
    samples.forEach(function(i){
      // Hover temporarily restores everyone's colours. Leaving the legend
      // restores the persistent multi-selection and its grey background.
      var muted=hovered===null && selected.size>0 && !selected.has(i);
      var size=i===hovered?5:(selected.has(i)?4:original[i].size);
      update['marker.color'].push(muted?'#BDBDBD':original[i].color);
      update['marker.size'].push(muted?2:size);
      update['marker.opacity'].push(muted?0.18:original[i].opacity);
      update['marker.line.color'].push(muted?'#BDBDBD':original[i].lineColor);
      update['marker.line.width'].push(muted?0:original[i].lineWidth);
      update.visible.push(true);
    });
    Promise.resolve(Plotly.restyle(el,update,samples)).then(done,done);
    function done(){pending=false;if(dirty) render();}
  }
  function legendEntry(target){
    if(!target || !target.closest) return null;
    var entry=target.closest('.legend .traces');
    return entry && el.contains(entry)?entry:null;
  }
  function entryIndex(entry){
    // Plotly's legend binds an array of legend items to each .traces group.
    var datum=entry && entry.__data__;
    var item=Array.isArray(datum)?datum[0]:datum;
    var trace=item && item.trace;
    if(!trace) return null;
    for(var j=0;j<samples.length;j++){
      var i=samples[j],candidate=el.data[i];
      if(candidate.legendgroup===trace.legendgroup && candidate.name===trace.name) return i;
    }
    return null;
  }
  // Native legend hover has no public plotly_legendhover event; use delegated
  // DOM mouse events, which also work after Plotly rebuilds the legend.
  el.addEventListener('mouseover',function(event){
    var entry=legendEntry(event.target);
    if(!entry || (event.relatedTarget && entry.contains(event.relatedTarget))) return;
    var index=entryIndex(entry);
    if(index!==null && hovered!==index){hovered=index;render();}
  });
  el.addEventListener('mouseout',function(event){
    var entry=legendEntry(event.target);
    if(!entry || (event.relatedTarget && entry.contains(event.relatedTarget))) return;
    if(hovered!==null){hovered=null;render();}
  });
  el.on('plotly_legendclick',function(event){
    var index=event.curveNumber;
    if(samples.indexOf(index)<0) return false;
    hovered=null; // Clicking immediately applies persistent selection/fading.
    if(selected.has(index)) selected.delete(index); else selected.add(index);
    render();
    return false;
  });
  el.on('plotly_legenddoubleclick',function(){return false;});
  el.on('plotly_buttonclicked',function(event){
    if(event.button && event.button.name==='nmds_show_all'){
      selected.clear();hovered=null;render();
    }
  });
  el.on('plotly_restyle',function(){
    if(restoring) return;
    var hidden=axes.filter(function(i){
      return el.data[i].visible===false || el.data[i].visible==='legendonly';
    });
    if(hidden.length){
      restoring=true;
      Promise.resolve(Plotly.restyle(el,{visible:true},hidden))
        .then(function(){restoring=false;},function(){restoring=false;});
    }
  });
}",data=list(group=group))
  plotly::config(widget,displaylogo=FALSE,scrollZoom=TRUE,responsive=TRUE,
    modeBarButtonsToRemove=c("pan3d","orbitRotation","resetCameraLastSave3d"),
    toImageButtonOptions=list(format="png",filename=paste0("NMDS_3D_",group),scale=2))
}

find_saved_nmds_3d <- function(host_data,mc,input_md5,figure_dir) {
  root <- dirname(figure_dir)
  folders <- list.dirs(root,recursive=FALSE,full.names=TRUE)
  files <- file.path(folders,"05_NMDS_3D_coordinates.rds")
  files <- files[file.exists(files)]
  if (!length(files)) return(NULL)
  files <- files[order(file.info(files)$mtime,decreasing=TRUE)]
  ids <- as.character(host_data$Repeat_group_ID)
  for (file in files) {
    cached <- tryCatch(readRDS(file),error=function(e) NULL)
    if (is.null(cached) || !identical(cached$MC,mc) ||
        !setequal(as.character(cached$input_md5),as.character(input_md5))) next
    d <- cached$host$data
    needed <- c("Repeat_group_ID","Host_type","TS",mc,"NMDS1","NMDS2","NMDS3")
    if (is.null(d) || !all(needed %in% names(d)) ||
        !identical(as.character(d$Repeat_group_ID),ids)) next
    same <- vapply(c("Host_type","TS",mc),function(column) {
      if (column %in% mc) isTRUE(all.equal(d[[column]],host_data[[column]],tolerance=1e-12)) else
        identical(as.character(d[[column]]),as.character(host_data[[column]]))
    },logical(1))
    xyz <- d[,c("NMDS1","NMDS2","NMDS3"),drop=FALSE]
    fit <- cached$host$fit
    if (!all(same) || !all(is.finite(as.matrix(xyz))) ||
        !identical(as.integer(fit$ndim),3L) || length(fit$stress)!=1L || !is.finite(fit$stress)) next
    ordination <- cached$host
    ordination$data <- cbind(host_data,xyz)
    ordination$cache_source <- file
    message("Reusing saved 3D NMDS coordinates: ",file)
    return(ordination)
  }
  NULL
}


export_interactive_nmds_3d <- function(host_data, host_nmds_2d, host_result, time_result,
  mc, seed, workers, trymax, figure_dir, selfcontained, reuse_saved=TRUE) {
  ordination <- if (reuse_saved) find_saved_nmds_3d(host_data,mc,
    unique(host_result$Input_MD5),figure_dir) else NULL
  if (is.null(ordination)) {
    message("Starting 3D NMDS (",workers," workers; up to ",trymax," random starts).")
    ordination <- fit_interactive_nmds_3d(host_data,mc,seed,workers,trymax)
  }
  saveRDS(list(host=ordination,MC=mc,input_md5=unique(host_result$Input_MD5),
    seed=seed,trymax=trymax), file.path(figure_dir,"05_NMDS_3D_coordinates.rds"))
  utils::write.csv(ordination$data,file.path(figure_dir,"05_NMDS_3D_coordinates.csv"),
    row.names=FALSE,fileEncoding="UTF-8")
  diagnostics <- data.frame(Dimensions=c(2L,3L),N_units=nrow(host_data),
    Stress=c(host_nmds_2d$fit$stress,ordination$fit$stress),
    Random_starts=c(host_nmds_2d$fit$tries,ordination$fit$tries),
    Repeated_best=c(host_nmds_2d$fit$converged,ordination$fit$converged),
    Stop_code=c(host_nmds_2d$fit$icause,ordination$fit$icause))
  utils::write.csv(diagnostics,file.path(figure_dir,"05_NMDS_2D_vs_3D.csv"),row.names=FALSE)
  if (ordination$fit$stress>host_nmds_2d$fit$stress)
    message("3D stress exceeds 2D stress: inspect convergence and consider more starts.")
  widgets <- list(
    "05_host_NMDS_3D"=make_interactive_nmds_3d(ordination,"Host_type",pv_hosts,host_result,mc),
    "06_time_NMDS_3D"=make_interactive_nmds_3d(ordination,"TS",pv_ts,time_result,mc))
  for (name in names(widgets)) {
    htmlwidgets::saveWidget(widgets[[name]],file=file.path(figure_dir,paste0(name,".html")),
      selfcontained=selfcontained,background="white",title="")
  }
  methods <- c("3D NMDS interactive figures; no figure titles and no ellipses.",
    "Both HOST and TS HTML: hover a legend entry to enlarge its points without fading any groups. Click to add/remove groups from a persistent multi-selection; unselected groups become grey (size 2, opacity 0.18). Show all resets selection.",
    "Only marker styles change; coordinates and camera stay fixed. Leaving a hovered legend entry restores the clicked selection. Without a selection all groups keep their original colours and size 3.",
    "Host and TS colours use exactly the same 3D coordinates. Axis units have equal spatial scale.",
    "Axis lines and tick labels are anchored to fixed 3D positions, rather than automatically switching box edges. Axis ranges stay fixed when groups are hidden.",
    "Turntable rotation keeps the vertical direction upright; orthographic projection avoids perspective size differences. Markers use size 3 and full opacity.",
    if (!is.null(ordination$cache_source)) paste("Reused saved 3D coordinates:",ordination$cache_source) else "A new 3D NMDS solution was fitted.",
    "One point is one Repeat_group mean. Hover includes original-denominator MC relative abundances.",
    "Distance: Bray-Curtis after square-root transformation of eight MC proportions; no renormalisation among the eight MCs.",
    "The 3D solution is fitted independently of the 2D solution; NMDS1 and NMDS2 may also change.",
    "PERMANOVA R-squared and BH P are read from the completed results, not recalculated from the plotted axes.",
    "Host PERMANOVA adjusts for TS; TS PERMANOVA adjusts for host. These ordination coordinates are unadjusted.",
    "Compare 2D/3D stress and repeated-best diagnostics in 05_NMDS_2D_vs_3D.csv.",
    "Repeated_best counts near-identical best solutions found across starts; Stop_code records the optimiser stopping condition.",
    if (selfcontained) "HTML files embed their resources and can be shared individually." else
      "Pandoc is unavailable: share each HTML together with its matching _files directory. These files work offline.",
    sprintf("PERMDISP P (raw / Lingoes, unadjusted), host: %.4g / %.4g; TS: %.4g / %.4g.",
      host_result$P_disp_raw,host_result$P_disp_lingoes,time_result$P_disp_raw,time_result$P_disp_lingoes))
  writeLines(enc2utf8(methods),file.path(figure_dir,"05_NMDS_3D_README.txt"),useBytes=TRUE)
  message("Interactive HTML saved: 05_host_NMDS_3D.html and 06_time_NMDS_3D.html")
  invisible(ordination)
}


save_plot <- function(name, plot) {
  plots[[name]] <<- plot
  ggsave(file.path(figure_dir, paste0(name, ".png")), plot,
         width = 12, height = 9.6, dpi = 300, bg = "white")
}
if (dir.exists(figure_dir)) stop("输出目录已存在，请稍后再运行。")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

# ===== 3. 宿主、采样批次：同一套坐标；颜色区分宿主或TS；统一圆点、细白边和细线椭圆 =====
host_nmds <- make_nmds(host_data)
host_colours <- pv_subset(pv_hosts, levels(host_data$Host_type))
save_plot("01_host_final", group_plot(host_nmds, "Host_type", "Selected MC composition by host",
  "Host", host_result, "Repeat_group means", colours = host_colours,
  letter_display = letter_displays[["Host_type"]]))
save_plot("02_time_final", group_plot(host_nmds, "TS", "Selected MC composition by sampling batch",
  "Sampling batch", time_result, "Repeat_group means",
  letter_display = letter_displays[["TS"]]))

# ===== 4. 干湿季：使用季节检验保留的Repeat_group重新计算坐标 =====
season_nmds <- make_nmds(season_data)
save_plot("03_season_final", group_plot(season_nmds, "Season", "Selected MC composition in dry and wet seasons",
  "Season", season_result, "Repeat_group means", colours = c(Dry = "#D55E00", Wet = "#0072B2"), use_bh = FALSE))

# ===== 5. 环境：只画原始主分析BH p≤0.05的指标，一项一张 =====
# 环境是连续变量：点的颜色表示环境数值，灰色椭圆的线型区分不同宿主。
significant_env <- env_results[env_results$Status == "ok" & is.finite(env_results$P_BH_raw) &
                               env_results$P_BH_raw <= 0.05, , drop = FALSE]
env_labels <- c(Salinity_psu = "Salinity (psu)", SST_C = "Temperature (deg C)", pH = "pH",
  DO_mg_L = "Dissolved oxygen (mg/L)", Chla_ug_L = "Chlorophyll a (ug/L)",
  DIN_mgN_L = "DIN (mg N/L)", PSM_mg_L = "Suspended matter (mg/L)", Turbidity_NTU = "Turbidity (NTU)")
env_nmds <- NULL
if (nrow(significant_env) == 0) {
  message("没有环境指标通过BH校正，仅输出宿主、时间、干湿季3张图。")
} else {
  env_ordinations <- list()
  for (variable in significant_env$Factor) {
    result <- significant_env[significant_env$Factor == variable, , drop = FALSE]
    env_data_i <- if(!is.null(season_env_saved$environment_inputs))
      as.data.frame(season_env_saved$environment_inputs[[variable]]) else env_data
    # Use exactly the observations used by this variable's saved test.
    stopifnot(nrow(env_data_i)==result$N_units,all(is.finite(env_data_i[[variable]])))
    env_nmds <- make_nmds(env_data_i)
    env_ordinations[[variable]] <- env_nmds
    ellipse <- ellipse_data(env_nmds$data,"Host_type")
    data <- env_nmds$data
    data$Env_value <- data[[variable]]
    caption <- paste0("Host-adjusted PERMANOVA: raw R2 = ", sprintf("%.3f", result$R2_raw),
      "; BH p = ", p_text(result$P_BH_raw), " (", stars(result$P_BH_raw), "). Lingoes BH p = ",
      p_text(result$P_BH_lingoes), ".\nBH adjustment uses all 8 environmental tests. Colour = measured environmental value; grey ",
      round(100 * ellipse_level), "% data ellipses = hosts (line type and point shape).\nNMDS coordinates are not host-adjusted; ellipses are not significance tests.",
      sprintf("\nStress = %.3f; n = %d Repeat_group means across %d batches.", env_nmds$fit$stress, nrow(data), length(unique(data$TS))))
    save_plot(paste0("04_environment_", variable, "_final"),
      ggplot(data, aes(NMDS1, NMDS2, colour = Env_value, shape = Host_type)) +
        stat_ellipse(data = ellipse$data, aes(NMDS1, NMDS2, group = Host_type, linetype = Host_type), inherit.aes = FALSE,
                     type = "norm", level = ellipse_level, colour = "grey50", linewidth = 0.6,
                     show.legend = TRUE) +
        geom_point(size = 3.6, alpha = 1, stroke = 0.9) + coord_equal() +
        scale_colour_gradientn(colours=pv_continuous, limits=range(env_data_i[[variable]], na.rm=TRUE),
          guide=guide_colourbar(barheight=grid::unit(55,"mm"), barwidth=grid::unit(6,"mm"),
                                ticks.colour="black", frame.colour="black")) +
        scale_shape_manual(values = pv_subset(pv_host_shapes, levels(env_data_i$Host_type)),
                           limits = levels(env_data_i$Host_type), drop = FALSE) +
        scale_linetype_manual(values = rep(1:6, length.out = nlevels(env_data_i$Host_type)),
                              limits = levels(env_data_i$Host_type), drop = FALSE) +
        labs(title = NULL,
          subtitle = sprintf("%d selected MCs; stress = %.3f; %d Repeat_group means across %d batches", length(mc_cols), env_nmds$fit$stress,
                             nrow(data), length(unique(data$TS))), colour = unname(env_labels[variable]), shape = "Host", linetype = "Host",
          caption = wrap_caption(caption)) + plot_theme)
  }
}

# ===== 6. Interactive 3D NMDS: same MC data and distance definition, new k=3 fit =====
if (make_interactive_3d) {
  interactive_nmds_3d <- export_interactive_nmds_3d(host_data,host_nmds,host_result,time_result,
    mc_cols,random_seed,nmds_3d_workers,nmds_3d_trymax,figure_dir,html_3d_selfcontained,
    reuse_saved=nmds_3d_reuse_saved)
}

# ===== 7. 保存二维PDF和说明 =====
# Cairo supports system fonts (including Arial); base pdf() uses its own font database.
# If Cairo is unavailable, use the built-in sans family only for PDF export.
use_cairo_pdf <- isTRUE(capabilities("cairo"))
pdf_device <- if (use_cairo_pdf) grDevices::cairo_pdf else grDevices::pdf
pdf_device(file.path(figure_dir, "final_results.pdf"), width = 12, height = 9.6,
           onefile = TRUE, family = if (use_cairo_pdf) "Arial" else "sans")
pdf_device_id <- grDevices::dev.cur()
tryCatch({
  for (plot in plots) {
    pdf_plot <- if (use_cairo_pdf) plot else
      plot + theme(text = element_text(family = "sans"))
    print(pdf_plot)
  }
}, finally = {
  if (pdf_device_id %in% grDevices::dev.list()) grDevices::dev.off(pdf_device_id)
})
guide <- c("本目录只包含最终展示图：宿主、采样批次、干湿季和显著环境指标。", paste("分析范围：", analysis_scope),
  "整体R²和p来自02/03。本脚本只重算NMDS作图坐标；宿主/TS总览不显示两两比较字母。",
  "宿主/TS总览使用统一圆点、细白边；颜色分别表示宿主/批次。椭圆线宽0.4，透明度0.65，无填充。",
  "总览仅显示stress、n、R²和BH p；完整英文图注见overview_figure_captions.txt。",
  "椭圆为95%数据分布椭圆，使用正态近似；不是组中心的95%置信区间，也不保证恰好圈住95%的点。",
  "椭圆大小、重叠程度及有无椭圆，均不能替代PERMANOVA或两两检验。",
  "不足4点或坐标共线的组不画椭圆，保留散点；遗漏组在图注/本说明中列出。",
  ellipse_notes,
  "环境图：点颜色表示连续环境值，点形状和灰色椭圆的线型区分宿主。",
  "环境指标仍以全部8项检验的原始BH p≤0.05筛选，不在筛选后重新校正。",
  paste("本次显著环境指标：", if (nrow(significant_env)) paste(significant_env$Factor, collapse = ", ") else "无"),
  "环境是批次层面的探索性关联；Repeat_group行共享批次环境记录，时间混杂未被消除。",
  "不同图使用各自检验的数据子集和分析单位，请勿直接按R²跨图排名。",
  "https://ggplot2.tidyverse.org/reference/stat_ellipse.html")
if (show_pairwise_letters) {
  letter_guide <- do.call(rbind, lapply(letter_displays, function(x) x$table))
  utils::write.csv(letter_guide, file.path(figure_dir, "pairwise_legend_letters.csv"),
                    row.names = FALSE, fileEncoding = "UTF-8")
  for (effect in names(letter_displays)) {
    display <- letter_displays[[effect]]
    if (length(display$excluded)) guide <- c(guide,
      paste(effect, "未标注字母的组：", paste(display$excluded, collapse = "；")))
  }
}
writeLines(enc2utf8(unname(overview_captions)), file.path(figure_dir, "overview_figure_captions.txt"), useBytes=TRUE)
writeLines(enc2utf8(guide), file.path(figure_dir, "绘图说明.txt"), useBytes = TRUE)
saveRDS(list(host = host_nmds, season = season_nmds, environment = env_nmds,
             environment_by_variable = if(exists("env_ordinations")) env_ordinations else list(),
             MC = mc_cols, input_md5 = unique(host_result$Input_MD5)),
        file.path(result_dir, "plot_ordinations.rds"))
message("已输出", length(plots), "张最终PNG和同页数的final_results.pdf：", figure_dir)
if (make_interactive_3d) message("另有两张3D交互HTML、三维坐标和2D/3D stress对照，保存在同一目录。")
