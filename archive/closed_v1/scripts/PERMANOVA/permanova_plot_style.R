# Shared display settings only; no data transformation or statistical analysis.
# Continuous: viridis (purple-blue-green-yellow), linear numeric mapping.
pv_continuous <- viridisLite::viridis(256, option = "D")
pv_hosts <- c("Acropora sp."="#0072B2", "Montipora sp."="#D55E00",
  "Platygyra sp."="#009E73", "CCA"="#CC79A7", "Diadema sp."="#E69F00",
  "Drupella sp."="#332288", "Haliclona sp."="#44AA99",
  "Holothuria sp."="#882255", "Perna sp."="#56B4E9",
  "Sea anemone"="#999933", "Shrimp"="#661100", "Trochus sp."="#555555",
  "Seawater"="#005F73", "Sediment"="#9B5DE5")
pv_sample_types <- c(HOST="#0072B2",Seawater="#009E73",Sediment="#D55E00")
pv_ts <- setNames(c("#0072B2", "#D55E00", "#009E73", "#CC79A7",
  "#E69F00", "#332288", "#56B4E9", "#333333"), paste0("TS", 1:8))
pv_season <- c(Dry="#D55E00", Wet="#0072B2")
pv_mc <- c(MC6="#0072B2", MC507="#E69F00", MC21="#009E73", MC520="#D55E00",
  MC10="#56B4E9", MC15="#CC79A7", MC0="#555555", MC25="#332288")
pv_host_shapes <- setNames(c(16,17,15,18,3,4,8,0,1,2,5,6,7,9), names(pv_hosts))
pv_ts_shapes <- setNames(c(16,17,15,18,3,4,8,0), names(pv_ts))
pv_subset <- function(palette, groups) {
  missing <- setdiff(groups, names(palette))
  if (length(missing)) stop("Add display settings for: ", paste(missing, collapse=", "))
  palette[groups]
}
# Match the heatmap's Lab-interpolated scale before choosing contrasting text.
pv_text_colour <- function(values, limits, colours=pv_continuous) {
  hex <- scales::col_numeric(colours, domain=limits)(values)
  rgb <- t(grDevices::col2rgb(hex))/255
  linear <- ifelse(rgb <= 0.04045, rgb/12.92, ((rgb+0.055)/1.055)^2.4)
  luminance <- drop(linear %*% c(0.2126,0.7152,0.0722))
  ifelse(luminance > 0.179, "black", "white")
}
pv_theme <- function(base_size=14) {
  ggplot2::theme_classic(base_size=base_size, base_family="Arial") +
    ggplot2::theme(plot.title=ggplot2::element_blank(), plot.subtitle=ggplot2::element_blank(),
      axis.text=ggplot2::element_text(colour="black"),
      axis.line=ggplot2::element_line(colour="black", linewidth=0.5),
      axis.ticks=ggplot2::element_line(colour="black", linewidth=0.45),
      legend.key=ggplot2::element_rect(fill="white", colour=NA),
      legend.text=ggplot2::element_text(size=12),
      plot.caption=ggplot2::element_text(hjust=0, size=10, lineheight=1.15),
      plot.margin=ggplot2::margin(12,18,12,12))
}
