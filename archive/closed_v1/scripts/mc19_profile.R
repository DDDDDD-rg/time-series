# Explicit workflow profiles. No analysis or file writes on source().
mc19_profile <- function() {
  id <- getOption("mc19.pipeline", "TS8_env5")
  if (length(id) != 1L || !id %in% c("TS8_env5", "TS7_temperature", "TS7_env6"))
    stop("mc19.pipeline must be TS8_env5, TS7_temperature or TS7_env6.")
  temperature <- id %in% c("TS7_temperature", "TS7_env6")
  incremental <- identical(id, "TS7_temperature") # Historical five-indicator reuse profile.
  n_ts <- if (temperature) 7L else 8L
  fields <- c("Nitrate_mgN_L", "Nitrite_mgN_L", "Ammonia_mgN_L",
    if(!incremental) "DIN_mgN_L", "DIP_mgP_L", if(temperature) "SST_C")
  list(id=id, temperature=temperature, incremental=incremental, n_TS=n_ts, ts=paste0("TS", seq_len(n_ts)),
    fields=fields, n_env=length(fields), n_models=19L*length(fields),
    tag=switch(id,TS8_env5="env5_TS8",TS7_temperature="env5_TS7_temperature_incremental",TS7_env6="env6_TS7_full"),
    output_tag=switch(id,TS8_env5="",TS7_temperature="_TS7_temp",TS7_env6="_TS7_env6"))
}

mc19_method_version <- function(scope, kind="GLMM") {
  p <- mc19_profile()
  switch(kind, GLMM=paste0("MC19_",p$tag,"_",scope,"_Wald_v1"),
    PERMANOVA=paste0("MC19_",p$tag,"_",scope,"_individual_TSlabels_v1_vegan"),
    two_part=paste0("MC19_two_part_",p$tag,"_v1"),
    talk=paste0("MC19_talk_",p$tag,"_v1"), stop("Unknown analysis kind."))
}

# Explicit seven-TS selection precedes environment matching. No missing-data fallback.
mc19_select_samples <- function(metadata, community) {
  p <- mc19_profile()
  if (!p$temperature) return(list(metadata=metadata, community=community))
  if (anyNA(metadata$TS) || any(!as.character(metadata$TS) %in% paste0("TS",1:8)))
    stop("Metadata contains an invalid TS.")
  metadata <- metadata[as.character(metadata$TS) %in% p$ts,,drop=FALSE]
  community <- community[as.character(community$New_ID) %in% as.character(metadata$New_ID),,drop=FALSE]
  list(metadata=metadata, community=community)
}

mc19_temperature_source <- function(project_dir)
  file.path(project_dir,"data","water_quality","Bluff Island_PM11_marine_water_quality.csv")

mc19_verify_temperature <- function(env, project_dir) {
  if (!mc19_profile()$temperature) return(invisible(NULL))
  path <- mc19_temperature_source(project_dir)
  raw <- utils::read.csv(path,check.names=FALSE,stringsAsFactors=FALSE,fileEncoding="UTF-8-BOM")
  column <- grep("^Temperature",names(raw),value=TRUE)
  if (length(column)!=1L || !all(c("Station","Dates","Depth") %in% names(raw)))
    stop("PM11 source columns are missing or ambiguous.")
  dates <- as.Date(raw$Dates,format="%Y/%m/%d")
  take <- raw$Station=="PM11" & raw$Depth=="Surface Water" & dates %in% as.Date(env$Water_quality_date)
  source <- raw[take,,drop=FALSE]; source_dates <- dates[take]
  if (nrow(source)!=nrow(env) || anyDuplicated(source_dates)) stop("PM11 temperature date/depth match failed.")
  value <- suppressWarnings(as.numeric(source[[column]][match(as.Date(env$Water_quality_date),source_dates)]))
  if (any(!is.finite(value)) || any(abs(value-env$SST_C)>1e-10))
    stop("Environment_TS SST_C differs from the PM11 surface-water source.")
  invisible(data.frame(TS=env$TS,Water_quality_date=env$Water_quality_date,
    Station="PM11",Depth="Surface Water",SST_C=value,Source=normalizePath(path,winslash="/")))
}
