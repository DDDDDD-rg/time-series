# Shared input rules for the current 8-TS, five-nutrient analysis. Definitions only.
sys.source(file.path(getOption("mc19.project_dir", "C:/Users/xzy/Desktop/CUHK_Study/GNBF6010_Research_Projects/Research/Msc main project-time series"), "scripts", "mc19_profile.R"), envir=environment())
mc19_ts_levels <- mc19_profile()$ts
mc19_env_fields <- mc19_profile()$fields
mc19_primary_env_fields <- c("Nitrate_mgN_L", "Nitrite_mgN_L", "Ammonia_mgN_L", "DIP_mgP_L")

mc19_prepare_environment <- function(env) {
  required <- c("TS", "Sampling_month", "Water_quality_date", "Season",
    paste0(mc19_primary_env_fields, "_raw"))
  if (!all(required %in% names(env))) stop("Environment_TS缺少字段：", paste(setdiff(required, names(env)), collapse=", "))
  env$TS <- as.character(env$TS)
  if (mc19_profile()$temperature) env <- env[env$TS %in% mc19_ts_levels,,drop=FALSE]
  if (anyNA(env$TS) || anyDuplicated(env$TS) || !setequal(env$TS, mc19_ts_levels))
    stop("Environment_TS必须每个TS一行，完整包含所选TS：",paste(mc19_ts_levels,collapse=", "))
  env <- env[match(mc19_ts_levels, env$TS), , drop=FALSE]
  missing <- character()
  for (field in mc19_primary_env_fields) {
    raw <- trimws(as.character(env[[paste0(field, "_raw")]]))
    absent <- is.na(raw) | !nzchar(raw)
    censored <- !absent & grepl("^<", raw)
    value <- suppressWarnings(as.numeric(sub("^<\\s*", "", raw)))
    invalid <- !absent & (!is.finite(value) | value < 0 | (censored & value <= 0))
    if (any(invalid)) stop(field, "原始报告值无效：", paste(env$TS[invalid], collapse=", "))
    value[censored] <- value[censored] / 2
    value[absent] <- NA_real_
    env[[field]] <- value
    env[[paste0(field, "_raw")]] <- raw
    env[[paste0(field, "_censored")]] <- censored
    if (any(absent)) missing <- c(missing, paste0(env$TS[absent], ":", field, "_raw"))
  }
  # All three components are expressed as mg N/L. Missing is never a zero.
  n_fields <- c("Nitrate_mgN_L", "Nitrite_mgN_L", "Ammonia_mgN_L")
  env$DIN_mgN_L <- rowSums(as.matrix(env[n_fields]))
  env$DIN_mgN_L_raw <- NA_character_ # Calculated DIN is not a directly reported measurement.
  env$DIN_mgN_L_censored <- rowSums(as.matrix(env[paste0(n_fields, "_censored")])) > 0
  if (mc19_profile()$temperature) {
    if (!"SST_C_raw" %in% names(env)) stop("Environment_TS lacks SST_C_raw.")
    raw <- trimws(as.character(env$SST_C_raw))
    value <- suppressWarnings(as.numeric(raw))
    if (anyNA(raw) || any(!nzchar(raw)) || any(!is.finite(value)) || any(value < -2 | value > 45))
      stop("Invalid or missing temperature in the selected TS; use observed degrees C.")
    env$SST_C <- value
    env$SST_C_raw <- raw
    env$SST_C_censored <- rep(FALSE,nrow(env))
  }
  date <- as.Date(env$Water_quality_date)
  absent_date <- is.na(date)
  if (any(absent_date)) missing <- c(missing, paste0(env$TS[absent_date], ":Water_quality_date"))
  if (length(missing)) stop(if(mc19_profile()$temperature) "7TS温度流程输入尚未填齐：" else "8TS输入尚未填齐：", paste(missing, collapse="; "),
    "。请填写Environment_TS的原始报告值和实际水样日期；程序不会回退为7TS。", call.=FALSE)
  if (anyNA(env$Sampling_month) || anyNA(env$Season) ||
      any(format(date, "%Y-%m") != as.character(env$Sampling_month)))
    stop("Environment_TS的实际水样日期与Sampling_month不一致。")
  env$Water_quality_date <- format(date, "%Y-%m-%d")
  env$Env5_complete <- as.integer(complete.cases(env[setdiff(mc19_env_fields,"SST_C")]))
  env$Environment_complete <- as.integer(complete.cases(env[mc19_env_fields]))
  env
}

mc19_match_environment <- function(metadata, env) {
  metadata$TS <- as.character(metadata$TS)
  index <- match(metadata$TS, env$TS)
  if (anyNA(index)) stop("Metadata包含无法匹配的TS。")
  for (field in c("Sampling_month", "Season")) {
    if (!field %in% names(metadata) || anyNA(metadata[[field]]) ||
        any(as.character(metadata[[field]]) != as.character(env[[field]][index])))
      stop("Metadata与Environment_TS的", field, "不一致。")
  }
  # Environment_TS raw reports are authoritative; no duplicate manual updates.
  for (field in c(mc19_env_fields, "Water_quality_date", "Env5_complete", "Environment_complete")) metadata[[field]] <- env[[field]][index]
  metadata
}
