# MC19 / coral or ecosystem / 原始独立样本：数据读取和统计工具。
# 本文件只定义设置与函数；由03调用，不读取旧01或Repeat_group流程。
project_dir <- "C:/Users/xzy/Desktop/CUHK_Study/GNBF6010_Research_Projects/Research/Msc main project-time series"
sys.source(file.path(project_dir, "scripts", "mc19_run_io.R"), envir=environment())
input_file <- file.path(project_dir, "data", "PERMANOVA_input_MC19_env7_20261026.xlsx")
analysis_scope <- getOption("mc19.analysis_scope", "coral")
stopifnot(length(analysis_scope) == 1L, analysis_scope %in% c("coral", "ecosystem"))
analysis_version <- paste0("MC19_env7_", analysis_scope, "_individual_TSlabels_v3_vegan")
compatible_analysis_versions <- c(analysis_version, if (analysis_scope == "coral")
  c("MC19_env7_coral_individual_TSlabels_v2_vegan", "MC19_env7_coral_individual_TSlabels_v1"))
scope_mask <- function(metadata, scope=analysis_scope) {
  niche <- tolower(trimws(as.character(metadata$Host_niche)))
  host <- tolower(trimws(as.character(metadata$Host_type)))
  biological <- !is.na(niche) & nzchar(niche) & !niche %in% c("env", "seawater", "sediment") &
    !is.na(host) & !host %in% c("seawater", "sediment")
  if (scope == "coral") biological & niche == "coral" else biological
}
mc19_scope_pointer <- function() paste0("00_latest_PERMANOVA_MC19_env7_", analysis_scope, ".txt")
mc_cols <- c("MC41", "MC95", "MC6", "MC507", "MC21", "MC451", "MC98", "MC2s1",
             "MC43", "MC25", "MC77", "MC376", "MC35", "MC472", "MC357", "MC80",
             "MC0", "MC10", "MC15")
env_cols <- c("Salinity_psu", "SST_C", "pH", "Nitrate_mgN_L", "Nitrite_mgN_L",
              "Ammonia_mgN_L", "DIP_mgP_L")
env_labels <- c(Salinity_psu="盐度", SST_C="水温", pH="pH", Nitrate_mgN_L="硝酸盐氮",
                Nitrite_mgN_L="亚硝酸盐氮", Ammonia_mgN_L="氨氮", DIP_mgP_L="DIP")
ts_levels <- paste0("TS", 1:8)

require_mc19_packages <- function() {
  packages <- c("readxl", "vegan")
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly=TRUE)]
  if (length(missing)) stop("请先安装R包：", paste(missing, collapse=", "))
}

read_mc19_samples <- function(path=input_file) {
  require_mc19_packages()
  if (!file.exists(path)) stop("输入文件不存在：", path)
  community <- as.data.frame(readxl::read_excel(path, sheet="Community"))
  metadata <- as.data.frame(readxl::read_excel(path, sheet="Metadata"))
  environment <- as.data.frame(readxl::read_excel(path, sheet="Environment_TS"))
  if (!identical(names(community), c("New_ID", mc_cols))) {
    stop("Community必须按约定顺序包含New_ID和19个MC；MC2S1使用原始拼写MC2s1。")
  }
  required <- c("New_ID", "Sample_ID", "Host_type", "Host_niche", "TS", "Sampling_month",
                "Season", "Is_coral", "Env7_complete", env_cols)
  if (!all(required %in% names(metadata))) stop("Metadata缺少必要字段。")
  if (!all(c("TS", "Sampling_month", "Season", "Water_quality_date", env_cols,
             paste0(env_cols, "_raw")) %in% names(environment))) stop("Environment_TS缺少必要字段。")
  stopifnot(!anyNA(community$New_ID), !anyNA(metadata$New_ID),
            !anyDuplicated(community$New_ID), !anyDuplicated(metadata$New_ID),
            setequal(community$New_ID, metadata$New_ID),
            !anyNA(metadata$TS), all(metadata$TS %in% ts_levels),
            !anyNA(environment$TS), !anyDuplicated(environment$TS))
  community <- community[match(metadata$New_ID, community$New_ID), , drop=FALSE]
  stopifnot(identical(community$New_ID, metadata$New_ID),
            all(vapply(community[mc_cols], is.numeric, logical(1))),
            all(vapply(metadata[env_cols], is.numeric, logical(1))))
  y_all <- as.matrix(community[mc_cols])
  if (anyNA(y_all) || any(!is.finite(y_all)) || any(y_all < 0 | y_all > 1)) {
    stop("MC相对丰度必须为有限的0–1数值；不重新归一化。")
  }
  # 每个样本沿用其TS对应的环境值；不对MC或环境行求分组均值。
  env_index <- match(metadata$TS, environment$TS)
  if (anyNA(env_index)) stop("有样本的TS在Environment_TS中找不到。")
  for (variable in env_cols) {
    a <- metadata[[variable]]; b <- environment[[variable]][env_index]
    same <- (is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & abs(a-b) <= 1e-12)
    if (!all(same)) stop(variable, "：Metadata与Environment_TS不一致。")
    raw <- trimws(as.character(environment[[paste0(variable, "_raw")]]))
    censored <- !is.na(raw) & grepl("^<", raw)
    parsed <- suppressWarnings(as.numeric(sub("^<\\s*", "", raw)))
    parsed[censored] <- parsed[censored] / 2
    actual <- environment[[variable]]
    same_raw <- (is.na(actual) & is.na(parsed)) |
      (!is.na(actual) & !is.na(parsed) & abs(actual-parsed) <= 1e-12)
    if (!all(same_raw)) stop(variable, "：数值与原始报告值/L/2规则不一致。")
  }
  if (!all(metadata$Sampling_month == environment$Sampling_month[env_index]) ||
      !all(metadata$Season == environment$Season[env_index])) stop("TS月份/季节匹配不一致。")
  is_coral <- !is.na(metadata$Host_niche) & metadata$Host_niche == "coral"
  stopifnot(!anyNA(metadata$Is_coral), all(metadata$Is_coral == as.integer(is_coral)))
  complete <- complete.cases(metadata[env_cols]) &
    apply(as.matrix(metadata[env_cols]), 1L, function(x) all(is.finite(x)))
  stopifnot(!anyNA(metadata$Env7_complete), all(metadata$Env7_complete == as.integer(complete)))
  in_scope <- scope_mask(metadata)
  keep <- in_scope & complete
  selected <- metadata[keep, , drop=FALSE]
  selected <- selected[, setdiff(names(selected), c("Repeat_group_ID", "Repeat_group_n", "Same_sample_n")), drop=FALSE]
  y <- y_all[keep, , drop=FALSE]
  if (nrow(selected) <= 2L || any(rowSums(y) == 0)) stop(analysis_scope, "样本不足，或存在19个MC全零的样本。")
  rownames(y) <- as.character(selected$New_ID)
  batches <- environment[environment$TS %in% selected$TS, , drop=FALSE]
  batches <- batches[order(match(batches$TS, ts_levels)), , drop=FALSE]
  if (nrow(batches) < 3L || nrow(batches) > 8L) stop("本版仅支持3–8个有效TS的完整枚举。")
  row_batch <- match(selected$TS, batches$TS)
  coverage <- data.frame(TS=ts_levels,
    N_scope_samples=as.integer(table(factor(metadata$TS[in_scope], levels=ts_levels))),
    N_environment_complete=as.integer(table(factor(selected$TS, levels=ts_levels))))
  coverage$N_excluded_missing_environment <- coverage$N_scope_samples - coverage$N_environment_complete
  list(scope=analysis_scope, metadata=selected, community=y, batches=batches, row_batch=row_batch,
       coverage=coverage, input_md5=unname(tools::md5sum(path)))
}

# Compatibility name retained for external callers; selection follows analysis_scope.
read_coral_samples <- read_mc19_samples

enumerate_ts_orders <- function(n) {
  stopifnot(n >= 3L, n <= 8L)
  recurse <- function(x) {
    if (length(x) == 1L) return(matrix(x, nrow=1L))
    do.call(rbind, lapply(seq_along(x), function(i) cbind(x[i], recurse(x[-i]))))
  }
  answer <- recurse(seq_len(n))
  storage.mode(answer) <- "integer"
  stopifnot(nrow(answer) == factorial(n), identical(as.integer(answer[1L, ]), seq_len(n)),
            !anyDuplicated(as.data.frame(answer)))
  answer
}

# 新输出只在顶层保留主结果表，其余分析文件集中在supporting_files。
# 读取时兼容旧版平铺目录；已有结果无需移动或重新分析。
mc19_support_dir <- function(result_dir) file.path(result_dir, "supporting_files")

mc19_support_file <- function(result_dir, name) {
  current <- file.path(mc19_support_dir(result_dir), name)
  legacy <- file.path(result_dir, name)
  if (file.exists(current)) return(current)
  if (file.exists(legacy)) return(legacy)
  current
}

make_coral_distance <- function(y) {
  vegan::vegdist(sqrt(y), method="bray")
}

fit_environment_vegan <- function(distance, env_value) {
  model_data <- data.frame(Env_value=env_value)
  # 每一种TS标签分配均由vegan拟合；这里只读取统计量，不让vegan置换样本行。
  fit <- vegan::adonis2(distance ~ Env_value, data=model_data,
                       permutations=0, by=NULL, parallel=1)
  required_rows <- c("Model", "Residual", "Total")
  required_columns <- c("Df", "SumOfSqs", "R2", "F")
  if (!all(required_rows %in% rownames(fit)) || !all(required_columns %in% names(fit))) {
    stop("vegan::adonis2返回的统计表缺少预期行或列。")
  }
  statistics <- c(SS_total=fit["Total", "SumOfSqs"],
    SS_environment=fit["Model", "SumOfSqs"], SS_residual=fit["Residual", "SumOfSqs"],
    Df_environment=fit["Model", "Df"], Df_residual=fit["Residual", "Df"],
    F_raw=fit["Model", "F"], R2_raw=fit["Model", "R2"])
  if (any(!is.finite(statistics)) || statistics["SS_total"] <= 0 ||
      statistics["SS_residual"] <= 0 || statistics["Df_environment"] != 1L ||
      statistics["Df_residual"] != length(env_value)-2L) {
    stop("vegan::adonis2未返回有效的单环境因子统计量。")
  }
  statistics
}

fit_one_environment <- function(variable, input, distance, orders) {
  values <- input$batches[[variable]]
  if (any(!is.finite(values)) || length(unique(values)) < 2L) stop(variable, "无有效变化。")
  # 所有指标共用同一套TS排列。先给TS分配值，再扩展到原始样本行。
  assigned <- matrix(values[as.vector(orders)], nrow=nrow(orders), ncol=ncol(orders))
  f <- numeric(nrow(orders))
  observed <- NULL
  for (i in seq_len(nrow(orders))) {
    env_value <- assigned[i, input$row_batch]
    statistics <- fit_environment_vegan(distance, env_value)
    f[i] <- unname(statistics["F_raw"])
    if (i == 1L) observed <- statistics
    if (i %% 500L == 0L) message("  ", variable, "：已拟合 ", i, "/", nrow(orders), " 种TS分配。")
  }
  # 第一行是实际分配；保留全部排列和重复环境值的重数，P值规则与v1完全相同。
  tolerance <- 1e-10 * max(1, abs(f[1L]))
  p <- sum(f >= f[1L]-tolerance) / length(f)
  raw <- as.character(input$batches[[paste0(variable, "_raw")]])
  result <- data.frame(Scope=analysis_scope, Factor=variable, Indicator=unname(env_labels[variable]),
    N_samples=nrow(input$community), N_MC=ncol(input$community), N_TS=nrow(input$batches),
    TS_used=paste(input$batches$TS, collapse=";"),
    SS_total=unname(observed["SS_total"]), SS_environment=unname(observed["SS_environment"]),
    SS_residual=unname(observed["SS_residual"]),
    Df_environment=unname(observed["Df_environment"]), Df_residual=unname(observed["Df_residual"]),
    F_raw=unname(observed["F_raw"]), R2_raw=unname(observed["R2_raw"]), P_raw=p,
    N_TS_orders_including_observed=nrow(orders),
    N_unique_value_assignments=nrow(unique(as.data.frame(assigned))),
    N_censored_TS=sum(!is.na(raw) & grepl("^<", raw)),
    Analysis_unit="Original individual sample", Model=paste("BrayCurtis ~", variable),
    Permutation_design="Whole TS environmental labels; exhaustive; observed included",
    Statistics_engine="vegan::adonis2(permutations=0)",
    Analysis_version=analysis_version, Input_MD5=input$input_md5, Status="ok",
    stringsAsFactors=FALSE)
  list(result=result, F_per_order=f)
}

write_csv_excel <- function(data, path) {
  con <- file(path, open="wb")
  on.exit(close(con), add=TRUE)
  writeBin(as.raw(c(239,187,191)), con)
  text <- utils::capture.output(utils::write.table(data, sep=",", row.names=FALSE,
    col.names=TRUE, na="", qmethod="double"))
  writeBin(charToRaw(enc2utf8(paste0(paste(text, collapse="\n"), "\n"))), con)
}
