# 七TS普通联合GLMM：19个MC，各自同时考虑四项营养盐和温度，不加DIN。
# 默认只读检查：options(mc19.joint.stage = "check")；拟合改为 "fit"。
# 范围：options(mc19.analysis_scope = "coral")；另一个范围为 "ecosystem"。
# 从 scripts.Rproj 打开后：source("GLMM/GLMM_MC19_joint_TS7.R")。
# 不调用八TS主入口；不拟合贝叶斯模型，不筛选环境，不修改旧结果。

# [1] 固定本轮目标：TS1–7、19MC、NO3/NO2/TAN/DIP和温度。
joint_TS7_config <- function() {
  override <- getOption("mc19.joint.project_dir", "")
  starts <- getwd()
  for (f in sys.frames()) if (exists("ofile", envir=f, inherits=FALSE)) {
    path <- get("ofile", envir=f, inherits=FALSE)
    if (is.character(path) && length(path)==1L) starts <- c(starts, dirname(path))
  }
  cli <- grep("^--file=", commandArgs(FALSE), value=TRUE)
  if (length(cli)) starts <- c(starts, dirname(sub("^--file=", "", cli[1])))
  candidates <- override[nzchar(override)]
  for (start in starts) {
    current <- normalizePath(start, winslash="/", mustWork=FALSE)
    for (level in 1:6) {
      candidates <- c(candidates, current, file.path(current, "Research", "Msc main project-time series"))
      current <- dirname(current)
    }
  }
  valid <- vapply(candidates, function(p) file.exists(file.path(p, "data", "MC19_TS8_env5.xlsx")) &&
    file.exists(file.path(p, "data", "mc", "MC_abundance_absolute.tsv")), logical(1))
  if (!any(valid)) stop("找不到研究目录。请打开 scripts.Rproj，或设置 mc19.joint.project_dir。")
  project <- normalizePath(candidates[which(valid)[1]], winslash="/", mustWork=TRUE)
  scope <- getOption("mc19.analysis_scope", "coral")
  stage <- getOption("mc19.joint.stage", "check")
  stopifnot(length(scope)==1L, scope %in% c("coral", "ecosystem"),
    length(stage)==1L, stage %in% c("check", "fit"))
  list(project=project, scope=scope, stage=stage, ts=paste0("TS", 1:7),
    workbook=file.path(project, "data", "MC19_TS8_env5.xlsx"),
    counts=file.path(project, "data", "mc", "MC_abundance_absolute.tsv"),
    temperature_source=file.path(project,"data","water_quality","Bluff Island_PM11_marine_water_quality.csv"),
    mc=c("MC41", "MC95", "MC6", "MC507", "MC21", "MC451", "MC98", "MC2s1", "MC43",
         "MC25", "MC77", "MC376", "MC35", "MC472", "MC357", "MC80", "MC0", "MC10", "MC15"),
    env=c("Nitrate_mgN_L", "Nitrite_mgN_L", "Ammonia_mgN_L", "DIP_mgP_L", "SST_C"),
    z=c("Nitrate_z", "Nitrite_z", "Ammonia_z", "DIP_z", "Temperature_z"),
    run_prefix="joint_glmm_TS7_env5",
    expected_n=c(coral=98L, ecosystem=310L)[[scope]], version="MC19_joint_TS7_env5_v1")
}

joint_TS7_require <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly=TRUE)]
  if (length(missing)) stop("缺少R包：", paste(missing, collapse=", "), "；请自行安装后再运行。")
}
joint_TS7_fields <- function(data, fields, label) {
  if (!all(fields %in% names(data))) stop(label, "缺少字段：", paste(setdiff(fields, names(data)), collapse=", "))
}
joint_TS7_ids <- function(data, label) {
  data$New_ID <- trimws(as.character(data$New_ID))
  if (anyNA(data$New_ID) || any(!nzchar(data$New_ID)) || anyDuplicated(data$New_ID))
    stop(label, "的New_ID缺失或重复。")
  data
}

# [2] 从原始报告值解析营养盐。保留<L标记，数值计算沿用L/2。
joint_TS7_report <- function(x) {
  raw <- trimws(as.character(x))
  absent <- is.na(raw) | !nzchar(raw)
  if (any(absent)) stop("七TS营养盐原始值有缺失；不自动删TS或样本。")
  censored <- grepl("^<", raw)
  numeric_text <- sub("^<\\s*", "", raw)
  valid <- grepl("^(?:[0-9]+(?:\\.[0-9]*)?|\\.[0-9]+)(?:[eE][+-]?[0-9]+)?$", numeric_text, perl=TRUE)
  value <- suppressWarnings(as.numeric(numeric_text))
  if (any(!valid | !is.finite(value) | value<0 | (censored & value<=0)))
    stop("营养盐报告值格式无效：", paste(raw[!valid | !is.finite(value) | value<0 | (censored & value<=0)], collapse=", "))
  value[censored] <- value[censored]/2
  list(raw=raw, value=value, censored=censored)
}

# 温度是实际摄氏度，不使用营养盐L/2规则；允许合法的零值/负摄氏度。
joint_TS7_temperature <- function(x) {
  raw <- trimws(as.character(x))
  value <- suppressWarnings(as.numeric(raw))
  if (anyNA(raw) || any(!nzchar(raw)) || any(!is.finite(value)) || any(grepl("[<>]",raw)))
    stop("七TS温度必须为实际有限数值；不接受缺失、<L或近似替代。")
  list(raw=raw,value=value,censored=rep(FALSE,length(raw)))
}

joint_TS7_load <- function(cfg) {
  joint_TS7_require("readxl")
  paths <- c(cfg$workbook,cfg$counts,cfg$temperature_source)
  if (!all(file.exists(paths))) stop("缺少输入：",paste(paths[!file.exists(paths)],collapse=", "))
  fingerprint <- data.frame(path=paths,md5=unname(tools::md5sum(paths)))
  read_sheet <- function(sheet) as.data.frame(readxl::read_excel(cfg$workbook, sheet=sheet))
  meta <- read_sheet("Metadata"); community <- read_sheet("Community"); env <- read_sheet("Environment_TS")
  joint_TS7_fields(meta, c("New_ID", "Sample_ID", "TS", "Sampling_month", "Host_type", "Host_niche"), "Metadata")
  joint_TS7_fields(community, c("New_ID", cfg$mc), "Community")
  joint_TS7_fields(env, c("TS", "Sampling_month", "Water_quality_date", paste0(cfg$env, "_raw")), "Environment_TS")
  meta <- joint_TS7_ids(meta, "Metadata"); community <- joint_TS7_ids(community, "Community")
  if (!setequal(meta$New_ID, community$New_ID)) stop("Metadata与Community的样本集合不一致。")
  env$TS <- as.character(env$TS); meta$TS <- as.character(meta$TS)
  if (anyNA(meta$TS) || any(!meta$TS %in% paste0("TS",1:8))) stop("Metadata的TS标签缺失或超出TS1–TS8。")
  if (anyNA(env$TS) || anyDuplicated(env$TS) || !all(cfg$ts %in% env$TS)) stop("环境TS缺失或重复。")
  # 显式选择七TS后才检查环境完整性；TS8缺测不会被当作零或自动补齐。
  env <- env[match(cfg$ts, env$TS), , drop=FALSE]
  for (field in cfg$env) {
    report <- if (field=="SST_C") joint_TS7_temperature(env[[paste0(field,"_raw")]]) else
      joint_TS7_report(env[[paste0(field, "_raw")]])
    env[[field]] <- report$value
    env[[paste0(field, "_raw")]] <- report$raw
    env[[paste0(field, "_censored")]] <- report$censored
  }
  date <- as.Date(env$Water_quality_date)
  if (anyNA(date) || anyNA(env$Sampling_month) ||
      any(format(date, "%Y-%m") != as.character(env$Sampling_month))) stop("七TS水样日期与采样月份不匹配。")
  env$Water_quality_date <- format(date, "%Y-%m-%d")
  # 水温与相同日期的PM11表层原CSV逐行核对，正式输出保存来源快照。
  water <- utils::read.csv(cfg$temperature_source,check.names=FALSE,fileEncoding="UTF-8-BOM")
  joint_TS7_fields(water,c("Station","Dates","Depth"),"PM11水质CSV")
  temperature_field <- grep("^Temperature",names(water),value=TRUE)
  if (length(temperature_field)!=1L) stop("PM11 CSV温度字段缺失或不唯一。")
  water <- water[which(water$Station=="PM11" & water$Depth=="Surface Water"),,drop=FALSE]
  source_dates <- as.Date(water$Dates,format="%Y/%m/%d")
  source_rows <- vapply(date,function(day) {
    hit <- which(source_dates==day)
    if (length(hit)!=1L) stop("PM11表层同日温度记录缺失或重复。")
    hit
  },integer(1))
  source_temperature <- joint_TS7_temperature(water[[temperature_field]][source_rows])$value
  if (any(abs(source_temperature-env$SST_C)>1e-10)) stop("工作簿温度与PM11原CSV不一致。")
  temperature_source <- data.frame(TS=env$TS,Water_quality_date=env$Water_quality_date,
    Station=water$Station[source_rows],Depth=water$Depth[source_rows],SST_C=source_temperature)
  niche <- tolower(trimws(as.character(meta$Host_niche)))
  host <- tolower(trimws(as.character(meta$Host_type)))
  biological <- !is.na(niche) & nzchar(niche) & !niche %in% c("env", "seawater", "sediment") &
    !is.na(host) & nzchar(host) & !host %in% c("seawater", "sediment")
  in_scope <- biological & if (cfg$scope=="coral") niche=="coral" else TRUE
  keep <- !is.na(meta$TS) & meta$TS %in% cfg$ts & in_scope
  inclusion <- meta[, c("New_ID", "Sample_ID", "TS", "Host_type", "Host_niche")]
  inclusion$included <- keep
  inclusion$reason <- ifelse(!meta$TS %in% cfg$ts, "outside_TS1_TS7", ifelse(!in_scope, "outside_scope", "included"))
  selected <- meta[keep, c("New_ID", "Sample_ID", "TS", "Sampling_month", "Host_type", "Host_niche"), drop=FALSE]
  if (nrow(selected)!=cfg$expected_n || !setequal(selected$TS, cfg$ts)) stop("七TS样本数或TS集合改变，请先复核输入。")
  index <- match(selected$TS, env$TS)
  if (anyNA(selected$Sampling_month) || any(as.character(selected$Sampling_month)!=as.character(env$Sampling_month[index])))
    stop("样本与环境的采样月份不一致。")

  # [3] 分母重新加总原始175行，包含Unassigned；不把19MC重新闭合。
  raw <- utils::read.delim(cfg$counts, check.names=FALSE, stringsAsFactors=FALSE, quote="", comment.char="")
  joint_TS7_fields(raw, c("MC", selected$New_ID), "原始计数表")
  if (nrow(raw)!=175L || anyNA(raw$MC) || any(!nzchar(trimws(raw$MC))) || anyDuplicated(raw$MC) ||
      !all(c(cfg$mc, "Unassigned") %in% raw$MC)) stop("原始MC行集合改变或缺少目标MC。")
  counts <- as.matrix(raw[, selected$New_ID, drop=FALSE])
  rownames(counts) <- raw$MC
  if (!is.numeric(counts) || anyNA(counts) || any(!is.finite(counts)) ||
      any(counts<0) || any(counts!=floor(counts))) stop("Reads必须为非负整数。")
  totals <- colSums(counts)
  if (any(totals<=0)) stop("有样本Ruegeria总reads为零；请核查，不自动删除。")
  selected_counts <- counts[match(cfg$mc, raw$MC), , drop=FALSE]
  book <- as.matrix(community[match(selected$New_ID, community$New_ID), cfg$mc, drop=FALSE])
  proportions <- t(sweep(selected_counts, 2L, totals, "/"))
  if (!is.numeric(book) || anyNA(book) || max(abs(book-proportions))>1e-10)
    stop("原始reads重算的19MC比例与工作簿不一致。")
  selected$Ruegeria_reads_total <- unname(totals)

  # [4] 标准化只用七个TS，各TS等权；生物样本不求平均。
  scaling <- data.frame(environment=cfg$env, term=cfg$z, center=NA_real_, SD=NA_real_)
  for (j in seq_along(cfg$env)) {
    field <- cfg$env[j]; mu <- mean(env[[field]]); s <- stats::sd(env[[field]])
    if (!is.finite(s) || s<=0) stop(field, "没有可估计的变化。")
    scaling$center[j] <- mu; scaling$SD[j] <- s
    env[[cfg$z[j]]] <- (env[[field]]-mu)/s
    selected[[field]] <- env[[field]][index]
    selected[[cfg$z[j]]] <- env[[cfg$z[j]]][index]
  }
  selected$TS <- factor(selected$TS, levels=cfg$ts)

  # [5] 在TS层检查秩、相关性和VIF；只提示共线性，不自动删变量。
  X <- stats::model.matrix(stats::reformulate(cfg$z), data=env)
  design_rank <- qr(X)$rank
  if (design_rank<ncol(X)) stop("联合环境设计矩阵不满秩，不能分开估计五项系数。")
  correlation <- stats::cor(env[cfg$z])
  vif <- diag(solve(correlation))
  collinearity <- data.frame(environment=cfg$env, term=cfg$z, VIF=unname(vif),
    condition_number=base::kappa(X, exact=TRUE), design_rank=design_rank, design_columns=ncol(X), n_TS=7L)
  if (any(vif>5) || collinearity$condition_number[1]>30)
    warning("环境指标相关性较强：请结合VIF、系数区间和后续稳定性检查解释。不会自动删项。", call.=FALSE)
  if (!identical(unname(tools::md5sum(paths)),fingerprint$md5)) stop("读取过程中输入发生变化。")
  list(samples=selected, counts=selected_counts, environment_TS=env, scaling=scaling,
    collinearity=collinearity, correlation=correlation, inclusion=inclusion, fingerprints=fingerprint,
    temperature_source=temperature_source)
}

joint_TS7_capture <- function(fun) {
  warnings <- character()
  value <- tryCatch(withCallingHandlers(fun(), warning=function(w) {
    warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
  }), error=function(e)e)
  list(value=if (inherits(value,"error")) NULL else value, warnings=unique(warnings),
    error=if (inherits(value,"error")) conditionMessage(value) else "")
}

# [6] 一个MC只拟合一个模型；同时估计五个环境斜率。
joint_TS7_fit_one <- function(mc, d, cfg) {
  model_formula <- cbind(MC_reads, Other_reads) ~
    Nitrate_z + Nitrite_z + Ammonia_z + DIP_z + Temperature_z + (1 | TS)
  started <- proc.time()[[3]]
  attempt <- joint_TS7_capture(function() glmmTMB::glmmTMB(
    formula=model_formula, data=d, family=glmmTMB::betabinomial(link="logit"),
    ziformula=~0, dispformula=~1, REML=FALSE, na.action=stats::na.fail,
    control=glmmTMB::glmmTMBControl(parallel=1L,
      optCtrl=list(iter.max=10000L, eval.max=10000L))))
  fit <- attempt$value
  model_row <- data.frame(MC_ID=mc, n_samples=nrow(d), n_TS=nlevels(d$TS),
    zero_count=sum(d$MC_reads==0), convergence_code=NA_integer_, pdHess=FALSE,
    max_abs_gradient=NA_real_, TS_SD_logit=NA_real_, phi=NA_real_, logLik=NA_real_, AIC=NA_real_,
    numerical_ok=FALSE, review_flags="", diagnostic_status="residual_diagnostics_not_run",
    warnings=paste(attempt$warnings,collapse=" | "), error=attempt$error, elapsed_seconds=NA_real_)
  result <- data.frame(MC_ID=mc, environment=cfg$env, term=cfg$z, analysis_scope=cfg$scope,
    beta_per_TS_SD=NA_real_, SE=NA_real_, CI95_low=NA_real_, CI95_high=NA_real_,
    odds_ratio_per_TS_SD=NA_real_, OR95_low=NA_real_, OR95_high=NA_real_,
    Wald_z=NA_real_, P_value=NA_real_, q_value=NA_real_, numerical_ok=FALSE,
    n_samples=nrow(d), n_TS=nlevels(d$TS), q_family_size=length(cfg$mc)*length(cfg$env),
    inference_status="unavailable_numerical_issue", stringsAsFactors=FALSE)
  if (!is.null(fit)) {
    extracted <- joint_TS7_capture(function() {
      tab <- summary(fit)$coefficients$cond
      stopifnot(all(cfg$z %in% rownames(tab)))
      mr <- model_row; rr <- result
      mr$convergence_code <- fit$fit$convergence
      mr$pdHess <- isTRUE(fit$sdr$pdHess)
      gradient <- fit$sdr$gradient.fixed
      if (length(gradient) && all(is.finite(gradient))) mr$max_abs_gradient <- max(abs(gradient))
      mr$TS_SD_logit <- unname(attr(glmmTMB::VarCorr(fit)$cond$TS,"stddev")[1])
      mr$phi <- exp(unname(glmmTMB::fixef(fit)$disp["(Intercept)"]))
      mr$logLik <- as.numeric(stats::logLik(fit)); mr$AIC <- stats::AIC(fit)
      rr$beta_per_TS_SD <- tab[cfg$z,"Estimate"]; rr$SE <- tab[cfg$z,"Std. Error"]
      mr$numerical_ok <- isTRUE(mr$convergence_code==0L) && mr$pdHess &&
        all(is.finite(c(tab[,'Estimate'],tab[,'Std. Error'],mr$logLik,mr$phi,mr$TS_SD_logit))) &&
        all(tab[,'Std. Error']>0) && mr$phi>0 && mr$TS_SD_logit>=0
      flags <- c(if (is.finite(mr$max_abs_gradient) && mr$max_abs_gradient>.01) "gradient_gt_0.01",
        if (any(abs(tab[,'Estimate'])>10)) "large_coefficient",
        if (is.finite(mr$TS_SD_logit) && mr$TS_SD_logit<1e-4) "TS_variance_near_zero")
      mr$review_flags <- paste(flags,collapse=" | "); rr$numerical_ok <- mr$numerical_ok
      if (mr$numerical_ok) {
        rr$CI95_low <- rr$beta_per_TS_SD-stats::qnorm(.975)*rr$SE
        rr$CI95_high <- rr$beta_per_TS_SD+stats::qnorm(.975)*rr$SE
        rr$odds_ratio_per_TS_SD <- exp(rr$beta_per_TS_SD)
        rr$OR95_low <- exp(rr$CI95_low); rr$OR95_high <- exp(rr$CI95_high)
        rr$Wald_z <- tab[cfg$z,"z value"]; rr$P_value <- tab[cfg$z,"Pr(>|z|)"]
        if (any(!is.finite(rr$P_value)) || any(rr$P_value<0 | rr$P_value>1)) stop("无效Wald P值。")
        rr$inference_status <- "exploratory_joint_Wald_7TS"
      }
      list(model=mr,result=rr)
    })
    model_row$warnings <- paste(unique(c(attempt$warnings,extracted$warnings)),collapse=" | ")
    if (!nzchar(extracted$error)) {
      warning_text <- model_row$warnings
      model_row <- extracted$value$model; model_row$warnings <- warning_text
      result <- extracted$value$result
    } else model_row$error <- extracted$error
  }
  model_row$elapsed_seconds <- proc.time()[[3]]-started
  list(MC_ID=mc, fit=fit, model_summary=model_row, coefficients=result,
    formula=paste(deparse(model_formula),collapse=" "))
}

joint_TS7_new_run <- function(cfg) {
  root <- file.path(cfg$project,"results")
  for (i in seq_len(9999L)) {
    path <- file.path(root,sprintf("%s_%s_run%03d",cfg$run_prefix,cfg$scope,i))
    if (dir.exists(path)) next
    if (dir.create(path,showWarnings=FALSE)) return(path)
  }
  stop("无法创建独立运行目录。")
}
joint_TS7_write <- function(x,path) utils::write.csv(x,path,row.names=FALSE,fileEncoding="UTF-8",na="")

# [7] check阶段不拟合、不写结果；fit阶段保存新目录，绝不覆盖旧分析。
joint_TS7_main <- function(cfg=joint_TS7_config()) {
  input <- joint_TS7_load(cfg)
  n_models <- length(cfg$mc); n_coefficients <- n_models*length(cfg$env)
  cat("\n范围：",cfg$scope,"；样本：",nrow(input$samples),"；TS：7；联合模型：",n_models,
    "；环境系数：",n_coefficients,"（四营养盐＋温度，不加DIN）。\n",sep="")
  print(input$scaling,row.names=FALSE); print(input$collinearity,row.names=FALSE)
  if (cfg$stage=="check") {
    cat("\n输入检查完成；未拟合模型，未创建结果目录。正式运行请设置 options(mc19.joint.stage='fit')。\n")
    return(invisible(input))
  }
  joint_TS7_require("glmmTMB")
  run <- joint_TS7_new_run(cfg); support <- file.path(run,"supporting_files")
  dir.create(support); dir.create(file.path(support,"models")); dir.create(file.path(support,"source_inputs"))
  final_status <- "incomplete"
  on.exit({
    writeLines(final_status,file.path(run,"run_status.txt"))
    writeLines(capture.output(utils::sessionInfo()),file.path(support,"sessionInfo.txt"))
  },add=TRUE)
  for (path in input$fingerprints$path) if (!file.copy(path,file.path(support,"source_inputs",basename(path))))
    stop("无法保存输入快照。")
  joint_TS7_write(input$samples,file.path(support,"sample_metadata.csv"))
  for (name in c("environment_TS","scaling","collinearity","inclusion","fingerprints","temperature_source"))
    joint_TS7_write(input[[name]],file.path(support,paste0("input_",name,".csv")))
  joint_TS7_write(data.frame(term=rownames(input$correlation),input$correlation,check.names=FALSE),
    file.path(support,"environment_correlation.csv"))
  source_path <- file.path(cfg$project,"scripts","GLMM","GLMM_MC19_joint_TS7.R")
  if (!file.copy(source_path,file.path(support,"GLMM_MC19_joint_TS7.R"))) stop("无法保存脚本快照。")
  all_results <- model_summaries <- analysis_inputs <- vector("list",length(cfg$mc))
  for (i in seq_along(cfg$mc)) {
    mc <- cfg$mc[i]; d <- input$samples
    d$MC_ID <- mc; d$MC_reads <- unname(input$counts[mc,])
    d$Other_reads <- d$Ruegeria_reads_total-d$MC_reads
    d$MC_proportion <- d$MC_reads/d$Ruegeria_reads_total
    message("联合模型 ",i,"/",length(cfg$mc),"：",mc)
    entry <- joint_TS7_fit_one(mc,d,cfg)
    saveRDS(entry,file.path(support,"models",paste0(mc,".rds")))
    all_results[[i]] <- entry$coefficients; model_summaries[[i]] <- entry$model_summary; analysis_inputs[[i]] <- d
    # 逐MC保存进度；即使后续中断，也保留已完成模型和失败原因。
    joint_TS7_write(do.call(rbind,all_results[seq_len(i)]),file.path(support,"coefficients_progress.csv"))
  }
  results <- do.call(rbind,all_results); models <- do.call(rbind,model_summaries)
  # [8] 一个范围内19MC×5条件系数统一BH；不可估项保留NA，计划族仍为95。
  stopifnot(nrow(results)==n_coefficients,nrow(models)==n_models)
  results$q_value <- stats::p.adjust(results$P_value,method="BH",n=n_coefficients)
  joint_TS7_write(results,file.path(run,"00_coefficients.csv"))
  joint_TS7_write(models,file.path(run,"01_model_checks.csv"))
  joint_TS7_write(do.call(rbind,analysis_inputs),file.path(support,"analysis_input.csv"))
  saveRDS(list(config=cfg,input_fingerprints=input$fingerprints,scaling=input$scaling,
    n_models=n_models,n_coefficients=n_coefficients,formula=entry$formula,diagnostics="residual_diagnostics_not_run"),
    file.path(support,"manifest.rds"))
  if (!identical(unname(tools::md5sum(input$fingerprints$path)),input$fingerprints$md5)) stop("运行期间原始输入发生变化。")
  final_status <- if (all(models$numerical_ok)) "fits_complete_residual_diagnostics_pending" else
    "fits_complete_some_models_need_review"
  cat("\n已保存：",run,"\n数值通过：",sum(models$numerical_ok),"/19。残差诊断及跨TS稳定性尚未执行。\n",sep="")
  invisible(run)
}

if (!isTRUE(getOption("mc19.joint.define_only",FALSE))) joint_TS7_main()
