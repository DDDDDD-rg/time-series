# 从scripts.Rproj运行。只验证输入和受控失败分支；不会拟合真实模型。
(function() {
  entry <- "GLMM/GLMM_MC19_joint_TS7.R"
  if (!file.exists(entry)) stop("请先在scripts.Rproj所在目录运行本检查。")
  old_options <- options(mc19.joint.define_only=TRUE,mc19.joint.stage="check",
    mc19.analysis_scope=getOption("mc19.analysis_scope","coral"))
  on.exit(options(old_options),add=TRUE)
  engine <- new.env(parent=globalenv())
  source(entry,local=engine,encoding="UTF-8")
  engine$joint_TS7_require(c("readxl","glmmTMB"))
  # 把拟合函数暂时替换为只会报错的哨兵：本测试绝不调用真实拟合。
  real_fit <- get("glmmTMB",envir=asNamespace("glmmTMB"))
  fit_calls <- 0L
  fit_arguments <- NULL
  fake_fit <- function(...) {fit_calls <<- fit_calls+1L;fit_arguments <<- list(...);stop("FIT_TRIPWIRE")}
  assignInNamespace("glmmTMB",fake_fit,ns="glmmTMB")
  on.exit(assignInNamespace("glmmTMB",real_fit,ns="glmmTMB"),add=TRUE)
  report <- engine$joint_TS7_report(c("<0.002","0.008","< 0.005","0"))
  stopifnot(identical(report$value,c(.001,.008,.0025,0)),
    identical(report$censored,c(TRUE,FALSE,TRUE,FALSE)))
  for (x in list(NA_character_,"","<0","abc","-0.002"))
    stopifnot(inherits(tryCatch(engine$joint_TS7_report(x),error=identity),"error"))
  temp <- engine$joint_TS7_temperature(c("25.8","0","-1.5"))
  stopifnot(identical(temp$value,c(25.8,0,-1.5)),!any(temp$censored))
  for(x in list(NA_character_,"","<16.7","NaN","Inf","warm"))
    stopifnot(inherits(tryCatch(engine$joint_TS7_temperature(x),error=identity),"error"))
  for (scope in c("coral","ecosystem")) {
    options(mc19.analysis_scope=scope)
    cfg <- engine$joint_TS7_config()
    before <- list.files(file.path(cfg$project,"results"),pattern="^joint_glmm_TS7_")
    printed <- capture.output(input <- suppressWarnings(engine$joint_TS7_main(cfg)))
    stopifnot(fit_calls==0L,nrow(input$samples)==cfg$expected_n,
      identical(levels(input$samples$TS),paste0("TS",1:7)),
      nrow(input$environment_TS)==7L,nrow(input$counts)==19L,
      identical(rownames(input$counts),cfg$mc),all(input$collinearity$design_rank==6L),
      length(cfg$env)==5L,!"DIN_mgN_L" %in% cfg$env,
      identical(cfg$run_prefix,"joint_glmm_TS7_env5"),
      identical(input$temperature_source$SST_C,c(25.8,31.4,28.6,16.7,24.6,26.6,25.2)),
      all(input$temperature_source$Station=="PM11"),all(input$temperature_source$Depth=="Surface Water"))
    stopifnot(max(abs(colMeans(input$environment_TS[cfg$z])))<1e-12,
      max(abs(vapply(input$environment_TS[cfg$z],sd,numeric(1))-1))<1e-12)
    stopifnot(identical(before,list.files(file.path(cfg$project,"results"),pattern="^joint_glmm_TS7_")))
    folder <- if (scope=="coral") "glmm_coral_v1" else "glmm_all_hosts_v1"
    old <- utils::read.csv(file.path(cfg$project,"results",folder,"supporting_files","analysis_input.csv"),
      stringsAsFactors=FALSE,check.names=FALSE)
    old <- old[old$MC_ID=="MC15" & old$environment=="Nitrite_mgN_L",,drop=FALSE]
    stopifnot(nrow(old)==cfg$expected_n,!anyDuplicated(old$New_ID))
    old <- old[match(input$samples$New_ID,as.character(old$New_ID)),,drop=FALSE]
    stopifnot(identical(as.character(old$TS),as.character(input$samples$TS)),
      all(old$Ruegeria_reads_total==input$samples$Ruegeria_reads_total),
      all(old$MC_reads==unname(input$counts["MC15",])))
    for (field in setdiff(cfg$env,"SST_C")) stopifnot(all(abs(old[[field]]-input$samples[[field]])<1e-12))
    reference <- utils::read.csv(file.path(cfg$project,"results",paste0("MC19_GLMM_TS7_env6_",scope,"_run001"),
      "supporting_files","input_environment_TS.csv"),stringsAsFactors=FALSE)
    reference <- reference[match(input$environment_TS$TS,reference$TS),,drop=FALSE]
    stopifnot(all(abs(reference$SST_C-input$environment_TS$SST_C)<1e-12))
    cat(scope,": original counts, denominator, environment and sample membership match saved 7TS input.\n")
  }
  d <- input$samples;d$MC_reads <- unname(input$counts["MC15",]);d$Other_reads <- d$Ruegeria_reads_total-d$MC_reads
  failed <- engine$joint_TS7_fit_one("MC15",d,cfg)
  stopifnot(fit_calls==1L,is.null(failed$fit),failed$model_summary$error=="FIT_TRIPWIRE",
    nrow(failed$coefficients)==5L,all(is.na(failed$coefficients$P_value)),
    all(failed$coefficients$q_family_size==95L),!failed$model_summary$numerical_ok,
    all(cfg$z %in% all.vars(fit_arguments$formula)),!"DIN_mgN_L" %in% all.vars(fit_arguments$formula),
    identical(fit_arguments$family$family,"betabinomial"),identical(fit_arguments$family$link,"logit"))
  # 故意改变一条原CSV温度，确认来源不一致时check停止，真实文件保持原样。
  wrong <- utils::read.csv(cfg$temperature_source,check.names=FALSE,fileEncoding="UTF-8-BOM")
  temp_field <- grep("^Temperature",names(wrong),value=TRUE)
  hit <- which(wrong$Station=="PM11" & wrong$Depth=="Surface Water" &
    as.Date(wrong$Dates,format="%Y/%m/%d")==as.Date("2024-04-15"))
  stopifnot(length(hit)==1L)
  wrong[[temp_field]][hit] <- as.numeric(wrong[[temp_field]][hit])+1
  wrong_path <- tempfile(fileext=".csv");on.exit(unlink(wrong_path),add=TRUE)
  utils::write.csv(wrong,wrong_path,row.names=FALSE,fileEncoding="UTF-8")
  bad_cfg <- cfg;bad_cfg$temperature_source <- wrong_path
  rejected <- tryCatch(engine$joint_TS7_load(bad_cfg),error=identity)
  stopifnot(inherits(rejected,"error"),grepl("PM11",conditionMessage(rejected)),fit_calls==1L)
  cat("PASS: parsing, censoring, sample alignment, preserved denominator, TS standardization inputs, no-output check mode and failed-fit record. Real model fits: 0.\n")
})()
