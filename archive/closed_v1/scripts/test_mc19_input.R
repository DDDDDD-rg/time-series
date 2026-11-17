# Input-only regression checks. No model fitting, permutation tests or figures.
local({
  old_options <- options(mc19.analysis_scope=getOption("mc19.analysis_scope"))
  on.exit(options(old_options),add=TRUE)
  project <- "C:/Users/xzy/Desktop/CUHK_Study/GNBF6010_Research_Projects/Research/Msc main project-time series"
  e <- new.env(parent=globalenv())
  for (f in c("glmm_MC19_pilot_helpers.R", "glmm_MC19_wald_helpers.R"))
    sys.source(file.path(project, "scripts", "GLMM", f), e)
  fail <- function(expr, pattern) {
    error <- tryCatch({force(expr); NULL}, error=function(x) x)
    stopifnot(inherits(error, "error"), grepl(pattern, conditionMessage(error), fixed=TRUE))
  }
  cfg <- e$mc19_wald_config(project, "coral")
  stopifnot(cfg$expected_TS==8L, cfg$n_models==95L, cfg$q_family_size==95L,
    length(cfg$environments)==5L, "DIN_mgN_L" %in% cfg$environments)
  env <- as.data.frame(readxl::read_excel(cfg$workbook, sheet="Environment_TS"))
  raw8 <- as.character(env$TS)=="TS8"
  # The delivered template intentionally has missing TS8; check its clear error.
  if (all(is.na(env$Nitrate_mgN_L_raw[raw8]))) {
    fail(e$mc19_wald_load_input(cfg), "8TS输入尚未填齐")
  }
  fixture <- getOption("mc19.test.input_fixture", "")
  if (nzchar(fixture)) {
    # This private fixture has synthetic TS8 environmental values only.
    for (scope in c("coral", "ecosystem")) {
      cfg <- e$mc19_wald_config(project, scope)
      cfg$workbook <- fixture
      input <- e$mc19_wald_load_input(cfg)
      n <- if (scope=="coral") 116L else 356L
      stopifnot(length(input$blocks)==95L, nrow(input$scaling)==5L,
        all(input$scaling$n_TS==8L), all(input$scaling$n_samples==n),
        all(vapply(input$blocks,nrow,integer(1))==n),
        all(vapply(input$blocks,function(x) nlevels(x$TS),integer(1))==8L))
      din <- input$environment_TS$DIN_mgN_L
      stopifnot(max(abs(din-c(.016,.018,.033,.132,.041,.187,.108,.036)))<1e-12)
      stopifnot(identical(input$environment_TS$DIN_mgN_L_censored,
        c(TRUE,TRUE,TRUE,FALSE,TRUE,FALSE,FALSE,TRUE)))
      options(mc19.analysis_scope=scope)
      p <- new.env(parent=globalenv())
      sys.source(file.path(project,"scripts","PERMANOVA","permanova_MC19_data_settings.R"),p)
      perma <- p$read_mc19_samples(fixture)
      stopifnot(nrow(perma$community)==n, nrow(perma$batches)==8L,
        length(p$env_cols)==5L, identical(as.character(perma$metadata$New_ID),input$blocks[[1]]$New_ID),
        max(abs(perma$batches$DIN_mgN_L-din))<1e-12)
    }
    orders <- p$enumerate_ts_orders(8L)
    stopifnot(nrow(orders)==40320L, identical(as.integer(orders[1,]),1:8))
  }
  # Validate guards independently of the on-disk template's current fill state.
  cfg <- e$mc19_wald_config(project,"coral")
  good <- env
  good$Water_quality_date <- as.character(as.Date(good$Water_quality_date))
  good$Water_quality_date[raw8] <- "2026-02-15"
  for (field in e$mc19_primary_env_fields) good[[paste0(field,"_raw")]][raw8] <- "0.01"
  missing <- good; missing$Ammonia_mgN_L_raw[raw8] <- NA_character_
  fail(e$mc19_prepare_environment(missing), "TS8:Ammonia_mgN_L_raw")
  duplicate <- good; duplicate$TS[8] <- "TS7"
  fail(e$mc19_prepare_environment(duplicate), "每个TS一行")
  invalid <- good; invalid$Nitrite_mgN_L_raw[8] <- "unknown"
  fail(e$mc19_prepare_environment(invalid), "原始报告值无效")
  wrong_month <- good; wrong_month$Water_quality_date[8] <- "2026-03-15"
  fail(e$mc19_prepare_environment(wrong_month), "Sampling_month不一致")
  cat("PASS: 8TS input checks, DIN/L2, missing/date guards, scope matching and 95-task preparation; no model fitted.\n")
})
