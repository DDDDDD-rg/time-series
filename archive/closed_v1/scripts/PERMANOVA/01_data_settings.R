# 历史MC8设置：保留供旧脚本参考；当前MC19入口不读取本文件。
# 当前设置与读取器见permanova_MC19_data_settings.R；入口为run_permanova.R。
# 数据和设置：02、03会自动读取本文件，不必先手动运行。
# 只分析输入表里的8个候选MC；不代表全部Ruegeria群落。

# ===== 1. 通常只需要修改这三个设置 =====
project_dir <- "C:/Users/xzy/Desktop/CUHK_Study/GNBF6010_Research_Projects/Research/Msc main project-time series"
analysis_scope <- getOption("mc.analysis_scope", "all_hosts")  # all_hosts / coral / sample_type
analysis_version <- "MC8_env_inclusive_v1"
n_permutations <- 9999

input_file <- file.path(project_dir, "data", "PERMANOVA_input_MC8_20261011.xlsx")
latest_result_pointer <- file.path(project_dir, "results", "00_latest_PERMANOVA_MC8.txt")
default_results_root <- if (file.exists(latest_result_pointer))
  readLines(latest_result_pointer, n = 1L, warn = FALSE, encoding = "UTF-8") else
  file.path(project_dir, "results", "PERMANOVA_MC8")
# Standalone tests must not overwrite the previous biological-host-only run.
if (is.null(getOption("mc.results_root"))) {
  manifest <- file.path(default_results_root, "00_analysis_version.txt")
  compatible <- file.exists(manifest) && identical(
    readLines(manifest, n=1L, warn=FALSE), analysis_version)
  if (!compatible) default_results_root <- file.path(project_dir,"results","PERMANOVA_MC8_env_inclusive_pending")
}
results_root <- getOption("mc.results_root", default_results_root)
output_dir <- file.path(results_root, analysis_scope)
random_seed <- 20261004

# 不自动安装包；缺包时明确列出。
packages <- c("readxl", "dplyr", "vegan", "permute")
missing_packages <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("请先安装这些R包：", paste(missing_packages, collapse = ", "))
}
suppressPackageStartupMessages(library(dplyr))
if (utils::packageVersion("dplyr") < "1.1.1") stop("请使用dplyr 1.1.1或更新版本。")
stopifnot(analysis_scope %in% c("all_hosts", "coral", "sample_type"), file.exists(input_file))

# ===== 2. 读取并核对输入 =====
community <- readxl::read_excel(input_file, sheet = "Community")
metadata <- readxl::read_excel(input_file, sheet = "Metadata")
mc_cols <- setdiff(names(community), "New_ID")
input_md5 <- unname(tools::md5sum(input_file))
env_cols <- c("Salinity_psu", "SST_C", "pH", "DO_mg_L", "Chla_ug_L",
              "DIN_mgN_L", "DIP_mgP_L", "PSM_mg_L", "Turbidity_NTU")
ts_levels <- paste0("TS", 1:8)

stopifnot(
  !anyNA(community$New_ID), !anyNA(metadata$New_ID),
  anyDuplicated(community$New_ID) == 0, anyDuplicated(metadata$New_ID) == 0,
  setequal(community$New_ID, metadata$New_ID),
  all(vapply(community[mc_cols], is.numeric, logical(1))),
  all(vapply(metadata[env_cols], is.numeric, logical(1))),
  !anyNA(community[mc_cols]),
  all(is.finite(as.matrix(community[mc_cols]))),
  all(as.matrix(community[mc_cols]) >= 0)
)
required_meta <- c("Repeat_group_ID", "Host_type", "Host_niche", "TS", "Season", "Sampling_month")
stopifnot(!anyNA(metadata[required_meta]), all(metadata$TS %in% ts_levels))

# 一个TS只能对应一套月份、季节和环境值；不能静默取多条记录的第一条。
batch_data <- metadata %>%
  select(TS, Sampling_month, Season, all_of(env_cols), DIP_raw) %>%
  distinct() %>%
  arrange(match(TS, ts_levels))
stopifnot(anyDuplicated(batch_data$TS) == 0)
sample_month <- as.integer(substr(batch_data$Sampling_month, 6, 7))
expected_season <- ifelse(sample_month %in% 4:9, "Wet", "Dry")
stopifnot(!anyNA(sample_month), all(batch_data$Season == expected_season))

# 同一重复组的宿主、生态类别和TS必须一致。
group_keys <- metadata %>% select(Repeat_group_ID, Host_type, Host_niche, TS) %>% distinct()
stopifnot(anyDuplicated(group_keys$Repeat_group_ID) == 0)

# ===== 3. 合并重复组：先对原始相对丰度取均值 =====
# 保留用户设定的前缀分组；不把组号解释成已确认的跨批次个体ID。
group_data <- metadata %>%
  left_join(community, by = "New_ID", relationship = "one-to-one") %>%
  group_by(Repeat_group_ID, Host_type, Host_niche, TS) %>%
  summarise(across(all_of(mc_cols), mean), n_records = n(), .groups = "drop") %>%
  left_join(batch_data, by = "TS", relationship = "many-to-one")
stopifnot(all(group_data$Host_type[group_data$Host_niche == "env"] %in% c("Seawater","Sediment")))
group_data$Sample_type <- factor(ifelse(group_data$Host_niche != "env", "HOST",
  as.character(group_data$Host_type)), levels=c("HOST","Seawater","Sediment"))
# Preserve all observed groups for the broad three-category comparison.
all_group_data <- group_data

if (analysis_scope == "coral") group_data <- filter(group_data, Host_niche == "coral")
host_counts <- group_data %>% count(Host_type, name = "N_repeat_groups")
excluded_hosts <- filter(host_counts, N_repeat_groups < 2)
# Keep the previous >=2-group rule for individual host/type tests; the pooled
# HOST/Seawater/Sediment analysis retains even singleton biological host types.
analysis_data <- if (analysis_scope == "sample_type") group_data else
  group_data %>% filter(!Host_type %in% excluded_hosts$Host_type)
analysis_data$Host_type <- factor(analysis_data$Host_type)
analysis_data$TS <- droplevels(factor(analysis_data$TS, levels = ts_levels))
analysis_data$Season <- factor(analysis_data$Season, levels = c("Dry", "Wet"))
stopifnot(nrow(analysis_data) > 0, !anyNA(analysis_data[mc_cols]))
message("分析范围：", analysis_scope, "；分析单位：", nrow(analysis_data),
        "；宿主数：", nlevels(analysis_data$Host_type))

# ===== 4. 小工具：保证两个入口使用同一种距离和结果格式 =====
# CSV加上UTF-8 BOM，让Excel直接打开时正确识别中文。
# 只补编码标记，write.csv生成的列、数值和缺失值保持不变。
write_csv_excel <- function(x, file, row.names = TRUE) {
  utils::write.csv(x, file, row.names = row.names, fileEncoding = "UTF-8")
  bytes <- readBin(file, what = "raw", n = file.info(file)$size)
  writeBin(c(as.raw(c(239, 187, 191)), bytes), file)
}

make_distance <- function(data) {
  mc_matrix <- as.matrix(data[, mc_cols, drop = FALSE])
  if (anyNA(mc_matrix) || any(!is.finite(mc_matrix)) || any(mc_matrix < 0)) {
    stop("MC矩阵含缺失、非有限值或负值。")
  }
  if (any(rowSums(mc_matrix) == 0)) stop("所选MC全部为0的分析单位需先单独处理。")
  vegan::vegdist(sqrt(mc_matrix), method = "bray")
}

add_bh <- function(p, family_size = length(p)) {
  adjusted <- rep(NA_real_, length(p))
  usable <- is.finite(p) & p >= 0 & p <= 1
  if (any(usable)) adjusted[usable] <- p.adjust(p[usable], "BH", n = family_size)
  adjusted
}

add_stars <- function(p) {
  case_when(is.na(p) ~ NA_character_, p <= 0.001 ~ "***",
            p <= 0.01 ~ "**", p <= 0.05 ~ "*", TRUE ~ "ns")
}

# 每次只提取正在检验的那一项；其他项的p值不能套用同一置换解释。
# 原始Bray为主要结果，Lingoes为敏感性检查；诊断失败不覆盖已得到的PERMANOVA。
test_factor <- function(data, rhs, effect, permutations, dispersion_group = NULL) {
  result <- data.frame(
    Factor = effect, N_units = nrow(data), N_TS = nlevels(droplevels(factor(data$TS))),
    N_hosts = length(unique(data$Host_type)), N_MC = length(mc_cols),
    MC_set = paste(mc_cols, collapse = ";"), Input_MD5 = input_md5,
    Analysis_version = analysis_version,
    TS_used = paste(unique(as.character(data$TS)), collapse = ";"),
    Hosts_used = paste(sort(unique(as.character(data$Host_type))), collapse = ";"),
    N_perm = nrow(permutations),
    F_raw = NA_real_, R2_raw = NA_real_, P_raw = NA_real_,
    F_lingoes = NA_real_, R2_lingoes = NA_real_, P_lingoes = NA_real_,
    P_disp_raw = NA_real_, P_disp_lingoes = NA_real_,
    Status = "ok", Dispersion_note = "连续变量：不做分类PERMDISP", Warnings = ""
  )
  warning_log <- character(0)
  withCallingHandlers({
    tryCatch({
      distance <- make_distance(data)
      if ("Env_value" %in% rhs && length(unique(data$Env_value)) < 2) {
        stop("环境指标没有变化，无法检验。")
      }
      model_formula <- reformulate(rhs, response = "distance")
      design <- model.matrix(reformulate(rhs), data = data)
      if (qr(design)$rank < ncol(design) || nrow(design) <= ncol(design)) {
        stop("模型共线或无残差自由度；请减少解释变量。")
      }
      stopifnot(ncol(permutations) == nrow(data))
      raw <- vegan::adonis2(model_formula, data = data, permutations = permutations, by = "margin")
      ling <- vegan::adonis2(model_formula, data = data, permutations = permutations, by = "margin", add = "lingoes")
      result[c("F_raw", "R2_raw", "P_raw")] <- as.list(as.numeric(raw[effect, c("F", "R2", "Pr(>F)")]))
      result[c("F_lingoes", "R2_lingoes", "P_lingoes")] <- as.list(as.numeric(ling[effect, c("F", "R2", "Pr(>F)")]))
      if (any(!is.finite(unlist(result[c("F_raw", "R2_raw", "P_raw", "F_lingoes", "R2_lingoes", "P_lingoes")])))) {
        result$Status <- "需检查非有限统计量"
      }
      if (!is.null(dispersion_group)) {
        result$Dispersion_note <- "分类变量诊断未完成"
        tryCatch({
          group <- droplevels(factor(data[[dispersion_group]]))
          bd_raw <- vegan::betadisper(distance, group, type = "median", bias.adjust = TRUE)
          bd_ling <- vegan::betadisper(distance, group, type = "median", bias.adjust = TRUE, add = "lingoes")
          pd_raw <- vegan::permutest(bd_raw, permutations = permutations)
          pd_ling <- vegan::permutest(bd_ling, permutations = permutations)
          result$P_disp_raw <- pd_raw$tab["Groups", "Pr(>F)"]
          result$P_disp_lingoes <- pd_ling$tab["Groups", "Pr(>F)"]
          disp_values <- c(pd_raw$tab["Groups", "F"], pd_ling$tab["Groups", "F"],
                           result$P_disp_raw, result$P_disp_lingoes)
          result$Dispersion_note <- if (any(!is.finite(disp_values))) "分散度统计量需检查" else
            if (min(result$P_disp_raw, result$P_disp_lingoes) < 0.05) "未校正分散度p<0.05：谨慎解释" else
              "未检出分散度差异，不等于证明相等"
        }, error = function(e) {
          result$Dispersion_note <<- paste("诊断失败：", conditionMessage(e))
        })
      }
    }, error = function(e) { result$Status <<- paste("error:", conditionMessage(e)) })
  }, warning = function(w) {
    warning_log <<- c(warning_log, conditionMessage(w))
    tryInvokeRestart("muffleWarning")
  })
  result$Warnings <- paste(unique(warning_log), collapse = " | ")
  result
}
