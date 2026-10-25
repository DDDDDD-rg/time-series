# 干湿季 / 环境：保留不完整宿主，置换批次标签。
# 在RStudio中打开本文件，运行全部代码即可；自动读取01的数据和设置。
# 本版的检验假设：参与检验的TS，其季节/环境标签在零假设下可交换。
# 香港观测时间序列并非随机分配实验，因此结果解释为探索性关联。

project_dir <- "C:/Users/xzy/Desktop/CUHK_Study/GNBF6010_Research_Projects/Research/Msc main project-time series"
source(file.path(project_dir, "scripts", "01_data_settings.R"), encoding = "UTF-8", local = TRUE)
# 单独保存新结果，避免覆盖旧的完整面板分析。测试时可临时指定此选项。
output_dir <- getOption("mc.batch_output_dir",
                       file.path(output_dir, "repeat_groups"))

# ===== 1. 保留每一个Repeat_group：这里不再取平均 =====
# 01已经在Repeat_group内完成第一轮均值；本函数只筛选和排序。
# 同一Host_type、同一TS可以有多行，每行保留自己的Repeat_group_ID和MC组成。
# all_hosts与coral使用相同规则；缺采的宿主×TS不补0。
make_repeat_data <- function(data, target_ts) {
  selected <- data %>%
    filter(as.character(TS) %in% target_ts) %>%
    arrange(TS, Host_type, Repeat_group_ID) %>%
    droplevels()
  stopifnot(nrow(selected) > 0, !anyNA(selected$Repeat_group_ID),
            anyDuplicated(selected$Repeat_group_ID) == 0)
  selected
}

# ===== 2. 列出批次标签的所有不同分配 =====
# 目前最多8个TS：8!=40320，可完整枚举。超过8批时需另定置换方案。
# 同样的值反复出现时去除重复分配；例如4湿4干只有70种不同分配。
# 返回的每行是一套批次值，不是用来重排MC样本行的索引。
make_label_assignments <- function(values) {
  if (length(values) < 2 || length(values) > 8) stop("本版枚举适用于2至8个批次。")
  if (any(!is.finite(values)) || length(unique(values)) < 2) {
    stop("批次值必须完整且有变化。")
  }
  index <- as.matrix(permute::allPerms(length(values),
    control = permute::how(observed = TRUE, maxperm = 50000)))
  assignments <- t(apply(index, 1, function(i) values[i]))
  unique(assignments)
}

# ===== 3. 只计算宿主调整后的F和R2，不在这里置换样本行 =====
# 用Host_type的类别编码扣除宿主背景。
# 统计量与adonis2(distance ~ Host_type + Batch_value, by="margin")的Batch_value行相同。
# 这里只用dbrda计算统计量；真正的p值来自下一节的批次标签置换。
fit_batch_statistics <- function(data, distance, batch_value, correction = FALSE) {
  data$Batch_value <- batch_value
  design <- model.matrix(~ Host_type + Batch_value, data)
  residual_df <- nrow(data) - ncol(design)
  if (qr(design)$rank < ncol(design) || residual_df < 1) {
    stop("批次指标与宿主混杂，或模型没有残差自由度。")
  }
  # 标准宿主类别编码，调整效果与Condition(Host_type)相同。
  # 避免当前vegan在负约束轴置换中计算分类质心时报错。
  host_background <- model.matrix(~ Host_type, data)[, -1, drop = FALSE]
  fit <- vegan::dbrda(distance ~ Batch_value + Condition(host_background),
                      data = data, add = correction)
  effect_ss <- fit$CCA$tot.chi
  residual_ss <- fit$CA$tot.chi
  effect_df <- fit$CCA$qrank
  F_value <- (effect_ss / effect_df) / (residual_ss / residual_df)
  R2_value <- effect_ss / fit$tot.chi
  if (any(!is.finite(c(F_value, R2_value)))) stop("F或R2不是有限值。")
  c(F = F_value, R2 = R2_value)
}

# ===== 4. 批次标签置换：MC行、宿主标签、缺采格局始终固定 =====
test_batch_labels <- function(data, factor_name, season = FALSE) {
  ts_used <- ts_levels[ts_levels %in% as.character(data$TS)]
  batch_values <- batch_data[[factor_name]][match(ts_used, batch_data$TS)]
  if (season) batch_values <- as.numeric(batch_values == "Wet") # 干=0，湿=1
  row_batch <- match(as.character(data$TS), ts_used)
  distance <- make_distance(data)
  observed_values <- batch_values[row_batch]
  observed_raw <- fit_batch_statistics(data, distance, observed_values)
  observed_lingoes <- fit_batch_statistics(data, distance, observed_values, "lingoes")

  assignments <- make_label_assignments(batch_values)
  is_observed <- rowSums(sweep(assignments, 2, batch_values, "!=")) == 0
  stopifnot(sum(is_observed) == 1)
  permuted_values <- assignments[!is_observed, , drop = FALSE]
  F_null_raw <- F_null_lingoes <- numeric(nrow(permuted_values))
  for (i in seq_len(nrow(permuted_values))) {
    # 先给每个TS一个新的值，再让该TS的所有宿主共享这个值。
    shuffled_values <- permuted_values[i, row_batch]
    F_null_raw[i] <- fit_batch_statistics(data, distance, shuffled_values)["F"]
    F_null_lingoes[i] <- fit_batch_statistics(data, distance, shuffled_values, "lingoes")["F"]
  }
  # 加1计入原始分配。本版完整枚举，所以这是全部允许标签分配中的比例。
  permutation_p <- function(null, observed) {
    tolerance <- 1e-8 * max(1, abs(observed))
    (1 + sum(null >= observed - tolerance)) / (1 + length(null))
  }
  result <- data.frame(
    Factor = factor_name, N_units = nrow(data), N_hosts = nlevels(data$Host_type),
    N_MC = length(mc_cols), MC_set = paste(mc_cols, collapse = ";"), Input_MD5 = input_md5,
    Analysis_version = analysis_version,
    N_TS = length(ts_used), TS_used = paste(ts_used, collapse = ";"),
    Hosts_used = paste(levels(data$Host_type), collapse = ";"),
    N_perm = nrow(permuted_values), N_label_assignments = nrow(assignments),
    F_raw = unname(observed_raw["F"]), R2_raw = unname(observed_raw["R2"]),
    P_raw = permutation_p(F_null_raw, observed_raw["F"]),
    F_lingoes = unname(observed_lingoes["F"]), R2_lingoes = unname(observed_lingoes["R2"]),
    P_lingoes = permutation_p(F_null_lingoes, observed_lingoes["F"]),
    P_disp_raw = NA_real_, P_disp_lingoes = NA_real_, Status = "ok",
    Dispersion_note = "连续变量：不人为分组做PERMDISP", Warnings = "",
    Analysis_unit = "Repeat_group均值",
    Permutation_design = "MC与宿主固定；完整枚举TS层面的指标标签分配",
    Scope = analysis_scope
  )

  # 季节是分类变量：另做批次标签置换的分散度诊断，不套用旧的行置换p值。
  # 这是按季节汇总的辅助诊断，未单独扣除宿主对分散度的影响。
  if (season) {
    dispersion_F <- function(value, correction = FALSE) {
      bd <- vegan::betadisper(distance, factor(value), type = "median",
                             bias.adjust = TRUE, add = correction)
      unname(anova(bd)["Groups", "F value"])
    }
    d_raw <- dispersion_F(observed_values)
    d_ling <- dispersion_F(observed_values, "lingoes")
    null_raw <- null_ling <- numeric(nrow(permuted_values))
    for (i in seq_len(nrow(permuted_values))) {
      values <- permuted_values[i, row_batch]
      null_raw[i] <- dispersion_F(values)
      null_ling[i] <- dispersion_F(values, "lingoes")
    }
    stopifnot(all(is.finite(c(d_raw, d_ling, null_raw, null_ling))))
    result$P_disp_raw <- permutation_p(null_raw, d_raw)
    result$P_disp_lingoes <- permutation_p(null_ling, d_ling)
    result$Dispersion_note <- if (min(result$P_disp_raw, result$P_disp_lingoes) <= 0.05)
      "季节汇总分散度有差异：谨慎解释；此诊断未调整宿主" else
      "季节汇总未检出分散度差异；此诊断未调整宿主，不证明相等"
  }
  result
}

# ===== 5. 干湿季：有季节标签的全部8个批次 =====
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
season_ts <- levels(analysis_data$TS)
season_data <- make_repeat_data(analysis_data, season_ts)
message("季节分析：", nlevels(season_data$Host_type), "类宿主，", nrow(season_data), "个Repeat_group分析单位。")
season_result <- test_batch_labels(season_data, "Season", season = TRUE)
season_result$Significance <- add_stars(season_result$P_raw)
print(season_result[, c("Factor", "N_units", "N_hosts", "N_TS", "R2_raw", "P_raw",
                       "P_lingoes", "Significance", "N_label_assignments")])
write_csv_excel(season_result, file.path(output_dir, "season_results.csv"), row.names = FALSE)
write_csv_excel(season_data, file.path(output_dir, "season_analysis_input.csv"), row.names = FALSE)

# ===== 6. 环境：每项指标按其有效批次检验；当前8项均为TS1–TS7 =====
# 当前TS1–TS7有这8项完整数据；TS8缺测，不补值。
# DIP有检测限记录，只保留在覆盖表，不把<0.002写成0或精确0.002。
environment_factors <- setdiff(env_cols, "DIP_mgP_L")
# Each variable uses its own available batches; no intersection-based data loss.
environment_inputs <- setNames(vector("list",length(environment_factors)),environment_factors)
environment_results <- vector("list",length(environment_factors))
for (i in seq_along(environment_factors)) {
  variable <- environment_factors[i]
  available_ts <- intersect(season_ts,batch_data$TS[is.finite(batch_data[[variable]])])
  if (length(available_ts)<6L) stop(variable,": fewer than 6 observed batches; review this design first.")
  data_i <- make_repeat_data(analysis_data,available_ts)
  environment_inputs[[variable]] <- data_i
  message("Environment ",i,"/",length(environment_factors),": ",variable,
    "; n=",nrow(data_i),"; batches=",length(available_ts))
  environment_results[[i]] <- test_batch_labels(data_i,variable)
}
# Legacy DO preview interface; other figures read environment_inputs[[variable]].
environment_data <- environment_inputs[["DO_mg_L"]]
environment_ts <- unique(as.character(environment_data$TS))
environment_results <- bind_rows(environment_results)
environment_results$P_BH_raw <- add_bh(environment_results$P_raw, length(environment_factors))
environment_results$P_BH_lingoes <- add_bh(environment_results$P_lingoes, length(environment_factors))
environment_results$Significance <- add_stars(environment_results$P_BH_raw)
print(environment_results[, c("Factor", "N_units", "N_hosts", "N_TS", "R2_raw", "P_raw",
                             "P_BH_raw", "P_BH_lingoes", "Significance")])
write_csv_excel(environment_results, file.path(output_dir, "environment_results.csv"), row.names = FALSE)
write_csv_excel(environment_data, file.path(output_dir, "environment_analysis_input.csv"), row.names = FALSE)
for (variable in environment_factors) write_csv_excel(environment_inputs[[variable]],
  file.path(output_dir,paste0("environment_input_",variable,".csv")),row.names=FALSE)

# ===== 7. 保存覆盖情况和作图输入 =====
environment_coverage <- batch_data
for (variable in environment_factors) environment_coverage[[paste0("Used_",variable)]] <-
  batch_data$TS %in% as.character(environment_inputs[[variable]]$TS)
environment_coverage$DIP_note <- ifelse(is.na(batch_data$DIP_mgP_L),
  ifelse(is.na(batch_data$DIP_raw), "缺测", "低于检测限，无精确浓度"), "有精确值")
write_csv_excel(environment_coverage, file.path(output_dir, "environment_coverage.csv"), row.names = FALSE)
correlation_data <- batch_data %>% filter(TS %in% season_ts)
write_csv_excel(cor(as.matrix(correlation_data[environment_factors]), use="pairwise.complete.obs"),
                file.path(output_dir, "environment_correlations.csv"))
# 覆盖表的数字是Repeat_group数量；0表示该宿主该TS未采到。
host_coverage <- as.data.frame.matrix(table(season_data$Host_type, season_data$TS))
write_csv_excel(host_coverage, file.path(output_dir, "host_batch_coverage.csv"))
saveRDS(list(season = season_result, environment = environment_results,
             season_input = season_data, environment_input = environment_data,
             environment_inputs = environment_inputs, analysis_version = analysis_version),
        file.path(output_dir, "season_environment_analysis.rds"))
writeLines(capture.output(sessionInfo()), file.path(output_dir, "season_environment_session.txt"))
message("输出目录：", output_dir)
message("p值依赖TS标签可交换；宿主调整不消除时间自相关或其他随时间变化的混杂。")
message("每项环境检验分别调整宿主，不是8项环境指标彼此调整后的独立作用。")
message("每行保留一个Repeat_group；同一Host_type×TS内的不同组没有再平均。")
