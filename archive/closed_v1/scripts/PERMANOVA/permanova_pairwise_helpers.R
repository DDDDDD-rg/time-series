# Shared pairwise functions only; sourcing this file never runs a test.
prepare_pair_data <- function(data, effect, groups) {
  block <- if (effect != "TS") "TS" else "Host_type"
  stopifnot(effect %in% c("Host_type", "Sample_type", "TS"), length(groups) == 2,
            !anyDuplicated(groups), all(groups %in% as.character(data[[effect]])))
  selected <- as.data.frame(data[as.character(data[[effect]]) %in% groups, , drop = FALSE])
  counts <- table(droplevels(factor(selected[[block]])),
                  factor(as.character(selected[[effect]]), levels = groups))
  common <- rownames(counts)[rowSums(counts > 0) == 2]
  coverage <- data.frame(Block_level = rownames(counts),
                        N1 = as.integer(counts[, 1]), N2 = as.integer(counts[, 2]),
                        Used_in_pair = rownames(counts) %in% common)
  selected <- selected[as.character(selected[[block]]) %in% common, , drop = FALSE]
  selected <- droplevels(selected)
  selected <- selected[order(selected[[block]], selected[[effect]],
                             selected$Repeat_group_ID), , drop = FALSE]
  stopifnot(!anyNA(selected$Repeat_group_ID),
            anyDuplicated(selected$Repeat_group_ID) == 0)
  rownames(selected) <- as.character(selected$Repeat_group_ID)
  list(data = selected, block = block, common = common, coverage = coverage)
}

# ===== 3. 每对先建完整结果行：跳过或报错也保留这一行 =====
new_pair_result <- function(effect, groups, pair_id) {
  data.frame(
    Pair_ID = pair_id, Factor = effect, Group1 = groups[1], Group2 = groups[2],
    N_MC = length(mc_cols), MC_set = paste(mc_cols, collapse = ";"), Input_MD5 = input_md5,
    Analysis_version = analysis_version,
    N1 = NA_integer_, N2 = NA_integer_, N_units = NA_integer_,
    N_hosts = NA_integer_, N_TS = NA_integer_, TS_used = "", Hosts_used = "",
    Block_factor = if (effect != "TS") "TS" else "Host_type",
    N_common_blocks = NA_integer_, Common_blocks = "", Residual_df = NA_integer_,
    N_label_assignments = NA_real_, N_perm = NA_integer_,
    F_raw = NA_real_, R2_raw = NA_real_, P_raw = NA_real_,
    F_lingoes = NA_real_, R2_lingoes = NA_real_, P_lingoes = NA_real_,
    P_disp_raw = NA_real_, P_disp_lingoes = NA_real_,
    Status = "未检验", Warnings = "", Dispersion_note = "",
    Analysis_unit = "Repeat_group均值", Scope = analysis_scope,
    Permutation_design = if (effect != "TS")
      "共同TS内置换Repeat_group；模型调整TS" else
      "共同Host_type内置换Repeat_group；模型调整Host_type",
    stringsAsFactors = FALSE
  )
}

# ===== 4. 检验一对：沿用01的test_factor，保留raw和Lingoes =====
test_one_pair <- function(data, effect, groups, pair_id, seed) {
  result <- new_pair_result(effect, groups, pair_id)
  coverage <- data.frame(Block_level = character(), N1 = integer(), N2 = integer(),
                         Used_in_pair = logical())
  warning_log <- character()
  result <- tryCatch(
    withCallingHandlers({
      prepared <- prepare_pair_data(data, effect, groups)
      pair_data <- prepared$data
      coverage <- prepared$coverage
      group_sizes <- table(factor(as.character(pair_data[[effect]]), levels = groups))
      result$N1 <- as.integer(group_sizes[1])
      result$N2 <- as.integer(group_sizes[2])
      result$N_units <- nrow(pair_data)
      result$N_hosts <- nlevels(pair_data$Host_type)
      result$N_TS <- nlevels(pair_data$TS)
      result$TS_used <- paste(levels(pair_data$TS), collapse = ";")
      result$Hosts_used <- paste(levels(pair_data$Host_type), collapse = ";")
      result$N_common_blocks <- length(prepared$common)
      result$Common_blocks <- paste(prepared$common, collapse = ";")
      if (length(prepared$common) < 2) {
        result$Status <- paste0("跳过：共同", prepared$block, "少于2")
      } else if (any(group_sizes < 2)) {
        result$Status <- "跳过：至少一组的Repeat_group少于2"
      } else {
        rhs <- c(effect, prepared$block)
        design <- model.matrix(reformulate(rhs), data = pair_data)
        result$Residual_df <- nrow(pair_data) - qr(design)$rank
        if (qr(design)$rank < ncol(design) || result$Residual_df < 1) {
          result$Status <- "跳过：模型共线或无残差自由度"
        } else {
          # 标签分配数供判断置换空间；不把9999当作所有对的实际置换数。
          used_counts <- coverage[coverage$Used_in_pair, , drop = FALSE]
          result$N_label_assignments <- prod(choose(used_counts$N1 + used_counts$N2,
                                                    used_counts$N1))
          set.seed(seed)
          permutations <- permute::shuffleSet(
            nrow(pair_data),
            control = permute::how(nperm = n_permutations,
                                  blocks = pair_data[[prepared$block]])
          )
          stopifnot(nrow(permutations) > 0, ncol(permutations) == nrow(pair_data))
          tested <- test_factor(
            data = pair_data, rhs = rhs, effect = effect,
            permutations = permutations, dispersion_group = effect
          )
          for (column in names(tested)) result[[column]] <- tested[[column]]
          result$N_hosts <- nlevels(pair_data$Host_type)
        }
      }
      result
    }, warning = function(w) {
      warning_log <<- c(warning_log, conditionMessage(w))
      tryInvokeRestart("muffleWarning")
    }),
    error = function(e) {
      result$Status <- paste("error:", conditionMessage(e))
      result
    }
  )
  result$Warnings <- paste(unique(c(result$Warnings[nzchar(result$Warnings)],
                                    warning_log)), collapse = " | ")
  coverage$Pair_ID <- rep(pair_id, nrow(coverage))
  coverage$Factor <- rep(effect, nrow(coverage))
  coverage$Group1 <- rep(groups[1], nrow(coverage))
  coverage$Group2 <- rep(groups[2], nrow(coverage))
  coverage$Block_factor <- rep(result$Block_factor, nrow(coverage))
  list(result = result, coverage = coverage)
}

# ===== 5. 一个因素的全部比较：按计划的组合数做BH =====
# 全样本的14类有91对，coral有3对；8个TS有28对；三类样本有3对。
# 跳过的组合保留NA，BH的family_size仍使用计划组合数。
run_pairwise_family <- function(data, effect, seed_offset = 0L) {
  group_levels <- levels(droplevels(data[[effect]]))
  if (length(group_levels) < 2) stop(effect, "少于2组，无法两两比较。")
  pairs <- combn(group_levels, 2)
  family_size <- ncol(pairs)
  pair_results <- pair_coverage <- vector("list", family_size)
  for (i in seq_len(family_size)) {
    message(effect, " ", i, "/", family_size, ": ",
            paste(pairs[, i], collapse = " vs "))
    tested <- test_one_pair(data, effect, pairs[, i], i,
                            random_seed + seed_offset + i - 1L)
    pair_results[[i]] <- tested$result
    pair_coverage[[i]] <- tested$coverage
  }
  results <- bind_rows(pair_results)
  stopifnot(nrow(results) == family_size, !anyDuplicated(results$Pair_ID))
  valid_raw <- with(results, is.finite(F_raw) & is.finite(R2_raw) &
                      is.finite(P_raw) & P_raw >= 0 & P_raw <= 1)
  valid_lingoes <- with(results, is.finite(F_lingoes) & is.finite(R2_lingoes) &
                          is.finite(P_lingoes) & P_lingoes >= 0 & P_lingoes <= 1)
  results$BH_family_size <- family_size
  results$P_BH_raw <- add_bh(ifelse(valid_raw, results$P_raw, NA_real_), family_size)
  results$P_BH_lingoes <- add_bh(ifelse(valid_lingoes, results$P_lingoes, NA_real_),
                               family_size)
  results$Significance <- add_stars(results$P_BH_raw)

  # 离散度是另一组诊断问题，单独BH，不与组成差异的p混在一起。
  results$P_disp_BH_raw <- add_bh(results$P_disp_raw, family_size)
  results$P_disp_BH_lingoes <- add_bh(results$P_disp_lingoes, family_size)
  results$Dispersion_note_BH <- vapply(seq_len(nrow(results)), function(i) {
    p <- c(results$P_disp_BH_raw[i], results$P_disp_BH_lingoes[i])
    if (!all(is.finite(p))) return("离散度诊断未完整完成，查看Status/Warnings")
    if (min(p) < 0.05) return("BH校正后离散度差异：结合组内变异解释组成检验")
    "未检出BH校正后的离散度差异"
  }, character(1))
  list(results = results, coverage = bind_rows(pair_coverage),
       levels = group_levels)
}

# ===== 6. 整理对称矩阵：校正p和R2；对角线留NA =====
make_pair_matrix <- function(results, group_levels, column) {
  matrix_result <- matrix(NA_real_, length(group_levels), length(group_levels),
                           dimnames = list(group_levels, group_levels))
  for (i in seq_len(nrow(results))) {
    a <- results$Group1[i]
    b <- results$Group2[i]
    matrix_result[a, b] <- matrix_result[b, a] <- results[[column]][i]
  }
  matrix_result
}

