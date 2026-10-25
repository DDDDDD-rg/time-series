# 将BH校正后的两两PERMANOVA结果转成字母显示；只做逻辑转换，不做检验。
# 使用insert-and-absorb：显著的两组不能共享字母，未检出差异的两组共享字母。
absorb_letter_columns <- function(columns) {
  columns <- columns[, colSums(columns) > 0, drop = FALSE]
  signatures <- apply(columns, 2, function(x) paste(which(x), collapse = ","))
  columns <- columns[, !duplicated(signatures), drop = FALSE]
  keep <- rep(TRUE, ncol(columns))
  if (ncol(columns) > 1) {
    for (j in seq_len(ncol(columns))) {
      others <- setdiff(seq_len(ncol(columns)), j)
      keep[j] <- !any(vapply(others, function(k) {
        all(!columns[, j] | columns[, k])
      }, logical(1)))
    }
  }
  columns[, keep, drop = FALSE]
}

letters_from_p_matrix <- function(p_matrix, alpha = 0.05) {
  stopifnot(is.matrix(p_matrix), nrow(p_matrix) == ncol(p_matrix),
            nrow(p_matrix) >= 2, is.finite(alpha), alpha > 0, alpha < 1,
            all(is.finite(p_matrix)), all(p_matrix >= 0 & p_matrix <= 1),
            isTRUE(all.equal(p_matrix, t(p_matrix), check.attributes = FALSE)))
  n <- nrow(p_matrix)
  significant <- p_matrix <= alpha
  diag(significant) <- FALSE
  columns <- matrix(TRUE, n, 1, dimnames = list(rownames(p_matrix), NULL))
  significant_pairs <- which(significant & upper.tri(significant), arr.ind = TRUE)
  for (r in seq_len(nrow(significant_pairs))) {
    a <- significant_pairs[r, 1]
    b <- significant_pairs[r, 2]
    split_columns <- list()
    for (j in seq_len(ncol(columns))) {
      current <- columns[, j]
      if (current[a] && current[b]) {
        first <- second <- current
        first[a] <- FALSE
        second[b] <- FALSE
        split_columns <- c(split_columns, list(first, second))
      } else {
        split_columns <- c(split_columns, list(current))
      }
    }
    columns <- absorb_letter_columns(do.call(cbind, split_columns))
  }
  column_order <- order(apply(columns, 2, function(x) which(x)[1]),
                         apply(columns, 2, function(x) paste(which(x), collapse = ",")))
  columns <- columns[, column_order, drop = FALSE]
  symbols <- c(letters, LETTERS, paste0(".", letters), paste0(".", LETTERS))
  if (ncol(columns) > length(symbols)) stop("字母列数过多，请改用两两p矩阵显示。")
  colnames(columns) <- symbols[seq_len(ncol(columns))]
  shared <- tcrossprod(columns * 1) > 0
  off_diagonal <- row(p_matrix) != col(p_matrix)
  stopifnot(all(shared[off_diagonal] == !significant[off_diagonal]))
  labels <- apply(columns, 1, function(x) paste(colnames(columns)[x], collapse = ""))
  list(letters = labels, matrix = columns)
}

# 有缺失检验时找最大的完整可比子集；并列按原图例顺序取第一个。
# 12宿主/8TS的组合规模很小。未入选的组标n/a；不把缺失p当作不显著。
largest_complete_subset <- function(available) {
  n <- nrow(available)
  if (all(available)) return(seq_len(n))
  if (n > 16) stop("缺失比较的自动子集选择目前最多支持16组。")
  if (n <= 2) return(integer())
  for (size in seq.int(n - 1L, 2L)) {
    subsets <- combn(n, size, simplify = FALSE)
    for (index in subsets) {
      if (all(available[index, index, drop = FALSE])) return(index)
    }
  }
  integer()
}

derive_pairwise_letters <- function(results, group_levels, alpha = 0.05) {
  required <- c("Factor", "Group1", "Group2", "P_BH_raw", "BH_family_size", "Status")
  stopifnot(all(required %in% names(results)), length(group_levels) >= 2,
            !anyDuplicated(group_levels), length(unique(results$Factor)) == 1)
  expected <- choose(length(group_levels), 2)
  if (nrow(results) != expected || any(results$BH_family_size != expected)) {
    stop("两两结果的计划比较数与图中组数不一致。")
  }
  p <- matrix(NA_real_, length(group_levels), length(group_levels),
              dimnames = list(group_levels, group_levels))
  diag(p) <- 1
  pair_keys <- character(nrow(results))
  for (i in seq_len(nrow(results))) {
    a <- match(results$Group1[i], group_levels)
    b <- match(results$Group2[i], group_levels)
    if (anyNA(c(a, b)) || a == b) stop("两两结果含未知组或组内比较。")
    pair_keys[i] <- paste(sort(c(a, b)), collapse = ":")
    value <- results$P_BH_raw[i]
    if (is.finite(value) && value >= 0 && value <= 1) p[a, b] <- p[b, a] <- value
  }
  if (anyDuplicated(pair_keys)) stop("两两结果有重复组合。")
  available <- is.finite(p)
  eligible <- largest_complete_subset(available)
  labels <- setNames(rep(NA_character_, length(group_levels)), group_levels)
  letter_matrix <- NULL
  if (length(eligible) >= 2) {
    display <- letters_from_p_matrix(p[eligible, eligible, drop = FALSE], alpha)
    labels[eligible] <- display$letters
    letter_matrix <- display$matrix
  }
  factor_name <- unique(results$Factor)
  table <- data.frame(
    Factor = factor_name, Group = group_levels, Letters = unname(labels),
    Letter_status = ifelse(is.na(labels), "n/a：完整比较子集之外", "已标注"),
    Alpha = alpha, P_column = "P_BH_raw", BH_family_size = expected,
    stringsAsFactors = FALSE
  )
  missing_pairs <- results[!is.finite(results$P_BH_raw),
                          c("Group1", "Group2", "Status"), drop = FALSE]
  list(table = table, matrix = letter_matrix, excluded = group_levels[is.na(labels)],
       missing_pairs = missing_pairs, family_size = expected, alpha = alpha)
}
