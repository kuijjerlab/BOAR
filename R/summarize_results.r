summarize_clustering_results <- function(result, modality) {
  if (inherits(result, "ConsensusPartitionList")) {
    run_list <- result@list
  } else if (is.list(result)) {
    run_list <- result
  } else {
    stop("`result` must be a ConsensusPartitionList or a named list.")
  }

  run_names <- names(run_list)
  if (is.null(run_names) || length(run_names) == 0) stop("No runs found.")
  if (missing(modality) || is.null(modality) || !nzchar(modality)) stop("`modality` must be a non-empty string.")

  all_stats <- lapply(run_names, function(run_name) {
    cp <- run_list[[run_name]]

    if (!inherits(cp, "ConsensusPartition")) {
      stop("Element `", run_name, "` is not a `ConsensusPartition` object.")
    }

    parts <- strsplit(run_name, ":", fixed = TRUE)[[1]]
    top_method <- if (length(parts) >= 1) parts[1] else NA_character_
    partition_method <- if (length(parts) >= 2) parts[2] else NA_character_

    stats <- as.data.frame(cola::get_stats(cp), stringsAsFactors = FALSE)
    if (nrow(stats) == 0) return(NULL)

    best_k_obj <- tryCatch(
      suppressMessages(cola::suggest_best_k(cp)),
      error = function(e) NA_integer_
    )

    best_k <- if (length(best_k_obj) > 0) best_k_obj[1] else NA_integer_

    if (is.na(best_k)) {
      best_row <- stats[which.max(stats$concordance), , drop = FALSE]
      best_k <- best_row$k[1]
    } else {
      best_row <- stats[stats$k == best_k, , drop = FALSE]

      if (nrow(best_row) == 0) {
        best_row <- stats[which.max(stats$concordance), , drop = FALSE]
        best_k <- best_row$k[1]
      }
    }

    optional_vals <- attr(best_k_obj, "optional")

    optional_k <- if (!is.null(optional_vals) && length(optional_vals) > 0) {
      paste(optional_vals, collapse = ",")
    } else {
      NA_character_
    }

    best_row$modality <- modality
    best_row$feature_selection <- top_method
    best_row$clustering_method <- partition_method

    mat <- cola::get_matrix(cp)
    best_row$Features <- nrow(mat)
    best_row$Samples <- ncol(mat)
    best_row$top_n <- cp@top_n
    best_row$k <- best_k
    best_row$Optional_k <- optional_k

    best_row
  })

  all_stats <- Filter(Negate(is.null), all_stats)
  if (length(all_stats) == 0) stop("No stats could be extracted from `result`.")

  out <- do.call(rbind, all_stats)
  rownames(out) <- NULL

  wanted_first <- c(
    "modality",
    "feature_selection",
    "clustering_method",
    "top_n",
    "k",
    "Optional_k",
    "Samples",
    "Features",
    "1-PAC",
    "mean_silhouette",
    "concordance"
  )

  existing_first <- intersect(wanted_first, colnames(out))
  remaining <- setdiff(colnames(out), existing_first)

  out[, c(existing_first, remaining), drop = FALSE]
}

select_best_clustering_combination <- function(
  summary_df,
  weights = c(
    one_minus_pac = 1,
    mean_silhouette = 1,
    concordance = 1
  )
) {
  df <- summary_df

  if ("1-PAC" %in% colnames(df) && !("one_minus_pac" %in% colnames(df))) {
    df$one_minus_pac <- df[["1-PAC"]]
  }

  required_cols <- c("one_minus_pac", "mean_silhouette", "concordance")
  missing_cols <- setdiff(required_cols, colnames(df))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  agg <- aggregate(
    cbind(one_minus_pac, mean_silhouette, concordance) ~
      feature_selection + clustering_method,
    data = df,
    FUN = mean,
    na.rm = TRUE
  )

  normalize <- function(x) {
    rng <- range(x, na.rm = TRUE)
    if (diff(rng) == 0) return(rep(0, length(x)))
    (x - rng[1]) / diff(rng)
  }

  agg$norm_pac <- normalize(agg$one_minus_pac)
  agg$norm_sil <- normalize(agg$mean_silhouette)
  agg$norm_conc <- normalize(agg$concordance)

  agg$score <-
    weights["one_minus_pac"] * agg$norm_pac +
    weights["mean_silhouette"] * agg$norm_sil +
    weights["concordance"] * agg$norm_conc

  best <- agg[which.max(agg$score), , drop = FALSE]

  list(
    best = best,
    ranked = agg[order(agg$score, decreasing = TRUE), ]
  )
}