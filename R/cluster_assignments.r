save_cluster_assignments <- function(
  cp,
  file,
  k = NULL,
  use_best_k = TRUE,
  append_k_to_filename = FALSE
) {
  if (!inherits(cp, "ConsensusPartition")) stop("`cp` must be a `ConsensusPartition` object.")
  if (missing(file) || is.null(file) || !nzchar(file)) stop("`file` must be a non-empty string.")

  if (is.null(k) && isTRUE(use_best_k)) {
    k <- tryCatch(
      suppressWarnings(suppressMessages(cola::suggest_best_k(cp))),
      error = function(e) NA_integer_
    )
  }

  if (is.null(k) || length(k) == 0 || is.na(k)) {
    possible_k <- tryCatch(cp@k, error = function(e) NULL)

    if (!is.null(possible_k) && length(possible_k) > 0) {
      k <- max(possible_k, na.rm = TRUE)
      message("Using fallback k = ", k)
    } else {
      message("Skipping cluster save (no valid k)")
      return(NULL)
    }
  }

  classes_df <- cola::get_classes(cp, k = k)
  if (is.null(classes_df) || nrow(classes_df) == 0) stop("No classes returned for k = ", k)
  if (!("class" %in% colnames(classes_df))) stop("`get_classes()` did not return a `class` column.")

  if (append_k_to_filename) {
    ext <- tools::file_ext(file)

    if (nzchar(ext)) {
      base <- sub(paste0("\\.", ext, "$"), "", file)
      file <- paste0(base, "_k", k, ".", ext)
    } else {
      file <- paste0(file, "_k", k)
    }
  }

  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  con <- file(file, open = "wt")
  on.exit(close(con), add = TRUE)

  split_samples <- split(rownames(classes_df), classes_df$class)

  for (cl in names(split_samples)) {
    cat(cl, "\n", file = con)
    cat(split_samples[[cl]], sep = "\n", file = con)
    cat("\n", file = con)
  }

  invisible(list(file = file, k = k, classes = classes_df))
}

read_cluster_txt_simple <- function(path, cancer_type, modality) {
  if (!file.exists(path)) return(NULL)

  lines <- trimws(readLines(path, warn = FALSE))
  lines <- lines[nchar(lines) > 0]

  cur_cluster <- NA_character_
  out <- list()
  idx <- 0L

  for (ln in lines) {
    if (grepl("^\\d+$", ln)) {
      cur_cluster <- ln
      next
    }

    if (!is.na(cur_cluster)) {
      idx <- idx + 1L
      out[[idx]] <- data.frame(
        cancer_type = cancer_type,
        modality = modality,
        cluster = cur_cluster,
        sample = ln,
        stringsAsFactors = FALSE
      )
    }
  }

  if (idx == 0L) return(NULL)
  dplyr::bind_rows(out)
}

plot_cluster_balance_from_txt <- function(
  cancer_types,
  modalities,
  feature_selection,
  clustering_method,
  cluster_root = "clusters",
  modality_dirs = NULL,
  output_prefix = NULL,
  cancer_order = NULL
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package `ggplot2` must be installed.")
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Package `dplyr` must be installed.")
  if (!requireNamespace("tidyr", quietly = TRUE)) stop("Package `tidyr` must be installed.")
  if (!requireNamespace("purrr", quietly = TRUE)) stop("Package `purrr` must be installed.")
  if (missing(cancer_types) || length(cancer_types) == 0) stop("`cancer_types` must be provided.")
  if (missing(modalities) || length(modalities) < 2) stop("Provide at least 2 modalities.")

  samples_df <- purrr::map_dfr(cancer_types, function(ct) {
    purrr::map_dfr(modalities, function(mod) {
      f <- find_cluster_file_for_analysis(
        cancer_type = ct,
        modality = mod,
        feature_selection = feature_selection,
        clustering_method = clustering_method,
        cluster_root = cluster_root,
        modality_dirs = modality_dirs
      )

      if (is.na(f)) {
        message("Missing cluster file for ", ct, " / ", mod)
        return(NULL)
      }

      read_cluster_txt_simple(f, ct, mod)
    })
  })

  if (nrow(samples_df) == 0) stop("No cluster assignments could be read.")

  counts_df <- samples_df |>
    dplyr::count(cancer_type, modality, cluster, name = "n")

  all_clusters <- sort(unique(counts_df$cluster), method = "radix")
  if (is.null(cancer_order)) cancer_order <- unique(counts_df$cancer_type)

  counts_plot <- counts_df |>
    dplyr::mutate(
      cancer_type = factor(cancer_type, levels = cancer_order),
      modality = factor(modality, levels = modalities),
      cluster = factor(cluster, levels = all_clusters)
    ) |>
    tidyr::complete(cancer_type, modality, cluster, fill = list(n = 0))

  p_sizes <- ggplot2::ggplot(
    counts_plot,
    ggplot2::aes(x = modality, y = n, fill = cluster)
  ) +
    ggplot2::geom_col(width = 0.8) +
    ggplot2::facet_wrap(~ cancer_type, ncol = min(11, length(unique(counts_plot$cancer_type)))) +
    ggplot2::labs(
      x = NULL,
      y = NULL,
      fill = "Cluster",
      title = NULL
    ) +
    ggplot2::theme_minimal(base_size = 16) +
    ggplot2::theme(
      panel.grid.major.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(size = 16, angle = 45, hjust = 1),
      axis.text.y = ggplot2::element_text(size = 16)
    )

  cluster_balance <- counts_df |>
    dplyr::group_by(cancer_type, modality) |>
    dplyr::summarise(
      total = sum(n),
      k = sum(n > 0),
      entropy = {
        p <- n[n > 0] / sum(n)
        -sum(p * log(p))
      },
      balance = ifelse(k > 1, entropy / log(k), 0),
      .groups = "drop"
    )

  print(cluster_balance)

  if (!is.null(output_prefix)) {
    ggplot2::ggsave(
      paste0(output_prefix, "_cluster_size.png"),
      p_sizes,
      width = 22,
      height = 7,
      dpi = 300
    )
  }

  invisible(list(
    sample_assignments = samples_df,
    counts = counts_df,
    balance = cluster_balance,
    size_plot = p_sizes
  ))
}

compare_clusterings_ari <- function(
  cancer_types,
  feature_selection,
  clustering_method,
  modality_1,
  modality_2,
  cluster_root = "clusters",
  modality_dirs = NULL
) {
  if (!requireNamespace("mclust", quietly = TRUE)) stop("Package `mclust` must be installed.")

  read_cluster_txt <- function(file) {
    lines <- trimws(readLines(file, warn = FALSE))
    lines <- lines[lines != ""]

    current_cluster <- NA_character_
    sample_ids <- character(0)
    cluster_ids <- character(0)

    for (ln in lines) {
      if (grepl("^[0-9]+$", ln)) {
        current_cluster <- ln
        next
      }

      if (is.na(current_cluster)) stop("Sample found before any cluster label in file: ", file)

      sample_ids <- c(sample_ids, ln)
      cluster_ids <- c(cluster_ids, current_cluster)
    }

    if (length(sample_ids) == 0) return(NULL)
    stats::setNames(cluster_ids, sample_ids)
  }

  results <- lapply(cancer_types, function(ct) {
    file_1 <- find_cluster_file_for_analysis(
      cancer_type = ct,
      modality = modality_1,
      feature_selection = feature_selection,
      clustering_method = clustering_method,
      cluster_root = cluster_root,
      modality_dirs = modality_dirs
    )

    file_2 <- find_cluster_file_for_analysis(
      cancer_type = ct,
      modality = modality_2,
      feature_selection = feature_selection,
      clustering_method = clustering_method,
      cluster_root = cluster_root,
      modality_dirs = modality_dirs
    )

    if (is.na(file_1) || is.na(file_2)) {
      return(data.frame(cancer = ct, ARI = NA_real_, stringsAsFactors = FALSE))
    }

    cl_1 <- read_cluster_txt(file_1)
    cl_2 <- read_cluster_txt(file_2)

    if (is.null(cl_1) || is.null(cl_2)) {
      return(data.frame(cancer = ct, ARI = NA_real_, stringsAsFactors = FALSE))
    }

    common_samples <- intersect(names(cl_1), names(cl_2))

    if (length(common_samples) == 0) {
      return(data.frame(cancer = ct, ARI = NA_real_, stringsAsFactors = FALSE))
    }

    ari <- mclust::adjustedRandIndex(
      as.vector(cl_1[common_samples]),
      as.vector(cl_2[common_samples])
    )

    data.frame(cancer = ct, ARI = ari, stringsAsFactors = FALSE)
  })

  out <- do.call(rbind, results)
  rownames(out) <- NULL
  out
}

compare_all_clusterings_ari <- function(
  cancer_types,
  feature_selection,
  clustering_method,
  modality_dirs = NULL,
  cluster_root = "clusters",
  output_plot = NULL,
  output_table = NULL,
  width = 4.6,
  height = 4,
  dpi = 300
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package `ggplot2` must be installed.")
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Package `dplyr` must be installed.")

  available_modalities <- if (!is.null(modality_dirs)) {
    names(modality_dirs)
  } else {
    c("Edge Weights", "Indegrees", "Expression")
  }

  comparisons <- list(
    c("Edge Weights", "Indegrees"),
    c("Edge Weights", "Expression"),
    c("Indegrees", "Expression")
  )

  comparison_names <- c(
    "Edge weights vs indegrees",
    "Edge weights vs expression",
    "Indegrees vs expression"
  )

  comparison_results <- list()

  for (i in seq_along(comparisons)) {
    mods <- comparisons[[i]]
    if (!all(mods %in% available_modalities)) next

    ari_df <- compare_clusterings_ari(
      cancer_types = cancer_types,
      feature_selection = feature_selection,
      clustering_method = clustering_method,
      modality_1 = mods[1],
      modality_2 = mods[2],
      cluster_root = cluster_root,
      modality_dirs = modality_dirs
    )

    ari_df$comparison <- comparison_names[i]
    comparison_results[[comparison_names[i]]] <- ari_df
  }

  if (length(comparison_results) == 0) stop("No modality pairs are available for ARI comparison.")

  combined_df <- dplyr::bind_rows(comparison_results)

  cancer_order <- combined_df |>
    dplyr::filter(!is.na(ARI)) |>
    dplyr::group_by(cancer) |>
    dplyr::summarise(mean_ARI = mean(ARI, na.rm = TRUE), .groups = "drop") |>
    dplyr::arrange(mean_ARI) |>
    dplyr::pull(cancer)

  if (length(cancer_order) == 0) stop("No valid ARI values available.")

  comparison_levels <- comparison_names[
    comparison_names %in% unique(combined_df$comparison)
  ]

  plot_df <- combined_df |>
    dplyr::filter(!is.na(ARI)) |>
    dplyr::mutate(
      cancer = factor(cancer, levels = cancer_order),
      comparison = factor(comparison, levels = comparison_levels)
    )

  mean_df <- plot_df |>
    dplyr::group_by(comparison) |>
    dplyr::summarise(mean_ARI = mean(ARI, na.rm = TRUE), .groups = "drop")

  cat("\nMean ARI by comparison:\n")
  print(mean_df |> dplyr::mutate(mean_ARI = round(mean_ARI, 3)))

  comparison_colors <- c(
    "Edge weights vs indegrees" = "#8522af",
    "Edge weights vs expression" = "#abab32",
    "Indegrees vs expression" = "#2696af"
  )

  comparison_sizes <- c(
    "Edge weights vs indegrees" = 2.2,
    "Edge weights vs expression" = 2.2,
    "Indegrees vs expression" = 2.2
  )

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = ARI, y = cancer)) +
    ggplot2::geom_point(
      ggplot2::aes(color = comparison, size = comparison),
      shape = 18,
      alpha = 0.85,
      position = ggplot2::position_dodge(width = 0.55)
    ) +
    ggplot2::geom_vline(
      data = mean_df,
      ggplot2::aes(xintercept = mean_ARI, color = comparison),
      linewidth = 0.55,
      linetype = "dashed",
      alpha = 0.85,
      show.legend = FALSE
    ) +
    ggplot2::scale_color_manual(name = NULL, values = comparison_colors) +
    ggplot2::scale_size_manual(name = NULL, values = comparison_sizes) +
    ggplot2::guides(
      color = ggplot2::guide_legend(
        title = NULL,
        override.aes = list(
          shape = 18,
          linetype = 0,
          linewidth = 0,
          size = 2.5,
          alpha = 0.9
        )
      ),
      size = "none"
    ) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_bw(base_size = 8) +
    ggplot2::theme(
      text = ggplot2::element_text(family = "sans"),
      axis.text.x = ggplot2::element_text(size = 7, angle = 90, vjust = 0.5, hjust = 1),
      axis.text.y = ggplot2::element_text(size = 7),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.title = ggplot2::element_text(size = 7),
      legend.text = ggplot2::element_text(size = 7),
      legend.key.width = grid::unit(0.7, "lines"),
      legend.key.height = grid::unit(0.7, "lines"),
      panel.grid.major = ggplot2::element_line(linewidth = 0.4, linetype = "dashed"),
      panel.grid.minor = ggplot2::element_blank()
    )

  if (!is.null(output_plot)) {
    dir.create(dirname(output_plot), recursive = TRUE, showWarnings = FALSE)
    ggplot2::ggsave(output_plot, p, width = width, height = height, dpi = dpi)
  }

  if (!is.null(output_table)) {
    dir.create(dirname(output_table), recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(combined_df, output_table, row.names = FALSE)
  }

  invisible(list(
    table = combined_df,
    plot = p,
    means = mean_df
  ))
}