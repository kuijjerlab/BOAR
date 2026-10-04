plot_clustering_summary <- function(
  summary_df,
  fs_order = c("ATC", "CV", "MAD", "SD"),
  modality_order = c("Edge Weights", "Indegrees", "Expression"),
  point_size = 4.5
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package `ggplot2` must be installed.")
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Package `dplyr` must be installed.")
  if (!requireNamespace("tidyr", quietly = TRUE)) stop("Package `tidyr` must be installed.")

  df <- as.data.frame(summary_df)
  if ("one_minus_pac" %in% colnames(df) && !("1-PAC" %in% colnames(df))) df[["1-PAC"]] <- df$one_minus_pac

  required_cols <- c("modality", "feature_selection", "clustering_method", "1-PAC", "mean_silhouette", "concordance")
  missing_cols <- setdiff(required_cols, colnames(df))
  if (length(missing_cols) > 0) stop("Missing required columns: ", paste(missing_cols, collapse = ", "))

  averaged <- df |>
    dplyr::mutate(
      feature_selection = trimws(as.character(feature_selection)),
      clustering_method = as.character(clustering_method)
    ) |>
    dplyr::filter(!is.na(feature_selection), feature_selection != "", feature_selection != "0") |>
    dplyr::group_by(modality, feature_selection, clustering_method) |>
    dplyr::summarise(
      one_minus_pac = mean(as.numeric(`1-PAC`), na.rm = TRUE),
      mean_silhouette = mean(as.numeric(mean_silhouette), na.rm = TRUE),
      concordance = mean(as.numeric(concordance), na.rm = TRUE),
      .groups = "drop"
    )

  long_all <- averaged |>
    tidyr::pivot_longer(
      cols = c(one_minus_pac, mean_silhouette, concordance),
      names_to = "metric",
      values_to = "score"
    ) |>
    dplyr::mutate(
      metric = dplyr::recode(
        metric,
        one_minus_pac = "1-PAC",
        mean_silhouette = "Mean Silhouette",
        concordance = "Concordance"
      ),
      metric = factor(metric, levels = c("1-PAC", "Mean Silhouette", "Concordance")),
      modality = factor(modality, levels = modality_order)
    )

  top5_data <- long_all |>
    dplyr::group_by(modality, feature_selection, clustering_method) |>
    dplyr::summarise(avg_score = mean(score, na.rm = TRUE), .groups = "drop") |>
    dplyr::group_by(modality) |>
    dplyr::arrange(dplyr::desc(avg_score), .by_group = TRUE) |>
    dplyr::slice_head(n = 5) |>
    dplyr::mutate(
      rank = dplyr::row_number(),
      label = paste0(rank, ". ", feature_selection, " + ", clustering_method, "  (", round(avg_score, 3), ")"),
      x_point = 2.30,
      x_text = 2.39,
      y = dplyr::case_when(
        rank == 1 ~ 0.305,
        rank == 2 ~ 0.25,
        rank == 3 ~ 0.195,
        rank == 4 ~ 0.14,
        rank == 5 ~ 0.085
      )
    ) |>
    dplyr::ungroup()

  top5_box <- top5_data |>
    dplyr::distinct(modality) |>
    dplyr::mutate(xmin = 2.20, xmax = 3.15, ymin = 0.045, ymax = 0.335)

  long_all <- long_all |>
    dplyr::mutate(
      feature_selection = factor(feature_selection, levels = fs_order),
      clustering_method = factor(clustering_method)
    )

  top5_data <- top5_data |>
    dplyr::mutate(
      feature_selection = factor(feature_selection, levels = fs_order),
      clustering_method = factor(clustering_method, levels = levels(long_all$clustering_method))
    )

  fs_colors <- c("ATC" = "#1f77b4", "CV" = "#ff7f0e", "MAD" = "#2ca02c", "SD" = "#d62728")
  shape_values <- c(15, 16, 8, 18, 17)
  cluster_methods <- levels(long_all$clustering_method)

  if (length(cluster_methods) > length(shape_values)) {
    stop("More clustering methods than available plotting shapes.")
  }

  shape_map <- stats::setNames(shape_values[seq_along(cluster_methods)], cluster_methods)
  fs_offsets <- c("ATC" = -0.18, "CV" = -0.06, "MAD" = 0.06, "SD" = 0.18)

  long_all <- long_all |>
    dplyr::mutate(
      metric_num = as.numeric(metric),
      fs_offset = fs_offsets[as.character(feature_selection)],
      xplot = metric_num + fs_offset
    )

  p <- ggplot2::ggplot(
    long_all,
    ggplot2::aes(x = xplot, y = score, colour = feature_selection, shape = clustering_method)
  ) +
    ggplot2::geom_vline(
      xintercept = c(1.5, 2.5),
      linetype = "dotted",
      linewidth = 0.8,
      colour = "gray50"
    ) +
    ggplot2::geom_point(size = point_size, stroke = 1.5, alpha = 0.95) +
    ggplot2::geom_rect(
      data = top5_box,
      ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      inherit.aes = FALSE,
      fill = "white",
      colour = "black",
      linewidth = 0.25,
      alpha = 0.9
    ) +
    ggplot2::geom_point(
      data = top5_data,
      ggplot2::aes(x = x_point, y = y, colour = feature_selection, shape = clustering_method),
      inherit.aes = FALSE,
      size = 3.1,
      stroke = 1.1,
      alpha = 0.95
    ) +
    ggplot2::geom_text(
      data = top5_data,
      ggplot2::aes(x = x_text, y = y, label = label),
      inherit.aes = FALSE,
      hjust = 0,
      size = 4.3
    ) +
    ggplot2::facet_wrap(~ modality, nrow = 1) +
    ggplot2::scale_x_continuous(
      breaks = c(1, 2, 3),
      labels = c("1-PAC", "Mean Silhouette", "Concordance"),
      limits = c(0.55, 3.45),
      expand = ggplot2::expansion(mult = c(0.02, 0.02))
    ) +
    ggplot2::scale_colour_manual(values = fs_colors, name = "Feature selection", drop = FALSE) +
    ggplot2::scale_shape_manual(values = shape_map, name = "Clustering method") +
    ggplot2::scale_y_continuous(limits = c(0, 1.08), breaks = seq(0, 1, by = 0.2)) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_bw(base_size = 18) +
    ggplot2::theme(
      strip.text = ggplot2::element_text(size = 20),
      strip.background = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(size = 17, angle = 30, hjust = 1),
      axis.text.y = ggplot2::element_text(size = 16),
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      legend.title = ggplot2::element_text(size = 18),
      legend.text = ggplot2::element_text(size = 16),
      legend.position = "bottom",
      legend.justification = "center",
      legend.box.just = "center",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      panel.border = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(
        colour = "lightgray",
        linewidth = 0.25
      ),
      panel.spacing = grid::unit(2.2, "lines"),
      plot.margin = ggplot2::margin(10, 20, 10, 10)
    ) +
    ggplot2::guides(
      colour = ggplot2::guide_legend(
        override.aes = list(shape = 1)
      ),
      shape = ggplot2::guide_legend()
    )

  attr(p, "averaged_data") <- averaged
  attr(p, "long_data") <- long_all
  attr(p, "top5") <- top5_data
  p
}

plot_best_combination_metrics <- function(
  summary_df,
  cancer_types,
  feature_selection,
  clustering_method,
  modalities = NULL,
  cluster_root = "clusters",
  modality_dirs = NULL,
  output_file = NULL,
  width = 4.2,
  height = 2.2,
  dpi = 300
) {
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Package `dplyr` must be installed.")
  if (!requireNamespace("tidyr", quietly = TRUE)) stop("Package `tidyr` must be installed.")
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package `ggplot2` must be installed.")

  best_df <- summary_df |>
    dplyr::filter(
      feature_selection == !!feature_selection,
      clustering_method == !!clustering_method
    )

  if (nrow(best_df) == 0) {
    stop("No summary rows found for ", feature_selection, " + ", clustering_method)
  }

  if (is.null(modalities)) modalities <- unique(as.character(best_df$modality))
  best_df <- best_df |> dplyr::filter(modality %in% modalities)

  metric_df <- best_df |>
    dplyr::select(cancer_type, modality, `1-PAC`, mean_silhouette, concordance) |>
    tidyr::pivot_longer(
      cols = c(`1-PAC`, mean_silhouette, concordance),
      names_to = "Metric",
      values_to = "Value"
    ) |>
    dplyr::mutate(
      Metric = dplyr::recode(
        Metric,
        mean_silhouette = "Mean Silhouette",
        concordance = "Concordance"
      )
    )

  balance_res <- plot_cluster_balance_from_txt(
    cancer_types = cancer_types,
    modalities = modalities,
    feature_selection = feature_selection,
    clustering_method = clustering_method,
    cluster_root = cluster_root,
    modality_dirs = modality_dirs,
    output_prefix = NULL
  )

  balance_long <- balance_res$balance |>
    dplyr::transmute(
      cancer_type = cancer_type,
      modality = modality,
      Metric = "Normalized Entropy",
      Value = balance
    )

  combined <- dplyr::bind_rows(metric_df, balance_long)

  modality_order <- c("Edge Weights", "Indegrees", "Expression")
  modality_order <- modality_order[modality_order %in% modalities]

  combined <- combined |>
    dplyr::mutate(
      modality = factor(modality, levels = modality_order),
      Metric = factor(
        Metric,
        levels = c("1-PAC", "Mean Silhouette", "Concordance", "Normalized Entropy")
      )
    )

  p <- ggplot2::ggplot(
    combined,
    ggplot2::aes(x = Metric, y = Value, fill = modality)
  ) +
    ggplot2::geom_violin(
      position = ggplot2::position_dodge(width = 0.8),
      alpha = 0.7,
      width = 0.8,
      trim = TRUE
    ) +
    ggplot2::geom_point(
      ggplot2::aes(color = modality),
      position = ggplot2::position_jitterdodge(
        jitter.width = 0.12,
        jitter.height = 0,
        dodge.width = 0.8
      ),
      size = 0.5,
      alpha = 0.85,
      show.legend = FALSE
    ) +
    ggplot2::scale_fill_manual(
      values = c(
        "Indegrees" = "#619CFF",
        "Edge Weights" = "#F8766D",
        "Expression" = "#4DAF4A"
      )
    ) +
    ggplot2::scale_color_manual(
      values = c(
        "Indegrees" = "#619CFF",
        "Edge Weights" = "#F8766D",
        "Expression" = "#4DAF4A"
      )
    ) +
    ggplot2::scale_y_continuous(limits = c(0, NA)) +
    ggplot2::theme_bw(base_size = 8) +
    ggplot2::labs(x = NULL, y = NULL, fill = NULL) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 20, hjust = 1),
      legend.position = "bottom",
      legend.text = ggplot2::element_text(size = 5),
      legend.key.size = grid::unit(0.3, "cm")
    )

  if (!is.null(output_file)) {
    ggplot2::ggsave(
      filename = output_file,
      plot = p,
      width = width,
      height = height,
      dpi = dpi
    )
  }

  invisible(list(
    feature_selection = feature_selection,
    clustering_method = clustering_method,
    data = combined,
    plot = p,
    balance = balance_res$balance
  ))
}