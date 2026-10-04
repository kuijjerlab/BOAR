plot_k_choices_from_summary <- function(
  cancer_types,
  summary_dir,
  feature_selection,
  clustering_method,
  output_file = NULL,
  cancer_order = NULL,
  width = 5.5,
  height = 4.5,
  dpi = 300
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package `ggplot2` must be installed.")
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Package `dplyr` must be installed.")

  all_results <- lapply(cancer_types, function(ct) {
    summary_file <- file.path(summary_dir, paste0(ct, "_combined.csv"))
    if (!file.exists(summary_file)) stop("Missing summary file: ", summary_file)

    df <- utils::read.csv(summary_file, stringsAsFactors = FALSE, check.names = FALSE)

    if (!all(c("feature_selection", "clustering_method", "modality") %in% names(df))) {
      stop("Summary file ", summary_file, " is missing required columns.")
    }

    k_col <- c("best_k", "Best_k", "suggested_k", "k")[c("best_k", "Best_k", "suggested_k", "k") %in% names(df)][1]
    if (is.na(k_col)) stop("No best-k column found in ", summary_file)

    df <- df[df$feature_selection == feature_selection & df$clustering_method == clustering_method, , drop = FALSE]
    if (nrow(df) == 0) return(NULL)

    df$cancer_type <- ct
    df$best_k <- as.numeric(df[[k_col]])
    df[, c("cancer_type", "modality", "best_k"), drop = FALSE]
  })

  dat <- dplyr::bind_rows(all_results)
  dat <- dat[!is.na(dat$best_k), , drop = FALSE]

  if (nrow(dat) == 0) stop("No best-k data available for plotting.")

  dat$modality <- dplyr::case_when(
    dat$modality %in% c("edge_weights", "Edge Weights") ~ "Edge Weights",
    dat$modality %in% c("indegrees", "Indegrees") ~ "Indegrees",
    dat$modality %in% c("expression", "Expression") ~ "Expression",
    TRUE ~ dat$modality
  )

  modality_levels <- c("Edge Weights", "Indegrees", "Expression")
  dat$modality <- factor(dat$modality, levels = modality_levels[modality_levels %in% unique(dat$modality)])

  if (is.null(cancer_order)) cancer_order <- sort(unique(dat$cancer_type))
  dat$cancer_type <- factor(dat$cancer_type, levels = rev(cancer_order))

  p <- ggplot2::ggplot(dat, ggplot2::aes(x = best_k, y = cancer_type, color = modality)) +
    ggplot2::geom_point(size = 1.8, alpha = 0.9) +
    ggplot2::facet_wrap(~ modality, nrow = 1, scales = "fixed") +
    ggplot2::scale_color_manual(values = c(
      "Edge Weights" = "#F8766D",
      "Indegrees" = "#619CFF",
      "Expression" = "#4DAF4A"
    ), drop = FALSE) +
    ggplot2::scale_x_continuous(breaks = 2:6, limits = c(2, 6)) +
    ggplot2::labs(x = NULL, y = NULL, color = NULL) +
    ggplot2::theme_bw(base_size = 8) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      legend.position = "none",
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(colour = "gray75", linewidth = 0.3, linetype = "dashed"),
      panel.grid.major.x = ggplot2::element_line(colour = "gray90", linewidth = 0.2),
      strip.background = ggplot2::element_rect(fill = "white"),
      strip.text = ggplot2::element_text(face = "bold")
    )

  if (!is.null(output_file)) ggplot2::ggsave(output_file, p, width = width, height = height, dpi = dpi)
  p
}

plot_k_choices <- function(
  k_summary_df,
  cancer_order = NULL,
  data_type_order = c("Edge Weights", "Indegrees", "Expression"),
  feature_selection = NULL,
  clustering_method = NULL,
  width = 3.5,
  height = 4,
  output_file = NULL
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package `ggplot2` must be installed.")

  required_cols <- c(
    "CANCER_TYPE",
    "DATA_TYPE",
    "FEATURE_SELECTION",
    "CLUSTERING_METHOD",
    "Best_k"
  )

  missing_cols <- setdiff(required_cols, colnames(k_summary_df))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  df <- k_summary_df

  if (!is.null(feature_selection)) {
    df <- df[df$FEATURE_SELECTION %in% feature_selection, , drop = FALSE]
  }

  if (!is.null(clustering_method)) {
    df <- df[df$CLUSTERING_METHOD %in% clustering_method, , drop = FALSE]
  }

  if (nrow(df) == 0) stop("No rows left after filtering.")

  df$Best_k <- suppressWarnings(as.numeric(df$Best_k))

  df$DATA_TYPE <- dplyr::case_when(
    df$DATA_TYPE == "edge_weights" ~ "Edge Weights",
    df$DATA_TYPE == "indegrees" ~ "Indegrees",
    df$DATA_TYPE == "expression" ~ "Expression",
    TRUE ~ as.character(df$DATA_TYPE)
  )

  if (is.null(cancer_order)) cancer_order <- unique(as.character(df$CANCER_TYPE))

  df$DATA_TYPE <- factor(df$DATA_TYPE, levels = data_type_order)
  df$CANCER_TYPE <- factor(df$CANCER_TYPE, levels = rev(cancer_order))

  df <- df[
    !is.na(df$CANCER_TYPE) &
      !is.na(df$DATA_TYPE) &
      !is.na(df$Best_k),
    ,
    drop = FALSE
  ]

  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(x = Best_k, y = CANCER_TYPE, color = DATA_TYPE)
  ) +
    ggplot2::geom_point(
      size = 1.8,
      alpha = 0.9,
      position = ggplot2::position_dodge(width = 0.5)
    ) +
    ggplot2::scale_color_manual(
      values = c(
        "Edge Weights" = "#F8766D",
        "Indegrees" = "#619CFF",
        "Expression" = "#4DAF4A"
      )
    ) +
    ggplot2::theme_bw(base_size = 8) +
    ggplot2::scale_x_continuous(breaks = 2:6, limits = c(2, 6)) +
    ggplot2::labs(x = NULL, y = NULL, color = NULL) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      legend.position = "top",
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(
        colour = "gray75",
        linewidth = 0.3,
        linetype = "dashed"
      ),
      panel.grid.major.x = ggplot2::element_line(
        colour = "gray90",
        linewidth = 0.2
      )
    )

  if (!is.null(output_file)) {
    ggplot2::ggsave(
      filename = output_file,
      plot = p,
      width = width,
      height = height,
      dpi = 300
    )
  }

  p
}