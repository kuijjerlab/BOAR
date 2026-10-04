clinical_associations_one_cancer <- function(
  cancer_type,
  modality,
  feature_selection,
  clustering_method,
  clinical_data,
  sample_id_fun = identity,
  exclude_cols = character(0),
  alpha = 0.05,
  verbose = TRUE,
  cluster_root = "clusters",
  modality_dirs = NULL
) {
  vmsg <- function(...) if (isTRUE(verbose)) message(...)
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Package `dplyr` must be installed.")

  vmsg("\n=== Processing: ", cancer_type, " / ", modality, " ===")

  classes_file <- find_cluster_file_for_analysis(
    cancer_type = cancer_type,
    modality = modality,
    feature_selection = feature_selection,
    clustering_method = clustering_method,
    cluster_root = cluster_root,
    modality_dirs = modality_dirs
  )

  if (is.na(classes_file)) {
    vmsg("    (skip: no matching cluster file found)")

    return(list(
      results_df = data.frame(),
      sig_df = data.frame(),
      sig_vars = character(0),
      present_vars = character(0),
      merged_df = data.frame()
    ))
  }

  classes_df <- tryCatch(
    read_classes_txt(classes_file),
    error = function(e) NULL
  )

  if (is.null(classes_df) || nrow(classes_df) == 0) {
    vmsg("    (skip: could not parse classes)")

    return(list(
      results_df = data.frame(),
      sig_df = data.frame(),
      sig_vars = character(0),
      present_vars = character(0),
      merged_df = data.frame()
    ))
  }

  classes_df$clinical_id <- as.character(sample_id_fun(classes_df$sample))
  classes_df$class <- as.character(classes_df$class)

  ids <- unique(classes_df$clinical_id)

  if (length(ids) < 5) {
    vmsg("    (skip: <5 samples)")

    return(list(
      results_df = data.frame(),
      sig_df = data.frame(),
      sig_vars = character(0),
      present_vars = character(0),
      merged_df = data.frame()
    ))
  }

  k <- nlevels(factor(classes_df$class))

  if (k < 2) {
    vmsg("    (skip: <2 classes)")

    return(list(
      results_df = data.frame(),
      sig_df = data.frame(),
      sig_vars = character(0),
      present_vars = character(0),
      merged_df = data.frame()
    ))
  }

  clinical <- as.data.frame(clinical_data, stringsAsFactors = FALSE)

  if (nrow(clinical) == 0) {
    stop("Clinical data are empty for ", cancer_type, ".")
  }

  required_cols <- c("sample_id", "survival_time", "survival_event")
  missing_cols <- setdiff(required_cols, colnames(clinical))

  if (length(missing_cols) > 0) {
    stop(
      "Clinical data for ", cancer_type, " are missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  clinical$clinical_id <- as.character(clinical$sample_id)

  merged_df <- merge(
    classes_df, clinical,
    by = "clinical_id", all.x = TRUE, sort = FALSE
  )

  standard_cols <- c("sample_id", "survival_time", "survival_event", "clinical_id")
  clinical_vars <- setdiff(colnames(clinical), c(exclude_cols, standard_cols))
  drop_cols <- unique(c(exclude_cols, standard_cols, "sample", "class"))

  categorical_cols <- names(merged_df)[sapply(merged_df, is.character)]
  categorical_cols <- setdiff(categorical_cols, drop_cols)

  numerical_cols <- names(merged_df)[sapply(merged_df, is.numeric)]
  numerical_cols <- setdiff(numerical_cols, drop_cols)

  results <- list()

  for (col in categorical_cols) {
    x <- merged_df[[col]]
    y <- merged_df$class

    ok <- !is.na(x) & !is.na(y)
    tab <- table(x[ok], y[ok])

    if (sum(tab) == 0 || nrow(tab) < 2 || ncol(tab) < 2) next
    if (min(rowSums(tab)) < 3 || min(colSums(tab)) < 3) next

    test <- tryCatch(chisq.test(tab), error = function(e) NULL)

    if (!is.null(test)) {
      results[[col]] <- data.frame(
        variable = col,
        type = "categorical",
        test = "Chi-squared",
        statistic = unname(ifelse(!is.null(test$statistic), test$statistic, NA)),
        p_value = test$p.value,
        stringsAsFactors = FALSE
      )
    }
  }

  for (col in numerical_cols) {
    x <- merged_df[[col]]
    y <- merged_df$class

    ok <- !is.na(x) & !is.na(y)

    if (length(unique(y[ok])) < 2 || length(unique(x[ok])) < 2) next
    if (min(table(y[ok])) < 3) next

    df_ok <- merged_df[ok, , drop = FALSE]
    df_ok$class <- as.factor(df_ok$class)

    fmla <- as.formula(paste(col, "~ class"))

    lm_fit <- tryCatch(lm(fmla, data = df_ok), error = function(e) NULL)
    if (is.null(lm_fit)) next

    res_norm <- tryCatch(
      shapiro.test(residuals(lm_fit))$p.value >= 0.05,
      error = function(e) FALSE
    )

    if (isTRUE(res_norm)) {
      aov_res <- tryCatch(aov(fmla, data = df_ok), error = function(e) NULL)

      if (!is.null(aov_res)) {
        pval <- tryCatch(
          summary(aov_res)[[1]][["Pr(>F)"]][1],
          error = function(e) NA_real_
        )

        results[[col]] <- data.frame(
          variable = col,
          type = "numeric",
          test = "ANOVA",
          statistic = NA,
          p_value = pval,
          stringsAsFactors = FALSE
        )
      }
    } else {
      kw <- tryCatch(kruskal.test(fmla, data = df_ok), error = function(e) NULL)

      if (!is.null(kw)) {
        results[[col]] <- data.frame(
          variable = col,
          type = "numeric",
          test = "Kruskal-Wallis",
          statistic = unname(kw$statistic),
          p_value = kw$p.value,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  if (length(results) == 0) {
    return(list(
      results_df = data.frame(),
      sig_df = data.frame(),
      sig_vars = character(0),
      present_vars = clinical_vars,
      merged_df = merged_df
    ))
  }

  results_df <- dplyr::bind_rows(results) |>
    dplyr::group_by(type) |>
    dplyr::mutate(adj_p = p.adjust(p_value, method = "BH")) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      cancer_type = cancer_type,
      modality = modality,
      feature_selection = feature_selection,
      clustering_method = clustering_method
    )

  sig_df <- results_df |>
    dplyr::filter(is.finite(adj_p) & adj_p <= alpha) |>
    dplyr::arrange(type, adj_p)

  list(
    results_df = results_df,
    sig_df = sig_df,
    sig_vars = unique(sig_df$variable),
    present_vars = clinical_vars,
    merged_df = merged_df
  )
}

clinical_associations_matrix <- function(
  cancer_types,
  modality,
  feature_selection,
  clustering_method,
  clinical_data,
  sample_id_fun = identity,
  exclude_cols = character(0),
  alpha = 0.05,
  verbose = TRUE,
  cluster_root = "clusters",
  modality_dirs = NULL,
  availability_map = NULL,
  output_file = NULL,
  filter_significant = TRUE
) {
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Package `dplyr` must be installed.")

  if (is.null(availability_map)) {
    availability_map <- build_clinical_availability_map(
      cancer_types = cancer_types,
      clinical_data = clinical_data,
      exclude_cols = exclude_cols,
      verbose = verbose
    )
  }

  clinical_vars_by_ct <- availability_map$clinical_vars_by_ct
  all_clinical_vars <- availability_map$all_clinical_vars

  sig_list <- setNames(vector("list", length(cancer_types)), cancer_types)
  results_list <- setNames(vector("list", length(cancer_types)), cancer_types)

  for (ct in cancer_types) {
    res <- tryCatch(
      clinical_associations_one_cancer(
        cancer_type = ct,
        modality = modality,
        feature_selection = feature_selection,
        clustering_method = clustering_method,
        clinical_data = clinical_data[[ct]],
        sample_id_fun = sample_id_fun,
        exclude_cols = exclude_cols,
        alpha = alpha,
        verbose = verbose,
        cluster_root = cluster_root,
        modality_dirs = modality_dirs
      ),
      error = function(e) {
        if (isTRUE(verbose)) message("Error for ", ct, ": ", e$message)

        list(
          results_df = data.frame(),
          sig_df = data.frame(),
          sig_vars = character(0),
          merged_df = data.frame()
        )
      }
    )

    sig_list[[ct]] <- res$sig_vars
    results_list[[ct]] <- res$results_df
  }

  presence_mat <- matrix(
    NA_integer_,
    nrow = length(all_clinical_vars),
    ncol = length(cancer_types),
    dimnames = list(all_clinical_vars, cancer_types)
  )

  for (ct in cancer_types) {
    clinically_present <- clinical_vars_by_ct[[ct]]
    presence_mat[clinically_present, ct] <- 0L

    sigs <- intersect(sig_list[[ct]], clinically_present)
    if (length(sigs) > 0) presence_mat[sigs, ct] <- 1L
  }

  presence_df <- data.frame(
    variable = rownames(presence_mat),
    presence_mat,
    check.names = FALSE
  )

  rownames(presence_df) <- NULL

  if (isTRUE(filter_significant)) {
    cancer_cols <- setdiff(names(presence_df), "variable")
    keep <- rowSums(presence_df[cancer_cols] == 1L, na.rm = TRUE) > 0
    presence_df <- presence_df[keep, , drop = FALSE]
  }

  all_results_df <- dplyr::bind_rows(results_list)

  if (!is.null(output_file)) {
    dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(presence_df, output_file, row.names = FALSE)
  }

  list(
    presence_df = presence_df,
    all_results_df = all_results_df,
    sig_list = sig_list,
    clinical_vars_by_ct = clinical_vars_by_ct
  )
}

plot_clinical_heatmap <- function(
  expression_df = NULL,
  indegrees_df = NULL,
  edge_weights_df = NULL,
  cancer_types = NULL,
  output_file = NULL,
  width = 6,
  height = 4.5,
  dpi = 300
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package `ggplot2` must be installed.")
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Package `dplyr` must be installed.")
  if (!requireNamespace("tidyr", quietly = TRUE)) stop("Package `tidyr` must be installed.")

  prep_df <- function(df, label) {
    if (is.null(df)) return(NULL)

    df <- as.data.frame(df, check.names = FALSE)

    if (!("variable" %in% names(df))) {
      stop(label, " input must contain a `variable` column.")
    }

    num_cols <- setdiff(names(df), c("variable", "Explained"))

    for (nm in num_cols) {
      df[[nm]] <- suppressWarnings(as.numeric(as.character(df[[nm]])))
    }

    df
  }

  edge_df <- prep_df(edge_weights_df, "Edge Weights")
  ind_df <- prep_df(indegrees_df, "Indegrees")
  expr_df <- prep_df(expression_df, "Expression")

  if (is.null(edge_df) && is.null(ind_df) && is.null(expr_df)) {
    stop("Provide at least one of expression_df, indegrees_df, or edge_weights_df.")
  }

  get_cancer_cols <- function(df) {
    if (is.null(df)) return(character(0))
    setdiff(names(df), c("variable", "Explained"))
  }

  all_types <- Reduce(
    union,
    list(
      get_cancer_cols(edge_df),
      get_cancer_cols(ind_df),
      get_cancer_cols(expr_df)
    )
  )

  if (!is.null(cancer_types)) {
    all_types <- intersect(all_types, cancer_types)
  } else {
    all_types <- sort(all_types)
  }

  if (length(all_types) == 0) stop("No cancer-type columns available after filtering.")

  add_missing_cols <- function(df) {
    if (is.null(df)) return(NULL)

    existing <- get_cancer_cols(df)

    for (nm in setdiff(all_types, existing)) {
      df[[nm]] <- NA_real_
    }

    keep <- c("variable", if ("Explained" %in% names(df)) "Explained", all_types)
    df[, keep, drop = FALSE]
  }

  edge_df <- add_missing_cols(edge_df)
  ind_df <- add_missing_cols(ind_df)
  expr_df <- add_missing_cols(expr_df)

  if (!is.null(edge_df)) edge_df$.in_edg <- rep(TRUE, nrow(edge_df))
  if (!is.null(ind_df)) ind_df$.in_ind <- rep(TRUE, nrow(ind_df))
  if (!is.null(expr_df)) expr_df$.in_gex <- rep(TRUE, nrow(expr_df))

  pivot_modality <- function(df, value_name) {
    if (is.null(df) || nrow(df) == 0) return(NULL)

    tidyr::pivot_longer(
      df,
      cols = dplyr::all_of(all_types),
      names_to = "cancer_type",
      values_to = value_name
    )
  }

  long_edg <- pivot_modality(edge_df, "edge_present")
  long_ind <- pivot_modality(ind_df, "indg_present")
  long_gex <- pivot_modality(expr_df, "gex_present")

  variables <- c(
    if (!is.null(edge_df)) edge_df$variable else character(0),
    if (!is.null(ind_df)) ind_df$variable else character(0),
    if (!is.null(expr_df)) expr_df$variable else character(0)
  )

  row_levels <- sort(unique(variables))
  label_list <- list()

  for (df in list(edge_df, ind_df, expr_df)) {
    if (!is.null(df) && "Explained" %in% names(df)) {
      label_list[[length(label_list) + 1L]] <- df[, c("variable", "Explained"), drop = FALSE]
    }
  }

  if (length(label_list) > 0) {
    label_df <- dplyr::bind_rows(label_list) |>
      dplyr::distinct(variable, .keep_all = TRUE)
  } else {
    label_df <- data.frame(
      variable = row_levels,
      Explained = row_levels,
      stringsAsFactors = FALSE
    )
  }

  missing_label <- is.na(label_df$Explained) | label_df$Explained == ""
  label_df$Explained[missing_label] <- label_df$variable[missing_label]

  long <- tidyr::expand_grid(variable = row_levels, cancer_type = all_types)

  if (!is.null(long_edg)) {
    edge_keep <- long_edg |>
      dplyr::select(variable, cancer_type, edge_present, dplyr::any_of(".in_edg"))

    long <- dplyr::left_join(
      long, edge_keep,
      by = c("variable", "cancer_type")
    )
  }

  if (!is.null(long_ind)) {
    ind_keep <- long_ind |>
      dplyr::select(variable, cancer_type, indg_present, dplyr::any_of(".in_ind"))

    long <- dplyr::left_join(
      long, ind_keep,
      by = c("variable", "cancer_type")
    )
  }

  if (!is.null(long_gex)) {
    expr_keep <- long_gex |>
      dplyr::select(variable, cancer_type, gex_present, dplyr::any_of(".in_gex"))

    long <- dplyr::left_join(
      long, expr_keep,
      by = c("variable", "cancer_type")
    )
  }

  long <- dplyr::left_join(long, label_df, by = "variable")

  if (!(".in_edg" %in% names(long))) long$.in_edg <- FALSE
  if (!(".in_ind" %in% names(long))) long$.in_ind <- FALSE
  if (!(".in_gex" %in% names(long))) long$.in_gex <- FALSE

  if (!("edge_present" %in% names(long))) long$edge_present <- NA_real_
  if (!("indg_present" %in% names(long))) long$indg_present <- NA_real_
  if (!("gex_present" %in% names(long))) long$gex_present <- NA_real_

  long <- long |>
    dplyr::mutate(
      .in_edg = tidyr::replace_na(.in_edg, FALSE),
      .in_ind = tidyr::replace_na(.in_ind, FALSE),
      .in_gex = tidyr::replace_na(.in_gex, FALSE),
      edge_value = dplyr::if_else(.in_edg, edge_present, 0),
      indg_value = dplyr::if_else(.in_ind, indg_present, 0),
      gex_value = dplyr::if_else(.in_gex, gex_present, 0),
      status = dplyr::case_when(
        is.na(edge_value) | is.na(indg_value) | is.na(gex_value) ~ "missing",
        edge_value == 1 & indg_value == 0 ~ "edge_only",
        edge_value == 0 & indg_value == 1 ~ "indegree_only",
        edge_value == 1 & indg_value == 1 ~ "both",
        edge_value == 0 & indg_value == 0 ~ "none"
      ),
      gex_mark = dplyr::if_else(gex_value == 1, "X", "")
    )

  explained_levels <- sort(unique(long$Explained))

  long <- long |>
    dplyr::mutate(
      variable_label = factor(Explained, levels = rev(explained_levels)),
      cancer_type = factor(cancer_type, levels = all_types),
      status = factor(
        status,
        levels = c("edge_only", "indegree_only", "both", "none", "missing")
      )
    )

  status_map <- c(
    edge_only = "Edge Weights only",
    indegree_only = "Indegrees only",
    both = "Both",
    none = "None",
    missing = "Missing"
  )

  present_status <- names(status_map)[
    vapply(
      names(status_map),
      function(x) any(as.character(long$status) == x, na.rm = TRUE),
      logical(1)
    )
  ]

  legend_items <- unname(status_map[present_status])

  if (any(long$gex_value == 1, na.rm = TRUE)) {
    legend_items <- c(legend_items, "Gene expression")
  }

  legend_df <- data.frame(
    cancer_type = factor(rep(all_types[1], length(legend_items)), levels = all_types),
    variable_label = factor(rep(explained_levels[1], length(legend_items)), levels = rev(explained_levels)),
    legend_item = factor(legend_items, levels = legend_items)
  )

  shape_values <- setNames(rep(22, length(legend_items)), legend_items)

  if ("Gene expression" %in% legend_items) {
    shape_values["Gene expression"] <- 4
  }

  legend_fill <- c(
    "Edge Weights only" = "#F8766D",
    "Indegrees only" = "#619CFF",
    "Both" = "#c85cb4",
    "None" = "#e3e3e3",
    "Missing" = "#fff0f0",
    "Gene expression" = "white"
  )

  legend_colour <- c(
    "Edge Weights only" = "#F8766D",
    "Indegrees only" = "#619CFF",
    "Both" = "#c85cb4",
    "None" = "#bdbdbd",
    "Missing" = "#e3cccc",
    "Gene expression" = "black"
  )

  p <- ggplot2::ggplot(
    long,
    ggplot2::aes(x = cancer_type, y = variable_label, fill = status)
  ) +
    ggplot2::geom_tile(color = "white", linewidth = 0.4) +
    ggplot2::geom_text(
      ggplot2::aes(label = gex_mark),
      size = 1.8,
      color = "black",
      fontface = "plain"
    ) +
    ggplot2::geom_point(
      data = legend_df,
      ggplot2::aes(x = cancer_type, y = variable_label, shape = legend_item),
      inherit.aes = FALSE,
      alpha = 0
    ) +
    ggplot2::scale_fill_manual(
      values = c(
        edge_only = "#F8766D",
        indegree_only = "#619CFF",
        both = "#c85cb4",
        none = "#e3e3e3",
        missing = "#fff0f0"
      ),
      guide = "none"
    ) +
    ggplot2::scale_shape_manual(
      name = "Significant association in",
      values = shape_values,
      breaks = legend_items,
      labels = legend_items
    ) +
    ggplot2::scale_x_discrete(drop = FALSE) +
    ggplot2::scale_y_discrete(drop = FALSE) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::guides(
      shape = ggplot2::guide_legend(
        nrow = 2,
        byrow = TRUE,
        override.aes = list(
          alpha = 1,
          size = 2.6,
          stroke = 0.45,
          colour = unname(legend_colour[legend_items]),
          fill = unname(legend_fill[legend_items])
        )
      )
    ) +
    ggplot2::theme_bw(base_size = 8) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = 7, angle = 90, vjust = 0.5, hjust = 1),
      axis.text.y = ggplot2::element_text(size = 6.5),
      legend.title = ggplot2::element_text(size = 6, face = "plain"),
      legend.text = ggplot2::element_text(size = 6, face = "plain"),
      legend.key.height = grid::unit(0.2, "cm"),
      legend.key.width = grid::unit(0.2, "cm"),
      legend.spacing.y = grid::unit(0.03, "cm"),
      legend.key.spacing.x = grid::unit(0.03, "cm"),
      legend.position = "bottom",
      panel.grid = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank()
    )

  if (!is.null(output_file)) {
    dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)

    ggplot2::ggsave(
      filename = output_file,
      plot = p,
      width = width,
      height = height,
      dpi = dpi
    )
  }

  invisible(p)
}