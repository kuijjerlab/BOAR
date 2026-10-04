clinical_survival_analysis <- function(
  cancer_type,
  modality,
  feature_selection,
  clustering_method,
  clinical_data,
  sample_id_fun = identity,
  cluster_root = "clusters",
  modality_dirs = NULL,
  output_file = NULL,
  pvalue_file = "results/survival/survival_pvalues.csv",
  verbose = TRUE) {
  vmsg <- function(...) if (isTRUE(verbose)) message(...)

  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("Package `survival` must be installed.")
  }

  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package `dplyr` must be installed.")
  }

  vmsg("\n=== Survival: ", cancer_type, " / ", modality, " ===")

  classes_file <- find_cluster_file_for_analysis(
    cancer_type = cancer_type,
    modality = modality,
    feature_selection = feature_selection,
    clustering_method = clustering_method,
    cluster_root = cluster_root,
    modality_dirs = modality_dirs)

  if (length(classes_file) != 1 || is.na(classes_file) || !file.exists(classes_file)) {
    vmsg("    (skip: no matching cluster file found)")
    return(list(
      p_value = NA_real_,
      merged_df = NULL,
      output_file = NULL,
      pvalue_file = pvalue_file))
  }

  classes_df <- read_classes_txt(classes_file)

  if (is.null(classes_df) || nrow(classes_df) == 0) {
    vmsg("    (skip: could not parse classes)")
    return(list(
      p_value = NA_real_,
      merged_df = NULL,
      output_file = NULL,
      pvalue_file = pvalue_file))
  }

  classes_df <- classes_df |>
    dplyr::distinct(sample, .keep_all = TRUE)

  classes_df$clinical_id <- as.character(sample_id_fun(classes_df$sample))
  classes_df$class <- as.character(classes_df$class)

  classes_df <- classes_df |>
    dplyr::distinct(clinical_id, .keep_all = TRUE)

  clinical <- as.data.frame(clinical_data, stringsAsFactors = FALSE)

  if (nrow(clinical) == 0) {
    stop("Clinical data are empty for ", cancer_type, ".")
  }

  required_cols <- c("sample_id", "survival_time", "survival_event")
  missing_cols <- setdiff(required_cols, colnames(clinical))

  if (length(missing_cols) > 0) {
    stop(
      "Missing required clinical columns for ",
      cancer_type,
      ": ",
      paste(missing_cols, collapse = ", "))
  }

  clinical$clinical_id <- as.character(clinical$sample_id)

  clinical <- clinical |>
    dplyr::distinct(clinical_id, .keep_all = TRUE)

  merged_df <- merge(
    classes_df,
    clinical,
    by = "clinical_id",
    all.x = TRUE,
    sort = FALSE)

  merged_df$time <- suppressWarnings(as.numeric(merged_df$survival_time))
  merged_df$event <- suppressWarnings(as.integer(merged_df$survival_event))

  invalid_events <- unique(
    merged_df$event[
      !is.na(merged_df$event) &
        !(merged_df$event %in% c(0L, 1L))])

  if (length(invalid_events) > 0) {
    stop(
      "`survival_event` must contain only 0, 1, or NA. Found: ",
      paste(invalid_events, collapse = ", "))
  }

  merged_df <- merged_df[
    !is.na(merged_df$time) &
      !is.na(merged_df$event) &
      !is.na(merged_df$class),
    ,
    drop = FALSE]

  merged_df <- merged_df[
    merged_df$time > 0,
    ,
    drop = FALSE]

  if (nrow(merged_df) < 2) {
    vmsg("    (skip: fewer than 2 samples)")
    return(list(
      p_value = NA_real_,
      merged_df = merged_df,
      output_file = NULL,
      pvalue_file = pvalue_file))
  }

  if (length(unique(merged_df$class)) < 2) {
    vmsg("    (skip: fewer than 2 clusters)")
    return(list(
      p_value = NA_real_,
      merged_df = merged_df,
      output_file = NULL,
      pvalue_file = pvalue_file))
  }

  if (sum(merged_df$event, na.rm = TRUE) == 0) {
    vmsg("    (skip: no observed events)")
    return(list(
      p_value = NA_real_,
      merged_df = merged_df,
      output_file = NULL,
      pvalue_file = pvalue_file))
  }

  surv_diff <- survival::survdiff(
    survival::Surv(time, event) ~ class,
    data = merged_df)

  pval <- stats::pchisq(
    surv_diff$chisq,
    df = length(surv_diff$n) - 1,
    lower.tail = FALSE)

  vmsg("    p-value: ", signif(pval, 4))

  surv_fit <- survival::survfit(
    survival::Surv(time, event) ~ class,
    data = merged_df)

  if (!is.null(output_file)) {
    grDevices::png(
      output_file,
      width = 1800,
      height = 1500,
      res = 300)

    plot(
      surv_fit,
      xlab = "Time",
      ylab = "Survival probability",
      mark.time = TRUE)

    legend(
      "bottomleft",
      legend = names(surv_fit$strata),
      lty = seq_along(surv_fit$strata),
      bty = "n")

    title(main = paste0(
      cancer_type,
      " - ",
      modality,
      "\np = ",
      signif(pval, 3)))

    grDevices::dev.off()
  }

  modality_out <- tolower(gsub(" ", "_", modality))

  pval_df <- data.frame(
    modality = modality_out,
    cancer_type = cancer_type,
    p_value = pval,
    stringsAsFactors = FALSE)

  if (file.exists(pvalue_file)) {
    existing <- utils::read.csv(pvalue_file, stringsAsFactors = FALSE)

    if ("adj_p_value" %in% names(existing)) {
      existing$adj_p_value <- NULL
    }

    existing <- existing[
      !(existing$cancer_type == cancer_type & existing$modality == modality_out),
      ,
      drop = FALSE]

    survival_pvalues <- rbind(existing, pval_df)
  } else {
    survival_pvalues <- pval_df
  }

  survival_pvalues$adj_p_value <- stats::p.adjust(
    survival_pvalues$p_value,
    method = "BH")

  utils::write.csv(survival_pvalues, pvalue_file, row.names = FALSE)

  vmsg("    Saved survival p-values to: ", pvalue_file)

  invisible(list(
    p_value = pval,
    merged_df = merged_df,
    survival_fit = surv_fit,
    output_file = output_file,
    pvalue_file = pvalue_file))
}

plot_survival_summary <- function(
  input_file = "results/survival/survival_pvalues.csv",
  output_file = "results/survival/survival_summary.png",
  width = 3.5,
  height = 2.5,
  dpi = 300) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package `ggplot2` must be installed.")
  }

  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package `dplyr` must be installed.")
  }

  if (!file.exists(input_file)) {
    stop("Input file not found: ", input_file)
  }

  df <- utils::read.csv(input_file, stringsAsFactors = FALSE)

  required_cols <- c("modality", "cancer_type", "p_value")
  missing_cols <- setdiff(required_cols, names(df))

  if (length(missing_cols) > 0) {
    stop(
      "Input file is missing required columns: ",
      paste(missing_cols, collapse = ", "))
  }

  long_df <- df |>
    dplyr::mutate(
      p_value = suppressWarnings(as.numeric(p_value)),
      p_value = ifelse(p_value == 0, 1e-16, p_value)) |>
    dplyr::filter(
      !is.na(p_value),
      p_value > 0) |>
    dplyr::mutate(
      neglogp = -log10(p_value),
      signif = p_value < 0.05,
      modality = factor(
        modality,
        levels = c("expression", "indegrees", "edge_weights"),
        labels = c("Expression", "Indegrees", "Edge Weights")),
      plot_colour = dplyr::case_when(
        modality == "Expression" & signif ~ "Expression",
        modality == "Indegrees" & signif ~ "Indegrees",
        modality == "Edge Weights" & signif ~ "Edge Weights",
        modality == "Expression" & !signif ~ "Expression_ns",
        modality == "Indegrees" & !signif ~ "Indegrees_ns",
        modality == "Edge Weights" & !signif ~ "Edge Weights_ns")) |>
    dplyr::group_by(cancer_type) |>
    dplyr::mutate(max_sig = max(neglogp, na.rm = TRUE)) |>
    dplyr::ungroup() |>
    dplyr::mutate(cancer_type = reorder(cancer_type, -max_sig))

  if (nrow(long_df) == 0) {
    stop("No valid survival p-values to plot.")
  }

  p <- ggplot2::ggplot(
    long_df,
    ggplot2::aes(
      x = cancer_type,
      y = neglogp,
      color = plot_colour)) +
    ggplot2::geom_point(size = 1.5) +
    ggplot2::scale_color_manual(
      values = c(
        "Expression" = "#4DAF4A",
        "Indegrees" = "#619CFF",
        "Edge Weights" = "#F8766D",
        "Expression_ns" = "#CFE8CF",
        "Indegrees_ns" = "#D6E5FF",
        "Edge Weights_ns" = "#FAD4CF"),
      breaks = c(
        "Edge Weights",
        "Indegrees",
        "Expression",
        "Edge Weights_ns",
        "Indegrees_ns",
        "Expression_ns"),
      labels = c(
        "Edge weights, p < 0.05",
        "Indegrees, p < 0.05",
        "Expression, p < 0.05",
        "Edge weights, p ≥ 0.05",
        "Indegrees, p ≥ 0.05",
        "Expression, p ≥ 0.05"),
      name = NULL) +
    ggplot2::guides(
      color = ggplot2::guide_legend(
        nrow = 2,
        byrow = TRUE)) +
    ggplot2::labs(
      x = NULL,
      y = "-log10(p-value)") +
    ggplot2::theme_bw(base_size = 8) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = 45,
        hjust = 1,
        size = 5),
      axis.text.y = ggplot2::element_text(size = 6),
      axis.title.y = ggplot2::element_text(size = 7),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = 6),
      legend.key.height = grid::unit(0.15, "cm"),
      legend.key.width = grid::unit(0.15, "cm"),
      legend.spacing.x = grid::unit(0.05, "cm"),
      legend.spacing.y = grid::unit(0.03, "cm"),
      legend.key.spacing.x = grid::unit(0.18, "cm"))

  if (!is.null(output_file)) {
    ggplot2::ggsave(
      filename = output_file,
      plot = p,
      width = width,
      height = height,
      dpi = dpi)
  }

  invisible(list(
    plot = p,
    data = long_df))
}