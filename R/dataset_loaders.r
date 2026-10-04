load_modalities <- function(config, ct) {
  matrices <- list()

  if ("Expression" %in% config$modalities) {
    matrices$Expression <- read_expression(
      path = config$paths$expression,
      cancer_type = ct,
      file_pattern = config$expression$file_pattern,
      gene_col = config$expression$gene_col,
      drop_cols = config$expression$drop_cols,
      remove_normals = config$expression$remove_normals,
      normal_pattern = config$expression$normal_pattern
    )
  }

  if ("Indegrees" %in% config$modalities) {
    matrices$Indegrees <- read_indegrees(
      path = config$paths$indegrees,
      cancer_type = ct,
      file_pattern = config$indegrees$file_pattern,
      object_name = config$indegrees$object_name
    )
  }

  if ("Edge Weights" %in% config$modalities) {
    matrices$`Edge Weights` <- read_edge_weights(
      path = config$paths$edge_weights,
      cancer_type = ct,
      file_pattern = config$edge_weights$file_pattern,
      object_name = config$edge_weights$object_name
    )
  }

  validation_args <- list()

  if ("Expression" %in% names(matrices)) validation_args$expr <- matrices$Expression
  if ("Indegrees" %in% names(matrices)) validation_args$indegrees <- matrices$Indegrees
  if ("Edge Weights" %in% names(matrices)) validation_args$edge_weights <- matrices$`Edge Weights`

  validated <- do.call(validate_modalities, validation_args)

  if (!is.null(validated$expr)) matrices$Expression <- validated$expr
  if (!is.null(validated$indegrees)) matrices$Indegrees <- validated$indegrees
  if (!is.null(validated$edge_weights)) matrices$`Edge Weights` <- validated$edge_weights

  matrices
}

load_clinical_data <- function(config, cancer_types) {
  clinical_data <- setNames(vector("list", length(cancer_types)), cancer_types)

  for (ct in cancer_types) {
    clinical_file <- file.path(
      config$paths$clinical,
      sprintf(config$clinical$file, ct)
    )

    if (!file.exists(clinical_file)) {
      stop("Clinical file not found: ", clinical_file)
    }

    clinical <- utils::read.delim(
      clinical_file,
      sep = "\t",
      header = TRUE,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )

    required_cols <- c(
      "sample_id",
      "survival_time",
      "survival_event"
    )

    missing_cols <- setdiff(
      required_cols,
      colnames(clinical)
    )

    if (length(missing_cols) > 0) {
      stop(
        "Clinical file is missing required columns: ",
        paste(missing_cols, collapse = ", ")
      )
    }

    if (anyDuplicated(clinical$sample_id)) {
      stop("Clinical file contains duplicated sample IDs: ", clinical_file)
    }

    clinical_data[[ct]] <- clinical
  }

  clinical_data
}