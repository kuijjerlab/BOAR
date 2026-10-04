read_expression <- function(
  path,
  cancer_type,
  file_pattern = "%s_expression.tsv",
  gene_col = "gene_name",
  drop_cols = "gene_id",
  remove_normals = TRUE,
  normal_pattern = "-11A-",
  min_samples = 3
) {
  file_name <- sprintf(file_pattern, cancer_type)
  file_path <- file.path(path, file_name)

  if (!file.exists(file_path)) stop("Expression file not found: ", file_path)

  expr <- read.table(
    file_path,
    sep = "\t",
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  if (!(gene_col %in% colnames(expr))) {
    stop("Gene column `", gene_col, "` not found in expression file: ", file_path)
  }

  expr[[gene_col]] <- make.unique(as.character(expr[[gene_col]]))
  rownames(expr) <- expr[[gene_col]]

  cols_to_drop <- intersect(unique(c(gene_col, drop_cols)), colnames(expr))
  expr <- expr[, setdiff(colnames(expr), cols_to_drop), drop = FALSE]

  expr[] <- lapply(expr, as.numeric)
  expr_mat <- as.matrix(expr)

  if (nrow(expr_mat) == 0 || ncol(expr_mat) == 0) {
    stop("Expression matrix is empty after conversion.")
  }

  keep_rows <- apply(expr_mat, 1, sd, na.rm = TRUE) > 0
  keep_rows[is.na(keep_rows)] <- FALSE
  expr_mat <- expr_mat[keep_rows, , drop = FALSE]

  colnames(expr_mat) <- gsub("\\.", "-", colnames(expr_mat))

  if (isTRUE(remove_normals) && !is.null(normal_pattern)) {
    expr_mat <- expr_mat[, !grepl(normal_pattern, colnames(expr_mat)), drop = FALSE]
  }

  if (ncol(expr_mat) < min_samples) {
    stop("Too few samples in expression matrix for ", cancer_type, ".")
  }

  expr_mat
}

read_indegrees <- function(
  path,
  cancer_type,
  file_pattern = "indegree_%s.RData",
  object_name = "ind"
) {
  file_name <- sprintf(file_pattern, cancer_type)
  file_path <- file.path(path, file_name)

  if (!file.exists(file_path)) stop("Indegree file not found: ", file_path)

  env <- new.env()
  load(file_path, envir = env)

  if (!(object_name %in% ls(env))) {
    stop("Object `", object_name, "` not found in: ", file_path)
  }

  obj <- env[[object_name]]
  if (inherits(obj, "tbl_df")) obj <- as.data.frame(obj)

  if (is.data.frame(obj) && !is.numeric(obj[[1]])) {
    rownames(obj) <- make.unique(as.character(obj[[1]]))
    obj <- obj[, -1, drop = FALSE]
  }

  mat <- as.matrix(obj)
  storage.mode(mat) <- "numeric"

  if (nrow(mat) == 0 || ncol(mat) == 0) {
    stop("Indegree matrix is empty after conversion.")
  }

  mat
}

read_edge_weights <- function(
  path,
  cancer_type,
  file_pattern = "%s_edge_weights.RData",
  object_name = NULL
) {
  file_name <- sprintf(file_pattern, cancer_type)
  file_path <- file.path(path, file_name)

  if (!file.exists(file_path)) stop("Edge-weight file not found: ", file_path)

  env <- new.env()
  load(file_path, envir = env)
  objs <- ls(env)

  if (length(objs) == 0) stop("No objects found in: ", file_path)

  if (is.null(object_name)) {
    object_name <- objs[1]
    message("Using object: ", object_name)
  }

  if (!(object_name %in% objs)) {
    stop("Object `", object_name, "` not found in: ", file_path)
  }

  obj <- env[[object_name]]
  if (inherits(obj, "tbl_df")) obj <- as.data.frame(obj)

  if (is.data.frame(obj) && !is.numeric(obj[[1]])) {
    rownames(obj) <- make.unique(as.character(obj[[1]]))
    obj <- obj[, -1, drop = FALSE]
  }

  mat <- as.matrix(obj)
  storage.mode(mat) <- "numeric"

  if (nrow(mat) == 0 || ncol(mat) == 0) {
    stop("Edge-weight matrix is empty after conversion.")
  }

  mat
}

validate_modalities <- function(
  expr = NULL,
  indegrees = NULL,
  edge_weights = NULL
) {
  if (is.null(expr) && is.null(indegrees) && is.null(edge_weights)) {
    stop("At least one of `expr`, `indegrees`, or `edge_weights` must be provided.")
  }

  check_numeric_matrix <- function(x, name) {
    if (is.null(x)) return(NULL)

    if (!is.matrix(x) && !is.data.frame(x)) {
      stop("`", name, "` must be a matrix or data frame.")
    }

    x <- as.matrix(x)
    storage.mode(x) <- "numeric"

    if (nrow(x) == 0 || ncol(x) == 0) stop("`", name, "` is empty.")
    if (is.null(rownames(x))) stop("`", name, "` must have row names.")
    if (anyDuplicated(rownames(x)) > 0) stop("`", name, "` has duplicated row names.")
    if (is.null(colnames(x))) stop("`", name, "` must have column names.")

    x
  }

  expr <- check_numeric_matrix(expr, "expr")
  indegrees <- check_numeric_matrix(indegrees, "indegrees")
  edge_weights <- check_numeric_matrix(edge_weights, "edge_weights")

  list(
    expr = expr,
    indegrees = indegrees,
    edge_weights = edge_weights
  )
}