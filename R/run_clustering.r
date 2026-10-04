run_clustering <- function(
  mat,
  cancer_type,
  top_value_method = "SD",
  partition_method = "kmeans",
  cores = 1,
  scale_rows = FALSE,
  report_dir = NULL,
  make_report = FALSE,
  save_clusters = FALSE,
  cluster_dir = NULL,
  cluster_k = NULL,
  seed = 12345,
  cluster_prefix = NULL
) {
  set.seed(seed)

  if (!requireNamespace("cola", quietly = TRUE)) {
    stop("Package `cola` must be installed to use `run_clustering()`.")
  }

  if (!is.matrix(mat) && !is.data.frame(mat)) {
    stop("`mat` must be a matrix or data frame.")
  }

  mat <- as.matrix(mat)

  if (nrow(mat) == 0 || ncol(mat) == 0) stop("`mat` is empty.")
  if (is.null(rownames(mat))) stop("`mat` must have row names.")
  if (is.null(colnames(mat))) stop("`mat` must have column names.")

  if (is.null(cluster_prefix)) {
    cluster_prefix <- paste0(
      "clusters_",
      top_value_method,
      "_",
      partition_method,
      "_",
      cancer_type
    )
  }

  mat <- cola::adjust_matrix(mat)

  rl <- cola::consensus_partition(
    mat,
    top_value_method = top_value_method,
    partition_method = partition_method,
    cores = cores
  )

  if (make_report) {
    if (is.null(report_dir)) {
      stop("`report_dir` must be provided when `make_report = TRUE`.")
    }

    dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
    cola::cola_report(rl, output_dir = report_dir)
  }

  if (save_clusters) {
    if (is.null(cluster_dir) || !nzchar(cluster_dir)) {
      stop("`cluster_dir` must be provided when `save_clusters = TRUE`.")
    }

    dir.create(cluster_dir, recursive = TRUE, showWarnings = FALSE)
    out_file <- file.path(cluster_dir, paste0(cluster_prefix, ".txt"))

    save_cluster_assignments(
      cp = rl,
      file = out_file,
      k = cluster_k,
      use_best_k = is.null(cluster_k)
    )
  }

  rl
}