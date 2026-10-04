read_exclude_cols <- function(path) {
  if (is.null(path)) return(character(0))
  if (!file.exists(path)) stop("Exclusion file not found: ", path)

  x <- readLines(path, warn = FALSE)
  x <- paste(x, collapse = " ")
  x <- gsub('"', "", x)
  x <- unlist(strsplit(x, ","))
  x <- trimws(x)

  x[x != ""]
}

read_classes_txt <- function(path) {
  if (!file.exists(path)) stop("Classes file not found: ", path)

  lines <- trimws(readLines(path, warn = FALSE))
  lines <- lines[nchar(lines) > 0]

  cur_class <- NULL
  out <- list()
  idx <- 0L

  for (ln in lines) {
    if (grepl("^\\d+$", ln)) {
      cur_class <- ln
      next
    }

    if (is.null(cur_class)) stop("Classes file must start with a class label.")

    idx <- idx + 1L
    out[[idx]] <- data.frame(
      sample = ln,
      class = cur_class,
      stringsAsFactors = FALSE
    )
  }

  if (idx == 0L) stop("No sample assignments found in: ", path)

  out <- do.call(rbind, out)
  out <- dplyr::distinct(out, sample, .keep_all = TRUE)
  out
}

find_cluster_file_for_analysis <- function(
  cancer_type,
  modality,
  feature_selection,
  clustering_method,
  cluster_root = "clusters",
  modality_dirs = NULL
) {
  if (!is.null(modality_dirs)) {
    if (is.null(names(modality_dirs)) || !(modality %in% names(modality_dirs))) {
      stop("`modality_dirs` must be a named list including: ", modality)
    }

    search_dir <- modality_dirs[[modality]]

    if (dir.exists(file.path(search_dir, cancer_type))) {
      search_dir <- file.path(search_dir, cancer_type)
    }
  } else {
    search_dir <- file.path(cluster_root, cancer_type, modality)
  }

  if (!dir.exists(search_dir)) return(NA_character_)

  expected_name <- paste0(
    "clusters_",
    feature_selection, "_",
    clustering_method, "_",
    cancer_type,
    ".txt"
  )

  expected_path <- file.path(search_dir, expected_name)
  if (file.exists(expected_path)) return(expected_path)

  available_files <- list.files(search_dir, pattern = "\\.txt$", full.names = FALSE)
  if (length(available_files) == 0) return(NA_character_)

  stop(
    "Cluster file not found.\n\n",
    "Expected:\n  ", expected_name, "\n\n",
    "In directory:\n  ", search_dir, "\n\n",
    "Available files:\n  ",
    paste(available_files, collapse = "\n  ")
  )
}

build_clinical_availability_map <- function(
  cancer_types,
  clinical_data,
  exclude_cols = character(0),
  verbose = TRUE
) {
  if (is.null(names(clinical_data))) stop("`clinical_data` must be a named list.")

  missing_cohorts <- setdiff(cancer_types, names(clinical_data))

  if (length(missing_cohorts) > 0) {
    stop("Clinical data missing for: ", paste(missing_cohorts, collapse = ", "))
  }

  clinical_by_ct <- setNames(vector("list", length(cancer_types)), cancer_types)
  clinical_vars_by_ct <- setNames(vector("list", length(cancer_types)), cancer_types)

  required_cols <- c("sample_id", "survival_time", "survival_event")

  for (ct in cancer_types) {
    if (isTRUE(verbose)) message("Loading clinical data for ", ct)

    clinical <- as.data.frame(
      clinical_data[[ct]],
      stringsAsFactors = FALSE
    )

    if (nrow(clinical) == 0) stop("No clinical data for ", ct)

    missing_cols <- setdiff(required_cols, colnames(clinical))

    if (length(missing_cols) > 0) {
      stop(
        "Clinical data for ", ct, " are missing required columns: ",
        paste(missing_cols, collapse = ", ")
      )
    }

    clinical_by_ct[[ct]] <- clinical

    clinical_vars_by_ct[[ct]] <- setdiff(
      colnames(clinical),
      c(required_cols, exclude_cols)
    )
  }

  all_clinical_vars <- sort(unique(unlist(clinical_vars_by_ct)))

  list(
    clinical_by_ct = clinical_by_ct,
    clinical_vars_by_ct = clinical_vars_by_ct,
    all_clinical_vars = all_clinical_vars
  )
}