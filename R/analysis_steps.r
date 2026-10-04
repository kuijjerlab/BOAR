prepare_analysis <- function(config) {
  modalities <- config$modalities
  results_dir <- config$paths$results

  dirs <- list(
    clusters = file.path(results_dir, "clusters"),
    reports = file.path(results_dir, "reports"),
    summary = file.path(results_dir, "summary"),
    figures = file.path(results_dir, "figures"),
    clinical = file.path(results_dir, "clinical"),
    survival = file.path(results_dir, "survival")
  )

  invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))

  modality_dirs <- setNames(
    lapply(modalities, function(x) file.path(dirs$clusters, x)),
    modalities
  )

  invisible(lapply(modality_dirs, dir.create, recursive = TRUE, showWarnings = FALSE))

  list(
    cancer_types = config$cancer_types,
    modalities = modalities,
    dirs = dirs,
    modality_dirs = modality_dirs,
    results_dir = results_dir
  )
}

run_clustering_for_cancer <- function(config, setup, ct) {
  modalities <- setup$modalities
  dirs <- setup$dirs
  modality_dirs <- setup$modality_dirs

  message("\n========================================")
  message("CANCER TYPE: ", ct)
  message("========================================")

  matrices <- load_modalities(config, ct)

  if (!all(modalities %in% names(matrices))) {
    missing_modalities <- setdiff(modalities, names(matrices))
    stop("Missing loaded modalities for ", ct, ": ", paste(missing_modalities, collapse = ", "))
  }

  for (modality in modalities) {
    message(modality, " samples: ", ncol(matrices[[modality]]))
  }

  all_summaries <- list()

  for (modality in modalities) {
    message("\n========================================")
    message(ct, " / MODALITY: ", modality)
    message("========================================")

    mat <- matrices[[modality]]
    modality_results <- list()

    for (fs in config$feature_methods) {
      for (cm in config$clustering_methods) {
        message("\n", ct, " / ", modality, " / ", fs, " / ", cm)

        report_dir <- file.path(
          dirs$reports,
          ct,
          gsub(" ", "_", modality),
          paste0(fs, "_", cm)
        )

        rl <- run_clustering(
          mat = mat,
          cancer_type = ct,
          top_value_method = fs,
          partition_method = cm,
          cores = config$cores,
          report_dir = report_dir,
          make_report = config$clustering$make_report,
          save_clusters = config$clustering$save_clusters,
          cluster_dir = modality_dirs[[modality]],
          cluster_k = NULL,
          seed = config$seed
        )

        modality_results[[paste(fs, cm, sep = ":")]] <- rl
      }
    }

    modality_summary <- summarize_clustering_results(
      result = modality_results,
      modality = modality
    )

    modality_summary$cancer_type <- ct
    all_summaries[[modality]] <- modality_summary

    rm(modality_results, modality_summary, mat)
    gc()
  }

  cancer_summary <- dplyr::bind_rows(all_summaries)

  cancer_summary <- cancer_summary[
    ,
    c("cancer_type", setdiff(names(cancer_summary), "cancer_type")),
    drop = FALSE
  ]

  summary_file <- file.path(dirs$summary, paste0(ct, "_combined.csv"))
  utils::write.csv(cancer_summary, summary_file, row.names = FALSE)

  message("Saved ", ct, " summary: ", summary_file)

  rm(matrices)
  gc()

  cancer_summary
}

combine_clustering_results <- function(config, setup, cancer_summaries) {
  dirs <- setup$dirs
  modalities <- setup$modalities

  combined_summary <- dplyr::bind_rows(cancer_summaries)

  combined_file <- file.path(dirs$summary, "combined_clustering_summary.csv")
  utils::write.csv(combined_summary, combined_file, row.names = FALSE)

  if (isTRUE(config$analysis$clustering_summary)) {
    p_clustering <- plot_clustering_summary(
      summary_df = combined_summary,
      modality_order = modalities
    )

    ggplot2::ggsave(
      filename = file.path(dirs$figures, "clustering_summary.png"),
      plot = p_clustering,
      width = 20,
      height = 7,
      dpi = 300
    )
  }

  best_result <- select_best_clustering_combination(combined_summary)
  best_fs <- as.character(best_result$best$feature_selection[1])
  best_cm <- as.character(best_result$best$clustering_method[1])

  message("\n========================================")
  message("OVERALL BEST COMBINATION: ", best_fs, " + ", best_cm)
  message("========================================")

  utils::write.csv(
    best_result$ranked,
    file.path(dirs$summary, "ranked_clustering_combinations.csv"),
    row.names = FALSE
  )

  list(
    combined_summary = combined_summary,
    best_result = best_result,
    best_fs = best_fs,
    best_cm = best_cm
  )
}

run_downstream_analysis <- function(config, setup, clustering) {
  cancer_types <- setup$cancer_types
  modalities <- setup$modalities
  dirs <- setup$dirs
  modality_dirs <- setup$modality_dirs
  best_fs <- clustering$best_fs
  best_cm <- clustering$best_cm

  best_k_results <- NULL
  balance_results <- NULL
  ari_results <- NULL

  if (isTRUE(config$analysis$best_k)) {
    best_k_results <- plot_k_choices_from_summary(
      cancer_types = cancer_types,
      summary_dir = dirs$summary,
      feature_selection = best_fs,
      clustering_method = best_cm,
      output_file = file.path(dirs$figures, "best_k.png")
    )
  }

  if (isTRUE(config$analysis$cluster_balance)) {
    balance_results <- plot_cluster_balance_from_txt(
      cancer_types = cancer_types,
      modalities = modalities,
      feature_selection = best_fs,
      clustering_method = best_cm,
      modality_dirs = modality_dirs,
      output_prefix = file.path(dirs$figures, "cluster_balance")
    )

    utils::write.csv(
      balance_results$balance,
      file.path(dirs$summary, "cluster_balance.csv"),
      row.names = FALSE
    )
  }

  best_metric_results <- plot_best_combination_metrics(
    summary_df = clustering$combined_summary,
    cancer_types = cancer_types,
    feature_selection = best_fs,
    clustering_method = best_cm,
    modalities = modalities,
    modality_dirs = modality_dirs,
    output_file = file.path(dirs$figures, "best_combination_metrics.png")
  )

  if (isTRUE(config$analysis$ari) && length(modalities) > 1) {
    ari_results <- compare_all_clusterings_ari(
      cancer_types = cancer_types,
      feature_selection = best_fs,
      clustering_method = best_cm,
      modality_dirs = modality_dirs,
      output_plot = file.path(dirs$figures, "ARI.png"),
      output_table = file.path(dirs$summary, "ARI.csv")
    )
  }

  list(
    best_k = best_k_results,
    balance = balance_results,
    best_metrics = best_metric_results,
    ari = ari_results
  )
}

run_clinical_workflow <- function(config, setup, clustering) {
  if (!isTRUE(config$analysis$clinical) && !isTRUE(config$analysis$survival)) return(NULL)

  cancer_types <- setup$cancer_types
  modalities <- setup$modalities
  dirs <- setup$dirs
  modality_dirs <- setup$modality_dirs
  best_fs <- clustering$best_fs
  best_cm <- clustering$best_cm

  clinical_data <- load_clinical_data(config, cancer_types)
  sample_id_fun <- config$clinical$sample_id_fun

  clinical_results <- NULL
  survival_results <- NULL
  availability_map <- NULL

  if (isTRUE(config$analysis$clinical)) {
    clinical_exclude_cols <- read_exclude_cols(config$paths$exclude_cols)

    clinical_exclude_cols <- unique(c(
      clinical_exclude_cols,
      "sample_id",
      "survival_time",
      "survival_event"
    ))

    availability_map <- build_clinical_availability_map(
      cancer_types = cancer_types,
      clinical_data = clinical_data,
      exclude_cols = clinical_exclude_cols,
      verbose = TRUE
    )

    clinical_results <- list()

    for (modality in modalities) {
      modality_key <- gsub(" ", "_", tolower(modality))

      clinical_results[[modality]] <- clinical_associations_matrix(
        cancer_types = cancer_types,
        modality = modality,
        feature_selection = best_fs,
        clustering_method = best_cm,
        clinical_data = clinical_data,
        sample_id_fun = sample_id_fun,
        exclude_cols = clinical_exclude_cols,
        availability_map = availability_map,
        modality_dirs = modality_dirs,
        output_file = file.path(
          dirs$clinical,
          paste0("clinical_", modality_key, ".csv")
        ),
        filter_significant = TRUE
      )

      utils::write.csv(
        clinical_results[[modality]]$all_results_df,
        file.path(
          dirs$clinical,
          paste0("clinical_", modality_key, "_all_tests.csv")
        ),
        row.names = FALSE
      )
    }

    has_associations <- any(vapply(
      clinical_results,
      function(x) nrow(x$presence_df) > 0,
      logical(1)
    ))

    if (has_associations) {
      heatmap_args <- list(
        cancer_types = cancer_types,
        output_file = file.path(dirs$figures, "clinical_associations.png")
      )

      if ("Expression" %in% names(clinical_results)) {
        heatmap_args$expression_df <- clinical_results$Expression$presence_df
      }

      if ("Indegrees" %in% names(clinical_results)) {
        heatmap_args$indegrees_df <- clinical_results$Indegrees$presence_df
      }

      if ("Edge Weights" %in% names(clinical_results)) {
        heatmap_args$edge_weights_df <- clinical_results$`Edge Weights`$presence_df
      }

      do.call(plot_clinical_heatmap, heatmap_args)
    } else {
      message("No clinical associations found.")
    }
  }

  if (isTRUE(config$analysis$survival)) {
    survival_pvalue_file <- file.path(
      dirs$survival,
      "survival_pvalues.csv"
    )

    if (file.exists(survival_pvalue_file)) file.remove(survival_pvalue_file)

    survival_results <- setNames(vector("list", length(cancer_types)), cancer_types)

    for (ct in cancer_types) {
      survival_results[[ct]] <- list()

      for (modality in modalities) {
        modality_key <- gsub(" ", "_", tolower(modality))

        survival_results[[ct]][[modality]] <- clinical_survival_analysis(
          cancer_type = ct,
          modality = modality,
          feature_selection = best_fs,
          clustering_method = best_cm,
          clinical_data = clinical_data[[ct]],
          sample_id_fun = sample_id_fun,
          modality_dirs = modality_dirs,
          output_file = file.path(
            dirs$survival,
            paste0(ct, "_survival_", modality_key, ".png")
          ),
          pvalue_file = survival_pvalue_file,
          verbose = TRUE
        )
      }
    }

    if (file.exists(survival_pvalue_file)) {
      plot_survival_summary(
        input_file = survival_pvalue_file,
        output_file = file.path(dirs$figures, "survival_summary.png"),
        width = 3.5,
        height = 2.5,
        dpi = 300
      )
    }
  }

  list(
    clinical = clinical_results,
    survival = survival_results,
    availability = availability_map
  )
}

save_analysis_settings <- function(config, setup, clustering) {
  settings <- data.frame(
    dataset = config$dataset,
    cancer_types = paste(config$cancer_types, collapse = ","),
    best_feature_selection = clustering$best_fs,
    best_clustering_method = clustering$best_cm,
    seed = config$seed,
    cores = config$cores,
    stringsAsFactors = FALSE
  )

  utils::write.csv(
    settings,
    file.path(setup$results_dir, "analysis_settings.csv"),
    row.names = FALSE
  )

  invisible(settings)
}

render_analysis_report <- function(config, setup, clustering, template = "report/analysis_report.Rmd") {
  if (!requireNamespace("rmarkdown", quietly = TRUE)) stop("Package `rmarkdown` must be installed.")
  if (!file.exists(template)) stop("Report template not found: ", template)

  results_dir <- normalizePath(setup$results_dir, mustWork = TRUE)
  output_file <- file.path(results_dir, "analysis_report.html")

  rmarkdown::render(
    input = template,
    output_file = basename(output_file),
    output_dir = results_dir,
    params = list(
      results_dir = results_dir,
      dataset = config$dataset,
      cancer_types = setup$cancer_types,
      modalities = setup$modalities,
      feature_methods = config$feature_methods,
      clustering_methods = config$clustering_methods,
      best_fs = clustering$best_fs,
      best_cm = clustering$best_cm,
      seed = config$seed,
      cores = config$cores,
      run_time = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    ),
    envir = new.env(parent = globalenv()),
    quiet = TRUE
  )

  output_file
}