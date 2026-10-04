library(targets)

invisible(lapply(
  list.files("R", pattern = "\\.[Rr]$", full.names = TRUE),
  source))

source("examples/config_tcga.r")

tar_option_set(
  packages = c("dplyr", "tidyr", "ggplot2", "cola", "survival", "mclust"),
  seed = config$seed)

always <- tar_cue(mode = "always")

list(
  tar_target(config_target, config),

  tar_target(
    setup,
    prepare_analysis(config_target),
    cue = always),

  tar_target(
    cancer_type,
    config_target$cancer_types,
    iteration = "vector"),

  tar_target(
    cancer_summary,
    run_clustering_for_cancer(
      config = config_target,
      setup = setup,
      ct = cancer_type),
    pattern = map(cancer_type),
    iteration = "list"),

  tar_target(
    clustering,
    combine_clustering_results(
      config = config_target,
      setup = setup,
      cancer_summaries = cancer_summary),
    cue = always),

  tar_target(
    downstream,
    run_downstream_analysis(
      config = config_target,
      setup = setup,
      clustering = clustering),
    cue = always),

  tar_target(
    clinical,
    run_clinical_workflow(
      config = config_target,
      setup = setup,
      clustering = clustering),
    cue = always),

  tar_target(
    settings,
    save_analysis_settings(
      config = config_target,
      setup = setup,
      clustering = clustering),
    cue = always),

  tar_target(
    analysis_report,
    {
      downstream
      clinical
      settings
      render_analysis_report(
        config = config_target,
        setup = setup,
        clustering = clustering)
    },
    cue = always)
)