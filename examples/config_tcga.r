config <- list(
  dataset = "TCGA",
  cancer_types = c("ACC", "CHOL"),
  feature_methods = c("SD"),
  clustering_methods = c("kmeans"),
  modalities = c("Expression", "Indegrees"),
  seed = 12345,
  cores = 1,

  paths = list(
    expression = "examples/data/expression",
    indegrees = "examples/data/indegrees",
    edge_weights = NULL,
    clinical = "examples/data/clinical",
    exclude_cols = "examples/data/clinical/exclude_cols.txt",
    results = "results"
  ),

  expression = list(
    file_pattern = "tcga_%s_expression.tsv",
    gene_col = "gene_name",
    drop_cols = "gene_id",
    remove_normals = TRUE,
    normal_pattern = "-11A-"
  ),

  indegrees = list(
    file_pattern = "indegree_%s.RData",
    object_name = "ind"
  ),

  edge_weights = list(
    file_pattern = NULL,
    object_name = NULL
  ),

  clustering = list(
    make_report = TRUE,
    save_clusters = TRUE
  ),

  clinical = list(
    file = "TCGA_%s_clinical.tsv",
    sample_id_fun = function(x) substr(x, 1, 12)
  ),

  analysis = list(
    clustering_summary = TRUE,
    best_k = TRUE,
    cluster_balance = TRUE,
    ari = TRUE,
    clinical = TRUE,
    survival = TRUE
  )
)