# BOAR: Benchmarking Of Approaches for Regulatory-Network Clustering

BOAR is a reproducible R pipeline for clustering cancer patients using gene expression and patient-specific gene regulatory network representations.

The pipeline was developed as part of a benchmarking study comparing different approaches for patient stratification using gene expression, gene indegrees, and transcription factor–gene edge weights.

The workflow is implemented with the `targets` package and uses the `cola` package for consensus clustering.

## What the pipeline does

The main workflow includes:

* reading expression and network-derived data;
* running consensus clustering with different feature-selection and clustering methods;
* comparing clustering performance using 1-PAC, mean silhouette and concordance;
* selecting the best-performing feature-selection and clustering combination;
* comparing cluster assignments between data representations using adjusted Rand index;
* evaluating cluster balance;
* testing associations between clusters and clinical variables;
* performing survival analysis;
* generating summary tables and figures.

## Installation

The project uses `renv` to record package versions and reproduce the software environment.

After cloning the repository, install `renv` if necessary:

install.packages("renv")

Then restore the project environment:

renv::restore()

This installs the package versions recorded in `renv.lock`.

The pipeline was developed and tested using R 4.5.1.

The main R packages used by the workflow are:

* targets
* cola
* dplyr
* tidyr
* ggplot2
* purrr
* mclust
* survival

## Running the pipeline

The repository contains a TCGA-specific configuration file under `examples/`.

The TCGA example is currently configured to analyze ACC and CHOL using two patient representations: gene expression and gene indegrees.

All files required to run this example are provided under:

examples/data/tcga/

The example data are organized as:

examples/data/tcga/
├── expression/
│   ├── tcga_ACC_expression.tsv
│   └── tcga_CHOL_expression.tsv
├── indegrees/
│   ├── indegree_ACC.RData
│   └── indegree_CHOL.RData
└── clinical/
    ├── TCGA_ACC_clinical.tsv
    ├── TCGA_CHOL_clinical.tsv
    └── TCGA_clinical_exclude_cols.txt

To run the example pipeline from the project root:

targets::tar_make()

The results are written to the output directory defined in `config_tcga.R`.

## Input data

The pipeline can work with three patient representations:

* Expression: gene expression values for each sample.
* Indegrees: gene-level summaries of regulatory input from patient-specific gene regulatory networks.
* Edge weights: individual transcription factor–gene regulatory interactions from patient-specific networks.

Matrices are expected to contain features in rows and samples in columns.

Not every dataset needs to contain all three representations. The modalities used in an analysis are defined in the corresponding configuration file.

Different modalities do not need to contain exactly the same samples. Sample matching is performed separately where required by downstream analyses.

## Configuration

Dataset-specific configuration files define the cancer types or cohorts, modalities, feature-selection methods, clustering methods, input and output paths, file naming conventions, clustering settings, clinical data settings, optional downstream analyses, random seed, and number of cores.

The example configuration files can be adapted to new datasets as long as the input files follow the expected formats.

## Clustering

Consensus clustering is performed using combinations of feature-selection and clustering methods.

The current workflow supports the following feature-selection methods:

SD
CV
MAD
ATC

and clustering methods:

kmeans
skmeans
hclust
mclust
pam

Clustering performance is summarized using 1-PAC, mean silhouette and concordance.

The best number of clusters is determined using `cola::suggest_best_k()`.

The pipeline identifies the best overall feature-selection and clustering combination across the analyzed cancer types, and uses this combination for the downstream analyses.

## Cluster comparison

Cluster assignments obtained from different data representations can be compared using adjusted Rand index.

ARI calculations use the samples shared between the two compared modalities.

Cluster balance is evaluated using normalized entropy, where values closer to 1 indicate more balanced cluster sizes.

## Clinical and survival analyses

Cluster assignments can be compared with clinical variables using categorical or numerical association tests followed by Benjamini-Hochberg correction.

Categorical variables are tested using chi-squared tests when the required sample-count criteria are satisfied.

Numerical variables are tested using either ANOVA or Kruskal-Wallis tests depending on residual normality.

Survival differences between clusters are evaluated using a log-rank test and can be visualized using Kaplan-Meier curves.

Clinical data are deliberately kept separate from the core clustering functions so that the same workflow can be applied to different cohorts.

Each clinical file must contain at least the following three columns:

sample_id
survival_time
survival_event

where:

* `sample_id` contains the patient or sample identifier used to match the clinical data to the clustering results.
* `survival_time` contains the survival or follow-up time as a numeric value.
* `survival_event` indicates whether the event occurred and must be coded as `0` for censored/alive and `1` for event/deceased.

Additional clinical variables can be included in the same file and are used for downstream association testing. For example:

sample_id    survival_time    survival_event    age    stage       sex
TCGA-XX-0001 1240             0                 63     Stage II    female
TCGA-XX-0002 561              1                 71     Stage III   male

The standardized identifier and survival columns are not included in the clinical association tests.

### Preparing TCGA clinical data

The expected clinical format is based on the TCGA clinical data used in the analyses accompanying this project. TCGA clinical data can be downloaded in R using `TCGAbiolinks`:

library(TCGAbiolinks)

clinical <- GDCquery_clinic(
  project = "TCGA-CHOL",
  type = "clinical"
)

The resulting table contains patient identifiers and survival-related variables including `submitter_id`, `vital_status`, `days_to_death`, and `days_to_last_follow_up`.

These can be converted to the format expected by the pipeline as follows:

clinical$sample_id <- clinical$submitter_id

clinical$survival_time <- ifelse(
  !is.na(clinical$days_to_death),
  clinical$days_to_death,
  clinical$days_to_last_follow_up
)

clinical$survival_event <- ifelse(
  is.na(clinical$vital_status),
  NA_integer_,
  ifelse(clinical$vital_status == "Dead", 1L, 0L)
)

The remaining clinical variables can be kept in the same table so that they are available for association testing:

clinical <- clinical[
  ,
  c(
    "sample_id",
    "survival_time",
    "survival_event",
    setdiff(
      colnames(clinical),
      c("sample_id", "survival_time", "survival_event")
    )
  )
]

write.table(
  clinical,
  file = "TCGA_CHOL_clinical.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

Clinical data from non-TCGA cohorts can be used in the same workflow as long as they are converted to the same standardized format before running the analysis.

## Output

The results directory is defined in the configuration file.

The pipeline creates separate directories for cluster assignments, `cola` reports, summary tables, figures, clinical analyses, and survival analyses.

Typical outputs include:

* cluster assignment files;
* `cola` clustering reports;
* per-cancer clustering summaries;
* a combined clustering summary;
* ranked feature-selection and clustering combinations;
* best-k figures;
* cluster-balance results;
* adjusted Rand index results;
* clinical association tables and heatmaps;
* Kaplan-Meier plots;
* survival p-value summaries;
* analysis settings.

## Reproducibility

The example configurations use:

seed = 12345

The configured seed is used by the clustering workflow.

The seed, number of cores, and selected best-performing feature-selection and clustering combination are saved with the analysis results.

Package versions are recorded in `renv.lock`. The software environment can be restored with:

renv::restore()

Intermediate results are managed by `targets` in the `_targets/` directory. This directory is a local cache and does not need to be included in the public repository.

## Troubleshooting

### Results directory was deleted or changed

The `targets` cache is stored separately from the generated analysis results. Therefore, deleting the `results/` directory does not automatically invalidate previously completed clustering targets.

If `results/` has been deleted, renamed, or the output path has been changed, downstream analyses may report missing summary or cluster-assignment files even though `targets` considers the clustering step up to date.

In this case, invalidate the cancer-level clustering targets and rerun the pipeline:

targets::tar_invalidate(cancer_summary)
targets::tar_make()

## Citation

BOAR: Benchmarking Of Approaches for Regulatory-Network Clustering

If you use BOAR, please cite the associated manuscript once available.

Manuscript in preparation:
Panos Sarantidis, Tatiana Belova, and Marieke L. Kuijjer.
"A Pan-Cancer Benchmark of Clustering Strategies for Patient-Specific Gene Regulatory Networks."

## License

This project is distributed under the terms of the license provided in the `LICENSE` file.