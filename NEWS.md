# microeda 0.0.0.9000

* Added a minimal QC summary skeleton for per-sample, per-feature, taxonomy-rank,
  and metadata-completeness diagnostics.
* Added compact library-size and sparsity diagnostics to `microeda_qc()`.
* Added compact prevalence and filtering diagnostics to `microeda_qc()`.
* Added feature read dominance diagnostics to `microeda_qc()` summaries and
  reports.
* Added human-readable QC observations to `microeda_qc()`.
* Added `as_qc_summary()` for compact QC display tables.
* Added `as_qc_issues()` for tabular QC flags and observations.
* Added `microeda_qc_report()` as a minimal text QC report skeleton.
* Added `microeda_qc_write_report()` for writing text QC reports to files.
* Added `include_flags` and `include_observations` content controls to
  `microeda_qc_report()`.
* `microeda_qc_report()` now includes compact QC flag and observation details
  when requested.
* Added `microeda_qc_plot()` as a minimal base R library-size QC plotting
  helper.
* Added `type = "sparsity"` support to `microeda_qc_plot()`.
* Added `type = "feature_abundance"` support to `microeda_qc_plot()`.
* Added `type = "prevalence"` support to `microeda_qc_plot()`.
* Added `microeda_alpha_plot()` for base R alpha diversity metric barplots.
* Added grouped boxplot support to `microeda_alpha_plot()`.
* Added `microeda_alpha_report()` for compact plain-text alpha diversity
  summaries and group-test reports.
* Improved `microeda_alpha_pairwise_report()` readability with compact
  rstatix-style columns and significance labels.
* Added Wilcoxon statistics to `as_alpha_pairwise()` and
  `microeda_alpha_pairwise_report()`.
* Improved QC and alpha text-report labels and fixed-width table formatting.
* Improved beta comparison text-report formatting and print hints.
* Clarified generic recommendation presentation as broad screening notes.
* Reworked README examples around a coherent check, QC, alpha, and beta
  workflow.
* Improved QC and alpha print hints for human-readable report helpers.
* Added alpha report index filtering and cleaner missing pairwise significance
  labels.
* Added optional `vegan`-backed beta group testing that pairs PERMANOVA with
  dispersion diagnostics.
* Added `microeda_beta()` for initial Bray-Curtis beta diversity distances.
* Added `as_beta_dist()` to extract beta diversity distances from
  `microeda_beta` objects.
* Added `as_beta_matrix()` to extract beta diversity distances as a square
  matrix.
* Added a compact print method for `microeda_beta` objects.
* Added `as_beta_samples()` to extract sample IDs and optional groups from
  `microeda_beta` objects.
* Added `microeda_beta_plot()` for base R heatmaps of beta diversity distances.
* Added `microeda_beta_ordination()` and `as_beta_coordinates()` for initial
  PCoA ordination of beta diversity distances.
* Added `microeda_alpha_pairwise_report()` for compact grouped alpha diversity
  pairwise comparison reports.
* Added `method = "jaccard"` to `microeda_beta()` for binary presence/absence
  Jaccard distances.
* Added `method = "hellinger"` to `microeda_beta()` for
  Hellinger-transformed Euclidean distances.
* Added `microeda_beta_compare()` and `as_beta_compare_summary()` for compact
  multi-method beta diversity comparison.
* Added `as_beta_compare_distances()` for long-form pairwise distances from
  beta method comparisons.
* Added `as_beta_compare_group_summary()` for within- and between-group beta
  distance summaries.
* Added `microeda_beta_compare_report()` for compact text reports of beta
  method comparisons.
* Added `microeda_beta_compare_ordination()` and
  `as_beta_compare_coordinates()` for side-by-side PCoA coordinates across beta
  methods.
* Added `as_beta_compare_distance_correlations()` for descriptive correlations
  between beta distance methods.
* `microeda_beta_compare_report()` now includes descriptive distance-method
  correlations.
