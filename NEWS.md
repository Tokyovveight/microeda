# microeda 0.0.0.9000

* Added a minimal QC summary skeleton for per-sample, per-feature, taxonomy-rank,
  and metadata-completeness diagnostics.
* Added compact library-size and sparsity diagnostics to `microeda_qc()`.
* Added compact prevalence and filtering diagnostics to `microeda_qc()`.
* Added human-readable QC observations to `microeda_qc()`.
* Added `as_qc_summary()` for compact QC display tables.
* Added `as_qc_issues()` for tabular QC flags and observations.
* Added `microeda_qc_report()` as a minimal text QC report skeleton.
* Added `microeda_qc_write_report()` for writing text QC reports to files.
* Added `include_flags` and `include_observations` content controls to
  `microeda_qc_report()`.
* Added `microeda_qc_plot()` as a minimal base R library-size QC plotting
  helper.
* Added `type = "sparsity"` support to `microeda_qc_plot()`.
* Added `type = "feature_abundance"` support to `microeda_qc_plot()`.
* Added `type = "prevalence"` support to `microeda_qc_plot()`.
* Added `microeda_alpha_plot()` for base R alpha diversity metric barplots.
* Added grouped boxplot support to `microeda_alpha_plot()`.
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
