# microeda 0.0.0.9000

* Added a minimal QC summary skeleton for per-sample, per-feature, taxonomy-rank,
  and metadata-completeness diagnostics.
* Added compact library-size and sparsity diagnostics to `microeda_qc()`.
* Added compact prevalence and filtering diagnostics to `microeda_qc()`.
* Added human-readable QC observations to `microeda_qc()`.
* Added `as_qc_summary()` for compact QC display tables.
* Added `microeda_qc_report()` as a minimal text QC report skeleton.
* Added `include_flags` and `include_observations` content controls to
  `microeda_qc_report()`.
* Added `microeda_qc_plot()` as a minimal base R library-size QC plotting
  helper.
* Added `type = "sparsity"` support to `microeda_qc_plot()`.
* Added `type = "feature_abundance"` support to `microeda_qc_plot()`.
* Added `type = "prevalence"` support to `microeda_qc_plot()`.
