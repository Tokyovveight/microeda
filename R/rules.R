#' Return the built-in evidence rules
#'
#' Return the built-in evidence rule table used by [microeda_check()].
#'
#' The same rule map is mirrored as `inst/extdata/evidence_rules.yml` so it can
#' be inspected or edited outside R.
#'
#' @return A data frame containing rule metadata.
#' @export
microeda_rules <- function() {
  data.frame(
    rule_id = c(
      "compositional_caveat",
      "relative_input",
      "non_integer_input",
      "high_sparsity",
      "library_size_imbalance",
      "high_taxa_sample_ratio",
      "low_group_n",
      "metadata_missingness",
      "taxonomy_missingness",
      "permanova_dispersion",
      "da_method_concordance",
      "covariate_modeling"
    ),
    topic = c(
      "composition",
      "input",
      "input",
      "sparsity",
      "sequencing_depth",
      "dimensionality",
      "group_design",
      "metadata",
      "taxonomy",
      "beta_diversity",
      "differential_abundance",
      "differential_abundance"
    ),
    severity = c(
      "info",
      "warning",
      "warning",
      "warning",
      "warning",
      "info",
      "warning",
      "warning",
      "warning",
      "info",
      "info",
      "info"
    ),
    recommendation = c(
      "Treat abundance differences as compositional signals unless absolute abundance information is available.",
      "Do not feed relative or transformed abundance tables into methods that expect raw counts; keep labels explicit.",
      "Verify whether the table is transformed or normalized before running count-based DA methods such as ANCOM-BC2 or ALDEx2.",
      "Compare Bray-Curtis or Jaccard views with log-ratio-aware approaches such as rCLR or robust Aitchison ordination.",
      "Inspect rarefaction curves and run sensitivity checks before strong alpha or beta diversity interpretation.",
      "Prefer exploratory summaries and feature filtering before high-dimensional testing.",
      "Treat DA and group comparisons as exploratory when any group has a very small sample size.",
      "Resolve missing values in design-critical metadata columns before model-based analysis.",
      "Report missing taxonomy by rank and be cautious with rank-level claims when annotation is incomplete.",
      "If PERMANOVA is used, pair it with dispersion diagnostics such as betadisper.",
      "Compare ANCOM-BC2 and ALDEx2 results and emphasize taxa supported by multiple methods.",
      "When covariates such as batch, host, site, or time are present, prefer adjusted models over simple two-group tests."
    ),
    caveat = c(
      "Relative abundance does not directly measure absolute biomass.",
      "Relative tables can be valid for visualization, but not for every downstream statistical method.",
      "Non-integer values may be appropriate for plots or ordination, but count models rely on count-scale assumptions.",
      "CLR with pseudocounts can be sensitive when zeros dominate; use sensitivity analysis.",
      "Rarefaction can help diagnose uneven depth, but rarefied counts should not be the default input for all analyses.",
      "Filtering changes the estimand; report thresholds and keep raw data available.",
      "Low n reduces power and stability; avoid strong biological conclusions from single-method hits.",
      "Missing metadata can silently change sample inclusion in models.",
      "Poor annotation can make genus-level or species-level claims unstable.",
      "Distance-based tests can mix location and dispersion effects.",
      "DA methods can disagree substantially across realistic microbiome datasets.",
      "Adjusted models require careful design checks and enough samples per model term."
    ),
    evidence = c(
      "Gloor 2017; Quinn 2019; Lutz 2022",
      "Gloor 2017; Quinn 2019; ANCOM-BC2 2024; ALDEx2 2014",
      "Quinn 2019; ANCOM-BC2 2024; ALDEx2 2014",
      "Martino 2019; Quinn 2019; Zhang 2024; Chan 2024",
      "Willis 2019; McMurdie and Holmes 2014; Schloss 2023; Lin and Peddada 2020",
      "Lutz 2022; Zhou 2023",
      "Nearing 2022; Weiss 2017",
      "Mirzayi 2021; Zhou 2023",
      "Zhou 2023; McMurdie and Holmes 2013",
      "Anderson 2006; Warton 2012",
      "Nearing 2022; Calgaro 2023",
      "ANCOM-BC2 2024; Weiss 2017"
    ),
    stringsAsFactors = FALSE
  )
}

build_recommendations <- function(diagnostics, group = NULL) {
  rules <- microeda_rules()
  selected <- "compositional_caveat"

  if (isTRUE(diagnostics$count_type$looks_relative)) {
    selected <- c(selected, "relative_input")
  }

  if (!isTRUE(diagnostics$count_type$integerish)) {
    selected <- c(selected, "non_integer_input")
  }

  if (diagnostics$sparsity$zero_fraction >= 0.7 ||
      diagnostics$sparsity$low_prevalence_fraction >= 0.5) {
    selected <- c(selected, "high_sparsity")
  }

  if (!is.na(diagnostics$library_size$imbalance_ratio) &&
      diagnostics$library_size$imbalance_ratio >= 10) {
    selected <- c(selected, "library_size_imbalance")
  }

  if (!is.na(diagnostics$taxa_sample_ratio) &&
      diagnostics$taxa_sample_ratio >= 10) {
    selected <- c(selected, "high_taxa_sample_ratio")
  }

  if (!is.null(diagnostics$metadata$group) &&
      diagnostics$metadata$group$min_n < 5) {
    selected <- c(selected, "low_group_n")
  }

  if (isTRUE(diagnostics$metadata$provided) &&
      any(diagnostics$metadata$missing_fraction_by_column > 0)) {
    selected <- c(selected, "metadata_missingness")
  }

  if (isTRUE(diagnostics$taxonomy$provided) &&
      any(diagnostics$taxonomy$missing_fraction_by_rank >= 0.3)) {
    selected <- c(selected, "taxonomy_missingness")
  }

  if (!is.null(group)) {
    selected <- c(selected, "permanova_dispersion", "da_method_concordance")
  }

  if (isTRUE(diagnostics$metadata$provided) &&
      diagnostics$metadata$n_columns > 1 &&
      !is.null(group)) {
    selected <- c(selected, "covariate_modeling")
  }

  selected <- unique(selected)
  rules[match(selected, rules$rule_id), , drop = FALSE]
}
