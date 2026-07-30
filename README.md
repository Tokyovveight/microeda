# microeda

`microeda` is an R package for evidence-guided exploratory analysis of
microbiome count data.

The core idea is deliberately conservative:

> inspect the data first, describe statistical risks, then use method-specific
> summaries and caveats instead of pretending there is one universal best
> microbiome workflow.

## Workflow Scope

- Accept a `phyloseq` object or a plain count matrix plus metadata.
- Inspect input scale, library sizes, sparsity, prevalence, metadata
  completeness, group sizes, taxonomy completeness, and feature read dominance.
- Report alpha diversity summaries and exploratory group comparisons.
- Compare Bray-Curtis, Jaccard, and Hellinger beta diversity views.
- Keep broad screening notes with `rule_id`, `topic`, `severity`,
  `recommendation`, `caveat`, and `evidence`.

PERMANOVA-style beta group testing is available when the optional `vegan`
package is installed and is reported together with dispersion diagnostics.
ALDEx2- and ANCOM-BC2-based exploratory differential representation is
available when the corresponding optional package is installed. ANCOM-BC2
currently supports one explicit contrast. DESeq2, cross-method dispatch, and
formal method ranking are not implemented yet.

## First Workflow

Prepare a count table with samples as rows and features as columns, plus sample
metadata and optional taxonomy.

```r
counts <- matrix(
  c(
    10, 0, 0, 5,
    20, 0, 1, 0,
    0,  4, 0, 0,
    2,  3, 0, 1
  ),
  nrow = 4,
  byrow = TRUE
)
rownames(counts) <- paste0("S", 1:4)
colnames(counts) <- paste0("ASV", 1:4)

metadata <- data.frame(
  group = c("A", "A", "B", "B"),
  batch = c("x", "y", "x", "y"),
  row.names = rownames(counts)
)

taxonomy <- data.frame(
  Phylum = c("Firmicutes", "Firmicutes", "Bacteroidota", "Actinobacteriota"),
  Genus = c("Lactobacillus", "Streptococcus", "Bacteroides", "Bifidobacterium"),
  row.names = colnames(counts)
)
```

Start with compact input screening:

```r
check <- microeda_check(
  counts,
  metadata = metadata,
  taxonomy = taxonomy,
  group = "group",
  taxa_are_rows = FALSE
)

check
```

Then build the QC, alpha diversity, and beta diversity report layers:

```r
qc <- microeda_qc(
  counts,
  metadata = metadata,
  taxonomy = taxonomy,
  group = "group",
  taxa_are_rows = FALSE
)
cat(microeda_qc_report(qc))

alpha <- microeda_alpha(
  counts,
  metadata = metadata,
  group = "group",
  taxa_are_rows = FALSE
)
alpha_cmp <- microeda_alpha_compare(alpha, group = "group")
cat(microeda_alpha_report(alpha, alpha_compare = alpha_cmp))
cat(microeda_alpha_pairwise_report(alpha_cmp))

beta_cmp <- microeda_beta_compare(
  counts,
  metadata = metadata,
  group = "group",
  taxa_are_rows = FALSE
)
cat(microeda_beta_compare_report(beta_cmp))

beta_bray <- microeda_beta(
  counts,
  metadata = metadata,
  group = "group",
  taxa_are_rows = FALSE,
  method = "bray"
)
if (requireNamespace("vegan", quietly = TRUE)) {
  beta_test <- microeda_beta_test(beta_bray, permutations = 99, seed = 1)
  cat(microeda_beta_test_report(beta_test))
}

if (requireNamespace("ALDEx2", quietly = TRUE)) {
  da <- microeda_da(
    counts,
    metadata = metadata,
    taxonomy = taxonomy,
    group = "group",
    contrast = c("A", "B"),
    methods = "aldex2",
    tax_rank = "Genus",
    taxa_are_rows = FALSE,
    mc.samples = 16
  )
  as_da_results(da)
}
```

For `phyloseq`, pass the object directly:

```r
check <- microeda_check(ps, group = "Treatment")
qc <- microeda_qc(ps, group = "Treatment")
alpha <- microeda_alpha(ps, group = "Treatment")
```

## Human-Readable Reports

Use these helpers when you want compact console or text output:

```r
cat(microeda_qc_report(qc))
cat(microeda_alpha_report(alpha, alpha_compare = alpha_cmp))
cat(microeda_alpha_pairwise_report(alpha_cmp))
cat(microeda_beta_compare_report(beta_cmp))
if (exists("beta_test")) {
  cat(microeda_beta_test_report(beta_test))
}
if (exists("beta_cmp_test")) {
  cat(microeda_beta_compare_test_report(beta_cmp_test))
}
microeda_qc_write_report(qc, tempfile(fileext = ".txt"))
```

`microeda_qc_report()` summarizes samples, features, reads, sparsity, QC flags,
observations, and feature dominance. `microeda_alpha_report()` formats alpha
summaries and omnibus group tests by diversity index.
`microeda_alpha_pairwise_report()` formats pairwise Wilcoxon comparisons with
statistics and adjusted p-value labels. `microeda_beta_compare_report()`
summarizes distance methods, method correlations, grouped distance summaries,
and caveats. `microeda_beta_test_report()` reports exploratory PERMANOVA
results together with dispersion diagnostics and caveats.
`microeda_beta_compare_test_report()` shows the same paired diagnostics
side-by-side across beta distance methods without ranking methods.

## Machine-Readable Extractors

The `as_*()` helpers return data frames, matrices, or `dist` objects for
downstream analysis and custom reporting.

```r
as_qc_summary(qc)
as_qc_issues(qc)

as_alpha_table(alpha)
as_alpha_summary(alpha)
as_alpha_tests(alpha_cmp)
as_alpha_pairwise(alpha_cmp)

beta_bray <- microeda_beta(
  counts,
  metadata = metadata,
  group = "group",
  taxa_are_rows = FALSE,
  method = "bray"
)
as_beta_dist(beta_bray)
as_beta_matrix(beta_bray)
as_beta_samples(beta_bray)

as_beta_compare_summary(beta_cmp)
as_beta_compare_distances(beta_cmp)
as_beta_compare_group_summary(beta_cmp)
as_beta_compare_distance_correlations(beta_cmp)

if (exists("beta_test")) {
  as_beta_test_summary(beta_test)
}
if (exists("beta_cmp_test")) {
  as_beta_compare_test_summary(beta_cmp_test)
}
if (exists("da")) {
  as_da_results(da)
}
```

## Beta Group Testing

`microeda_beta_test()` accepts a grouped `microeda_beta` object and uses the
stored distance object for paired PERMANOVA and betadisper-style dispersion
diagnostics. This workflow requires the optional `vegan` package.

```r
if (requireNamespace("vegan", quietly = TRUE)) {
  beta_test <- microeda_beta_test(beta_bray, permutations = 999, seed = 1)
  as_beta_test_summary(beta_test)
  cat(microeda_beta_test_report(beta_test))
}
```

These tests are exploratory. PERMANOVA can be confounded by group dispersion
differences, so the dispersion diagnostics and caveats should be inspected
alongside the PERMANOVA table.

For side-by-side exploratory diagnostics across beta distance methods, pass an
existing grouped `microeda_beta_compare` object to
`microeda_beta_compare_test()`. This reuses the stored distance objects and does
not rank methods.

```r
if (requireNamespace("vegan", quietly = TRUE)) {
  beta_cmp_test <- microeda_beta_compare_test(
    beta_cmp,
    permutations = 999,
    seed = 1
  )
  as_beta_compare_test_summary(beta_cmp_test)
  cat(microeda_beta_compare_test_report(beta_cmp_test))
}
```

The alpha table includes classic indices (`observed`, `chao1`, `shannon`,
`simpson`, `inverse_simpson`) and Hill/effective-diversity equivalents
(`hill_q0`, `hill_q1`, `hill_q2`). In practice, `hill_q1 = exp(Shannon)` and
`hill_q2 = inverse Simpson`, which makes the values easier to interpret as
effective numbers of taxa.

## Differential Representation

`microeda_da()` exposes ALDEx2- and ANCOM-BC2-backed exploratory differential
representation workflows. It uses method-native p-value adjustment, returns
standardized results with `as_da_results()`, and preserves complete native
backend output in `da$raw_outputs`.

```r
if (requireNamespace("ALDEx2", quietly = TRUE)) {
  da <- microeda_da(
    counts,
    metadata = metadata,
    taxonomy = taxonomy,
    group = "group",
    contrast = "pairwise",
    methods = "aldex2",
    tax_rank = "Genus",
    taxa_are_rows = FALSE
  )
  as_da_results(da)
  as_da_summary(da)
  cat(microeda_da_report(da, top_n = 5))
  write_da_results(da, "microeda_da_results.csv")
}
```

This workflow does not rank methods or claim ground-truth differentially
abundant taxa. Public multi-method dispatch and DESeq2 are planned for later
slices.

ANCOM-BC2 is optional and can be installed with
`BiocManager::install("ANCOMBC")`. Its current backend accepts one explicit
contrast, does not filter or round counts, and reports the native natural-log
coefficient as `group2 - group1`. Native q-values are not adjusted again.

```r
if (requireNamespace("ANCOMBC", quietly = TRUE)) {
  da_ancombc2 <- microeda_da(
    counts,
    metadata = metadata,
    taxonomy = taxonomy,
    group = "group",
    contrast = c("A", "B"),
    methods = "ancombc2",
    tax_rank = "Genus",
    taxa_are_rows = FALSE,
    ancombc2_p_adj_method = "holm"
  )
  as_da_summary(da_ancombc2)
  cat(microeda_da_report(da_ancombc2, top_n = 5))
}
```

For matched samples, identify pairs explicitly. microeda validates each
contrast independently and never infers pairs from sample order. Native ALDEx2
effect values are oriented as `group2 - group1`.

```r
metadata_paired <- metadata
metadata_paired$pair_id <- c("P1", "P2", "P1", "P2")

if (requireNamespace("ALDEx2", quietly = TRUE)) {
  da_paired <- microeda_da(
    counts,
    metadata = metadata_paired,
    group = "group",
    contrast = c("A", "B"),
    taxa_are_rows = FALSE,
    paired.test = TRUE,
    pair_id = "pair_id"
  )
}
```

## Plots And Ordinations

The current plotting helpers use base R.

```r
microeda_qc_plot(qc, type = "library_size")
microeda_qc_plot(qc, type = "sparsity")
microeda_qc_plot(qc, type = "feature_abundance")
microeda_qc_plot(qc, type = "prevalence")

microeda_alpha_plot(alpha)
microeda_alpha_plot(alpha, type = "boxplot")
microeda_alpha_plot(alpha, metric = "shannon", type = "boxplot", group = "group")

microeda_beta_plot(beta_bray)

ord <- microeda_beta_ordination(beta_bray)
as_beta_coordinates(ord)

ord_cmp <- microeda_beta_compare_ordination(beta_cmp)
as_beta_compare_coordinates(ord_cmp)
```

PCoA coordinates are method-specific; axes from different distance methods are
intended for side-by-side inspection, not direct axis-by-axis equivalence.

## Broad Screening Notes

`as_recommendations(check)` extracts broad screening notes from
`microeda_check()`. These notes are caveats for initial review, not contextual
workflow recommendations, formal method ranking, or a substitute for the QC,
alpha, and beta reports.

```r
as_recommendations(check)
microeda_rules()
```

The package keeps the current evidence map in
`inst/extdata/evidence_rules.yml`. The R functions use a built-in version of
the same rules so the package does not require a YAML parser at runtime.
