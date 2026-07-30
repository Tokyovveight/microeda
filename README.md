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
ALDEx2-, ANCOM-BC2-, and DESeq2-based exploratory differential representation
is available when the corresponding optional package is installed. The methods
can be run separately or side-by-side over explicit or pairwise contrasts.
Formal method ranking is not implemented.

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

`microeda_da()` exposes ALDEx2-, ANCOM-BC2-, and DESeq2-backed exploratory
differential representation workflows. It uses method-native p-value
adjustment, returns standardized results with `as_da_results()`, and supports
controlled retention of method-native backend output in `da$raw_outputs`.

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
abundant taxa. Adjusted p-values remain method-native and scoped to each
method and contrast; microeda does not apply cross-method adjustment or create
consensus calls.

ANCOM-BC2 is optional and can be installed with
`BiocManager::install("ANCOMBC")`. Its current backend accepts one explicit
contrast at a time, does not filter or round counts, and reports the native
natural-log coefficient as `group2 - group1`. A pairwise request runs these
explicit primary analyses sequentially and does not use `res_pair`. Native
q-values are not adjusted again.

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

DESeq2 is optional and can be installed with
`BiocManager::install("DESeq2")`. It is provided as a sensitivity/comparison
view based on a negative-binomial count model, not as a compositional
correction. The default `poscounts` size-factor estimator supports sparse count
tables without adding a pseudocount. Standardized effects are native unshrunk
log2 fold changes for `group2 / group1`; native Cook's filtering, independent
filtering, BH adjustment, and missing values are preserved. Pairwise requests
fit each two-group comparison separately.

```r
if (requireNamespace("DESeq2", quietly = TRUE)) {
  da_deseq2 <- microeda_da(
    counts,
    metadata = metadata,
    taxonomy = taxonomy,
    group = "group",
    contrast = c("A", "B"),
    methods = "deseq2",
    tax_rank = "Genus",
    taxa_are_rows = FALSE
  )
  as_da_summary(da_deseq2)
  cat(microeda_da_report(da_deseq2, top_n = 5))
}
```

When all optional packages are installed, the same contrast plan can be
evaluated side-by-side. The standardized rows retain method-specific effect
semantics and follow method order, then contrast order. Agreement among methods
is not treated as a consensus or confirmation.

```r
if (requireNamespace("ALDEx2", quietly = TRUE) &&
    requireNamespace("ANCOMBC", quietly = TRUE) &&
    requireNamespace("DESeq2", quietly = TRUE)) {
  da_multi <- microeda_da(
    counts,
    metadata = metadata,
    taxonomy = taxonomy,
    group = "group",
    contrast = "pairwise",
    methods = c("aldex2", "ancombc2", "deseq2"),
    tax_rank = "Genus",
    taxa_are_rows = FALSE,
    mc.samples = 128,
    ancombc2_p_adj_method = "holm",
    raw_storage = "compact",
    progress = TRUE
  )
  cat(microeda_da_report(da_multi))
  results <- as_da_results(da_multi)
  summary <- as_da_summary(da_multi)
  write_da_results(da_multi, "microeda_da_results.csv")
}
```

Native-output retention is explicit:

```r
da_compact <- microeda_da(
  ps,
  group = "Location",
  contrast = c("Arsk", "Laishevo"),
  methods = c("aldex2", "ancombc2", "deseq2"),
  raw_storage = "compact",
  progress = TRUE
)
```

- `raw_storage = "full"` is the backward-compatible default and keeps complete
  fit objects for maximum auditability, at the largest returned-object size.
- `raw_storage = "compact"` keeps native result tables, mapping metadata,
  warnings/messages, and diagnostics without large fit/intermediate objects.
- `raw_storage = "none"` keeps standardized results, caveats, parameters, and
  timings only; `as_da_raw_output()` then explains how to rerun with retention.

Compact or no raw storage reduces the returned object size, but it does not
reduce peak memory required while a backend is running. Per-contrast timings
are always available in `da$params$timings`, and total elapsed time is in
`da$params$total_elapsed_seconds`. `progress = TRUE` emits method/contrast
start and finish messages without mixing them with captured backend warnings.
The `method_results` and `raw_outputs` paths reference the same retained raw
lists when the object is built; recursive `object.size()` reporting can count
both paths and overstate independently allocated memory.

With `raw_storage = "full"` or `"compact"`, raw outputs remain separate.
Explicit runs use
`da$raw_outputs$aldex2`, `da$raw_outputs$ancombc2`, or
`da$raw_outputs$deseq2`. For pairwise runs, ALDEx2 retains
`da$raw_outputs$aldex2$contrasts$A_vs_B`, while ANCOM-BC2 and DESeq2 use
`da$raw_outputs$ancombc2$A_vs_B` and `da$raw_outputs$deseq2$A_vs_B`.
In full mode, each ANCOM-BC2 or DESeq2 contrast contains its complete native
explicit result. Compact outputs are marked with
`raw_storage = "compact"` and retain native tables and diagnostics without the
largest fit objects.

Use the comparison helpers to inspect method-specific standardized values next
to one another without creating a consensus or ranking. Native effects retain
their different scales, so their magnitudes must not be compared directly.

```r
comparison <- as_da_comparison(
  da_multi,
  contrast = "A_vs_B"
)

cat(microeda_da_comparison_report(
  da_multi,
  contrast = "A_vs_B",
  alpha = 0.05,
  max_features = 15
))

write_da_comparison(
  da_multi,
  "da_comparison.csv",
  contrast = "A_vs_B"
)

raw_deseq2 <- as_da_raw_output(
  da_multi,
  method = "deseq2",
  contrast = "A_vs_B"
)
```

`as_da_results()` retains the long standardized table,
`as_da_summary()` provides method-by-contrast counts,
`microeda_da_report()` describes the full DA object, and
`write_da_results()` exports the long table. `as_da_comparison()` and
`write_da_comparison()` use only standardized rows; native objects remain
available through `as_da_raw_output()`.

For matched samples, identify pairs explicitly. microeda validates each
contrast independently and never infers pairs from sample order. Native ALDEx2
effect values are oriented as `group2 - group1`. Paired execution currently
requires ALDEx2 to be the only requested method because repeated-measures
ANCOM-BC2 and DESeq2 designs are not implemented.

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
