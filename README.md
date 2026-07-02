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

Compositional/log-ratio methods, PERMANOVA, dispersion tests, differential
abundance methods, and formal method ranking are not implemented yet.

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
microeda_qc_write_report(qc, tempfile(fileext = ".txt"))
```

`microeda_qc_report()` summarizes samples, features, reads, sparsity, QC flags,
observations, and feature dominance. `microeda_alpha_report()` formats alpha
summaries and omnibus group tests by diversity index.
`microeda_alpha_pairwise_report()` formats pairwise Wilcoxon comparisons with
statistics and adjusted p-value labels. `microeda_beta_compare_report()`
summarizes distance methods, method correlations, grouped distance summaries,
and caveats.

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
```

The alpha table includes classic indices (`observed`, `chao1`, `shannon`,
`simpson`, `inverse_simpson`) and Hill/effective-diversity equivalents
(`hill_q0`, `hill_q1`, `hill_q2`). In practice, `hill_q1 = exp(Shannon)` and
`hill_q2 = inverse Simpson`, which makes the values easier to interpret as
effective numbers of taxa.

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
