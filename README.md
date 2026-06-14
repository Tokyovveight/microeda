# microeda

`microeda` is an early R package skeleton for evidence-guided exploratory
analysis of microbiome count data.

The core idea is deliberately conservative:

> inspect the data first, describe the risks, then recommend method families
> with caveats and evidence instead of pretending there is one universal best
> microbiome workflow.

## MVP scope

- Accept a `phyloseq` object or a plain count matrix plus metadata.
- Summarize library sizes, sparsity, prevalence, metadata completeness, group
  sizes, and taxonomy completeness.
- Flag common microbiome EDA risks:
  - compositional input and relative-abundance interpretation
  - non-integer or transformed values used with count-based methods
  - high zero fraction and low prevalence
  - uneven sequencing depth
  - low group sizes
  - missing metadata or taxonomy ranks
- Return a recommendations table with `rule_id`, `topic`, `severity`,
  `recommendation`, `caveat`, and `evidence`.

## Example

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

report <- microeda_check(
  counts,
  metadata = metadata,
  group = "group",
  taxa_are_rows = FALSE
)

report
as_recommendations(report)
```

## Alpha diversity

```r
alpha <- microeda_alpha(
  counts,
  metadata = metadata,
  group = "group",
  taxa_are_rows = FALSE
)

alpha
as_alpha_table(alpha)
as_alpha_summary(alpha)
```

The alpha table includes classic indices (`observed`, `chao1`, `shannon`,
`simpson`, `inverse_simpson`) and Hill/effective-diversity equivalents
(`hill_q0`, `hill_q1`, `hill_q2`). In practice, `hill_q1 = exp(Shannon)` and
`hill_q2 = inverse Simpson`, which makes the values easier to interpret as
effective numbers of taxa.

Group comparisons are intentionally exploratory and include depth/coverage
diagnostics:

```r
alpha_cmp <- microeda_alpha_compare(alpha)

alpha_cmp
as_alpha_tests(alpha_cmp)
as_alpha_pairwise(alpha_cmp)
```

For `phyloseq`, pass the object directly:

```r
report <- microeda_check(ps, group = "Treatment")
alpha <- microeda_alpha(ps, group = "Treatment")
alpha_cmp <- microeda_alpha_compare(alpha)
```

## Evidence rules

The package keeps the current evidence map in
`inst/extdata/evidence_rules.yml`. The R functions use a built-in version of
the same rules so the MVP does not require a YAML parser at runtime.
