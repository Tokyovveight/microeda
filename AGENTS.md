# AGENTS.md

This repository is an R package called `microeda`.

## Product goal

`microeda` is an evidence-guided exploratory data analysis package for microbiome count data.

It should inspect microbiome data first, describe statistical risks, and recommend method families with caveats instead of pretending that there is one universal best workflow.

## Input scope

The package should support:

- `phyloseq` objects;
- plain count matrices;
- sample metadata;
- taxonomy tables when available.

## Coding rules

- Use small functions.
- One function should do one thing.
- Avoid monolithic scripts.
- Do not silently modify user input.
- Return structured objects, not only printed text.
- Keep heavy microbiome/statistics packages optional when possible.
- Add tests for every exported function.
- Run `devtools::document()`, `devtools::test()`, and `devtools::check()` before commits.

## Statistical rules

- Separate QC, alpha diversity, beta diversity, and differential abundance.
- Treat microbiome data as compositional unless absolute abundance information is supplied.
- Report sparsity, prevalence, and library size imbalance before statistical interpretation.
- Do not recommend rarefying as a universal normalization method.
- Use rarefaction curves primarily as QC/sensitivity diagnostics.
- PERMANOVA results should be paired with dispersion diagnostics such as `betadisper`.
- Differential abundance methods should be compared, not treated as a single ground truth.
- ANCOM-BC2 and ALDEx2 should be interpreted with their assumptions and caveats.

## Current MVP

The current MVP includes:

- `microeda_check()`
- `microeda_alpha()`
- `microeda_alpha_compare()`
- `as_recommendations()`
- `as_alpha_table()`
- `as_alpha_summary()`
- `as_alpha_tests()`
- `as_alpha_pairwise()`

## Next priority

The next priority is a QC report module:

- HTML report template;
- rarefaction curves;
- library size plots;
- sparsity/prevalence plots;
- metadata completeness plots.
