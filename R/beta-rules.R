beta_compare_rule_context <- function(x) {
  if (!inherits(x, "microeda_beta_compare")) {
    stop("`x` must be a microeda_beta_compare object.", call. = FALSE)
  }

  method_summary <- as_beta_compare_summary(x)
  correlations <- as_beta_compare_distance_correlations(x)

  group_summary <- beta_empty_group_summary()
  if (!is.null(x$group)) {
    group_summary <- as_beta_compare_group_summary(x)
  }

  summary <- data.frame(
    n_methods = length(x$methods),
    methods = paste(method_summary$method, collapse = ", "),
    n_samples = length(x$sample_ids),
    group = if (is.null(x$group)) NA_character_ else x$group,
    has_group = !is.null(x$group),
    has_group_summary = nrow(group_summary) > 0,
    has_distance_correlations = nrow(correlations) > 0,
    stringsAsFactors = FALSE
  )

  list(
    summary = summary,
    methods = method_summary[c("method", "n_samples", "n_pairs")],
    caveats = beta_compare_rule_caveats(
      methods = method_summary$method,
      has_distance_correlations = nrow(correlations) > 0
    )
  )
}

beta_compare_rule_caveats <- function(methods, has_distance_correlations) {
  rows <- list()

  if ("jaccard" %in% methods) {
    rows[[length(rows) + 1L]] <- beta_compare_rule_caveat(
      context_id = "jaccard_incidence",
      topic = "beta_method",
      method = "jaccard",
      message = paste(
        "Jaccard is incidence-based and complementary;",
        "it does not use abundance differences."
      )
    )
  }

  if ("hellinger" %in% methods) {
    rows[[length(rows) + 1L]] <- beta_compare_rule_caveat(
      context_id = "hellinger_not_log_ratio",
      topic = "beta_method",
      method = "hellinger",
      message = paste(
        "Hellinger is an abundance-profile ecological distance after",
        "square-root relative abundance transformation; it is not a",
        "log-ratio compositional distance."
      )
    )
  }

  if (isTRUE(has_distance_correlations)) {
    rows[[length(rows) + 1L]] <- beta_compare_rule_caveat(
      context_id = "distance_correlations_descriptive",
      topic = "beta_comparison",
      method = NA_character_,
      message = paste(
        "Distance correlations are descriptive and do not validate a method",
        "or identify the correct method."
      )
    )
  }

  if (length(methods) > 1L) {
    rows[[length(rows) + 1L]] <- beta_compare_rule_caveat(
      context_id = "pcoa_axes_method_specific",
      topic = "ordination",
      method = NA_character_,
      message = paste(
        "PCoA axes are method-specific and should not be compared across",
        "distance methods as identical constructs."
      )
    )
  }

  if (length(rows) == 0L) {
    return(beta_empty_rule_caveats())
  }

  out <- do.call(rbind, rows)
  row.names(out) <- NULL
  out
}

beta_compare_rule_caveat <- function(context_id, topic, method, message) {
  data.frame(
    context_id = context_id,
    topic = topic,
    method = method,
    severity = "info",
    message = message,
    stringsAsFactors = FALSE
  )
}

beta_empty_rule_caveats <- function() {
  data.frame(
    context_id = character(),
    topic = character(),
    method = character(),
    severity = character(),
    message = character(),
    stringsAsFactors = FALSE
  )
}
