#' Compute beta diversity distances
#'
#' `microeda_beta()` computes a minimal sample-by-sample beta diversity distance
#' object from microbiome counts. The first implementation supports
#' Bray-Curtis distances only.
#'
#' @inheritParams microeda_check
#' @param method Distance method. Only `"bray"` is currently supported.
#'
#' @details
#' Bray-Curtis distances are calculated directly as
#' `sum(abs(a - b)) / sum(a + b)`. Pairs where both samples have zero total
#' abundance are assigned distance `0`.
#'
#' @return A `microeda_beta` object with the distance object, method, sample
#'   IDs, optional group information, count-type diagnostics, and matched call.
#' @examples
#' counts <- matrix(c(1, 2, 0, 2, 1, 0), nrow = 2, byrow = TRUE)
#' rownames(counts) <- c("S1", "S2")
#' colnames(counts) <- paste0("ASV", 1:3)
#'
#' beta <- microeda_beta(counts, taxa_are_rows = FALSE)
#' beta$distance
#' @export
microeda_beta <- function(x,
                          metadata = NULL,
                          group = NULL,
                          taxa_are_rows = TRUE,
                          method = "bray") {
  method <- validate_beta_method(method)
  extracted <- microeda_extract(
    x = x,
    metadata = metadata,
    taxonomy = NULL,
    taxa_are_rows = taxa_are_rows
  )

  counts <- extracted$counts
  validate_beta_counts(counts)

  sample_ids <- rownames(counts)
  group_values <- resolve_beta_group(
    metadata = extracted$metadata,
    group = group,
    sample_ids = sample_ids
  )
  library_sizes <- rowSums(counts)

  structure(
    list(
      distance = beta_bray_distance(counts),
      method = method,
      sample_ids = sample_ids,
      group = group,
      group_values = group_values,
      count_type = diagnose_count_type(counts, library_sizes),
      call = match.call()
    ),
    class = "microeda_beta"
  )
}

validate_beta_method <- function(method) {
  supported_methods <- "bray"
  if (!is.character(method) || length(method) != 1 ||
      is.na(method) || !nzchar(method) || !method %in% supported_methods) {
    stop(
      "`method` must be one of: ",
      paste0("\"", supported_methods, "\"", collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  method
}

validate_beta_counts <- function(counts) {
  if (any(!is.finite(counts))) {
    stop("`counts` must contain only finite numeric values.", call. = FALSE)
  }

  if (any(counts < 0)) {
    stop("`counts` cannot contain negative values.", call. = FALSE)
  }

  invisible(counts)
}

resolve_beta_group <- function(metadata, group = NULL, sample_ids) {
  if (is.null(group)) {
    return(NULL)
  }

  if (!is.character(group) || length(group) != 1 ||
      is.na(group) || !nzchar(group)) {
    stop("`group` must be a single non-missing character string.", call. = FALSE)
  }

  if (is.null(metadata)) {
    stop("`metadata` is required when `group` is supplied.", call. = FALSE)
  }

  if (!group %in% colnames(metadata)) {
    stop("`group` is not a column in `metadata`.", call. = FALSE)
  }

  group_values <- metadata[[group]]
  names(group_values) <- sample_ids
  group_values
}

beta_bray_distance <- function(counts) {
  n_samples <- nrow(counts)
  distances <- matrix(
    0,
    nrow = n_samples,
    ncol = n_samples,
    dimnames = list(rownames(counts), rownames(counts))
  )

  if (n_samples > 1) {
    for (i in seq_len(n_samples - 1)) {
      for (j in seq.int(i + 1, n_samples)) {
        denominator <- sum(counts[i, ] + counts[j, ])
        distance <- 0
        if (denominator > 0) {
          distance <- sum(abs(counts[i, ] - counts[j, ])) / denominator
        }
        distances[i, j] <- distance
        distances[j, i] <- distance
      }
    }
  }

  stats::as.dist(distances)
}
