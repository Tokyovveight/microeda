#' Ordinate beta diversity distances
#'
#' `microeda_beta_ordination()` computes a minimal principal coordinates
#' analysis (PCoA) ordination from a `microeda_beta` distance object using
#' [stats::cmdscale()].
#'
#' @param x A `microeda_beta` object.
#' @param method Ordination method. Only `"pcoa"` is currently supported.
#' @param dimensions Number of ordination axes to return.
#'
#' @return A `microeda_beta_ordination` object with coordinates, eigenvalues,
#'   variance explained, method metadata, optional group information, and matched
#'   call.
#' @examples
#' counts <- matrix(
#'   c(
#'     1, 2, 0,
#'     2, 1, 0,
#'     0, 0, 3
#'   ),
#'   nrow = 3,
#'   byrow = TRUE
#' )
#' rownames(counts) <- c("S1", "S2", "S3")
#' colnames(counts) <- paste0("ASV", 1:3)
#'
#' beta <- microeda_beta(counts, taxa_are_rows = FALSE)
#' ord <- microeda_beta_ordination(beta)
#' as_beta_coordinates(ord)
#' @export
microeda_beta_ordination <- function(x,
                                     method = "pcoa",
                                     dimensions = 2) {
  if (!inherits(x, "microeda_beta")) {
    stop("`x` must be a microeda_beta object.", call. = FALSE)
  }

  method <- validate_beta_ordination_method(method)
  dimensions <- validate_beta_ordination_dimensions(
    dimensions = dimensions,
    n_samples = length(x$sample_ids)
  )

  ordination <- stats::cmdscale(as_beta_dist(x), k = dimensions, eig = TRUE)
  coordinates <- beta_ordination_coordinates(
    points = ordination$points,
    sample_ids = x$sample_ids,
    dimensions = dimensions,
    group = x$group,
    group_values = x$group_values
  )

  structure(
    list(
      coordinates = coordinates,
      eigenvalues = unname(ordination$eig),
      variance_explained = beta_ordination_variance(ordination$eig),
      method = method,
      distance_method = x$method,
      dimensions = dimensions,
      sample_ids = x$sample_ids,
      group = x$group,
      group_values = x$group_values,
      call = match.call()
    ),
    class = "microeda_beta_ordination"
  )
}

#' Extract beta ordination coordinates
#'
#' @param x A `microeda_beta_ordination` object.
#'
#' @return A data frame with one row per sample and one column per ordination
#'   axis.
#' @export
as_beta_coordinates <- function(x) {
  if (!inherits(x, "microeda_beta_ordination")) {
    stop("`x` must be a microeda_beta_ordination object.", call. = FALSE)
  }

  x$coordinates
}

#' Ordinate beta diversity method comparisons
#'
#' `microeda_beta_compare_ordination()` computes one ordination for each
#' beta diversity result stored in a `microeda_beta_compare` object.
#'
#' PCoA axes from different distance methods are method-specific coordinate
#' systems. Use these coordinates for side-by-side inspection, downstream
#' plotting, and reports; do not interpret `Axis1` from one distance method as
#' the same axis as `Axis1` from another distance method.
#'
#' @param x A `microeda_beta_compare` object.
#' @param method Ordination method. Only `"pcoa"` is currently supported.
#' @param dimensions Number of ordination axes to return.
#'
#' @return A `microeda_beta_compare_ordination` object containing named
#'   `microeda_beta_ordination` results, method metadata, sample IDs, optional
#'   group information, and matched call.
#' @examples
#' counts <- matrix(
#'   c(
#'     1, 2, 0,
#'     2, 1, 0,
#'     0, 0, 3
#'   ),
#'   nrow = 3,
#'   byrow = TRUE
#' )
#' rownames(counts) <- c("S1", "S2", "S3")
#' colnames(counts) <- paste0("ASV", 1:3)
#'
#' beta_cmp <- microeda_beta_compare(counts, taxa_are_rows = FALSE)
#' ord_cmp <- microeda_beta_compare_ordination(beta_cmp)
#' as_beta_compare_coordinates(ord_cmp)
#' @export
microeda_beta_compare_ordination <- function(x,
                                             method = "pcoa",
                                             dimensions = 2) {
  if (!inherits(x, "microeda_beta_compare")) {
    stop("`x` must be a microeda_beta_compare object.", call. = FALSE)
  }

  results <- lapply(x$methods, function(method_name) {
    microeda_beta_ordination(
      x$results[[method_name]],
      method = method,
      dimensions = dimensions
    )
  })
  names(results) <- x$methods
  first_result <- results[[1]]

  structure(
    list(
      results = results,
      methods = x$methods,
      ordination_method = first_result$method,
      dimensions = first_result$dimensions,
      sample_ids = x$sample_ids,
      group = x$group,
      group_values = x$group_values,
      call = match.call()
    ),
    class = "microeda_beta_compare_ordination"
  )
}

#' Extract beta comparison ordination coordinates
#'
#' @param x A `microeda_beta_compare_ordination` object.
#'
#' @return A data frame with one row per method and sample. The first columns
#'   are `method`, `sample_id`, and one column per ordination axis. A `group`
#'   column is included when group metadata is present.
#' @export
as_beta_compare_coordinates <- function(x) {
  if (!inherits(x, "microeda_beta_compare_ordination")) {
    stop(
      "`x` must be a microeda_beta_compare_ordination object.",
      call. = FALSE
    )
  }

  rows <- lapply(x$methods, function(method_name) {
    beta_compare_coordinate_rows(
      method = method_name,
      coordinates = as_beta_coordinates(x$results[[method_name]])
    )
  })
  out <- do.call(rbind, rows)
  row.names(out) <- NULL
  out
}

validate_beta_ordination_method <- function(method) {
  supported_methods <- "pcoa"
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

validate_beta_ordination_dimensions <- function(dimensions, n_samples) {
  if (!is.numeric(dimensions) || length(dimensions) != 1 ||
      is.na(dimensions) || !is.finite(dimensions) ||
      dimensions != floor(dimensions)) {
    stop(
      "`dimensions` must be a single finite whole number.",
      call. = FALSE
    )
  }

  dimensions <- as.integer(dimensions)
  if (dimensions < 1) {
    stop("`dimensions` must be at least 1.", call. = FALSE)
  }

  max_dimensions <- n_samples - 1L
  if (max_dimensions < 1L) {
    stop("PCoA requires at least 2 samples.", call. = FALSE)
  }

  if (dimensions > max_dimensions) {
    stop(
      "`dimensions` must be no larger than ",
      max_dimensions,
      " for ",
      n_samples,
      " samples.",
      call. = FALSE
    )
  }

  dimensions
}

beta_ordination_coordinates <- function(points,
                                        sample_ids,
                                        dimensions,
                                        group = NULL,
                                        group_values = NULL) {
  points <- as.matrix(points)
  if (ncol(points) != dimensions) {
    points <- matrix(points, ncol = dimensions)
  }

  colnames(points) <- paste0("Axis", seq_len(dimensions))
  coordinates <- data.frame(
    sample_id = sample_ids,
    points,
    row.names = sample_ids,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  if (!is.null(group)) {
    coordinates$group <- unname(as.character(group_values))
  }

  coordinates
}

beta_ordination_variance <- function(eigenvalues) {
  positive_eigenvalues <- pmax(eigenvalues, 0)
  total_positive <- sum(positive_eigenvalues)
  if (total_positive == 0) {
    return(rep(NA_real_, length(eigenvalues)))
  }

  unname(positive_eigenvalues / total_positive)
}

beta_compare_coordinate_rows <- function(method, coordinates) {
  data.frame(
    method = rep(method, nrow(coordinates)),
    coordinates,
    row.names = NULL,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}
