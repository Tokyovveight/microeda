#' Compute beta diversity distances
#'
#' `microeda_beta()` computes a minimal sample-by-sample beta diversity distance
#' object from microbiome counts. The current implementation supports
#' Bray-Curtis, binary Jaccard, and Hellinger-transformed Euclidean distances.
#'
#' @inheritParams microeda_check
#' @param method Distance method. One of `"bray"`, `"jaccard"`, or
#'   `"hellinger"`.
#'
#' @details
#' Bray-Curtis distances are calculated directly as
#' `sum(abs(a - b)) / sum(a + b)`. Pairs where both samples have zero total
#' abundance are assigned distance `0`.
#'
#' Jaccard distances use binary presence/absence and are calculated as
#' `1 - intersection / union`. Pairs where both samples have no detected
#' features are assigned distance `0`.
#'
#' Hellinger distances are calculated as Euclidean distances after transforming
#' each sample to square-root relative abundances with
#' `sqrt(counts / rowSums(counts))`. Samples with zero library size are
#' pragmatically transformed to all-zero vectors. Hellinger distance is not a
#' log-ratio or compositional method.
#'
#' @return A `microeda_beta` object with the distance object, method, sample
#'   IDs, optional group information, count-type diagnostics, and matched call.
#' @examples
#' counts <- matrix(c(1, 2, 0, 2, 1, 0), nrow = 2, byrow = TRUE)
#' rownames(counts) <- c("S1", "S2")
#' colnames(counts) <- paste0("ASV", 1:3)
#'
#' beta <- microeda_beta(counts, taxa_are_rows = FALSE)
#' beta_jaccard <- microeda_beta(counts, taxa_are_rows = FALSE, method = "jaccard")
#' beta_hellinger <- microeda_beta(counts, taxa_are_rows = FALSE, method = "hellinger")
#' as_beta_dist(beta)
#' as_beta_matrix(beta)
#' as_beta_samples(beta)
#' microeda_beta_plot(beta)
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
      distance = beta_distance(counts, method = method),
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

#' Extract beta diversity distances
#'
#' @param x A `microeda_beta` object.
#'
#' @return A [stats::dist()] object.
#' @export
as_beta_dist <- function(x) {
  if (!inherits(x, "microeda_beta")) {
    stop("`x` must be a microeda_beta object.", call. = FALSE)
  }

  x$distance
}

#' Extract beta diversity distances as a matrix
#'
#' @param x A `microeda_beta` object.
#'
#' @return A square numeric matrix of beta diversity distances.
#' @export
as_beta_matrix <- function(x) {
  as.matrix(as_beta_dist(x))
}

#' Extract beta diversity sample labels
#'
#' @param x A `microeda_beta` object.
#'
#' @return A data frame with `sample_id` and, when present, `group`.
#' @export
as_beta_samples <- function(x) {
  if (!inherits(x, "microeda_beta")) {
    stop("`x` must be a microeda_beta object.", call. = FALSE)
  }

  out <- data.frame(
    sample_id = x$sample_ids,
    stringsAsFactors = FALSE
  )

  if (!is.null(x$group)) {
    out$group <- unname(as.character(x$group_values))
  }

  out
}

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

#' Plot beta diversity distances
#'
#' `microeda_beta_plot()` draws a minimal base R heatmap of beta diversity
#' distances. The current implementation supports `type = "heatmap"` only.
#'
#' @param x A `microeda_beta` object.
#' @param type Plot type. Only `"heatmap"` is currently supported.
#' @param ... Additional arguments passed to [graphics::image()].
#'
#' @return The value returned by [graphics::image()], invisibly.
#' @examples
#' counts <- matrix(c(1, 2, 0, 2, 1, 0), nrow = 2, byrow = TRUE)
#' rownames(counts) <- c("S1", "S2")
#' colnames(counts) <- paste0("ASV", 1:3)
#'
#' beta <- microeda_beta(counts, taxa_are_rows = FALSE)
#' microeda_beta_plot(beta)
#' @export
microeda_beta_plot <- function(x, type = "heatmap", ...) {
  if (!inherits(x, "microeda_beta")) {
    stop("`x` must be a microeda_beta object.", call. = FALSE)
  }

  type <- validate_beta_plot_type(type)

  if (identical(type, "heatmap")) {
    return(plot_beta_heatmap(x, ...))
  }
}

validate_beta_method <- function(method) {
  supported_methods <- c("bray", "jaccard", "hellinger")
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

beta_distance <- function(counts, method) {
  if (identical(method, "bray")) {
    return(beta_bray_distance(counts))
  }

  if (identical(method, "jaccard")) {
    return(beta_jaccard_distance(counts))
  }

  beta_hellinger_distance(counts)
}

validate_beta_plot_type <- function(type) {
  supported_types <- "heatmap"
  if (!is.character(type) || length(type) != 1 ||
      is.na(type) || !nzchar(type) || !type %in% supported_types) {
    stop(
      "`type` must be one of: ",
      paste0("\"", supported_types, "\"", collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  type
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

beta_jaccard_distance <- function(counts) {
  presence <- counts > 0
  n_samples <- nrow(presence)
  distances <- matrix(
    0,
    nrow = n_samples,
    ncol = n_samples,
    dimnames = list(rownames(counts), rownames(counts))
  )

  if (n_samples > 1) {
    for (i in seq_len(n_samples - 1)) {
      for (j in seq.int(i + 1, n_samples)) {
        union <- sum(presence[i, ] | presence[j, ])
        distance <- 0
        if (union > 0) {
          intersection <- sum(presence[i, ] & presence[j, ])
          distance <- 1 - intersection / union
        }
        distances[i, j] <- distance
        distances[j, i] <- distance
      }
    }
  }

  stats::as.dist(distances)
}

beta_hellinger_distance <- function(counts) {
  library_sizes <- rowSums(counts)
  transformed <- counts
  transformed[] <- 0

  positive_libraries <- library_sizes > 0
  if (any(positive_libraries)) {
    relative_abundance <- sweep(
      counts[positive_libraries, , drop = FALSE],
      1,
      library_sizes[positive_libraries],
      "/"
    )
    transformed[positive_libraries, ] <- sqrt(relative_abundance)
  }

  stats::dist(transformed, method = "euclidean")
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

plot_beta_heatmap <- function(x, ...) {
  distances <- as_beta_matrix(x)
  sample_ids <- rownames(distances)
  positions <- seq_along(sample_ids)

  dots <- list(...)
  draw_axes <- TRUE
  if (!is.null(dots$axes)) {
    draw_axes <- isTRUE(dots$axes)
  }
  dots$axes <- FALSE
  if (is.null(dots$main)) {
    dots$main <- paste("Beta diversity distances:", x$method)
  }
  if (is.null(dots$xlab)) {
    dots$xlab <- "Sample"
  }
  if (is.null(dots$ylab)) {
    dots$ylab <- "Sample"
  }

  result <- do.call(
    graphics::image,
    c(list(x = positions, y = positions, z = distances), dots)
  )

  if (isTRUE(draw_axes)) {
    graphics::axis(1, at = positions, labels = sample_ids, las = 2)
    graphics::axis(2, at = positions, labels = sample_ids, las = 2)
    graphics::box()
  }

  invisible(result)
}
