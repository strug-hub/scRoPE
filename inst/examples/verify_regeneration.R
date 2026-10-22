# Source-only example verifier; no package namespace exports or dependencies.
# Fixed counts, design and RNG state remain exact. Only the named continuous
# fields below permit last-bit differences between numerical environments.
.verify_example_regeneration <- function(actual, expected, continuous, label) {
  attribute_set <- function(x) {
    a <- attributes(x)
    # Match identical()'s attribute-set semantics: attribute values and names
    # are exact, but ordinary data.frame replacement may reorder attributes.
    a[sort(names(a))]
  }
  compare <- function(a, b, path) {
    fail <- function(detail) stop(path, ": ", detail, call. = FALSE)
    if (!identical(typeof(a), typeof(b)) ||
        !identical(length(a), length(b)) ||
        !identical(attribute_set(a), attribute_set(b))) {
      fail("type, length or attributes differ")
    }
    if (path %in% continuous) {
      if (!is.double(a) || !all(is.finite(a)) || !all(is.finite(b)))
        fail("continuous values must be finite doubles")
      bound <- 64 * .Machine$double.eps * pmax(1, abs(a), abs(b))
      if (any(abs(a - b) > bound))
        fail("continuous difference exceeds the 64-epsilon scaled bound")
    } else if (is.list(a)) {
      for (i in seq_along(a)) {
        suffix <- if (is.null(names(a))) paste0("[[", i, "]]") else
          paste0("$", names(a)[i])
        compare(a[[i]], b[[i]], paste0(path, suffix))
      }
    } else if (!identical(a, b)) {
      fail("values must match exactly")
    }
  }
  compare(actual, expected, label)
  invisible(TRUE)
}

verify_synthetic_regeneration <- function(actual, expected) {
  .verify_example_regeneration(actual, expected,
    c("sample_data$pred$X1", "sample_data$offset"), "sample_data")
}

verify_toy24_regeneration <- function(actual, expected) {
  .verify_example_regeneration(actual, expected,
    c("toy_data$donor_truth$intercept_multiplier", "toy_data$donor_truth$donor_slope"),
    "toy_data")
}
