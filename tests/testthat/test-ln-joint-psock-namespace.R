ln_joint_normalise_psock_objective_ids <- function(x) {
  if (!is.list(x)) {
    return(x)
  }
  object_names <- names(x)
  for (index in seq_along(x)) {
    value <- x[[index]]
    field <- if (is.null(object_names)) "" else object_names[[index]]
    if (identical(field, "objective_instance") && is.list(value) &&
        identical(
          value$scheme,
          "process_local_exact_joint_ln_objective_context_v1"
        )) {
      value$id <- paste0("<run-local-id-gene-", value$gene_index, ">")
    }
    x[index] <- list(ln_joint_normalise_psock_objective_ids(value))
  }
  x
}

test_that("joint-LN dispatchers expose their fitter dependencies", {
  namespace <- asNamespace("scRoPE")
  expected <- c(
    fit_gene_unconstrained = "fit_gene_joint_ln_unconstrained",
    fit_gene_constrained = "fit_gene_joint_ln_constrained"
  )

  for (name in names(expected)) {
    dispatcher <- get(name, envir = namespace, inherits = FALSE)
    symbols <- all.names(
      body(dispatcher),
      functions = TRUE,
      unique = TRUE
    )
    expect_true(expected[[name]] %in% symbols)
    expect_false("getFromNamespace" %in% symbols)
  }
})

test_that("installed joint-LN PSOCK fits use the running package", {
  package_path <- normalizePath(
    getNamespaceInfo(asNamespace("scRoPE"), "path"),
    mustWork = FALSE
  )
  skip_if_not(
    file.exists(file.path(package_path, "Meta", "package.rds")),
    "requires an installed package"
  )
  skip_if(parallelly::availableCores() < 2L)

  contrast <- list(list(
    name = "psock_x",
    L = matrix(c(1, 0, 0), nrow = 1L),
    b = 0
  ))
  one <- ln_joint_fit_fixture(
    additional_tests = c("score", "lrt"),
    lrt_contrasts = contrast,
    ncore = 1L
  )$fit

  original_libraries <- .libPaths()
  on.exit(.libPaths(original_libraries), add = TRUE)
  package_library <- normalizePath(dirname(package_path), mustWork = FALSE)
  other_libraries <- original_libraries[
    normalizePath(original_libraries, mustWork = FALSE) != package_library
  ]
  alternate_package <- tryCatch(
    find.package("scRoPE", lib.loc = other_libraries),
    error = function(e) character(0)
  )
  if (length(alternate_package)) {
    # Reproduce a lib.loc install while an older default scRoPE is visible.
    .libPaths(other_libraries)
  }

  two <- ln_joint_fit_fixture(
    additional_tests = c("score", "lrt"),
    lrt_contrasts = contrast,
    ncore = 2L
  )$fit

  expect_equal(one$summary, two$summary, tolerance = 1e-8)
  expect_equal(one$overdispersion, two$overdispersion, tolerance = 1e-8)
  expect_identical(one$convergence, two$convergence)
  expect_equal(
    one$provenance$per_gene,
    two$provenance$per_gene,
    tolerance = 1e-8
  )
  expect_equal(one$lrt_details, two$lrt_details, tolerance = 1e-8)
  expect_equal(
    ln_joint_normalise_psock_objective_ids(one$diagnostics$raw),
    ln_joint_normalise_psock_objective_ids(two$diagnostics$raw),
    tolerance = 1e-8
  )
})
