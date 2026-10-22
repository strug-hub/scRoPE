test_that("bundled synthetic counts retain the declared teaching design", {
  data(sample_data)

  counts <- as.matrix(sample_data$count)
  expect_identical(dim(counts), c(8L, 1200L))
  expect_identical(rownames(counts), sprintf("G%02d", seq_len(8L)))
  expect_true(all(is.finite(counts)))
  expect_true(all(counts >= 0 & counts == floor(counts)))
  expect_false(anyDuplicated(colnames(counts)) > 0L)
  expect_identical(rownames(sample_data$pred), colnames(counts))

  donor <- as.character(sample_data$sid)
  expect_length(donor, ncol(counts))
  expect_length(unique(donor), 30L)
  expect_identical(rle(donor)$lengths, rep(40L, 30L))
  expect_true(all(vapply(split(sample_data$pred$X2, donor),
                        function(x) length(unique(x)) == 1L, logical(1))))
  treatment <- sample_data$pred$X2[match(unique(donor), donor)]
  expect_equal(sum(treatment == 0), 20L)
  expect_equal(sum(treatment == 1), 10L)
  expect_true(all(is.finite(sample_data$pred$X1)))
  expect_equal(sum(sample_data$pred$cc == "case"), 600L)
  expect_equal(sum(sample_data$pred$cc == "control"), 600L)
  expect_true(all(vapply(split(sample_data$pred$cc, donor), function(x) {
    sum(x == "case") == 20L && sum(x == "control") == 20L
  }, logical(1))))
  expect_identical(
    colnames(model.matrix(~ X1 + X2 + cc, data = sample_data$pred)),
    c("(Intercept)", "X1", "X2", "cccontrol")
  )
  expect_length(sample_data$offset, ncol(counts))
  expect_true(all(is.finite(sample_data$offset) & sample_data$offset > 0))
})

synthetic_example_test_functions <- function() {
  directory <- system.file("examples", package = "scRoPE")
  if (!file.exists(file.path(directory, "generate_synthetic.R"))) {
    directory <- testthat::test_path("..", "..", "inst", "examples")
  }
  functions <- new.env(parent = globalenv())
  sys.source(file.path(directory, "synthetic_parameters.R"), envir = functions)
  sys.source(file.path(directory, "generate_synthetic.R"), envir = functions)
  sys.source(file.path(directory, "verify_regeneration.R"), envir = functions)
  functions
}

test_that("the generator preserves exact counts and bounded continuous values", {
  data(sample_data)
  functions <- synthetic_example_test_functions()
  parameters <- functions$synthetic_example_parameters()
  generated <- functions$generate_synthetic_examples(parameters)
  replay <- functions$generate_synthetic_examples(parameters)

  expect_identical(parameters$seed, 20261020L)
  expect_identical(parameters$design_version, "synthetic-teaching-v1")
  expect_identical(generated$parameters, parameters)
  expect_identical(generated, replay)
  expect_true(functions$verify_synthetic_regeneration(generated$reference, sample_data))
  expect_identical(generated$truth$gene, parameters$genes)
  expect_identical(generated$truth$beta_X2, parameters$beta_X2)
  expect_identical(generated$truth$phi, parameters$phi)
  for (component in c("sid", "pred", "offset")) {
    expect_identical(generated$reference[[component]],
                     generated$misspecified[[component]])
  }
  expect_identical(dim(generated$misspecified$count),
                   dim(generated$reference$count))
  expect_identical(dimnames(generated$misspecified$count),
                   dimnames(generated$reference$count))
})

test_that("synthetic verification allows only bounded continuous differences", {
  data(sample_data)
  verify <- synthetic_example_test_functions()$verify_synthetic_regeneration
  changed <- sample_data
  changed$pred$X1[1] <- changed$pred$X1[1] +
    4 * .Machine$double.eps * max(1, abs(changed$pred$X1[1]))
  changed$offset[1] <- changed$offset[1] + 4 * .Machine$double.eps
  expect_true(verify(changed, sample_data))
  changed$pred$X1[1] <- sample_data$pred$X1[1] +
    128 * .Machine$double.eps * max(1, abs(sample_data$pred$X1[1]))
  expect_error(verify(changed, sample_data), "64-epsilon scaled bound")
  changed$pred$X1[1] <- sample_data$pred$X1[1] + 1e-8
  expect_error(verify(changed, sample_data), "continuous difference")
  changed <- sample_data
  changed$count[1, 1] <- changed$count[1, 1] + 1
  expect_error(verify(changed, sample_data), "sample_data\\$count:")
  changed <- sample_data
  changed$pred$X2[1] <- 1 - changed$pred$X2[1]
  expect_error(verify(changed, sample_data), "X2.*exactly")
  changed <- sample_data
  names(changed$offset) <- seq_along(changed$offset)
  expect_error(verify(changed, sample_data), "offset.*attributes")
  changed <- sample_data
  changed$pred$X1 <- as.integer(changed$pred$X1)
  expect_error(verify(changed, sample_data), "X1.*type")
  changed <- sample_data
  changed$offset[1] <- NA_real_
  expect_error(verify(changed, sample_data), "finite doubles")
})

test_that("the generator preserves RNG settings and seed with Inversion", {
  functions <- synthetic_example_test_functions()
  local({
    original_kind <- RNGkind()
    original_seed_exists <- exists(".Random.seed", envir = .GlobalEnv,
                                    inherits = FALSE)
    if (original_seed_exists) original_seed <- get(".Random.seed", .GlobalEnv)
    on.exit({
      do.call(RNGkind, as.list(original_kind))
      if (original_seed_exists) {
        assign(".Random.seed", original_seed, envir = .GlobalEnv)
      } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    }, add = TRUE)

    # This checks settings and .Random.seed with Inversion. R does not expose
    # the hidden Box-Muller normal cache for a generator to save and restore.
    RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection")
    set.seed(2718)
    caller_kind <- RNGkind()
    caller_seed <- get(".Random.seed", .GlobalEnv)
    functions$generate_synthetic_examples()
    expect_identical(RNGkind(), caller_kind)
    expect_identical(get(".Random.seed", .GlobalEnv), caller_seed)

    rm(".Random.seed", envir = .GlobalEnv)
    functions$generate_synthetic_examples()
    expect_identical(RNGkind(), caller_kind)
    expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
  })
})

test_that("the synthetic Seurat constructor reproduces the bundled contents", {
  skip_if_not_installed("Seurat")
  data(sample_data)
  data(sample_seurat)
  functions <- synthetic_example_test_functions()
  reconstructed <- functions$make_synthetic_seurat(sample_data)
  # Seurat can record version-specific object metadata. Compare the teaching
  # content and public conversion, not serialization bytes or software stamps.
  expect_identical(reconstructed[[]], sample_seurat[[]])
  convert <- function(object) scToNeb(
    object, assay = "RNA", id = "replicate", pred = c("X1", "X2", "cc"),
    offset = "offset", verbose = FALSE
  )
  expect_identical(convert(reconstructed), convert(sample_seurat))
})

test_that("matching synthetic Seurat data roundtrips through scToNeb", {
  skip_if_not_installed("Seurat")
  data(sample_data)
  data(sample_seurat)

  expect_s4_class(sample_seurat, "Seurat")
  expect_identical(colnames(sample_seurat), colnames(sample_data$count))
  expect_true(all(c("replicate", "celltype", "tech", "X1", "X2", "cc",
                    "offset") %in% colnames(sample_seurat[[]])))

  conversion_warnings <- character(0)
  converted <- withCallingHandlers(
    scToNeb(
      sample_seurat, assay = "RNA", id = "replicate",
      pred = c("X1", "X2", "cc"), offset = "offset", verbose = FALSE
    ),
    warning = function(w) {
      conversion_warnings <<- c(conversion_warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  # Capture lifecycle/deprecation warnings as well as ordinary warnings.
  expect_identical(conversion_warnings, character(0))
  expect_identical(as.matrix(converted$count), as.matrix(sample_data$count))
  expect_identical(converted$pred, sample_data$pred[c("X1", "X2", "cc")])
  expect_identical(converted$id, as.character(sample_data$sid))
  expect_identical(converted$offset, as.double(sample_data$offset))
  expect_identical(
    model.matrix(~ X1 + X2 + cc, data = converted$pred),
    model.matrix(~ X1 + X2 + cc, data = sample_data$pred)
  )
})
