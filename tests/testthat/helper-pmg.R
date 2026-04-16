pmg_test_fixture <- function() {
  prepare_fn <- getFromNamespace("hl_prepare_data", "scRoPE")
  posindy_fn <- getFromNamespace("call_posindy", "scRoPE")

  count <- matrix(
    c(0, 1, 3,
      2, 0, 1,
      1, 4, 0,
      3, 1, 2),
    nrow = 1,
    byrow = TRUE
  )
  rownames(count) <- "gene1"

  id <- rep(paste0("s", 1:4), each = 3)
  x <- rep(c(-1, 0, 1), 4)
  pred <- model.matrix(~x)

  prep <- prepare_fn(
    count = count,
    id = id,
    pred = pred,
    cpc = 0,
    mincp = 1,
    verbose = FALSE
  )
  gene_index <- prep$gid[[1]]
  posv <- posindy_fn(prep$count, gene_index - 1L, prep$ctx$nind)

  list(
    count = count,
    id = id,
    pred = pred,
    prep = prep,
    gene_index = gene_index,
    posv = posv
  )
}

pmg_interior_test_fixture <- function() {
  prepare_fn <- getFromNamespace("hl_prepare_data", "scRoPE")
  posindy_fn <- getFromNamespace("call_posindy", "scRoPE")

  count <- matrix(
    c(1, 4, 1, 5,
      4, 2, 9, 12,
      0, 3, 11, 14,
      1, 1, 7, 11,
      3, 2, 11, 17,
      2, 4, 6, 12),
    nrow = 1,
    byrow = TRUE
  )
  rownames(count) <- "gene1"

  id <- rep(paste0("s", 1:6), each = 4)
  x <- rep(c(-1.5, -0.5, 0.5, 1.5), 6)
  pred <- model.matrix(~x)

  prep <- prepare_fn(
    count = count,
    id = id,
    pred = pred,
    cpc = 0,
    mincp = 1,
    verbose = FALSE
  )
  gene_index <- prep$gid[[1]]
  posv <- posindy_fn(prep$count, gene_index - 1L, prep$ctx$nind)

  list(
    count = count,
    id = id,
    pred = pred,
    prep = prep,
    gene_index = gene_index,
    posv = posv
  )
}

pmg_subject_level_test_fixture <- function() {
  prepare_fn <- getFromNamespace("hl_prepare_data", "scRoPE")
  posindy_fn <- getFromNamespace("call_posindy", "scRoPE")

  count <- matrix(
    c(1, 4, 1, 5,
      4, 2, 9, 12,
      0, 3, 11, 14,
      1, 1, 7, 11,
      3, 2, 11, 17,
      2, 4, 6, 12),
    nrow = 1,
    byrow = TRUE
  )
  rownames(count) <- "gene1"

  id <- rep(paste0("s", 1:6), each = 4)
  z_subject <- rep(c(-1.5, -0.9, -0.2, 0.3, 0.8, 1.4), each = 4)
  pred <- model.matrix(~z_subject)

  prep <- prepare_fn(
    count = count,
    id = id,
    pred = pred,
    cpc = 0,
    mincp = 1,
    verbose = FALSE
  )
  gene_index <- prep$gid[[1]]
  posv <- posindy_fn(prep$count, gene_index - 1L, prep$ctx$nind)

  list(
    count = count,
    id = id,
    pred = pred,
    prep = prep,
    gene_index = gene_index,
    posv = posv
  )
}

pmg_truth_sim_fixture <- function(seed = 20260416,
                                  ng = 40L,
                                  cps = 20L,
                                  sig2 = 0.2,
                                  lambda_s = 0.5,
                                  effcell = 0.35,
                                  effsub = 0) {
  prepare_fn <- getFromNamespace("hl_prepare_data", "scRoPE")
  posindy_fn <- getFromNamespace("call_posindy", "scRoPE")

  set.seed(seed)

  ng <- as.integer(ng)
  cps <- as.integer(cps)
  n <- ng * cps

  x <- rnorm(n)
  z_subject <- rbinom(ng, 1, 0.5)
  z <- rep(z_subject, each = cps)
  id <- rep(paste0("s", sprintf("%02d", seq_len(ng))), each = cps)

  if (sig2 > 0) {
    shape_subj <- 1 / (exp(sig2) - 1)
    rate_subj <- 1 / (exp(sig2 / 2) * (exp(sig2) - 1))
    omega_subject <- rgamma(ng, shape = shape_subj, rate = rate_subj)
  } else {
    omega_subject <- rep(1, ng)
  }
  ref <- rep(omega_subject, each = cps)

  eta <- lambda_s + effcell * x + effsub * z + log(ref)
  y <- rpois(n, lambda = exp(eta))

  count <- Matrix::Matrix(matrix(y, nrow = 1), sparse = TRUE)
  rownames(count) <- "gene1"
  colnames(count) <- paste0("cell", seq_len(n))
  pred <- model.matrix(~x + z)

  prep <- prepare_fn(
    count = count,
    id = id,
    pred = pred,
    cpc = 0,
    mincp = 1,
    verbose = FALSE
  )
  gene_index <- prep$gid[[1]]
  posv <- posindy_fn(prep$count, gene_index - 1L, prep$ctx$nind)

  predictor_names <- colnames(pred)
  truth_beta_original <- c("(Intercept)" = lambda_s, x = effcell, z = effsub)
  truth_beta_scaled <- setNames(numeric(length(predictor_names)), predictor_names)
  for (jj in seq_along(predictor_names)) {
    nm <- predictor_names[[jj]]
    if (identical(nm, "(Intercept)")) {
      truth_beta_scaled[[jj]] <- lambda_s
    } else {
      truth_beta_scaled[[jj]] <- truth_beta_original[[nm]] * prep$design_sds[[jj]]
    }
  }

  list(
    count = count,
    id = id,
    pred = pred,
    prep = prep,
    gene_index = gene_index,
    posv = posv,
    truth = list(
      sig2 = sig2,
      lambda_s = lambda_s,
      effcell = effcell,
      effsub = effsub,
      beta_original = truth_beta_original,
      beta_scaled = truth_beta_scaled
    )
  )
}
