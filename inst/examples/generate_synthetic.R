# Reproducible teaching data; see data-raw/SYNTHETIC_DESIGN.md in the source.
# Source synthetic_parameters.R first. No fitting or selection occurs here.
# Restores RNGkind and .Random.seed. Base R cannot restore a pre-existing
# Box-Muller hidden normal cache; the declared example uses Inversion.
generate_synthetic_examples <- function(parameters = synthetic_example_parameters()) {
  old_kind <- RNGkind()
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    do.call(RNGkind, as.list(old_kind))
    if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  p <- parameters
  RNGversion(p$rng_version)
  do.call(RNGkind, as.list(p$rng_kind))
  set.seed(p$seed)
  ng <- length(p$genes)
  nd <- p$donors
  nc <- p$cells_per_donor
  stopifnot(nc %% 2L == 0L, p$control_donors + p$treatment_donors == nd)
  n <- nd * nc
  donors <- sprintf("D%02d", seq_len(nd))
  sid <- rep(donors, each = nc)
  cells <- sprintf("%s_C%02d", sid, rep(seq_len(nc), nd))
  X1 <- rnorm(n)
  cc <- factor(unlist(lapply(seq_len(nd), function(i)
    sample(rep(c("case", "control"), each = nc / 2L))), use.names = FALSE),
    levels = c("case", "control"))
  X2 <- rep(c(rep(0, p$control_donors), rep(1, p$treatment_donors)), each = nc)
  offset <- exp(runif(n, p$offset_log_range[1], p$offset_log_range[2]))
  donor_u <- matrix(runif(ng * nd), nrow = ng,
                    dimnames = list(p$genes, donors))
  count_u <- matrix(runif(ng * n), nrow = ng,
                    dimnames = list(p$genes, cells))
  pred <- data.frame(X1 = X1, X2 = X2, cc = cc, row.names = cells)
  eta <- outer(p$beta_intercept, rep(1, n)) + outer(p$beta_X1, X1) +
    outer(p$beta_X2, X2) + outer(p$beta_cccontrol, as.numeric(cc == "control"))
  fixed_mu <- sweep(exp(eta), 2, offset, "*")
  v <- p$donor_variance
  donor_effects <- list(
    reference = matrix(qgamma(donor_u, shape = 1 / v, rate = 1 / v), ng,
                       dimnames = dimnames(donor_u)),
    misspecified = matrix(qlnorm(donor_u, meanlog = -log1p(v) / 2,
                                sdlog = sqrt(log1p(v))), ng,
                          dimnames = dimnames(donor_u)))
  conditional_mean <- lapply(donor_effects, function(w)
    fixed_mu * w[, rep(seq_len(nd), each = nc), drop = FALSE])
  examples <- lapply(conditional_mean, function(mu) {
    counts <- matrix(qnbinom(count_u, size = 1 / p$phi, mu = mu), ng,
                     dimnames = list(p$genes, cells))
    list(count = methods::as(Matrix::Matrix(counts, sparse = TRUE), "dgCMatrix"),
         sid = sid, pred = pred, offset = offset)
  })
  truth <- data.frame(gene = p$genes, beta_intercept = p$beta_intercept,
                      beta_X1 = p$beta_X1, beta_X2 = p$beta_X2,
                      beta_cccontrol = p$beta_cccontrol, phi = p$phi)
  c(examples, list(truth = truth, parameters = p, donor_effects = donor_effects,
                   conditional_mean = conditional_mean))
}

make_synthetic_seurat <- function(example) {
  if (!requireNamespace("Seurat", quietly = TRUE))
    stop("Regenerating the Seurat teaching object requires the suggested package Seurat.")
  metadata <- example$pred
  metadata$replicate <- example$sid
  metadata$celltype <- example$pred$cc
  metadata$tech <- factor(ifelse(example$pred$X2 == 0, "synthetic_control",
                                 "synthetic_treatment"))
  metadata$offset <- example$offset
  # Use a conventional RNA Assay without changing global Seurat options.
  assay <- Seurat::CreateAssayObject(counts = example$count)
  Seurat::CreateSeuratObject(counts = assay, project = "scRoPE_synthetic",
                             meta.data = metadata)
}
