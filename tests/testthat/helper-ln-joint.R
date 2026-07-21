ln_joint_test_fixture <- function() {
  set.seed(20260720)
  n_subjects <- 16L
  cells_per_subject <- 24L
  n <- n_subjects * cells_per_subject
  id <- rep(seq_len(n_subjects), each = cells_per_subject)
  x <- rnorm(n)
  z_subject <- rbinom(n_subjects, 1, 0.5)
  z <- rep(z_subject, each = cells_per_subject)
  offset <- exp(runif(n, -0.25, 0.25))
  sigma2 <- 0.35
  phi <- 0.25
  omega <- rgamma(
    n_subjects,
    shape = 1 / (exp(sigma2) - 1),
    rate = 1 / (exp(sigma2 / 2) * (exp(sigma2) - 1))
  )
  cell_effect <- rgamma(n, shape = 1 / phi, rate = 1 / phi)
  mu <- offset * exp(-0.3 + 0.2 * x - 0.15 * z) *
    rep(omega, each = cells_per_subject) * cell_effect
  count <- Matrix::Matrix(matrix(rpois(n, mu), nrow = 1L), sparse = TRUE)
  rownames(count) <- "gene1"
  pred <- cbind(x = x, "(Intercept)" = 1, z = z)
  list(count = count, id = id, pred = pred, offset = offset)
}

ln_joint_prepare_fixture <- function(fx) {
  count <- as(fx$count, "dgCMatrix")
  nind <- ncol(count)
  centered <- center_m(as.matrix(fx$pred))
  offset_state <- cv_offset(as.double(fx$offset), 1, nind)
  id <- as.numeric(factor(
    as.character(fx$id),
    levels = unique(as.character(fx$id))
  ))
  k <- length(unique(id))
  fid <- as.integer(c(which(c(1, diff(id)) == 1), nind + 1L))
  cumsumy <- matrix(call_cumsumy(count, fid - 1L, k, 1L), nrow = 1L)
  posind <- list(which(cumsumy[1L, ] > 0))
  posv <- call_posindy(Matrix::t(count), 0L, nind)
  sds <- as.numeric(centered$sds)
  intcol <- which(sds == 0)
  public_scale <- sds
  public_scale[intcol] <- 1
  ctx <- list(
    pred = centered$pred,
    offset = offset_state$offset,
    fid = fid,
    cumsumy = cumsumy,
    posind = posind,
    nb = ncol(centered$pred),
    nind = nind,
    k = k,
    intcol = intcol,
    moffset = log(offset_state$mexpoffset),
    min = c(1e-4, 1e-4),
    max = c(10, 1000),
    allow_per_gene_switch = FALSE,
    use_betas = "betae",
    cutoff_cell = 20,
    cell_ind = integer(0),
    ncell = 0L,
    cv2 = offset_state$cv^2,
    kappa = 800,
    eps = 1e-6,
    opt = "lbfgs",
    cell_init = 1,
    output_re = FALSE,
    id = id,
    ln_estimator = "joint_ln",
    coef_internal_to_public = diag(
      1 / public_scale,
      nrow = length(public_scale)
    )
  )
  list(
    count = count,
    centered = centered,
    public_scale = public_scale,
    posv = posv,
    ctx = ctx
  )
}

ln_joint_fit_fixture <- function(output_re = FALSE,
                                 additional_tests = character(0),
                                 lrt_contrasts = NULL,
                                 max = c(10, 1000),
                                 ncore = 1L) {
  fx <- ln_joint_test_fixture()
  fit <- scrope(
    fx$count,
    fx$id,
    pred = fx$pred,
    offset = fx$offset,
    min = c(1e-4, 1e-4),
    max = max,
    verbose = FALSE,
    cpc = 0,
    mincp = 1,
    ncore = ncore,
    output_re = output_re,
    additional_tests = additional_tests,
    lrt_contrasts = lrt_contrasts,
    lrt_details = "lrt" %in% additional_tests,
    keep_diagnostics = TRUE,
    ln_estimator = "joint_ln"
  )
  list(fx = fx, fit = fit)
}
