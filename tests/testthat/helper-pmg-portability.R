# Independent reference for helper-pmg.R's balanced cell-predictor fixture.
# Every donor has the same x values and zero offsets, so the marginal likelihood
# factors into NB donor totals and multinomial counts conditional on each total.
# At the maximum, the NB mean is mean(totals); its shape is alpha = 1/expm1(s).
# The conditional slope solves E_b[x] = sum(x*y)/sum(y). The package's intercept
# is log(mean(totals)) - log(sum(exp(b*x))) - s/2 because E[donor effect]=exp(s/2).
pmg_balanced_cell_reference <- function(fx) {
  groups <- split(seq_along(fx$id), factor(fx$id, levels = unique(fx$id)))
  x <- as.numeric(fx$pred[, 2L])
  x_one <- x[groups[[1L]]]
  stopifnot(ncol(fx$pred) == 2L, fx$prep$ctx$intcol == 1L,
            all(fx$pred[, 1L] == 1), all(fx$prep$ctx$offset == 0),
            all(vapply(groups, function(i) identical(x[i], x_one), logical(1))))
  y <- as.numeric(fx$count)
  totals <- vapply(groups, function(i) sum(y[i]), numeric(1))
  mean_total <- mean(totals)
  target_x <- sum(y * x) / sum(y)
  slope_score <- function(b) {
    weights <- exp(b * x_one)
    sum(x_one * weights) / sum(weights) - target_x
  }
  slope <- uniroot(slope_score, c(-3, 3), tol = 1e-13)$root
  # For integer totals, digamma(t+alpha)-digamma(alpha) is this finite sum.
  # Avoid special-function subtraction in the independent scalar reference.
  shape_score <- function(alpha) {
    sum(vapply(totals, function(total) {
      sum(1 / (alpha + seq_len(total) - 1))
    }, numeric(1))) - length(totals) * log1p(mean_total / alpha)
  }
  shape <- uniroot(shape_score, c(0.1, 1000), tol = 1e-11)$root
  sigma <- log1p(1 / shape)
  intercept <- log(mean_total) - log(sum(exp(slope * x_one))) - sigma / 2
  theta <- c(intercept, slope * fx$prep$design_sds[2L], sigma)
  weights <- exp(slope * x_one)
  probabilities <- weights / sum(weights)
  conditional_information <- sum(y) * sum(probabilities * (x_one - target_x)^2)
  donor_slope_scores <- vapply(groups, function(i) {
    sum(y[i] * (x[i] - target_x))
  }, numeric(1))
  list(theta = theta, public_slope = slope, shape = shape,
       model_se = 1 / sqrt(conditional_information),
       robust_se = sqrt(sum(donor_slope_scores^2)) / conditional_information,
       shape_score = shape_score(shape), slope_score = slope_score(slope),
       mean_total = mean_total, totals = totals)
}

pmg_cell_endpoint_receipt <- function(raw, fx, legacy) {
  reference <- pmg_balanced_cell_reference(fx)
  theta <- c(raw$theta$beta, raw$theta$subVar)
  H <- raw$observed_info
  gradient <- as.numeric(raw$score)
  newton <- as.numeric(solve(H, gradient))
  # This independent expression includes the factorial constants omitted by the
  # native objective. Only differences within the same expression are compared.
  objective <- function(theta) {
    sigma <- theta[3L]
    alpha <- 1 / expm1(sigma)
    lambda <- alpha * exp(-sigma / 2)
    eta <- as.numeric(fx$prep$ctx$pred %*% theta[1:2] + fx$prep$ctx$offset)
    donor <- factor(fx$id, levels = unique(fx$id))
    total_mu <- as.numeric(rowsum(exp(eta), donor))
    -sum(lgamma(reference$totals + alpha) - lgamma(alpha) +
           alpha * log(lambda) -
           (reference$totals + alpha) * log(lambda + total_mu)) -
      sum(as.numeric(fx$count) * eta) + sum(lgamma(as.numeric(fx$count) + 1))
  }
  mean_predictor <- function(theta) {
    as.numeric(fx$prep$ctx$pred %*% theta[1:2] + fx$prep$ctx$offset) + theta[3L] / 2
  }
  list(
    reference = reference, theta = theta, gradient = gradient,
    observed_info = H, information_eigenvalues = eigen(H, symmetric = TRUE)$values,
    optimizer = raw$optimizer, optimizer_diagnostics = raw$diagnostics,
    newton = newton,
    max_scaled_newton = max(abs(newton) / pmax(1, abs(theta))),
    gradient_quadratic = as.numeric(crossprod(gradient, newton)),
    corrected_root_error = max(abs(theta - newton - reference$theta)),
    objective_delta = objective(theta) - objective(reference$theta),
    max_relative_prediction_error = max(abs(expm1(
      mean_predictor(theta) - mean_predictor(reference$theta)))),
    legacy = legacy,
    sigma_difference = abs(theta[3L] - legacy$overdispersion),
    # The historical endpoint is not an exact root. Account for its measured
    # root error and this endpoint's Newton correction, plus the independently
    # checked 1e-10 second-order/root-computation remainder.
    sigma_error_budget = abs(newton[3L]) +
      abs(legacy$overdispersion - reference$theta[3L]) + 1e-10
  )
}

pmg_cell_endpoint_is_precise <- function(receipt) {
  # These are test accuracy requirements, not runtime acceptance changes.
  # g'H^-1g <= 1e-8 bounds first-order error to 1e-4 information/SE units;
  # the scaled step bound is ten times tighter than the existing 1e-5
  # constrained-fit diagnostic threshold. The 1e-10 corrected-root check
  # independently guards against an inaccurate local quadratic approximation.
  isTRUE(all(is.finite(c(receipt$theta, receipt$gradient, receipt$newton)))) &&
    isTRUE(all(receipt$information_eigenvalues > 0)) &&
    isTRUE(receipt$optimizer_diagnostics$optimizer_success) &&
    isTRUE(receipt$max_scaled_newton <= 1e-6) &&
    isTRUE(receipt$gradient_quadratic >= 0 && receipt$gradient_quadratic <= 1e-8) &&
    isTRUE(receipt$corrected_root_error <= 1e-10) &&
    isTRUE(abs(receipt$objective_delta) <= 1e-8) &&
    isTRUE(receipt$max_relative_prediction_error <= 1e-6) &&
    isTRUE(receipt$sigma_difference <= receipt$sigma_error_budget)
}

pmg_write_cell_endpoint_receipt <- function(receipt, summary) {
  path <- Sys.getenv("SCROPE_TEST_DIAGNOSTICS_DIR",
                     unset = testthat::test_path("_diagnostics"))
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  record <- list(receipt = receipt, public_summary = summary,
                 session = utils::sessionInfo())
  saveRDS(record, file.path(path, "pgmm-cell-endpoint.rds"))
  # dput preserves more information than rounded failure messages.
  dput(record, file = file.path(path, "pgmm-cell-endpoint.txt"))
}
