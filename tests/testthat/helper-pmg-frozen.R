# Fixed outputs from the preserved scRoPE 0.2.5.9005 nebula(model="PMM")
# runtime for helper-pmg.R's deterministic fixtures. The retired fitting route
# must not be called by the candidate test suite. Independent analytic PGMM
# score/Hessian/covariance checks remain in test-pmg-core.R.
pmg_frozen_legacy <- function(which) {
  values <- list(
    core = list(
      summary = data.frame(
        "logFC_(Intercept)" = 0.405415102469892, logFC_x = 0,
        "se_(Intercept)" = 0.235755290754404, se_x = 0.288675135408261,
        check.names = FALSE),
      overdispersion = 1e-4, convergence = 1, algorithm = "PGMM"
    ),
    cell = list(
      summary = data.frame(logFC_x = 0.652327232781207,
                           se_x = 0.0883153649868509, p_x = 1.50926598070044e-13),
      overdispersion = 0.0526288721672751, convergence = 1, algorithm = "PGMM"
    ),
    subject = list(
      summary = data.frame(logFC_z_subject = 0.170415867217784,
                           se_z_subject = 0.113453753861086, p_z_subject = 0.13307815046668),
      overdispersion = 0.026113523099255, convergence = 1, algorithm = "PGMM"
    )
  )
  values[[which]]
}
