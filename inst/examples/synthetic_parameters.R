# Frozen before generating counts or inspecting fitted outputs, 2026-10-20.
# Teaching design only: no calibration experiment or paper simulation.
synthetic_example_parameters <- function() {
  list(
    design_version = "synthetic-teaching-v1",
    seed = 20261020L,
    rng_version = "4.0.0",
    rng_kind = c("Mersenne-Twister", "Inversion", "Rejection"),
    donors = 30L, cells_per_donor = 40L,
    control_donors = 20L, treatment_donors = 10L,
    genes = sprintf("G%02d", 1:8),
    beta_intercept = c(1.6, 2.0, 2.4, 1.8, 2.2, 1.6, 2.0, 2.4),
    beta_X1 = rep(0.2, 8),
    beta_X2 = c(-0.6, -0.3, 0, 0, 0, 0.3, 0.6, 0.9),
    beta_cccontrol = rep(-0.25, 8),
    phi = rep(0.25, 8), donor_variance = 0.12,
    offset_log_range = c(-0.25, 0.25),
    reference_distribution = "mean_one_gamma",
    misspecified_distribution = "mean_one_lognormal",
    target = "X2", null = 0, confidence_level = 0.95,
    model = "LN", ln_estimator = "joint_ln",
    ln_operational_fallback = "none", ncore = 1L
  )
}
