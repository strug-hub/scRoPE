# Run in a fresh R session after loading the intended scRoPE installation.
# All 24 genes; one fit; no files written, seed search, retries or filtering.
library(scRoPE)
data(toy_data)
fit_time <- system.time(
  fit <- scrope(
    toy_data$count, toy_data$id,
    pred = toy_data$pred, offset = toy_data$offset,
    model = "LN", ln_estimator = "joint_ln",
    ln_operational_fallback = "none",
    ncore = 1, cpc = 0, mincp = 1, verbose = FALSE
  )
)
.toy_helpers <- new.env()
sys.source(system.file("examples", "toy24", "extract_results.R", package = "scRoPE"),
           envir = .toy_helpers)
results <- .toy_helpers$extract_toy_results(toy_data, fit)
rm(.toy_helpers)
print(results[, c("gene", "truth", "estimate", "model_se", "robust_se",
                  "model_includes_truth", "robust_includes_truth", "availability")])
print(fit_time)
