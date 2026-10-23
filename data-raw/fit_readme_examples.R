# From the source root with the intended scRoPE library selected:
# Rscript data-raw/fit_readme_examples.R [optional-private-fit-receipt.rds]
# Exactly two fits, eight genes each. No seed/parameter search or filtering.
library(scRoPE)
source("inst/examples/synthetic_parameters.R")
source("inst/examples/generate_synthetic.R")
x <- generate_synthetic_examples()
p <- x$parameters
fits <- lapply(c("reference", "misspecified"), function(scenario) {
  d <- x[[scenario]]
  scrope(d$count, d$sid, pred = model.matrix(~ X1 + X2 + cc, d$pred),
         offset = d$offset, model = p$model, ln_estimator = p$ln_estimator,
         ln_operational_fallback = p$ln_operational_fallback, ncore = p$ncore,
         cpc = 0, mincp = 1, keep_diagnostics = TRUE)
})
names(fits) <- c("reference", "misspecified")
results <- do.call(rbind, lapply(names(fits), function(scenario) {
  f <- fits[[scenario]]
  s <- f$summary[match(p$genes, f$summary$gene), , drop = FALSE]
  a <- f$availability[match(p$genes, f$availability$gene), , drop = FALSE]
  stopifnot(identical(s$gene, p$genes), nrow(a) == length(p$genes))
  ans <- data.frame(scenario = scenario, gene = p$genes, truth = p$beta_X2,
    estimate = s$logFC_X2, se_model = s$se_X2, se_robust = s$se_robust_X2,
    p_model = s$p_naive_X2, p_robust = s$p_robust_X2,
    availability = a$availability, stringsAsFactors = FALSE)
  ans$se_ratio <- ans$se_robust / ans$se_model
  ans$estimation_error <- ans$estimate - ans$truth
  z <- qnorm((1 + p$confidence_level) / 2)
  ans$model_lower <- ans$estimate - z * ans$se_model
  ans$model_upper <- ans$estimate + z * ans$se_model
  ans$robust_lower <- ans$estimate - z * ans$se_robust
  ans$robust_upper <- ans$estimate + z * ans$se_robust
  ans
}))
rownames(results) <- NULL
write.csv(results, "inst/examples/illustration_results.csv", row.names = FALSE,
          na = "NA")
args <- commandArgs(trailingOnly = TRUE)
if (length(args)) saveRDS(list(parameters = p, fits = fits, results = results,
  session = sessionInfo(), loaded_package = find.package("scRoPE")), args[[1]])
print(results)
print(sessionInfo())
