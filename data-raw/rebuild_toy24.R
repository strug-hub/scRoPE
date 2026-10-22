# From the source root: Rscript data-raw/rebuild_toy24.R [OUTPUT_DIRECTORY]
# No fits, installs or changes to reference results. Default rebuilds data/.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 1L) stop("Usage: Rscript data-raw/rebuild_toy24.R [OUTPUT_DIRECTORY]")
out <- if (length(args)) args[[1L]] else "data"
if (!dir.exists(out)) stop("Output directory must already exist.")
e <- new.env()
sys.source("inst/examples/toy24/generate.R", e)
toy_data <- e$generate_toy24()
save(toy_data, file = file.path(out, "toy_data.rda"), compress = "xz", version = 2)
cat("Saved all 24 genes; fixed seed", toy_data$seed, "\n")
