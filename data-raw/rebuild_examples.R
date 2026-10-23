# From the package source root: Rscript data-raw/rebuild_examples.R
# Requires existing Matrix and Seurat installations; installs nothing.
source("inst/examples/synthetic_parameters.R")
source("inst/examples/generate_synthetic.R")
examples <- generate_synthetic_examples()
sample_data <- examples$reference
sample_seurat <- make_synthetic_seurat(sample_data)
save(sample_data, file = "data/sample_data.rda", compress = "xz", version = 2)
save(sample_seurat, file = "data/sample_seurat.rda", compress = "xz", version = 2)
write.csv(examples$truth, "inst/examples/generating_coefficients.csv", row.names = FALSE)
cat("Regenerated synthetic teaching examples; seed", examples$parameters$seed, "\n")
cat("R:", R.version.string, "\nMatrix:", as.character(packageVersion("Matrix")),
    "\nSeurat:", as.character(packageVersion("Seurat")),
    "\nSeuratObject:", as.character(packageVersion("SeuratObject")), "\n")
