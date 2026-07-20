#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
repo_root <- if (length(args)) normalizePath(args[[1L]], mustWork = TRUE) else
  normalizePath(getwd(), mustWork = TRUE)
helper <- file.path(repo_root, "tests", "testthat", "helper-hl-first-order.R")
if (!file.exists(helper)) stop("Missing fixture helper: ", helper)
source(helper, local = FALSE)

output <- file.path(repo_root, "temp", "hl_first_order_validation", "fixtures")
dir.create(output, recursive = TRUE, showWarnings = FALSE)
ids <- hl_fo_fixture_ids()
rows <- vector("list", length(ids))
for (ii in seq_along(ids)) {
  id <- ids[[ii]]
  fixture <- hl_fo_fixture(id, fresh = TRUE)
  path <- file.path(output, paste0(id, ".rds"))
  saveRDS(fixture, path, version = 3L)
  rows[[ii]] <- data.frame(
    id = id,
    file = basename(path),
    md5 = unname(tools::md5sum(path)),
    expected = paste(fixture$expected, collapse = "|"),
    requested_seed = if (!is.null(fixture[["requested_seed"]]))
      fixture[["requested_seed"]] else if (!is.null(fixture[["seed"]]))
        fixture[["seed"]] else paste(fixture[["seeds"]], collapse = "|"),
    actual_seed = if (!is.null(fixture[["actual_seed"]]))
      fixture[["actual_seed"]] else NA_character_,
    stringsAsFactors = FALSE
  )
}
manifest <- do.call(rbind, rows)
write.csv(manifest, file.path(output, "manifest.csv"), row.names = FALSE)
saveRDS(
  list(
    artifact = "corrected_first_order_hl_seeded_fixture_arrays",
    rng = c(
      kind = "L'Ecuyer-CMRG",
      normal.kind = "Inversion",
      sample.kind = "Rejection"
    ),
    fixture_ids = ids,
    manifest = manifest,
    package_version = unname(read.dcf(
      file.path(repo_root, "DESCRIPTION"), fields = "Version"
    )[[1L]])
  ),
  file.path(output, "manifest.rds"),
  version = 3L
)
print(manifest, row.names = FALSE)
cat("Saved fixture arrays to", output, "\n")
