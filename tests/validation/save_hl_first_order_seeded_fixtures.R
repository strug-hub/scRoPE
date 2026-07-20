#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
repo_root <- if (length(args)) normalizePath(args[[1L]], mustWork = TRUE) else
  normalizePath(getwd(), mustWork = TRUE)
helper <- file.path(repo_root, "tests", "testthat", "helper-hl-first-order.R")
if (!file.exists(helper)) stop("Missing fixture helper: ", helper)
source(helper, local = FALSE)

output <- file.path(repo_root, "temp", "hl_first_order_validation", "fixtures")
dir.create(output, recursive = TRUE, showWarnings = FALSE)
ids <- hl_fo_all_fixture_ids()
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
    contract_version = fixture$contract_version,
    spec_expected = paste(fixture$spec_expected, collapse = "|"),
    expected = paste(fixture$expected, collapse = "|"),
    expected_by_case = if (is.null(fixture$expected_by_case)) NA_character_ else
      paste(
        paste(names(fixture$expected_by_case), fixture$expected_by_case, sep = "="),
        collapse = "|"
      ),
    boundary_faces = if (is.null(fixture$boundary_contract$faces))
      NA_character_ else paste(fixture$boundary_contract$faces, collapse = "|"),
    selection_rule = if (is.null(fixture$selection$rule)) NA_character_ else
      fixture$selection$rule,
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
    artifact = "corrected_first_order_hl_seeded_fixture_arrays_v2",
    contract_version = hl_fo_fixture_contract_version(),
    rng = c(
      kind = "L'Ecuyer-CMRG",
      normal.kind = "Inversion",
      sample.kind = "Rejection"
    ),
    specification_fixture_ids = hl_fo_fixture_ids(),
    supplemental_fixture_ids = hl_fo_supplemental_fixture_ids(),
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
