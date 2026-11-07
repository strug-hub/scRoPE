library(testthat)
library(scRoPE)

if (identical(Sys.getenv("SCROPE_HL_HISTORY_DIAGNOSTICS"), "true")) {
  source(file.path("testthat", "helper-hl-history-diagnostic.R"))
  restore <- .hl_history_diagnostic_install(file.path("testthat", "_diagnostics"))
  tryCatch(test_check("scRoPE"), finally = restore())
} else {
  test_check("scRoPE")
}
