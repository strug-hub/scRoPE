# Canonical labels for adjusted-test failure reasons
failure_codes <- list(
  wald = list(
    missing_godambe = "missing_godambe",
    missing_blocks = "missing_blocks",
    singular_S = "singular_S",
    singular_G_ll = "singular_G_ll",
    singular_G_beta = "singular_G_beta",
    singular_cov = "singular_cov",
    undefined_stat = "undefined_stat",
    singular_godambe = "singular_godambe"
  ),
  lrt = list(
    negative_raw_lr = "negative_raw_lr",
    singular_Gpsi = "singular_Gpsi",
    singular_Hpsi = "singular_Hpsi",
    singular_Hlambda = "singular_Hlambda",
    nonconverged_constrained = "nonconverged_constrained",
    boundary_constrained = "boundary_constrained",
    nonpositive_curvature = "nonpositive_curvature",
    nonpositive_scale = "nonpositive_scale",
    undefined_scale = "undefined_scale",
    wald_fallback_failed = "wald_fallback_failed",
    missing_profile_lr = "missing_profile_lr",
    undefined_adjustment = "undefined_adjustment"
  ),
  score = list(
    singular_Spsi = "singular_Spsi",
    undefined_stat = "undefined_stat"
  )
)

# Helper to normalise failure reasons to the known vocabulary
normalise_failure <- function(reason, type = c("wald", "lrt", "score")) {
  type <- match.arg(type)
  if (is.null(reason)) {
    return(NA_character_)
  }

  reason_vec <- as.character(reason)
  original_na <- is.na(reason)
  reason_vec[original_na] <- NA_character_
  out <- rep(NA_character_, length(reason_vec))
  missing_idx <- is.na(reason_vec) | (!original_na & !nzchar(reason_vec))
  if (all(missing_idx)) {
    return(if (length(out) == 1) out[[1]] else out)
  }

  known_idx <- reason_vec %in% failure_codes[[type]]
  out[known_idx] <- reason_vec[known_idx]
  other_idx <- !(known_idx | missing_idx)
  out[other_idx] <- paste0("other_", reason_vec[other_idx])

  if (length(out) == 1) {
    return(out[[1]])
  }
  out
}
