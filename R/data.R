#' An example data set for testing scRoPE
#'
#' A dataset containing a count matrix, subject IDs, a data frame of predictors and scaling factors.
#'
#' @format A list of four objects:
#' \describe{
#'   \item{count}{A raw count matrix}
#'   \item{sid}{A vector of subject IDs}
#'   \item{pred}{A data frame of three predictors}
#'   \item{offset}{A vector of scaling factors}
#' }
#' @usage data(sample_data)
#' 
"sample_data"

#' Example Seurat object for scRoPE workflows
#'
#' A compact Seurat object for demonstrating `scToNeb()` and related helpers.
#'
#' @format A Seurat object with a count assay and per-cell metadata.
#' @usage data(sample_seurat)
#'
"sample_seurat"
