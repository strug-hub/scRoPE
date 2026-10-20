#' An example data set for testing scRoPE
#'
#' A dataset containing 10 genes and 6,176 cells from 30 subject identifiers,
#' with predictors X1, X2, and cc and scaling factors.
#'
#' @source Bundled unchanged from scRoPE's initial import and byte-identical to
#' the example in the preserved NEBULA 1.5.3 source archive.
#' \url{https://github.com/lhe17/nebula}
#' @details The inspected upstream documentation does not record the original
#' generation procedure or study source. Archive agreement does not uniquely
#' identify the original imported NEBULA version.
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
#' A Seurat object containing 1,000 genes and 1,000 cells from the panc8
#' pancreas example, for demonstrating `scToNeb()` and related helpers.
#'
#' @source Bundled unchanged from scRoPE's initial import and byte-identical to
#' NEBULA 1.5.3's example. NEBULA attributes it to Satija Lab's
#' \emph{panc8.SeuratData: Eight Pancreas Datasets Across Five Technologies}
#' (2019; cited package version 3.0.2), distributed through
#' \url{https://github.com/satijalab/seurat-data}.
#' @details The exact subset-selection procedure and dataset-specific
#' redistribution record have not been recovered.
#'
#' @format A Seurat object with a count assay and per-cell metadata.
#' @usage data(sample_seurat)
#'
"sample_seurat"
