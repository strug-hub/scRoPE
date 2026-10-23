nbresidual_beta_matrix <- function(summary, predictor_names) {
  if (!is.data.frame(summary) && !is.matrix(summary)) {
    stop("The fitted summary must be a data frame or matrix.")
  }

  nb <- length(predictor_names)
  summary_names <- colnames(summary)
  expected_names <- paste0("logFC_", predictor_names)
  if (anyDuplicated(expected_names)) {
    duplicated_names <- unique(expected_names[duplicated(expected_names)])
    stop(sprintf(
      "The residual design has duplicate coefficient names: %s.",
      paste(duplicated_names, collapse = ", ")
    ))
  }

  logfc_names <- if (is.null(summary_names)) {
    character(0)
  } else {
    summary_names[startsWith(summary_names, "logFC_")]
  }
  if (length(logfc_names)) {
    duplicated_names <- unique(logfc_names[duplicated(logfc_names)])
    if (length(duplicated_names)) {
      stop(sprintf(
        "The fitted summary has duplicate coefficient columns: %s.",
        paste(duplicated_names, collapse = ", ")
      ))
    }

    missing_names <- setdiff(expected_names, logfc_names)
    unexpected_names <- setdiff(logfc_names, expected_names)
    if (length(missing_names) || length(unexpected_names)) {
      details <- character(0)
      if (length(missing_names)) {
        details <- c(details, paste0(
          "missing: ", paste(missing_names, collapse = ", ")
        ))
      }
      if (length(unexpected_names)) {
        details <- c(details, paste0(
          "unexpected: ", paste(unexpected_names, collapse = ", ")
        ))
      }
      stop(sprintf(
        "The fitted summary coefficient columns do not match the residual design (%s).",
        paste(details, collapse = "; ")
      ))
    }
    coefficient_columns <- match(expected_names, summary_names)
  } else {
    # Compatibility with legacy fit objects whose leading coefficient block did
    # not use the logFC_<predictor> naming convention.
    if (ncol(summary) < nb) {
      stop(sprintf(
        "The legacy fitted summary has %d columns but the residual design requires %d coefficients.",
        ncol(summary), nb
      ))
    }
    coefficient_columns <- seq_len(nb)
  }

  beta <- summary[, coefficient_columns, drop = FALSE]
  numeric_columns <- vapply(seq_len(ncol(beta)), function(jj) {
    is.numeric(beta[, jj])
  }, logical(1))
  if (!all(numeric_columns)) {
    bad_names <- if (is.null(colnames(beta))) {
      as.character(which(!numeric_columns))
    } else {
      colnames(beta)[!numeric_columns]
    }
    stop(sprintf(
      "The fitted summary coefficient columns must be numeric: %s.",
      paste(bad_names, collapse = ", ")
    ))
  }
  as.matrix(beta)
}

nbresidual_hl_gene_rows <- function(summary, count) {
  gene_id <- summary[, "gene_id"]
  if (!is.numeric(gene_id) || anyNA(gene_id) ||
      any(!is.finite(gene_id) | gene_id < 1 | gene_id != floor(gene_id)) ||
      anyDuplicated(gene_id)) {
    stop("First-order HL residuals require unique positive integer gene_id values.")
  }
  rows <- gene_id
  count_names <- rownames(count)
  if (!is.null(count_names) && "gene" %in% colnames(summary)) {
    genes <- as.character(summary[, "gene"])
    if (anyNA(genes) || any(!nzchar(genes))) {
      stop("First-order HL residuals require nonmissing fitted gene names.")
    }
    original_order <- all(rows <= nrow(count)) &&
      identical(as.character(count_names[rows]), genes)
    if (!original_order) {
      if (anyDuplicated(count_names) || anyDuplicated(genes)) {
        stop("Cannot align first-order HL count rows with duplicated gene names; supply the original count row order.")
      }
      rows <- match(genes, count_names)
      if (anyNA(rows)) {
        stop("The count matrix is missing genes from the first-order HL fitted summary.")
      }
    }
  }
  if (any(rows > nrow(count))) {
    stop("First-order HL gene_id values exceed the count rows; supply original rows or matching gene names.")
  }
  rows
}

#' Extract Pearson residuals from sc-RoPE and compatible NEBULA fits
#'
#' @param nebulaSand An object of the result obtained from running the function scrope.
#' @param count A raw count matrix of the single-cell data. The rows are the genes, and the columns are the cells. The matrix can be a matrix object or a sparse dgCMatrix object. Supply the original cell order. For first-order HL fits, uniquely named gene rows may be subsetted or reordered if all fitted genes remain present; otherwise supply the original gene row order.
#' @param id A vector of subject IDs. The length should be the same as the number of columns of the count matrix.
#' @param pred The original design matrix used for fitting, with cells and donors in their original order. First-order HL coefficients use these original coordinates, including the original constant-column value. Joint-LN, PGMM and legacy objects retain their centered-public-design convention. If not specified, a unit intercept column is generated.
#' @param offset The original vector of positive scaling factors used for fitting, in the original cell order. If not specified, a vector of all ones is generated.
#' @param conditional A logical value. By default (FALSE), the function returns
#'   marginal Pearson residuals. If TRUE, the function will return conditional
#'   Pearson residuals. Conditional residuals are not available for `PGMM`
#'   fits; in that case the function warns and returns marginal residuals.
#' @details First-order HL fits are identified by their retained
#'   \code{provenance$estimator = "first_order_nested_hl"}; diagnostics need
#'   not have been retained. Keep this provenance when saving fit objects.
#'   Objects without estimator provenance use the historical centered-design
#'   convention, which is not suitable for first-order HL coefficients.
#' @return residuals: A matrix of Pearson residuals. The number of columns is the number of cells in the count matrix. The rows correspond to gene IDs reported in the result from scrope.
#' @return gene: Gene names corresponding to the row names of the count matrix.
#' @export
#' @examples
#' library(scRoPE)
#' data(sample_data)
#' pred = model.matrix(~X1+X2+cc,data=sample_data$pred)
#' counts = sample_data$count[1L, , drop = FALSE]
#' re = scrope(count=counts, id=sample_data$sid, pred=pred,
#'             model="LN", ln_estimator="joint_ln", ncore=1, verbose=FALSE)
#' resid = nbresidual(re,count=counts,id=sample_data$sid,pred=pred)
#'
nbresidual = function(nebulaSand, count, id, pred = NULL, offset = NULL, conditional = FALSE)
{
  estimator_metadata <- nebulaSand$provenance
  is_first_order_hl <- is.list(estimator_metadata) &&
    identical(estimator_metadata$estimator, "first_order_nested_hl")
  if (is_first_order_hl &&
      (!identical(estimator_metadata$route, "HL") ||
       !identical(estimator_metadata$effective_model, "HL"))) {
    stop("First-order HL residual metadata must identify the HL route and effective model.")
  }
  if (is.list(estimator_metadata) && identical(estimator_metadata$route, "HL") &&
      !is_first_order_hl) {
    stop("The HL residual estimator metadata does not declare supported prediction coordinates.")
  }
  is_joint_ln <- is.list(estimator_metadata) &&
    identical(estimator_metadata$estimator, "ln_joint_root")
  if (is_joint_ln) {
    parameterization <- estimator_metadata$parameterization
    if (!is.list(parameterization) ||
        !identical(parameterization$beta_public, "centered_public_design") ||
        !identical(parameterization$subject, "sigma2_natural") ||
        !identical(parameterization$cell_public, "Cell_equals_1_over_gamma")) {
      stop("Joint-LN residual metadata does not declare the supported prediction coordinates.")
    }
    if (isTRUE(conditional) &&
        !identical(estimator_metadata$random_effect,
                   "pml_mode_given_ln_joint_root")) {
      stop("Conditional joint-LN residuals require validated fixed-beta donor modes.")
    }
  }

  if((conditional==TRUE) & (nebulaSand$algorithm[1]=='PGMM'))
  {
    warning("The option conditional=TRUE does not support PGMM yet. Marginal Pearson residuals will be returned. Or try using NBLMM or NBGMM.")
    conditional <- FALSE
  }
  
  if((conditional==TRUE) & (is.null(nebulaSand$random_effect)))
  {stop("To obtain conditional Pearson residuals, please first run scrope with output_re=TRUE.")}
  
  ncell = ncol(count)
  ngene = nrow(nebulaSand$summary)
  
  if (is.null(pred)) {
    pred = matrix(1, ncol = 1, nrow = ncell)
    nb = 1
    predictor_names <- "(Intercept)"
  }else{
    pred = as.matrix(pred)
    nb = ncol(pred)
    predictor_names <- colnames(pred)
    if (is.null(predictor_names)) {
      predictor_names <- paste0("X", seq_len(nb))
    }
  }

  beta <- nbresidual_beta_matrix(nebulaSand$summary, predictor_names)

  if (is_first_order_hl) {
    if (nrow(pred) != ncell || !is.numeric(pred) || any(!is.finite(pred))) {
      stop("First-order HL residuals require a finite numeric design with one row per count column.")
    }
    if (length(id) != ncell || anyNA(id)) {
      stop("First-order HL residuals require one nonmissing subject ID per count column.")
    }
    # hl_fo_public_scrope_hl() already maps beta to the original input design.
    # Centering again or replacing a nonunit intercept would change the mean.
  } else {
    sds = apply(pred,2,sd)
    intercept = which(sds==0)
    if(nb>1)
    {pred[,-intercept] = scale(pred[,-intercept],scale=FALSE)}
    pred[,intercept] = 1
  }
  
  if (is.null(offset)) {
    offset = rep(1,ncell)
  }
  if (is_first_order_hl &&
      (!is.numeric(offset) || length(offset) != ncell ||
       any(!is.finite(offset) | offset <= 0))) {
    stop("First-order HL residual offsets must be finite positive scaling factors, one per cell.")
  }
  logoff = log(offset)
  
  if(conditional==TRUE){
    id = as.character(id)
    levels = unique(id)
    id = as.numeric(factor(id,levels=levels))
  }
  
  # nbres = matrix(NA,ngene,ncell)
  
  overdisp = as.matrix(nebulaSand$overdispersion)
  if (is_first_order_hl &&
      (!is.numeric(overdisp) || nrow(overdisp) != ngene || ncol(overdisp) != 2L ||
       !setequal(colnames(overdisp), c("Subject", "Cell")))) {
    stop("First-order HL overdispersion must have Subject and Cell values aligned with each fitted gene.")
  }
  if (is_first_order_hl) overdisp <- overdisp[, c("Subject", "Cell"), drop = FALSE]
  sigma2ind = overdisp[,1]
  sigma2cell = if (ncol(overdisp) >= 2) {
    overdisp[,2]
  } else {
    rep(0, ngene)
  }
  sigma2cell[which(nebulaSand$algorithm=='PGMM')] = 0
  
  if(conditional==FALSE)
  {
    ey = exp(pred%*%t(beta) + logoff + rep(sigma2ind/2,each=ncell))
    vary = ey + t(t(ey*ey)*(exp(sigma2ind)*(1+sigma2cell)-1))
  }else{
    if (is_first_order_hl &&
        (nrow(nebulaSand$random_effect) != ngene ||
         ncol(nebulaSand$random_effect) != length(unique(id)))) {
      stop("First-order HL random effects must align with the fitted genes and original donor order.")
    }
    ey = exp(pred%*%t(beta) + logoff + as.numeric(t(nebulaSand$random_effect[,id])))
    vary = ey + t(t(ey*ey)*sigma2cell)
  }
  gene_rows <- if (is_first_order_hl) {
    nbresidual_hl_gene_rows(nebulaSand$summary, count)
  } else nebulaSand$summary[, "gene_id"]
  obs = as.matrix(t(count[gene_rows, , drop = FALSE]))
  residuals = as.matrix(obs - ey)/sqrt(vary)

  residuals = t(residuals)
  rownames(residuals) = nebulaSand$summary[,'gene_id']
  colnames(residuals) = colnames(count)
  
  return(list(residuals=residuals,gene=rownames(count)[gene_rows]))
}
