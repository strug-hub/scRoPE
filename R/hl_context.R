hl_prepare_data <- function(
    count,
    id,
    pred = NULL,
    offset = NULL,
    min = c(1e-4, 1e-4),
    max = c(10, 1000),
    opt = "lbfgs",
    verbose = TRUE,
    cpc = 0.005,
    mincp = 5,
    output_re = FALSE,
    reml = 0,
    eps = 1e-6) {

  if (is.vector(count)) {
    count <- matrix(count, nrow = 1)
  }
  count <- as(count, "dgCMatrix")
  nind <- ncol(count)
  if (nind < 2) {
    stop("The count matrix must have >= 2 columns.")
  }
  if (is.null(rownames(count))) {
    rownames(count) <- paste0("Gene", seq_len(nrow(count)))
  }
  gname <- rownames(count)
  ngene <- nrow(count)

  if (is.null(pred)) {
    pred <- matrix(1, nrow = nind, ncol = 1)
    sds <- 1
    intcol <- 1
    predn <- "(Intercept)"
  } else {
    predn <- colnames(pred)
    if ((!is.null(predn)) && length(unique(predn)) < ncol(pred)) {
      stop("All columns of the design matrix must have unique names.")
    }
    predC <- center_m(as.matrix(pred))
    sds <- predC$sds
    if ((sum(sds == 0) > 1) || any(sds < 0)) {
      stop("Some predictors have zero variation or are zero vectors.")
    }
    if (sum(sds == 0) == 0) {
      stop("Design matrix must have an intercept column.")
    }
    intcol <- which(sds == 0)
    pred <- as.matrix(predC$pred)
    if (Matrix::rankMatrix(pred) < ncol(pred)) {
      warning("Predictors are collinear or rank-deficient.")
    }
    if (nrow(pred) != nind) {
      stop("Design matrix row != count matrix column.")
    }
  }
  nb <- ncol(pred)

  if (is.null(offset)) {
    of_re <- cv_offset(0, 0, nind)
  } else {
    if (length(offset) != nind) {
      stop("Offset length != #columns of count.")
    }
    of_re <- cv_offset(as.double(offset), 1, nind)
  }
  offset_vec <- of_re$offset
  moffset <- log(of_re$mexpoffset)
  cv2 <- of_re$cv * of_re$cv

  if (length(id) != nind) {
    stop("ID length != #columns of count.")
  }
  id <- as.character(id)
  levels <- unique(id)
  id_num <- as.numeric(factor(id, levels = levels))
  if (is.unsorted(id_num)) {
    stop("Cells for the same subject must be contiguous.")
  }
  k <- length(levels)
  fid <- which(c(1, diff(id_num)) == 1)
  fid <- as.integer(c(fid, nind + 1))
  mfs <- nind / k

  cumsumy <- call_cumsumy(count, fid - 1, k, ngene)
  if (ngene == 1) {
    cumsumy <- matrix(cumsumy, ncol = k)
  }
  posind <- lapply(seq_len(ngene), function(x) which(cumsumy[x, ] > 0))

  gid2 <- which(tabulate(count@i + 1L, ngene) >= mincp)
  gidc <- which((rowSums(cumsumy) / nind) > cpc)
  gid <- intersect(gid2, gidc)
  lgid <- length(gid)
  if (verbose) {
    cat("Remove", ngene - lgid, "genes with low expression.\n")
    cat("Analyzing", lgid, "genes with", k, "subjects and", nind, "cells.\n")
  }
  if (lgid == 0) {
    stop("No gene passed the filtering.")
  }
  count_t <- Matrix::t(count)

  cell_ind_raw <- get_cell(pred, fid - 1, nb, k)
  cell_ind <- which(cell_ind_raw > 0)
  ncell <- length(cell_ind)

  ctx <- list(
    pred = pred,
    offset = offset_vec,
    fid = fid,
    cumsumy = cumsumy,
    posind = posind,
    nb = nb,
    nind = nind,
    k = k,
    intcol = intcol,
    moffset = moffset,
    min = min,
    max = max,
    eps = eps,
    opt = opt,
    cell_init = 1,
    output_re = output_re,
    id = id_num,
    cv2 = cv2,
    cell_ind = as.integer(cell_ind),
    ncell = ncell,
    mfs = mfs,
    reml = reml
  )

  list(
    ctx = ctx,
    count = count_t,
    gid = gid,
    gene_names = gname[gid],
    design_sds = sds,
    predictor_names = predn
  )
}
