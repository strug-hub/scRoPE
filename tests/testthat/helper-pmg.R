pmg_test_fixture <- function() {
  prepare_fn <- getFromNamespace("hl_prepare_data", "scRoPE")
  posindy_fn <- getFromNamespace("call_posindy", "scRoPE")

  count <- matrix(
    c(0, 1, 3,
      2, 0, 1,
      1, 4, 0,
      3, 1, 2),
    nrow = 1,
    byrow = TRUE
  )
  rownames(count) <- "gene1"

  id <- rep(paste0("s", 1:4), each = 3)
  x <- rep(c(-1, 0, 1), 4)
  pred <- model.matrix(~x)

  prep <- prepare_fn(
    count = count,
    id = id,
    pred = pred,
    cpc = 0,
    mincp = 1,
    verbose = FALSE
  )
  gene_index <- prep$gid[[1]]
  posv <- posindy_fn(prep$count, gene_index - 1L, prep$ctx$nind)

  list(
    count = count,
    id = id,
    pred = pred,
    prep = prep,
    gene_index = gene_index,
    posv = posv
  )
}

pmg_interior_test_fixture <- function() {
  prepare_fn <- getFromNamespace("hl_prepare_data", "scRoPE")
  posindy_fn <- getFromNamespace("call_posindy", "scRoPE")

  count <- matrix(
    c(1, 4, 1, 5,
      4, 2, 9, 12,
      0, 3, 11, 14,
      1, 1, 7, 11,
      3, 2, 11, 17,
      2, 4, 6, 12),
    nrow = 1,
    byrow = TRUE
  )
  rownames(count) <- "gene1"

  id <- rep(paste0("s", 1:6), each = 4)
  x <- rep(c(-1.5, -0.5, 0.5, 1.5), 6)
  pred <- model.matrix(~x)

  prep <- prepare_fn(
    count = count,
    id = id,
    pred = pred,
    cpc = 0,
    mincp = 1,
    verbose = FALSE
  )
  gene_index <- prep$gid[[1]]
  posv <- posindy_fn(prep$count, gene_index - 1L, prep$ctx$nind)

  list(
    count = count,
    id = id,
    pred = pred,
    prep = prep,
    gene_index = gene_index,
    posv = posv
  )
}
