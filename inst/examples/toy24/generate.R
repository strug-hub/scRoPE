# Base R only. No files are changed by this function.
# Prespecified extension documented in README.md. First eight genes and streams
# reproduce the earlier eight-gene prototype exactly; append streams 9 through 24.
generate_toy24 <- function() {
  old_kind <- RNGkind()
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    do.call(RNGkind, as.list(old_kind))
    if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  })
  RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection")
  set.seed(20261021)
  n <- 200L; cells <- 40L
  id <- rep(seq_len(n), each = cells)
  x <- rep(rep(c(-0.5, 0.5), each = cells / 2L), n)
  genes <- data.frame(
    gene = c("reference_null_1", "reference_null_2", "reference_up", "reference_down",
             "omitted_null_1", "omitted_null_2", "omitted_up", "omitted_down"),
    scenario = rep(c("No omitted slopes", "Omitted donor slopes"), each = 4L),
    baseline_at_zero = rep(c(1, 2, 1, 2), 2L),
    truth = rep(c(0, 0, log(1.5), -log(1.5)), 2L),
    tau_squared = rep(c(0, 1), each = 4L), stringsAsFactors = FALSE)
  # Appended parameters fixed before seeing any new-gene outcomes.
  added <- data.frame(
    gene = c(paste0("reference_null_", 3:6),
             "reference_up_2", "reference_down_2", "reference_up_3", "reference_down_3",
             paste0("omitted_null_", 3:6),
             "omitted_up_2", "omitted_down_2", "omitted_up_3", "omitted_down_3"),
    scenario = rep(c("No omitted slopes", "Omitted donor slopes"), each = 8L),
    baseline_at_zero = rep(c(1, 2, 1, 2, 2, 1, 1, 2), 2L),
    truth = rep(c(0, 0, 0, 0, log(1.5), -log(1.5), log(1.5), -log(1.5)), 2L),
    tau_squared = rep(c(0, 1), each = 8L), stringsAsFactors = FALSE)
  genes <- rbind(genes, added)
  genes$mean_minus <- genes$baseline_at_zero * exp(-genes$truth / 2)
  genes$mean_plus <- genes$baseline_at_zero * exp(genes$truth / 2)
  genes$fold_change <- exp(genes$truth)
  cell_names <- sprintf("cell_%05d", seq_along(id))
  count <- matrix(0L, nrow = nrow(genes), ncol = length(id),
                  dimnames = list(genes$gene, cell_names))
  latents <- vector("list", nrow(genes))
  streams <- vector("list", nrow(genes))
  stream <- .Random.seed
  for (g in seq_len(nrow(genes))) {
    streams[[g]] <- stream
    assign(".Random.seed", stream, envir = .GlobalEnv)
    u <- rgamma(n, shape = 1 / 0.35, rate = 1 / 0.35)
    b <- rnorm(n, sd = sqrt(genes$tau_squared[g]))
    mu <- genes$baseline_at_zero[g] * u[id] *
      exp(genes$truth[g] * x + b[id] * x - genes$tau_squared[g] * x^2 / 2)
    count[g, ] <- as.integer(rnbinom(length(id), size = 4, mu = mu))
    latents[[g]] <- data.frame(gene = genes$gene[g], donor = seq_len(n),
                               intercept_multiplier = u, donor_slope = b)
    stream <- parallel::nextRNGStream(stream)
  }
  list(count = count, id = id, pred = cbind("(Intercept)" = 1, x = x),
       offset = rep(1, length(id)), gene_metadata = genes,
       cell_metadata = data.frame(cell = cell_names, donor = id, x = x),
       donor_truth = do.call(rbind, latents), seed = 20261021L,
       rng = c("L'Ecuyer-CMRG", "Inversion", "Rejection"), streams = streams)
}
