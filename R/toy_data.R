#' Synthetic omitted-slope quick-start data
#'
#' Twenty-four genes and 8,000 cells from 200 independent simulated donors,
#' with 20 cells in each of two conditions per donor. These are synthetic
#' teaching data, separate from paper simulations and biological analyses.
#'
#' @source Generated with fixed seed 20261021 by the installed script
#' \code{system.file("examples", "toy24", "generate.R", package="scRoPE")}.
#' The accompanying README records parameters, draw order and design history.
#' @details Each of four categories has six genes: reference null, reference
#' non-null, omitted-slope null and omitted-slope non-null. Conditional counts
#' are negative binomial with size 4 and mean
#' \eqn{a U \exp(\beta x + bx - \tau^2 x^2/2)}. Here U is mean-one Gamma
#' with variance 0.35, b is independent Normal with variance \eqn{\tau^2},
#' and x is -0.5 or 0.5. Reference genes have \eqn{\tau^2=0}; omitted-slope
#' genes have \eqn{\tau^2=1}. Baseline a is 1 or 2. The population mean is
#' \eqn{a\exp(\beta x)}, so beta is the generating population log-mean
#' contrast: zero or plus/minus log(1.5).
#'
#' JLN retains a donor random intercept and a common x coefficient, omitting
#' donor-specific slopes. Both uncertainty estimates share that coefficient.
#' Equality of the approximate JLN population target and the non-null generating
#' contrast under omission has not been established. Interval inclusion in this
#' draw is not estimated coverage; sandwich variance does not correct bias.
#' Eight genes were retained from an earlier prototype and 16 additional genes
#' were prespecified after viewing those eight. No seed search or gene selection
#' was used. RNG settings and seed are restored by the generator, except that
#' base R does not expose a caller's hidden Box-Muller cache.
#'
#' @format A list containing:
#' \describe{
#'   \item{count}{A 24 by 8,000 integer count matrix, genes by cells.}
#'   \item{id}{Integer donor IDs 1--200, in contiguous blocks of 40.}
#'   \item{pred}{An 8,000 by 2 matrix: intercept and condition x.}
#'   \item{offset}{8,000 unit multiplicative cell offsets.}
#'   \item{gene_metadata}{Gene names, scenarios, baseline means, generating
#'     log contrasts (truth), slope variances, condition means and fold changes.}
#'   \item{cell_metadata}{Cell names, donor IDs and x in count-column order.}
#'   \item{donor_truth}{All simulated donor multipliers and slopes for each gene.}
#'   \item{seed}{The fixed date-based seed 20261021.}
#'   \item{rng}{L'Ecuyer-CMRG, Inversion and Rejection.}
#'   \item{streams}{Initial RNG state for each gene's independent stream.}
#' }
#' @usage data(toy_data)
#' @examples
#' data(toy_data)
#' dim(toy_data$count)
#' head(toy_data$gene_metadata)
#' # See the installed examples/toy24/quick_start.R for the complete JLN fit.
"toy_data"
