# Saved results only: no generation, fitting, retries or gene selection.
# Rscript render.R RESULTS_CSV NEW_FIGURE_DIRECTORY
started <- proc.time()[["elapsed"]]
started_utc <- format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: Rscript render.R RESULTS_CSV NEW_FIGURE_DIRECTORY")
if (file.exists(args[[2L]]) || dir.exists(args[[2L]]))
  stop("Choose a new figure directory; existing outputs are never overwritten.")

r <- read.csv(args[[1L]], stringsAsFactors = FALSE, check.names = FALSE)
required <- c("gene", "scenario", "truth", "estimate", "model_lower", "model_upper",
              "robust_lower", "robust_upper", "availability", "reason")
if (!all(required %in% names(r)))
  stop("Missing result columns: ", paste(setdiff(required, names(r)), collapse = ", "))

null_suffixes <- paste0("null_", seq_len(6L))
nonnull_suffixes <- c("up", "down", "up_2", "down_2", "up_3", "down_3")
panels <- list(
  "Reference: null genes" = paste0("reference_", null_suffixes),
  "Reference: non-null genes" = paste0("reference_", nonnull_suffixes),
  "Omitted donor slopes: null genes" = paste0("omitted_", null_suffixes),
  "Omitted donor slopes: non-null genes" = paste0("omitted_", nonnull_suffixes))
gene_order <- unlist(panels, use.names = FALSE)
if (nrow(r) != length(gene_order) || anyNA(r$gene) || anyDuplicated(r$gene) ||
    !setequal(r$gene, gene_order))
  stop("Results must contain exactly the 24 prespecified genes, once each.")
r <- r[match(gene_order, r$gene), , drop = FALSE]
numeric_columns <- c("truth", "estimate", "model_lower", "model_upper",
                     "robust_lower", "robust_upper")
# read.csv represents an entirely missing column as logical. Preserve the
# all-unavailable case while still rejecting non-numeric saved results.
for (nm in numeric_columns)
  if (is.logical(r[[nm]]) && all(is.na(r[[nm]]))) r[[nm]] <- as.numeric(r[[nm]])
if (!all(vapply(r[numeric_columns], is.numeric, logical(1))))
  stop("Truth, estimates and interval endpoints must be numeric.")
if (any(!is.finite(r$truth))) stop("Every gene must retain its finite generating truth.")
if (anyNA(r$scenario) || any(r$scenario != ifelse(
    startsWith(r$gene, "reference_"), "No omitted slopes", "Omitted donor slopes")))
  stop("Scenario labels do not match the prespecified gene identities.")
for (method in c("model", "robust")) {
  lo <- r[[paste0(method, "_lower")]]
  hi <- r[[paste0(method, "_upper")]]
  if (any(is.finite(lo) & is.finite(hi) & lo > hi))
    stop("Reversed ", method, " interval endpoints.")
}
r$has_inference <- !is.na(r$availability) & r$availability == "inference"
if (any(r$has_inference & !apply(is.finite(as.matrix(r[numeric_columns])), 1L, all)))
  stop("An inference row has a non-finite estimate or interval endpoint.")
r$unavailable_label <- ""
r$unavailable_label[!r$has_inference] <- paste0("U", seq_len(sum(!r$has_inference)))

# Use one range from every saved truth, estimate and finite endpoint. No
# result-dependent gene selection or separate per-panel horizontal zooming.
all_values <- unlist(r[numeric_columns], use.names = FALSE)
data_range <- range(c(0, all_values[is.finite(all_values)]))
padding <- max(diff(data_range) * 0.065, 0.025)
x_limits <- data_range + c(-padding, padding)
x_breaks <- pretty(data_range, n = 5L)
x_breaks <- x_breaks[x_breaks >= x_limits[[1L]] & x_breaks <= x_limits[[2L]]]

suppressPackageStartupMessages(library(ggplot2))
colors <- c("Model-based" = "#27788E", "Donor-robust" = "#C25B37")
make_panel <- function(ids, title) {
  s <- r[match(ids, r$gene), , drop = FALSE]
  s$position <- rev(seq_len(nrow(s)))
  labels <- if (all(grepl("_null_", s$gene))) {
    paste("Null", sub("^.*_null_", "", s$gene))
  } else {
    paste0(rep(c("Up", "Down"), 3L), " ", rep(seq_len(3L), each = 2L),
           " (", formatC(exp(s$truth), digits = 3L, format = "fg"), "-fold)")
  }
  labels[!s$has_inference] <- paste0(labels[!s$has_inference], " [",
                                    s$unavailable_label[!s$has_inference], "]")
  d <- do.call(rbind, lapply(c("model", "robust"), function(method) {
    data.frame(estimate = s$estimate, lower = s[[paste0(method, "_lower")]],
      upper = s[[paste0(method, "_upper")]],
      y = s$position + if (method == "model") 0.12 else -0.12,
      method = if (method == "model") "Model-based" else "Donor-robust",
      has_inference = s$has_inference)
  }))
  # Non-inference genes retain their row, truth diamond and U label. Their
  # saved numerical values are not promoted to valid inference by this plot.
  d <- d[d$has_inference, , drop = FALSE]
  ggplot() +
    geom_vline(xintercept = 0, linetype = "dashed", color = "#B5BEC3") +
    geom_segment(data = d, aes(x = lower, xend = upper, y = y, yend = y,
                              color = method), linewidth = 0.9) +
    geom_point(data = d, aes(x = estimate, y = y, color = method), size = 2.2) +
    geom_point(data = s, aes(x = truth, y = position), shape = 18,
               size = 2.8, color = "#172D3A") +
    scale_color_manual(values = colors, limits = names(colors), drop = FALSE) +
    scale_x_continuous(limits = x_limits, breaks = x_breaks, expand = expansion(0)) +
    scale_y_continuous(breaks = rev(s$position), labels = rev(labels),
                       limits = c(0.55, nrow(s) + 0.45), expand = expansion(0)) +
    labs(title = title, x = "Log mean contrast and nominal 95% Wald interval", y = NULL) +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
      plot.title = element_text(face = "bold", size = 13, color = "#172D3A"),
      axis.text.y = element_text(size = 10, color = "#394A53"),
      axis.text.x = element_text(size = 9.5), axis.title.x = element_text(size = 10),
      legend.position = "none", plot.margin = margin(8, 18, 8, 12))
}
plots <- Map(make_panel, panels, names(panels))
caption <- c(
  sprintf("Diamonds mark population generating contrasts; paired dots share the same estimate. All %d genes are shown.", nrow(r)),
  "One prespecified draw per gene. Generating contrasts are not certified approximate-JLN population targets.",
  sprintf("No calibration, coverage, power or bias-correction claim. %d of %d genes provide inference.",
          sum(r$has_inference), nrow(r)))
if (any(!r$has_inference)) caption <- c(caption,
  "U labels retain unavailable genes; their availability and reasons are recorded in plotting-receipt.txt.")

draw_figure <- function() {
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(5L, 2L,
    heights = grid::unit(c(0.8, 3.7, 3.7, 0.45, 1.15), "null"))))
  on.exit(grid::popViewport())
  grid::pushViewport(grid::viewport(layout.pos.row = 1L, layout.pos.col = 1:2))
  grid::grid.text(sprintf("%d genes: one estimate per gene, two uncertainty calculations", nrow(r)),
    x = 0.025, y = 0.68, just = "left",
    gp = grid::gpar(fontsize = 18, fontface = "bold", col = "#172D3A"))
  grid::grid.text(paste0(
    "200 donors × 40 cells | Eight retained + 16 prespecified new genes | Seed 20261021\n",
    "Shared horizontal scale | Local scRoPE 0.3.0 candidate"),
    x = 0.025, y = 0.23, just = "left", gp = grid::gpar(fontsize = 10.5, lineheight = 1.15))
  grid::popViewport()
  for (i in seq_along(plots)) print(plots[[i]], newpage = FALSE,
    vp = grid::viewport(layout.pos.row = 2L + (i - 1L) %/% 2L,
                        layout.pos.col = 1L + (i - 1L) %% 2L))
  grid::pushViewport(grid::viewport(layout.pos.row = 4L, layout.pos.col = 1:2))
  for (i in seq_along(colors)) {
    x <- c(0.23, 0.43)[[i]]
    grid::grid.lines(x = c(x, x + 0.035), y = c(0.5, 0.5),
                     gp = grid::gpar(col = colors[[i]], lwd = 2))
    grid::grid.points(x = x + 0.0175, y = 0.5, pch = 16,
                     size = grid::unit(2.3, "mm"), gp = grid::gpar(col = colors[[i]]))
    grid::grid.text(names(colors)[[i]], x = x + 0.045, y = 0.5, just = "left",
                    gp = grid::gpar(fontsize = 11))
  }
  grid::grid.points(x = 0.665, y = 0.5, pch = 18,
                   size = grid::unit(2.8, "mm"), gp = grid::gpar(col = "#172D3A"))
  grid::grid.text("Generating truth", x = 0.685, y = 0.5, just = "left",
                  gp = grid::gpar(fontsize = 11))
  grid::popViewport()
  grid::pushViewport(grid::viewport(layout.pos.row = 5L, layout.pos.col = 1:2))
  grid::grid.text(paste(caption, collapse = "\n"), x = 0.025, y = 0.83,
    just = c("left", "top"), gp = grid::gpar(fontsize = 10.5, lineheight = 1.4))
  grid::popViewport()
}

if (!dir.create(args[[2L]], recursive = TRUE)) stop("Could not create figure directory.")
svg_started <- proc.time()[["elapsed"]]
svglite::svglite(file.path(args[[2L]], "twenty-four-gene-intervals.svg"),
                width = 13.4, height = 10.0, bg = "white")
tryCatch(draw_figure(), finally = grDevices::dev.off())
svg_seconds <- proc.time()[["elapsed"]] - svg_started
png_started <- proc.time()[["elapsed"]]
ragg::agg_png(file.path(args[[2L]], "twenty-four-gene-intervals.png"),
              width = 13.4, height = 10.0, units = "in", res = 150, background = "white")
tryCatch(draw_figure(), finally = grDevices::dev.off())
png_seconds <- proc.time()[["elapsed"]] - png_started
elapsed <- proc.time()[["elapsed"]] - started
receipt <- c("Saved-result plotting receipt; no generation or fitting.",
  paste("Started UTC:", started_utc),
  paste("Completed UTC:", format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")),
  paste("Input:", normalizePath(args[[1L]])),
  sprintf("Plotting elapsed seconds: %.6f", elapsed),
  "Scope: input reading, validation, library loading, and both figure renders/saves.",
  "Excludes R process startup and writing this receipt. Individual local timing; not a laptop benchmark.",
  sprintf("SVG render/save seconds: %.6f", svg_seconds),
  sprintf("PNG render/save seconds: %.6f", png_seconds),
  sprintf("Genes retained: %d; inference: %d; unavailable: %d", nrow(r),
          sum(r$has_inference), sum(!r$has_inference)),
  paste("Common horizontal limits:", paste(format(x_limits, digits = 16), collapse = ", ")),
  paste("R:", getRversion()),
  paste("Packages:", paste(vapply(c("ggplot2", "svglite", "ragg"), function(pkg)
    paste(pkg, packageVersion(pkg)), character(1)), collapse = "; ")))
if (any(!r$has_inference)) {
  u <- r[!r$has_inference, , drop = FALSE]
  receipt <- c(receipt, "Unavailable rows (none replaced or hidden):", vapply(seq_len(nrow(u)),
    function(i) paste(u$unavailable_label[[i]], u$gene[[i]],
      "availability:", u$availability[[i]], "reason:", u$reason[[i]]), character(1)))
}
writeLines(receipt, file.path(args[[2L]], "plotting-receipt.txt"))
cat(sprintf("Plotting elapsed: %.6f s (read, validation, both renders/saves; no fits).\n", elapsed))
