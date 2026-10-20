# From package source root: Rscript data-raw/render_readme.R
# Uses the complete fixed-run CSV. No fitting, filtering or gene selection.
library(ggplot2)
library(grid)
d <- read.csv("inst/examples/illustration_results.csv", stringsAsFactors = FALSE)
stopifnot(nrow(d) == 16L, identical(unique(d$gene), sprintf("G%02d", 1:8)))
scenario_names <- c(reference = "Gamma donor effects | reference",
                    misspecified = "Lognormal donor effects | misspecified")
d$scenario <- factor(d$scenario, levels = names(scenario_names), labels = scenario_names)
d$y <- 9 - match(d$gene, sprintf("G%02d", 1:8))
cols <- c("Model" = "#327A9A", "Donor-robust" = "#C55736")
base_theme <- theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 19, color = "#182C3A"),
        plot.subtitle = element_text(size = 12, margin = margin(b = 14)),
        strip.text = element_text(face = "bold", size = 12),
        legend.position = "bottom", legend.title = element_blank(),
        plot.caption = element_text(hjust = 0, lineheight = 1.2, color = "#46525C"),
        plot.margin = margin(18, 22, 16, 16))
save_figure <- function(p, name, width, height) {
  ggsave(file.path("man/figures", paste0(name, ".svg")), p,
         width = width, height = height, device = svglite::svglite, bg = "white")
  ggsave(file.path("man/figures", paste0(name, ".png")), p,
         width = width, height = height, dpi = 150, device = ragg::agg_png, bg = "white")
}
intervals <- rbind(
  transform(d, method = "Model", lower = model_lower, upper = model_upper,
            position = y + 0.13),
  transform(d, method = "Donor-robust", lower = robust_lower, upper = robust_upper,
            position = y - 0.13))
intervals$method <- factor(intervals$method, levels = names(cols))
limits <- range(c(intervals$lower, intervals$upper, d$truth), finite = TRUE) + c(-0.08, 0.08)
missing <- subset(d, availability != "inference")
available_count <- sum(d$availability == "inference")
availability_caption <- if (available_count == nrow(d)) {
  sprintf("all %d fits provided inference.", nrow(d))
} else {
  sprintf("%d of %d fits provided inference; unavailable fits are labeled.",
          available_count, nrow(d))
}
p <- ggplot() +
  geom_vline(xintercept = 0, color = "#CBD2D8", linetype = 2) +
  geom_segment(data = intervals, aes(x = lower, xend = upper, y = position,
    yend = position, color = method), linewidth = 1.0, na.rm = TRUE) +
  geom_point(data = subset(intervals, availability == "inference"),
             aes(x = estimate, y = position, color = method),
             size = 2.3, na.rm = TRUE) +
  geom_point(data = d, aes(x = truth, y = y, shape = "Generating coefficient"),
             size = 3.1, stroke = 1.0, color = "#182C3A") +
  geom_text(data = missing, aes(x = limits[2], y = y, label = "Unavailable"),
            hjust = 1, size = 3, color = "#7C2831") +
  facet_wrap(~scenario, nrow = 1) +
  scale_color_manual(values = cols) + scale_shape_manual(values = c("Generating coefficient" = 4)) +
  scale_x_continuous(limits = limits, breaks = seq(-1, 1.5, 0.5)) +
  scale_y_continuous(breaks = 1:8, labels = rev(sprintf("G%02d", 1:8))) +
  labs(title = "Same fitted effect, different uncertainty",
       subtitle = "Treatment coefficient X2: paired examples, all eight genes, 30 independent donors",
       x = "Log coefficient and nominal 95% Wald interval", y = NULL,
       caption = paste("Dots share one estimate per gene; crosses mark generating coefficients. No bias correction is applied.",
                       paste("One realization per scenario: interval inclusion is not a coverage estimate;", availability_caption), sep = "\n")) +
  base_theme + theme(panel.grid.major.y = element_blank())
save_figure(p, "synthetic-intervals", 12, 6.6)
q <- ggplot(d, aes(x = gene, y = se_ratio)) +
  geom_hline(yintercept = 1, color = "#7B8791", linetype = 2) +
  geom_segment(aes(xend = gene, y = 1, yend = se_ratio), color = "#C55736", linewidth = 0.8) +
  geom_point(color = "#C55736", size = 2.8, na.rm = TRUE) +
  geom_text(data = missing, aes(x = gene, y = 1, label = "Unavailable"),
            angle = 90, vjust = -0.5, size = 3, color = "#7C2831") +
  facet_wrap(~scenario, nrow = 1) +
  scale_y_continuous(limits = c(min(0.8, min(d$se_ratio, na.rm = TRUE) - 0.04),
                               max(1.32, max(d$se_ratio, na.rm = TRUE) + 0.04)),
                     breaks = seq(0.8, 1.4, 0.1)) +
  labs(title = "The adjustment can widen or narrow an interval",
       subtitle = "Robust SE / model SE; 1 means equal standard errors",
       x = NULL, y = "SE ratio",
       caption = paste("The distributions differ modestly, and their observed adjustments are similar.",
                       "Neither a larger ratio nor this single draw establishes better calibration or power.", sep = "\n")) + base_theme
save_figure(q, "synthetic-se-ratio", 12, 4.2)

# Design schematic: counts and labels, not fitted or selected outcomes.
draw_design <- function() {
  grid.newpage()
  grid.rect(gp = gpar(fill = "white", col = NA))
  grid.text("A transparent, paired teaching example", .04, .93, just = "left",
            gp = gpar(fontsize = 24, fontface = "bold", col = "#182C3A"))
  grid.text("Seed 20261020  |  8 genes  |  30 donors x 40 cells  |  no seed or gene selection",
            .04, .855, just = "left", gp = gpar(fontsize = 13, col = "#46525C"))
  box <- function(x, y, w, h, title, lines, fill) {
    grid.roundrect(x, y, w, h, r = unit(0.015, "npc"),
                   gp = gpar(fill = fill, col = "#CAD4DC"))
    grid.text(title, x, y + h * .26, gp = gpar(fontface = "bold", fontsize = 15, col = "#182C3A"))
    grid.text(lines, x, y - h * .12, gp = gpar(fontsize = 12, lineheight = 1.4, col = "#263E4D"))
  }
  box(.255, .67, .43, .25, "Shared study design",
      "20 control + 10 treated donors\nX2: donor treatment; X1 and cc: cell predictors\nPositive offsets; conditional NB dispersion = 0.25", "#EFF5F8")
  box(.745, .67, .43, .25, "Known generating effects",
      "X2 coefficients, G01 through G08\n-0.6, -0.3, 0, 0, 0, 0.3, 0.6, 0.9\nIndependent donors; artificial category labels", "#EFF5F8")
  grid.lines(c(.5, .5), c(.53, .48), arrow = arrow(length = unit(.10, "inches")),
             gp = gpar(col = "#647581", lwd = 1.2))
  box(.255, .34, .43, .24, "Gamma reference",
      "Mean-one gamma donor multiplier\nVariance = 0.12\nRepresentable by the JLN working family", "#E8F1F5")
  box(.745, .34, .43, .24, "Lognormal misspecification",
      "Mean-one lognormal donor multiplier\nSame variance = 0.12; different shape\nSame distribution in both treatment groups", "#FBEEE8")
  grid.text("Paired random draws and identical covariates/offsets; only donor-distribution shape changes.",
            .5, .145, gp = gpar(fontsize = 12, col = "#263E4D"))
  grid.text("Compare model and donor-robust uncertainty at the same JLN estimate. This is not a calibration study.",
            .5, .075, gp = gpar(fontsize = 12, fontface = "bold", col = "#182C3A"))
}
svglite::svglite("man/figures/synthetic-design.svg", width = 12, height = 5.6)
draw_design(); dev.off()
ragg::agg_png("man/figures/synthetic-design.png", width = 12, height = 5.6, units = "in", res = 150)
draw_design(); dev.off()
