#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) args[1L] else "dev/covariance-screen/results"
data <- read.csv(file.path(root, "pooled-performance.csv"))
draw <- function() {
  old <- par(mfrow = c(2, 3), mar = c(4.4, 4, 2.7, .7), oma = c(.5, .5, 2.5, .5),
    mgp = c(2.5, .7, 0), cex.lab = .95, cex.axis = .9)
  on.exit(par(old))
  colors <- c(known = "#245F8F", fitted = "#B14E22")
  for (beta in c(0, .5)) for (rho in c(-.5, 0, .5)) {
    plot(NA, xlim = c(-.57, .57), ylim = c(0, 1), xlab = "Assumed cross-error correlation",
      ylab = "95% interval coverage", xaxt = "n",
      main = sprintf("True rho = %s; beta = %s", rho, beta))
    axis(1, at = c(-.5, 0, .5))
    abline(h = .95, col = "gray40", lty = 2)
    abline(v = rho, col = "gray75", lty = 3)
    for (mode in c("known", "fitted")) {
      d <- data[data$mode == mode & data$policy == "minimum" & data$rho == rho & data$beta == beta, ]
      d <- d[order(d$assumed_rho), ]
      x <- d$assumed_rho + if (mode == "known") -.015 else .015
      lines(x, d$coverage_usable, col = colors[mode], lty = if (mode == "known") 1 else 2)
      arrows(x, d$coverage_mc_lower, x, d$coverage_mc_upper, angle = 90, code = 3,
        length = .035, col = colors[mode])
      points(x, d$coverage_usable, col = colors[mode], pch = if (mode == "known") 16 else 17)
    }
    if (beta == 0 && rho == -.5) legend("bottomleft", c("Known marginal parameters", "Fitted marginal models"),
      col = colors, pch = c(16, 17), lty = c(1, 2), bty = "n", cex = .73)
  }
  mtext("Cross-model covariance assumptions: basic simulation screen", outer = TRUE, cex = 1.15)
}
png(file.path(root, "coverage.png"), width = 1800, height = 1200, res = 180)
draw(); dev.off()
pdf(file.path(root, "coverage.pdf"), width = 10, height = 6.7, useDingbats = FALSE)
draw(); dev.off()
