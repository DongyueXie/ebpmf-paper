# Plotting functions for flashier fits ------------------------------------------------
library(ggplot2)
plot.factors <- function(res,
                         cell.types,
                         kset = NULL,
                         rm_kset = NULL,
                         max.pt.size = 2,
                         title = NULL,
                         nonnegative = FALSE) {
  # Sort loadings according to proportion of variance explained.
  if (is.null(kset)) {
    kset <- setdiff(order(res$pve, decreasing = TRUE), which(res$pve == 0))
  }
  kset = kset[!kset%in%rm_kset]

  if (is.null(res$cell.prescaling.factors)) {
    res$cell.prescaling.factors <- 1
  }

  # Re-normalize loadings so that factors are equally spread out.
  LL <- res$L.pm[,kset,drop=F]
  LL <- LL * res$cell.prescaling.factors
  LL <- t(t(LL) / apply(abs(LL), 2, max))

  # To make it easier to compare factors, flip them to make the largest
  #   loadings positive.
  flip <- 2 * (colSums(LL > 0.75) > colSums(LL < -0.75)) - 1
  LL <- t(t(LL) * flip)

  # Make the size of the point depend on how many of that type there are.
  sizes <- max.pt.size / sqrt(table(cell.types) / min(table(cell.types)))

  if (nonnegative) {
    ylims <- c(-0.05, 1.05)
  } else {
    ylims <- c(-1.05, 1.05)
  }

  df <- reshape2::melt(LL, value.name = "loading")
  df$cell.type <- rep(as.factor(cell.types), length(kset))
  ggplot(df, aes(x = Var2, y = loading, color = cell.type)) +
    geom_jitter(position = position_jitter(0.45),
                size = rep(sizes[cell.types], length(kset))) +
    labs(title = title, x = NULL) +
    lims(y = ylims)
}
