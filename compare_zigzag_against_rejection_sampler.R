library(hdtg)
library(MASS)
library(ggplot2)
set.seed(1918)


# Define a truncated Gaussian target
d <- 16
mode <- rep(.25, d) + .1 * rnorm(d) 
Sigma_mean <- diag(.5, d) + .5
df <- 2 * d
Sigma <- drop(rWishart(1, df, Sigma_mean / df))
Prec <- solve(Sigma)
lower_bd <- rep(0, d)
upper_bd <- rep(Inf, d)

# Run Hamiltonian zigzag
n_sample <- 10^6
burnin <- 100
x0 <- rep(1, d)
useNuts <- TRUE

hzz_samples <- hdtg::zigzagHMC(
  n_sample, burnin, mean = mode, prec = Prec,
  lowerBounds = lower_bd, upperBounds = upper_bd,
  init = x0, nutsFlg = useNuts, seed = 615
)

# Generate ground-truth with rejection sampling
n_ref_samples <- 10^7 
ref_samples <- mvrnorm(n_ref_samples, mode, Sigma) 
within_bdry <- sapply(
  1:n_ref_samples, 
  function(i) all(ref_samples[i, ] > lower_bd)
)
ref_samples <- ref_samples[within_bdry, ]


# Compare the two sets of samples in terms of their one-dimensional marginals
coord_index <- sample.int(d, size = 2)

make_comparison_hist <- function(
  hzz_samples, ref_samples, coord_index, legend_loc,
  plot_range_upper_bd = 4, ylim = c(0, 0.65)
) {
  hzz_marg_samples <- hzz_samples[, coord_index]
  ref_marg_samples <- ref_samples[, coord_index]
  
  breaks <- seq(0, plot_range_upper_bd, length.out = 41)
  colors <- list(
    hzz = rgb(0, 128, 255, 255 * .4, maxColorValue = 255), # rgb(1, 0, 0, 1/4),
    ref = rgb(255, 127, 0, 255 * .25, maxColorValue = 255) # rgb(0, 0, 1, 1/4)
  )
  hzz_marg_samples <- hzz_marg_samples[hzz_marg_samples < plot_range_upper_bd]
  ref_marg_samples <- ref_marg_samples[ref_marg_samples < plot_range_upper_bd]

  hist(
    hzz_marg_samples,
    col=colors$hzz,
    breaks=breaks, probability=T,
    xlim=c(0, plot_range_upper_bd),
    main = "",
    xlab = "",
    cex.axis = 1.1,
    cex.lab = 1.2
  )
  mtext(bquote(italic(x[.(coord_index)])), side=1, line=3, cex=1.4)
  
  hist(
    ref_marg_samples, 
    col=colors$ref,
    breaks=breaks, probability=T, 
    xlim=c(0, plot_range_upper_bd), 
    add=T
  )
  
  if (!is.null(legend_loc)) {
    legend(
      legend_loc[1], legend_loc[2],
      legend = c("Zigzag-NUTS", "Reference\n(rejection sampler)"), 
      cex = 1.3, fill = unlist(colors), bty = "n"
    )
  }
}

make_comparison_hist(
  hzz_samples, 
  ref_samples,
  coord_index = coord_index[1],
  legend_loc = c(1.7, 0.6)
)

make_comparison_hist(
  hzz_samples, 
  ref_samples,
  coord_index = coord_index[2],
  legend_loc = NULL
)

# Compare 2-dimensional marginals

# Source: https://stackoverflow.com/questions/28562288/how-to-use-the-hsl-hue-saturation-lightness-cylindric-color-model
hsl_to_rgb <- function(h, s, l, alpha) {
  # h = 0 - 360 (whole input degrees)
  # s = 0.0 - 1 (0 - 100%)
  # l = 0.0 - 1, (0 - 100%)
  h <- h / 360
  r <- g <- b <- 0.0
  if (s == 0) {
    r <- g <- b <- l
  } else {
    hue_to_rgb <- function(p, q, t) {
      if (t < 0) { t <- t + 1.0 }
      if (t > 1) { t <- t - 1.0 }
      if (t < 1/6) { return(p + (q - p) * 6.0 * t) }
      if (t < 1/2) { return(q) }
      if (t < 2/3) { return(p + ((q - p) * ((2/3) - t) * 6)) }
      return(p)
    }
    q <- ifelse(l < 0.5, l * (1.0 + s), l + s - (l*s))
    p <- 2.0 * l - q
    r <- hue_to_rgb(p, q, h + 1/3)
    g <- hue_to_rgb(p, q, h)
    b <- hue_to_rgb(p, q, h - 1/3)
  }
  return(rgb(r, g, b))
}
make_2d_hist <- function(samples, coord_index, hue, plot_range_upper_bd=4) {
  
  x <- samples[, coord_index[1]]
  y <- samples[, coord_index[2]]
  within_plot_range <- (x < plot_range_upper_bd) & (y < plot_range_upper_bd)
  x <- x[within_plot_range]
  y <- y[within_plot_range]
  data <- data.frame(x = x, y = y)
  
  ggplot(data, aes(x = x, y = y) ) +
    geom_bin2d(bins = 40, drop = F, aes(fill = ..density..)) +
    scale_fill_gradient(
      low = hsl_to_rgb(hue, .9, .05), 
      high = hsl_to_rgb(hue, .9, .7),
      guide = "none" 
    ) +
    labs(
      x = bquote(italic(x[.(coord_index[1])])), 
      y = bquote(italic(x[.(coord_index[2])]))
    ) +
    theme_minimal() +
    theme(
      aspect.ratio = 1,
      axis.text = element_text(size = 13),
      axis.title.y = element_text(size = 17, margin = margin(r = 10)),
      axis.title.x = element_text(size = 17, margin = margin(t = 10))
    )
  
}

make_2d_hist(hzz_samples, coord_index, hue = 210)
make_2d_hist(ref_samples, coord_index, hue = 30)
