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
  hzz_samples, ref_samples, coord_index, legend_loc, save_to_pdf = FALSE,
  plot_range_upper_bd = 4, ylim = c(0, 0.65)
) {
  if (save_to_pdf) {
    filename <- sprintf(
      "hamiltonian_zigzag_against_rejection_univar_marginal_comparison_along_coord_%d.pdf",
      coord_index
    )
    pdf(filename, width = 6, height = 4.5)
  }
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
  
  if (save_to_pdf) { dev.off() }
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
make_2d_hist <- function(
  samples, coord_index, color_limit = NULL, title = NULL, 
  plot_range_upper_bd = 4, save_to_pdf = FALSE, filename_prefix = NULL
) {
  
  x <- samples[, coord_index[1]]
  y <- samples[, coord_index[2]]
  within_plot_range <- (x < plot_range_upper_bd) & (y < plot_range_upper_bd)
  x <- x[within_plot_range]
  y <- y[within_plot_range]
  data <- data.frame(x = x, y = y)
  
  bivar_hist_plot <- 
    ggplot(data, aes(x = x, y = y) ) +
      geom_bin2d(bins = 40, drop = F, aes(fill = ..density..)) +
      scale_fill_viridis_c(limits = color_limit) +
      labs(
        x = bquote(italic(x[.(coord_index[1])])), 
        y = bquote(italic(x[.(coord_index[2])])),
        title = title
      ) +
      theme_minimal() +
      theme(
        aspect.ratio = 1,
        axis.text = element_text(size = 13),
        axis.title.y = element_text(size = 17, margin = margin(r = 10)),
        axis.title.x = element_text(size = 17, margin = margin(t = 10)),
        plot.title = element_text(size = 18, hjust = 0.5, vjust = 1)
      )
  
  if (save_to_pdf) {
    filename <- sprintf(
      "2d_marginal_along_coord_%d_and_%d.pdf", coord_index[1], coord_index[2]
    )
    if (!is.null(filename_prefix)) {
      filename <- paste(filename_prefix, filename, sep = "_")
    }
    ggsave(filename, width = 4, height = 4, units = "in")
  } else {
    bivar_hist_plot
  }
}

color_limit <- c(0, 0.0032)
make_2d_hist(hzz_samples, coord_index, color_limit, title = "Zigzag HMC")
make_2d_hist(ref_samples, coord_index, color_limit, title = "Reference (rejection sampler)")
