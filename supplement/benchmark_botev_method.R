library(TruncatedNormal)

warning(paste0(
  "The computation within this script will take at least days to complete and ", 
  "it has *not* been confirmed how much time it actually takes."
))


source(path.expand("supplement/helper.R"))
list2env(
  load_hiv_target_params(),
  envir=.GlobalEnv
)
print("Loaded target dist parameters.")
  
print(paste0(
  "Inverting precision matrix via cholesky to get covariance; this will take some ",
  "time, but should be in the order of minutes and not hours on typical machines."
))
tic <- proc.time()
sigma <- chol2inv(chol(Prec))
toc <- proc.time()
print(sprintf(
  "Inversion via cholesky took %.3g sec.", 
  toc['elapsed'] - tic['elapsed']
))

lower_bd <- rep(0, nrow(sigma))
lower_bd[orthant_ind == -1] <- -Inf
lower_bd[is.na(orthant_ind)] <- -Inf
upper_bd <- rep(0, nrow(sigma))
upper_bd[orthant_ind == 1] <- Inf
upper_bd[is.na(orthant_ind)] <- Inf

n_samples <- 1
set.seed(615)
print("Starting the sampler of Botev (2017).")
tic <- proc.time()
botev_samples <- rtmvnorm(
  n = n_samples, mu = mode, sigma = sigma, lb = lower_bd, ub = upper_bd
)
toc <- proc.time()
print(sprintf("The sampling took %.3g sec.", toc['elapsed'] - tic['elapsed'])) 