args <- commandArgs(trailingOnly = TRUE)
precond_type <- args[1]
sampler_seed <- as.integer(args[2])

library(hdtg)
library(coda)

# Define utility functions and variables
source(path.expand("supplement/helper.R"))
list2env(
  load_hiv_target_params(),
  envir=.GlobalEnv
)
path_to_stat_sample <- get_path_to_stationary_sample(sampler_seed)
print("Loaded target dist parameters.")


# Burn-in to reach stationarity
n_burnin <- 100L

pc_var <- 1 / min(compute_extreme_eigvals(Prec))
zigzag_base_step <- 0.1 * sqrt(pc_var)
x_stationary <- zigzagHMC(
  n = 1L, burnin = n_burnin - 1L, 
  mean = mode, prec = Prec, 
  lowerBounds = lowerBounds, upperBounds = upperBounds, 
  stepsize = zigzag_base_step, nutsFlg = T, seed = sampler_seed
)
x_stationary <- drop(x_stationary)
print("Burn-in finished.")


# Run samplers from stationarity for benchmarking
n_sample <- 1000L

## Set flags and scale parameters appropriately
if (precond_type == "none") {
  precondition <- FALSE
  
} else if (precond_type == "cond_var") {
  # Parameters of the target are scaled implicitly by hdtg
  precondition <- TRUE 
  
} else if (precond_type == "marg_var") {
  precondition <- FALSE
  marginal_sd <- sqrt(marginal_vars)
  Prec <- outer(marginal_sd, marginal_sd) * Prec 
  x_stationary <- x_stationary / marginal_sd
  mode <- mode / marginal_sd
} else {
  stop("Unknown precondition type.")
}

## Compute appropriate base integration time
if (precond_type == "cond_var") {
  eigvals <- compute_extreme_eigvals(cov2cor(Prec))
} else {
  eigvals <- compute_extreme_eigvals(Prec)
}
condnum <- max(eigvals) / min(eigvals)
pc_var <- 1 / min(eigvals)
zigzag_base_step <- 0.1 * sqrt(pc_var)

## Now benchmark
print("Starting the sampler benchmarking.")
invisible(gc()) # Run garbage collection before timing
start <- proc.time()
samples <- zigzagHMC(
  nSample = n_sample, burnin = 0, init = x_stationary,
  mean = mode, prec = Prec, 
  lowerBounds = lowerBounds, upperBounds = upperBounds, 
  nutsFlg = T, seed = sampler_seed, 
  precondition = precondition,
  stepsize = zigzag_base_step
)
runtime <- proc.time() - start
ess <- effectiveSize(samples)

## Save results
print("Saving the benchmark result.")
path_to_samples_result <- get_path_to_sampler_result(
  sampler_type = "zigzag_nuts", precond_type = precond_type, sampler_seed = sampler_seed,
  output_type = "samples"
)
saveRDS(samples, file = path_to_samples_result)

path_to_summary_result <- get_path_to_sampler_result(
  sampler_type = "zigzag_nuts", precond_type = precond_type, sampler_seed = sampler_seed,
  output_type = "summary"
)
saveRDS(
  list(samples = samples, runtime = runtime[["elapsed"]], ess = ess, 
       base_integ_time = zigzag_base_step, condnum = condnum), 
  file = path_to_summary_result
)
