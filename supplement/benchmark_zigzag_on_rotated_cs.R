args <- commandArgs(trailingOnly = TRUE)
if (args[1] == "null") {
  target_seed <- NULL
  print("Seed specified as NULL. Will use the un-rotated target.")
} else {
  target_seed <- as.integer(args[1])
}
sampler_seed <- as.integer(args[2])
rho <- as.numeric(args[3])
d <- as.integer(args[4])
stepsize_multiplier <- as.numeric(args[5])
sampling_algorithm <- as.character(args[6])
truncated <- FALSE
  
library(hdtg)
library(coda)

# Create covariance matrix with the rotated but otherwise same covariance matrix as the compound symmetric case.
source(path.expand("supplement/helper.R"))
list2env(
  load_rotated_cs_target_params(rho, d, target_seed, truncated),
  envir=.GlobalEnv
)

# Define utility variables.
nuts_flag <- (sampling_algorithm == "nuts")
target_type <- get_rotated_cs_target_name(truncated, target_seed)

if (sampling_algorithm == "markovian") {
  sampler_name <- "markovian_zigzag"
} else {
  sampler_name <- get_zigzag_sampler_name(nuts_flag, stepsize_multiplier)
}
path_to_samples_result <- get_path_to_sampler_result(
  sampler_type = sampler_name, sampler_seed = sampler_seed,
  target_type = target_type, rho = rho, d = d, output_type = "samples"
)
path_to_summary_result <- get_path_to_sampler_result(
  sampler_type = sampler_name, sampler_seed = sampler_seed,
  target_type = target_type, rho = rho, d = d, output_type = "summary"
)

# Exit if the simulation has previously been run.
if (file.exists(path_to_samples_result) && file.exists(path_to_summary_result)) {
  print(sprintf("The simulation result '%s' exists already. Quitting.", path_to_samples_result))
  quit(status = 0)
}

# Burn-in to reach stationarity
n_burnin <- 100L
zigzag_base_step <- .1 * sqrt(max_eigval)
if (truncated) {
  x_stationary <- zigzagHMC(
    nSample = 1L, burnin = n_burnin - 1L, 
    mean = mode, prec = Prec, 
    lowerBounds = lowerBounds, upperBounds = upperBounds, 
    stepsize = zigzag_base_step, nutsFlg = T, seed = sampler_seed
  )
  x_stationary <- drop(x_stationary)
} else {
  x_stationary <- MASS::mvrnorm(mu = mode, Sigma = Sigma)
}

# Benchmark Markovian zigzag
n_hzz_sample <- 25000L
n_mzz_sample <- 250000L 
zigzag_base_step <- stepsize_multiplier * sqrt(max_eigval)
  # In case of Markovian zigzag, it's just a time interval in between sample collection

print(sprintf(
  "Starting %s zigzag benchmarking on %d dimensional %s correlation %.2g target with stepsize %.2g and seed %d.", 
  sampling_algorithm, d, target_type, rho, stepsize_multiplier, as.integer(sampler_seed)
))
invisible(gc()) # Run garbage collection before timing
start <- proc.time()

if (sampling_algorithm == "markovian") {
  samples <- hdtg:::markovianZigzag(
    nSample = n_mzz_sample, burnin = 0, 
    init = x_stationary, stepsize = zigzag_base_step, 
    mean = mode, prec = Prec, 
    lowerBounds = lowerBounds, upperBounds = upperBounds, 
    seed = sampler_seed
  )
} else {
  samples <- zigzagHMC(
    nSample = n_hzz_sample, burnin = 0, 
    init = x_stationary, stepsize = zigzag_base_step,
    mean = mode, prec = Prec, 
    lowerBounds = lowerBounds, upperBounds = upperBounds, 
    nutsFlg = nuts_flag, seed = sampler_seed
  )
}
runtime <- proc.time() - start
print(sprintf("Finished %s zigzag benchmarking.", sampling_algorithm))

## Compute ESS
ess <- effectiveSize(samples)
pc_samples <- samples %*% pc
pc_comp_samples <- samples %*% pc_complement
pc_ess <- effectiveSize(pc_samples)
pc_comp_ess <- effectiveSize(pc_comp_samples)

## Save results
print("Saving the benchmark result.")
saveRDS(samples, file = path_to_samples_result)
saveRDS(
  list(
    ess = ess, pc_ess = pc_ess, pc_comp_ess = pc_comp_ess, 
    runtime = runtime[["elapsed"]], base_integ_time = zigzag_base_step
  ), 
  file = path_to_summary_result
)
