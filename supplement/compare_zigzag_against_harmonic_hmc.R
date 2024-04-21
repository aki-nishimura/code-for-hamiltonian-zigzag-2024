args <- commandArgs(trailingOnly = TRUE)
sampler_type <- args[1]
sampler_seed <- as.integer(args[2])
hhmc_integ_time <- args[3] # Ignored if sampler_type == "zigzag"

library(hdtg)
library(coda)

# Define utility functions and variables
source(path.expand("supplement/helper.R"))
list2env(
  load_hiv_target_params(),
  envir=.GlobalEnv
)
target_type <- "hiv"

path_to_stat_sample <- get_path_to_stationary_sample(
  sampler_seed, target_type, rho, d
)
print("Loaded target dist parameters.")


# Burn-in, using Hamiltonian zigzag, to reach stationarity
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


# Benchmark Hamiltonian zigzag or harmonic HMC of Pakman and Paninski (2014)

if (sampler_type == "zigzag") {
  
  n_sample <- 1500L
  
  ## Specify where to save the results.
  path_to_samples_result <- get_path_to_sampler_result(
    sampler_type = "zigzag_nuts", precond_type = "none", sampler_seed = sampler_seed,
    output_type = "samples", target_type = target_type, rho = rho, d = d
  )
  path_to_summary_result <- get_path_to_sampler_result(
    sampler_type = "zigzag_nuts", precond_type = "none", sampler_seed = sampler_seed,
    output_type = "summary", target_type = target_type, rho = rho, d = d
  )

  ## Pick base stepsize
  eigvals <- compute_extreme_eigvals(Prec)
  condnum <- max(eigvals) / min(eigvals)
  pc_var <- 1 / min(eigvals)
  zigzag_base_step <- 0.1 * sqrt(pc_var)
    
  ## Start sampling
  print("Starting the sampler benchmarking.")
  invisible(gc()) # Run garbage collection before timing
  start <- proc.time()
  samples <- zigzagHMC(
    nSample = n_sample, burnin = 0, init = x_stationary,
    mean = mode, prec = Prec, 
    lowerBounds = lowerBounds, upperBounds = upperBounds, 
    nutsFlg = T, seed = sampler_seed, 
    stepsize = zigzag_base_step
  )
  runtime <- proc.time() - start
  ess <- effectiveSize(samples)
  
  ## Save results
  print("Saving the zigzag NUTS benchmark result.")
  saveRDS(samples, file = path_to_samples_result)
  saveRDS(
    list(runtime = runtime[["elapsed"]], ess = ess, 
         base_integ_time = zigzag_base_step, condnum = condnum), 
    file = path_to_summary_result
  )
  
} else {
  
  n_sample <- 100L
  
  integ_time <- pi * switch(hhmc_integ_time, 
    short = { c(1 / 8, 1 / 4) }, 
    medium = { c(3 / 16, 3 / 8) },
    long = { c(1 / 4, 1 / 2) }
  )
  
  ## Specify where to save the results.
  sampler_type <- sprintf("harmonic_hmc_%s_integ_time", hhmc_integ_time)
  path_to_hhmc_samples_result <- get_path_to_sampler_result(
    sampler_type = sampler_type, precond_type = "none", sampler_seed = sampler_seed,
    output_type = "samples", target_type = target_type, rho = rho, d = d
  )
  path_to_hhmc_summary_result <- get_path_to_sampler_result(
    sampler_type = sampler_type, precond_type = "none", sampler_seed = sampler_seed,
    output_type = "summary", target_type = target_type, rho = rho, d = d
  )
  
  ## Set input parameters
  constraintDirec <- diag(orthant_ind)[!is.na(orthant_ind), ]
  g <- rep(0, sum(!is.na(orthant_ind)))
  chol_factor <- hdtg::cholesky(Prec) # chol(Prec)
  
  ## Start the sampler
  invisible(gc()) 
  start <- proc.time()
  hhmc_output <- harmonicHMC(
    nSample = n_sample, burnin = 0, init = x_stationary,
    mean = mode, choleskyFactor = chol_factor, precFlg = T,
    F = constraintDirec, g = g,
    time = integ_time, seed = sampler_seed,
    extraOutputs = c("numBounces")
  )
  harmonic_runtime <- proc.time() - start
  n_bounces <- hhmc_output$numBounces
  harmonic_samples <- hhmc_output$samples
  harmonic_pc_samples <- drop(harmonic_samples %*% pc)
  
  
  harmonic_ess <- effectiveSize(harmonic_samples)
  harmonic_pc_ess <- effectiveSize(harmonic_pc_samples)
  
  ## Save results
  print("Saving the harmonic HMC benchmark result.")
  saveRDS(harmonic_samples, file = path_to_hhmc_samples_result)
  saveRDS(
    list(ess = harmonic_ess, pc_ess = harmonic_pc_ess, n_bounces = n_bounces,
         runtime = harmonic_runtime[["elapsed"]], integ_time = integ_time),
    file = path_to_hhmc_summary_result
  )
}
