data_folder <- path.expand("supplement/data/")
output_folder <- path.expand("supplement/results/")
subsample_data <- FALSE
subsample_size <- 1024L
subsampling_seed <- 615

load_hiv_target_params <- function() {
  # Load parameters of truncated Gaussians
  Prec <- scan(paste0(data_folder, "precision_matrix.csv"), sep=',')
  Prec <- matrix(Prec, ncol = sqrt(length(Prec)), byrow = TRUE)
  orthant_ind <- scan(paste0(data_folder, "orthant_indicator.txt"))
  mode <- scan(paste0(data_folder, "mean.txt"))
  marginal_vars <- scan(paste0(data_folder, "marginal_variances.txt"))
  pc <- scan(paste0(data_folder, "principal_component.txt"))
  
  ## Subsample data as warranted for testing purpose
  if (subsample_data) {
    set.seed(subsampling_seed)
    subset_ind <- sample.int(dim(Prec)[1], size = subsample_size)
    Prec <- Prec[subset_ind, subset_ind]
    orthant_ind <- orthant_ind[subset_ind]
    mode <- mode[subset_ind]
    marginal_vars <- marginal_vars[subset_ind]
    pc <- pc[subset_ind]
  }
  
  ## Calculate lower and upper bounds corresponding to the orthant indicator
  lowerBounds <- ifelse(orthant_ind == 1, 0, -Inf)
  lowerBounds[is.na(lowerBounds)] <- -Inf
  upperBounds <- ifelse(orthant_ind == -1, 0, Inf)
  upperBounds[is.na(upperBounds)] <- Inf
  
  return(list(
    Prec = Prec, orthant_ind = orthant_ind, mode = mode, marginal_vars = marginal_vars,
    pc = pc, lowerBounds = lowerBounds, upperBounds = upperBounds,
    subsample_data = subsample_data, subsample_size = subsample_size
  ))
}

load_cs_target_params <- function(rho, d) {
  Sigma <- (1 - rho) * diag(d) + rho * array(1., dim = c(d, d))
  Prec <- 1 / (1 - rho) * (
    diag(d) - rho / (1 - rho + rho * d) * array(1., dim = c(d, d))
  ) 
  mode <- rep(0, d)
  lowerBounds <- rep(0, d)
  upperBounds <- rep(Inf, d)
  orthant_ind <- rep(1, d)
  pc <- rep(1, d) / sqrt(d)
  return(list(
    Prec = Prec, orthant_ind = orthant_ind, mode = mode,
    pc = pc, lowerBounds = lowerBounds, upperBounds = upperBounds
  ))
}

load_rotated_cs_target_params <- function(rho, d, seed, truncated = FALSE, pathological = FALSE) {
  if (pathological) {
    # "Pathologically" easy case of independent coordinates
    pc <- c(1, rep(0, d - 1))
    pc_complement <- diag(d)[, 2:d] 
  } else {
    if (is.null(seed)) {
      pc <- rep(1, d) / sqrt(d)
    } else {
      set.seed(seed)
      pc <- rnorm(d)
      pc <- pc / sqrt(sum(pc^2))
    }
    pc_complement <- diag(d)[, 2:d] 
    pc_complement <- pc_complement - pc %*% (t(pc) %*% pc_complement)
    pc_complement <- qr.Q(qr(pc_complement))
  }
  max_eigval <- (1 - rho) + rho * d
  min_eigval <- (1 - rho) 
  Sigma <- (max_eigval - min_eigval) * outer(pc, pc) + min_eigval * diag(d)
  Prec <- 1 / min_eigval * diag(d) + (1 / max_eigval - 1 / min_eigval) * outer(pc, pc)
  mode <- rep(0, d)
  if (truncated) {
    lowerBounds <- ifelse(sign(pc) <= 0, -Inf, 0)
    upperBounds <- ifelse(sign(pc) <= 0, 0, Inf)
  } else {
    lowerBounds <- rep(-Inf, d)
    upperBounds <- rep(Inf, d)
  }
  return(list(
    Sigma = Sigma, Prec = Prec, mode = mode, 
    lowerBounds = lowerBounds, upperBounds = upperBounds,
    pc = pc, pc_complement = pc_complement, 
    max_eigval = max_eigval, min_eigval = min_eigval
  ))
}

get_path_to_stationary_sample <- function(sampler_seed, target_type = "hiv", rho = 0, d = NA) {
  filename <- "stationary_sample"
  if (target_type == "cs") {
    filename <- paste(
      filename, 
      sprintf("%d_dim_cs_with_%d_corr", d, as.integer(100 * rho)), 
      sep = "_"
    )
  } else if (subsample_data) {
    filename <- paste(filename, sprintf("%d_dimensional", subsample_size), sep = "_")
  }
  filename <- paste(filename, sprintf("seed_%d", sampler_seed), sep = "_")
  filename <- paste0(filename, ".rds")
  path <- paste0(output_folder, filename)
  return(path)
}

get_path_to_sampler_result <- function(
    sampler_type, precond_type = "none", sampler_seed, output_type, 
    target_type = "hiv", d = NA, rho = 0
  ) {
  filename <- sampler_type
  if (precond_type != "none") {
    filename <- paste(filename, precond_type, "preconditioned", sep = "_")
  }
  if (target_type == "cs" || 
      substr(target_type, 1, 10) == "rotated_cs" || 
      substr(target_type, 1, 12) == "unrotated_cs") {
    filename <- paste(
      filename, 
      sprintf("%d_dim_%s_with_%d_corr", d, target_type, as.integer(100 * rho)), 
      sep = "_"
    )
  } else if (subsample_data) {
    filename <- paste(filename, sprintf("%d_dimensional", subsample_size), sep = "_")
  }
  filename <- paste(filename, output_type, sep = "_")
  filename <- paste(filename, sprintf("seed_%d", sampler_seed), sep = "_")
  filename <- paste0(filename, ".rds")
  path <- paste0(output_folder, filename)
  return(path)
}

get_zigzag_sampler_name <- function(nuts_flag, stepsize_multiplier) {
  if (nuts_flag) {
    sampler_name <- sprintf(
      'zigzag_nuts_%s_base_integ_time',
      ifelse(stepsize_multiplier == 0.1, 'shorter', 'longer')
    )
  } else {
    time_length <- to_verbal_expression(stepsize_multiplier)
    sampler_name <- sprintf('zigzag_hmc_%s_integ_time', time_length)
  }
  return(sampler_name)
}

to_verbal_expression <- function(stepsize_multiplier) {
  if (stepsize_multiplier <= 0.5) {
    time_length <- "shortest"
  } else if (stepsize_multiplier <= 1 / sqrt(2)) {
    time_length <- "shorter"
  } else if (stepsize_multiplier <= 1) {
    time_length <- "medium"
  } else if (stepsize_multiplier <= sqrt(2)) {
    time_length <- "longer"
  } else {
    time_length <- "longest"
  }
  return(time_length)
}

get_rotated_cs_target_name <- function(truncated, seed = 111) {
  target_type <- "rotated_cs"
  if (truncated) { target_type <- paste0(target_type, "_truncated") }
  if (is.null(seed) || is.na(seed)) {
    target_type <- paste0("un", target_type)
  } else {
    target_type <- paste(target_type, as.character(seed / 111), sep = "_")
  }
  return(target_type)
}

load_rotated_cs_result_summary <- function(sampler_type, sampler_seed, target_type, rho, d) {
  path_to_summary_result <- get_path_to_sampler_result(
    sampler_type = sampler_type, sampler_seed = sampler_seed,
    target_type = target_type, rho = rho, d = d, output_type = "summary"
  )
  summary <- readRDS(path_to_summary_result)
  if (is.null(summary$pc_ess)) {
    library(coda)
    list2env(
      load_rotated_cs_target_params(rho, d, target_seed, truncated),
      envir=environment()
    )
    samples <- readRDS(get_path_to_sampler_result(
      sampler_type = sampler_type, sampler_seed = sampler_seed,
      target_type = target_type, rho = rho, d = d, output_type = "samples"
    ))
    pc_samples <- samples %*% pc
    pc_comp_samples <- samples %*% pc_complement
    pc_ess <- effectiveSize(pc_samples)
    pc_comp_ess <- effectiveSize(pc_comp_samples)
    summary <- c(summary, list(pc_ess = pc_ess, pc_comp_ess = pc_comp_ess))
    saveRDS(summary, file = path_to_summary_result)
  }
  return(summary)
}

average_ess_over_seeds <- function(
    sampler_type, target_type, rho, d, target_seed, truncated, 
    normalize_by_time = TRUE, sec_per_time_unit = 60
  ) {
  target_type <- get_rotated_cs_target_name(truncated, target_seed)
  summary <- lapply(
    1:5, function (sampler_seed) {
      load_rotated_cs_result_summary(
        sampler_type = sampler_type, sampler_seed = sampler_seed,
        target_type = target_type, rho = rho, d = d
      )
    }
  )
  if (normalize_by_time) {
    ess <- lapply(
      summary, 
      function (per_seed) sec_per_time_unit * per_seed$ess / per_seed$runtime
    )
    pc_ess <- sapply(
      summary, 
      function (per_seed) sec_per_time_unit * per_seed$pc_ess / per_seed$runtime
    )
  } else {
    ess <- lapply(summary, function (per_seed) per_seed$ess)
    ps_ess <- sapply(summary, function (per_seed) per_seed$pc_ess)
  }
  ave_ess <- tibble(
    rho = rho,
    target_seed = target_seed,
    min = mean(sapply(ess, min)),
    pc = mean(pc_ess),
  )
  return(ave_ess)
}

format_numeric <- function(numeric_var, digits) {
  return(formatC(signif(numeric_var, digits), digits, format = "fg"))
}

print_ave_ess_table <- function(
    ess_table, flatten_across = "rho", signif_digits = 2
  ) {
  ess_table <- 
    ess_table %>%
    mutate(
      min = format_numeric(min, signif_digits), 
      pc = format_numeric(pc, signif_digits)
    ) 
  if ("stepsize_multiplier" %in% colnames(ess_table)) {
    ess_table <- ess_table %>% group_by(stepsize_multiplier) 
  }
  if (flatten_across == "rho") {
    ess_table %>% select(-target_seed)
    print(
      ess_table %>% 
        select(-target_seed) %>%
        pivot_wider(
          names_from = !!sym(flatten_across),
          values_from = c(min, pc),
          names_vary = "slowest"
        ) %>%
        kable(format = "latex") %>%
        kable_styling()
    )
  } else if (flatten_across == "target_seed") {
    print(
      ess_table %>% 
        select(-rho, -min) %>%
        pivot_wider(
          names_from = !!sym(flatten_across),
          values_from = pc,
          names_vary = "fastest"
        ) %>%
        kable(format = "latex") %>%
        kable_styling()
    )
  }
}

compute_extreme_eigvals <- function(symMatrix) {
  return(mgcv::slanczos(A = symMatrix, k = 1, kl = 1)[['values']])
}

