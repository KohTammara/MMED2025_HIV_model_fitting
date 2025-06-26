## Move over necessary functions and create 3rd model function with new c(t) and make progression rate constant (See george)
library(dplyr)
library(tibble)
library(deSolve)
library(ellipse)

# Define inverse logit and logit functions
logit <- function(x) log(x / (1 - x))
invlogit <- function(x) 1 / (1 + exp(-x))


# Model no intervention ---------------------------------------------------
SImod <- function(tt, yy, parms) with(c(parms, as.list(yy)), {
	I <- I1 + I2 + I3 + I4
	N <- I + S
	transmissionCoef <- Beta * exp(-alpha * I/N)
	g <- 4 * progRt
	deriv <- rep(NA, 7)
	deriv[1] <- birthRt * N - deathRt * S - transmissionCoef * S * I/N
	deriv[2] <- transmissionCoef * S * I/N - progRt * I1 - deathRt * I1
	deriv[3] <- g * I1 - g * I2 - deathRt * I2
	deriv[4] <- g * I2 - g * I3 - deathRt * I3
	deriv[5] <- g * I3 - g * I4 - deathRt * I4
	deriv[6] <- transmissionCoef * S * I/N
	deriv[7] <- g * I4
	return(list(deriv))
})

SI4control.zero <- function(time, state, parameters) {
	with(c(as.list(state), parameters), {
		N <- sum(state[1:5])
		I <- I1 + I2 + I3 + I4
		lambda <- Beta * exp(-alpha * I / N)
		g <- 4 * progRt
		control_effect <- 1
		dS <- birthRt * N - control_effect * lambda * S * I / N - deathRt * S
		progression <- g * c(I1, I2, I3, I4)
		dI <- c(control_effect * lambda * S * I / N, progression[1:3]) - progression - deathRt * c(I1, I2, I3, I4)
		dCI <- control_effect * lambda * S * I / N
		dCD <- progression[4]
		list(c(dS, dI, dCI, dCD))
	})
}


# Model original C(t) -----------------------------------------------------
SI4control <- function(time, state, parameters) {
	with(c(as.list(state), parameters), {
		N <- sum(state[1:5])
		I <- I1 + I2 + I3 + I4
		lambda <- Beta * exp(-alpha * I / N)
		g <- 4 * progRt
		control_effect <- min(1, 1 - cMax / (1 + exp(-cRate * (time - cHalf))))
		
		dS <- birthRt * N - control_effect * lambda * S * I / N - deathRt * S
		progression <- g * c(I1, I2, I3, I4)
		dI <- c(control_effect * lambda * S * I / N, progression[1:3]) - progression - deathRt * c(I1, I2, I3, I4)
		dCI <- control_effect * lambda * S * I / N
		dCD <- progression[4]
		
		list(c(dS, dI, dCI, dCD))
	})
}


# Modified C(t) -----------------------------------------------------------
SI4control_time <- function(time, state, parameters) {
	with(c(as.list(state), parameters), {
		N <- sum(state[1:5])
		I <- I1 + I2 + I3 + I4
		lambda <- Beta * exp(-alpha * I / N)
		g <- 4 * progRt
		control_effect <- pmin(1, cMax+(1-cMax)*exp(-(time-cStart)*cRate))
		dS <- birthRt * N - control_effect * lambda * S * I / N - deathRt * S
		progression <- g * c(I1, I2, I3, I4)
		dI <- c(control_effect * lambda * S * I / N, progression[1:3]) - progression - deathRt * c(I1, I2, I3, I4)
		dCI <- control_effect * lambda * S * I / N
		dCD <- progression[4]
		
		list(c(dS, dI, dCI, dCD))
	})
}


# Substitute parameter function -------------------------------------------
substitute_parameters <- function(estimated_params, fixed_params) {
	for (name in names(estimated_params)) {
		if (grepl('log_', name)) {
			param_name <- gsub('log_', '', name)
			fixed_params[[param_name]] <- exp(estimated_params[[name]])
		} else if (grepl('logit_', name))  {
			param_name <- gsub('logit_', '', name)
			fixed_params[[param_name]] <- 1/(1+exp(-estimated_params[[name]]))
		} else {
			fixed_params[[name]] <- estimated_params[[name]]
		}
	}
	return(fixed_params)
}


# Simulate epidemic and sample ------------------------------------------------------

simulate_epidemic <- function(initial_conditions, time_sequence, model_function, parameters) {
	sim_data <- as.data.frame(lsoda(initial_conditions, time_sequence, model_function, parms = parameters))
	sim_data$I <- rowSums(sim_data[, infection_states])
	sim_data$N <- rowSums(sim_data[, c('S', infection_states)])
	sim_data$P <- sim_data$I / sim_data$N
	return(sim_data)
}

sample_epidemic <- function(simulated_data, sample_times = seq(1985, 2024, by = 1), sample_sizes = rep(1000, length(sample_times))) {
	prevalence <- simulated_data$P[round(simulated_data$time, 4) %in% round(sample_times, 4)]
	num_positive <- rbinom(length(sample_sizes), sample_sizes, prevalence)
	lci <- mapply(function(x, n) binom.test(x, n)$conf.int[1], x = num_positive, n = sample_sizes)
	uci <- mapply(function(x, n) binom.test(x, n)$conf.int[2], x = num_positive, n = sample_sizes)
	
	data.frame(time = sample_times, numPos = num_positive, numSamp = sample_sizes, sampPrev = num_positive / sample_sizes, lci = lci, uci = uci)
}


# Create disease parameters -----------------------------------------------
create_disease_parameters <- function(Beta = 0.6 ## transmission coefficient when prevalence is 0 
																			, alpha = 3.5 ## for transmission coefficient: decline with prevalence
																			, progRt = (1/15) ## rate of of progression through each of the I classes, for 10 years total
																			, birthRt = .03 ## birth rate, 3% of people give birth per year
																			, deathRt = 1/60 ## 60 year natural life expectancy
																			, cMax = 0.7 # 1 - THIS is intervention effect
																			, cRate = 0.5
																			, cHalf = 1998
																			, cStart = 1990
){
	return(as.list(environment())) ## ARG
}


# Calculate negative log likelihood function ------------------------------
calculate_nll <- function(parameters, observed_data, model_function) {
	simulated_data <- simulate_epidemic(initial_conditions, time_sequence, model_function, parameters)
	matched_times <- simulated_data$time %in% observed_data$time
	if (any(is.na(simulated_data$P[matched_times])) || any(simulated_data$P[matched_times] < 0 | simulated_data$P[matched_times] > 1)) {
		warning(paste0("Invalid prevalence values (NA or outside [0,1]) detected"))
		return(1e6)
	}
	nll_values <- dbinom(observed_data$numPos, observed_data$numSamp, prob = simulated_data$P[matched_times], log = TRUE)
	return(-sum(nll_values))
}


# MLE function ------------------------------------------------------------
estimate_mle <- function(observed_data, model_function, initial_guess, fixed_params) {
	print(observed_data)
	objective_function <- function(fit_params) {
		updated_params <- substitute_parameters(fit_params, fixed_params)
		calculate_nll(updated_params, observed_data, model_function)
	}
	
	mle_fit_pre <- optim(par = initial_guess, fn = objective_function, method = "SANN", control = list(trace = 1, maxit = 150))
	mle_fit <- optim(par = mle_fit_pre$par, fn = objective_function, method = "Nelder-Mead", control = list(trace = 1, maxit = 1000), hessian = TRUE)
	
	print(mle_fit)
	estimated_params <- mle_fit$par
	covariance_matrix <- solve(mle_fit$hessian)
	print(covariance_matrix)
	standard_errors <- sqrt(diag(covariance_matrix))
	
	list(params = estimated_params, se = standard_errors, loglik = mle_fit$value, fisherInfMatrix = covariance_matrix)
}


# LSE function ------------------------------------------------------------
estimate_lse <- function(observed_data, model_function, initial_guess, fixed_params) {
	lse_objective <- function(fit_params) {
		updated_params <- substitute_parameters(fit_params, fixed_params)
		simulated_data <- simulate_epidemic(initial_conditions, time_sequence, model_function, updated_params)
		
		matched_times <- simulated_data$time %in% observed_data$time
		observed_prev <- observed_data$numPos / observed_data$numSamp
		predicted_prev <- simulated_data$P[matched_times]
		sum((observed_prev - predicted_prev)^2)
	}
	
	lse_fit <- optim(par = initial_guess, fn = lse_objective, method = "Nelder-Mead", control = list(trace = 1, maxit = 1000), hessian = TRUE)
	
	estimated_params <- lse_fit$par
	covariance_matrix <- solve(lse_fit$hessian)
	standard_errors <- sqrt(diag(covariance_matrix))
	
	list(params = estimated_params, se = standard_errors, ssr = lse_fit$value, fit = lse_fit)
}


# profile likelihood 
objXcMax_cRate <- function(cMax, cRate, fixed_params, observed_data, model_function) {
	fit_params <- c(logit_cMax = logit(cMax), log_cRate = log(cRate))
	updated_params <- substitute_parameters(fit_params, fixed_params)
	calculate_nll(updated_params, observed_data, model_function)
}


is_in_conf_region <- function(cMax_est, cRate_est, fixed_params, observed_data, model_function) {
	res <- 40
	grid_cMax <- seq(0.5, 0.95, length.out = res)
	grid_cRate <- seq(0.1, 2.5, length.out = res)
	prof_fun <- Vectorize(function(cM, cR) objXcMax_cRate(cM, cR, fixed_params, observed_data, model_function),
												vectorize.args = c("cM", "cR"))
	surface <- outer(grid_cMax, grid_cRate, prof_fun)
	min_ll <- min(surface, na.rm = TRUE)
	chi_cut <- min_ll + qchisq(0.95, df = 2) / 2
	est_ll <- objXcMax_cRate(cMax_est, cRate_est, fixed_params, observed_data, model_function)
	return(est_ll <= chi_cut)
}

# Wrapper function for each dataset
process_one <- function(observed_data, scenario_id, dataset_id, true_cMax, true_cRate) {
	init_guess_mle <- c(logit_cMax = logit(0.7), log_cRate = log(0.5))
	# init_guess_lse <- c(log_cMax = log(0.7), log_cRate = log(0.5))
	fixed_params <- create_disease_parameters(cMax = true_cMax, cRate = true_cRate)
	
	# Estimate MLE
	mle_result <- estimate_mle(observed_data, SI4control_time, init_guess_mle, fixed_params)
	# Estimate LSE
	# lse_result <- estimate_lse(observed_data, SI4control_time, init_guess_lse, fixed_params)
	# 
	# # MLE estimates and CIs
	# mle_cMax <- invlogit(mle_result$params["logit_cMax"])
	# mle_cRate <- exp(mle_result$params["log_cRate"])
	# mle_se_cMax <- mle_result$se["logit_cMax"] * (mle_cMax * (1 - mle_cMax))
	# mle_se_cRate <- mle_result$se["log_cRate"] * mle_cRate
	# mle_ci_cMax <- mle_cMax + c(-1.96, 1.96) * mle_se_cMax
	# mle_ci_cRate <- mle_cRate + c(-1.96, 1.96) * mle_se_cRate
	# 
	# # LSE estimates and CIs
	# lse_cMax <- exp(lse_result$params["log_cMax"])
	# lse_cRate <- exp(lse_result$params["log_cRate"])
	# lse_se_cMax <- lse_result$se["log_cMax"] * (lse_cMax * (1 - lse_cMax))
	# lse_se_cRate <- lse_result$se["log_cRate"] * lse_cRate
	# lse_ci_cMax <- lse_cMax + c(-1.96, 1.96) * lse_se_cMax
	# lse_ci_cRate <- lse_cRate + c(-1.96, 1.96) * lse_se_cRate
	# 
	# # Check if the estimated MLE points are within the 95% profile likelihood region
	# in_region_mle <- is_in_conf_region(mle_cMax, mle_cRate, fixed_params, observed_data, SI4control_time)
	# 
	# tibble(
		# scenario = scenario_id,
		# dataset = dataset_id,
		# mle_cMax = mle_cMax,
		# mle_cRate = mle_cRate,
		# mle_cMax_lci = mle_ci_cMax[1],
		# mle_cMax_uci = mle_ci_cMax[2],
		# mle_cRate_lci = mle_ci_cRate[1],
		# mle_cRate_uci = mle_ci_cRate[2],
		# lse_cMax = lse_cMax,
		# lse_cRate = lse_cRate,
		# lse_cMax_lci = lse_ci_cMax[1],
		# lse_cMax_uci = lse_ci_cMax[2],
		# lse_cRate_lci = lse_ci_cRate[1],
		# lse_cRate_uci = lse_ci_cRate[2],
		# true_cMax = true_cMax,
		# true_cRate = true_cRate,
		# mle_se_cMax = mle_se_cMax,
		# mle_se_cRate = mle_se_cRate,
		# lse_se_cMax = lse_se_cMax,
		# lse_se_cRate = lse_se_cRate,
		# mle_in_region = in_region_mle
	# )
	
	return (c(invlogit(mle_result$params["logit_cMax"]),exp(mle_result$params["log_cRate"])))
}
# data1 <- readRDS("/Users/tkoh/Documents/MMED2025/MMED2025_HIV_model_fitting/simPdata_list2.rds")
# process_one(data1, scenario_id = 1, c(seq(1:20)), 0.8, 1)

# Main loop for processing all datasets
analyze_all <- function(path1, path2, cMax1, cRate1, cMax2, cRate2, output_path) {
	data1 <- readRDS(path1)
	print(data1)
	# data2 <- readRDS(path2)
	
	results1 <- purrr::map2_dfr(data1, 1:20, ~ process_one(.x, scenario_id = 1, .y, cMax1, cRate1))
	print(results1)
	# results2 <- purrr::map2_dfr(data2, 1:20, ~ process_one(.x, scenario_id = 2, .y, cMax2, cRate2))
	
	# full_results <- bind_rows(results1, results2)
	# saveRDS(full_results, output_path)
	# return(full_results)
}

# Run analysis
# results <- analyze_all(
# 	path1 = "scenario1.rds",
# 	path2 = "scenario2.rds",
# 	cMax1 = 0.7,
# 	cRate1 = 0.5,
# 	cMax2 = 0.8,
# 	cRate2 = 1.0,
# 	output_path = "fitting_results_with_se_and_region.rds"
# )

# ============================================================
# Global Parameters and Initial Conditions
# ============================================================

initial_prevalence <- exp(-9.5)
time_sequence <- seq(1976, 2025, by = 1)
initial_conditions <- c(S = 1 - initial_prevalence, I1 = initial_prevalence, I2 = 0, I3 = 0, I4 = 0, CI = 0, CD = 0)
infection_states <- paste0('I', 1:4)





##### Test
# Example of running the analysis
results <- analyze_all(
	path1 = "/Users/tkoh/Documents/MMED2025/MMED2025_HIV_model_fitting/simPdata_list1.rds",
	path2 = "/Users/tkoh/Documents/MMED2025/MMED2025_HIV_model_fitting/simPdata_list2.rds",
	cMax1 = 0.7,
	cRate1 = 0.5,
	cMax2 = 0.7,
	cRate2 = 0.5,
	output_path = "fitting_results_with_se_and_region.rds"
)