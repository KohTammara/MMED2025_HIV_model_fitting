# HIV Model Fitting for Evaluating Interventions: A Simulation Exercise


# Things I noted:
# ProgRt being used different in the two ODE systems



# ============================================================
# Setup
# ============================================================

rm(list = ls())
graphics.off()

require(dplyr)
require(boot)
require(deSolve)
require(ellipse)
require(ggplot2)

# ============================================================
# Global Parameters and Initial Conditions
# ============================================================

initial_prevalence <- exp(-9.5)
time_sequence <- seq(1976, 2025, by = 1)
initial_conditions <- c(S = 1 - initial_prevalence, I1 = initial_prevalence, I2 = 0, I3 = 0, I4 = 0, CI = 0, CD = 0)
infection_states <- paste0('I', 1:4)

# ============================================================
# Disease Model Functions
# ============================================================

SImodel <- function(time, state, parameters) {
	with(c(parameters, as.list(state)), {
		I <- I1 + I2 + I3 + I4
		N <- I + S
		transmission_rate <- Beta * exp(-alpha * I / N)
		
		dS <- birthRt * N - deathRt * S - transmission_rate * S * I / N
		dI1 <- transmission_rate * S * I / N - progRt * I1 - deathRt * I1
		dI2 <- progRt * I1 - progRt * I2 - deathRt * I2
		dI3 <- progRt * I2 - progRt * I3 - deathRt * I3
		dI4 <- progRt * I3 - progRt * I4 - deathRt * I4
		dCI <- transmission_rate * S * I / N
		dCD <- progRt * I4
		
		list(c(dS, dI1, dI2, dI3, dI4, dCI, dCD))
	})
}

SI4control <- function(time, state, parameters) {
	with(c(as.list(state), parameters), {
		N <- sum(state[1:5])
		I <- I1 + I2 + I3 + I4
		lambda <- Beta * exp(-alpha * I / N)
		g <- 4 * progRt
		control_effect <- min(1, cMax+(1-cMax)*exp(-(time-cStart)*cRate))
		
		dS <- birthRt * N - control_effect * lambda * S * I / N - deathRt * S
		progression <- g * c(I1, I2, I3, I4)
		dI <- c(control_effect * lambda * S * I / N, progression[1:3]) - progression - deathRt * c(I1, I2, I3, I4)
		dCI <- control_effect * lambda * S * I / N
		dCD <- progression[4]
		
		list(c(dS, dI, dCI, dCD))
	})
}

create_disease_parameters <- function(Beta = 0.6 ## transmission coefficient when prevalence is 0 
																			, alpha = 3.5 ## for transmission coefficient: decline with prevalence
																			, progRt = (1/15) ## rate of of progression through each of the I classes, for 10 years total
																			, birthRt = .03 ## birth rate, 3% of people give birth per year
																			, deathRt = 1/60 ## 60 year natural life expectancy
																			, cMax = 0.3 # 1 - THIS is intervention effect
																			, cRate = 0.5
																			, cStart = 1998
){
	return(as.list(environment())) ## ARG
}

# ============================================================
# Simulation and Sampling Functions
# ============================================================

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

# sim1 <- simulate_epidemic(initial_conditions, time_sequence, SI4control, create_disease_parameters(cMax = 1.1))
# sample_epidemic(sim1)

# sim1 <- simulate_epidemic(initial_conditions, time_sequence, SI4control, create_disease_parameters(cMax = 1.1))
# sample_epidemic(sim1)

sim2 <- simulate_epidemic(initial_conditions, time_sequence, SI4control, create_disease_parameters(cMax = 0.9))
simdata2 <- sample_epidemic(sim2)



# ============================================================
# Estimation Functions
# ============================================================

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

substitute_parameters(c(log_cRate = 0.5, logit_cMax = log(0.2/0.8)), create_disease_parameters())

calculate_nll <- function(parameters, observed_data, model_function) {
	# browser()
	simulated_data <- simulate_epidemic(initial_conditions, time_sequence, model_function, parameters)
	matched_times <- simulated_data$time %in% observed_data$time
	if (any(is.na(simulated_data$P[matched_times])) || any(simulated_data$P[matched_times] < 0 | simulated_data$P[matched_times] > 1)) {
		warning(paste0("Invalid prevalence values (NA or outside [0,1]) detected"))
		# warning(paste0("Invalid prevalence values (NA or outside [0,1]) detected", parameters$cMax ))
		return(1e6)
	}
	nll_values <- dbinom(observed_data$numPos, observed_data$numSamp, prob = simulated_data$P[matched_times], log = TRUE)
	return(-sum(nll_values))
}



estimate_mle <- function(observed_data, model_function, initial_guess, fixed_params) {
	
	# browser()
	
	objective_function <- function(fit_params) {
		# browser()
		updated_params <- substitute_parameters(fit_params, fixed_params)
		calculate_nll(updated_params, observed_data, model_function)
	}
	
	mle_fit_pre <- optim(par = initial_guess, fn = objective_function, method = "SANN", control = list(trace = 1, maxit = 150))
	mle_fit <- optim(par = mle_fit_pre$par, fn = objective_function, method = "Nelder-Mead", control = list(trace = 1, maxit = 1000), hessian = TRUE)
	
	estimated_params <- mle_fit$par
	covariance_matrix <- solve(mle_fit$hessian)
	standard_errors <- sqrt(diag(covariance_matrix))
	
	list(params = estimated_params, se = standard_errors, loglik = mle_fit$value, fit = mle_fit, fisherInfMatrix = covariance_matrix)
}

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



logit <- function(x){
	log(x/(1-x)) 
}

invlogit <- function(x){
	1/(1+exp(-x)) 
}

# ============================================================
# Run Simulation, Sampling, Estimation and Grid Evaluation
# ============================================================

# gather_cmax <- gather_crate <- c()
# for (ii in 1:100){
# set.seed(45)

true_params <- create_disease_parameters()
simulated_data <- simulate_epidemic(initial_conditions, time_sequence, SI4control, true_params)
observed_data <- sample_epidemic(simulated_data)

plot(simulated_data$time, simulated_data$P, type = 'l', col = 'red', ylab = 'Prevalence', xlab = 'Year',  xlim = c(1960, 2050), ylim = c(0, max(simulated_data$P, observed_data$uci)))
points(observed_data$time, observed_data$sampPrev, col = 'blue', pch = 16)
# Add confidence intervals
arrows(observed_data$time, observed_data$lci, observed_data$time, observed_data$uci, angle = 90, code = 3, length = 0.05, col = 'blue')
legend('topleft', legend = c('Simulated Prevalence', 'Observed Prevalence', '95% CI'), col = c('red', 'blue', 'blue'), pch = c(NA, 16, NA), lty = c(1, NA, 1), bty = 'n')

initial_guess <- c(logit_cMax = log(0.6/0.4), log_cRate = log(0.4))
# read from rds so specify paths
data1 <- readRDS("/Users/tkoh/Documents/MMED2025/MMED2025_HIV_model_fitting/simPdata_list1.rds")
data2 <- readRDS("/Users/tkoh/Documents/MMED2025/MMED2025_HIV_model_fitting/simPdata_list2.rds")
data3 <- readRDS("/Users/tkoh/Documents/MMED2025/MMED2025_HIV_model_fitting/simPdata_list2.rds")
data4 <- readRDS("/Users/tkoh/Documents/MMED2025/MMED2025_HIV_model_fitting/simPdata_list2.rds")
# tibble to store results
datas <- list(data1, data2, data3, data4)
results_tbl <- tibble()

# 
# data <- scenarios[[ii]]
# 
# 
# mle_result <- estimate_mle(data, SI4control, initial_guess, true_params)
# 
# 
# 
# ml.val <- mle_result$value
# conf.cutoff <- ml.val + qchisq(.95,2)/2
# nll <- calculate_nll(true_params, data, SI4control)
# mle_within_region <- (nll <= conf.cutoff)

# 
# 
# mle_cRate <- exp(mle_result$params['log_cRate'])
# mle_cMax <- 1/(1+exp(-mle_result$params['logit_cMax']))
# 
# 
# mle_matrix <- mle_result$fisherInfMatrix
# mle_se <- mle_result$se
# mle_cMax_lci <- invlogit(1/(1+exp(-mle_result$params['logit_cMax']))-qnorm(0.975)*mle_matrix[1,1])
# mle_cMax_uci <- invlogit(1/(1+exp(-mle_result$params['logit_cMax']))+qnorm(0.975)*mle_matrix[1,1])
# mle_cRate_lci <- exp(mle_result$params['log_cRate']-qnorm(0.975)*mle_matrix[1,1])
# mle_cRate_uci <- exp(mle_result$params['log_cRate']*mle_matrix[2,2])
# 
# 
# lse_result <- estimate_lse(observed_data, SI4control, initial_guess, true_params)
# lse_cRate <- exp(lse_result$params['log_cRate'])
# lse_cMax <- 1/(1+exp(-lse_result$params['logit_cMax']))
# lse_se <- lse_result$se
# 
# results_tbl <- bind_rows(
# 	results_tbl,
# 	tibble(
# 		scenario = jj,
# 		dataset = ii,
# 		mle_cRate = mle_cRate,
# 		mle_cMax = mle_cMax,
# 		mle_matrix = mle_matrix,
# 		mle_within_region <- length(mle_within_region), # 0 = true
# 		mle_cMax_lci <- mle_cMax_lci,
# 		mle_cMax_uci <- mle_cMax_uci,
# 		mle_cRate_lci <- mle_cRate_lci,
# 		mle_cRate_uci <- mle_cRate_uci,
# 		
# 		lse_cRate = lse_cRate,
# 		lse_cMax = lse_cMax,
# 		lse_se = lse_se,
# 		true_cRate = true_params$cRate,
# 		true_cMax = true_params$cMax
# 	)








for (jj in 1:length(datas)){
	scenarios <- datas[[jj]]
	for (ii in 1:50){

			data <- scenarios[[ii]]


			mle_result <- estimate_mle(data, SI4control, initial_guess, true_params)



			ml.val <- mle_result$value
			conf.cutoff <- ml.val + qchisq(.95,2)/2
			nll <- calculate_nll(true_params, data, SI4control)
			mle_within_region <- (nll <= conf.cutoff)



			mle_cRate <- exp(mle_result$params['log_cRate'])
			mle_cMax <- 1/(1+exp(-mle_result$params['logit_cMax']))


			mle_matrix <- mle_result$fisherInfMatrix
			mle_se <- mle_result$se
			mle_cMax_lci <- invlogit(1/(1+exp(-mle_result$params['logit_cMax']))-qnorm(0.975)*mle_matrix[1,1])
			mle_cMax_uci <- invlogit(1/(1+exp(-mle_result$params['logit_cMax']))+qnorm(0.975)*mle_matrix[1,1])
			mle_cRate_lci <- exp(mle_result$params['log_cRate']-qnorm(0.975)*mle_matrix[1,1])
			mle_cRate_uci <- exp(mle_result$params['log_cRate']*mle_matrix[2,2])


			lse_result <- estimate_lse(observed_data, SI4control, initial_guess, true_params)
			lse_cRate <- exp(lse_result$params['log_cRate'])
			lse_cMax <- 1/(1+exp(-lse_result$params['logit_cMax']))
			lse_se <- lse_result$se

			results_tbl <- bind_rows(
				results_tbl,
				tibble(
					scenario = jj,
					dataset = ii,
					mle_cRate = mle_cRate,
					mle_cMax = mle_cMax,
					mle_matrix = mle_matrix,
					mle_within_region <- length(mle_within_region), # 0 = true
					mle_cMax_lci <- mle_cMax_lci,
					mle_cMax_uci <- mle_cMax_uci,
					mle_cRate_lci <- mle_cRate_lci,
					mle_cRate_uci <- mle_cRate_uci,

					lse_cRate = lse_cRate,
					lse_cMax = lse_cMax,
					lse_se = lse_se,
					true_cRate = true_params$cRate,
					true_cMax = true_params$cMax
				)
			)
	}
}
saveRDS(results_tbl, "/Users/tkoh/Documents/MMED2025/MMED2025_HIV_model_fitting/results.rds")
