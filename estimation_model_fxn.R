# ===============================================
# Model functions
# ===============================================


# BASELINE SI MODEL WITHOUT INTERVENTION ----------------------------------
SImod <- function(tt, yy, parms) with(c(parms, as.list(yy)), {
	I <- I1 + I2 + I3 + I4
	N <- I + S
	transmissionCoef <- Beta * exp(-alpha * I/N)
	deriv <- rep(NA, 7)
	deriv[1] <- birthRt * N - deathRt * S - transmissionCoef * S * I/N
	deriv[2] <- transmissionCoef * S * I/N - progRt * I1 - deathRt * I1
	deriv[3] <- progRt * I1 - progRt * I2 - deathRt * I2
	deriv[4] <- progRt * I2 - progRt * I3 - deathRt * I3
	deriv[5] <- progRt * I3 - progRt * I4 - deathRt * I4
	deriv[6] <- transmissionCoef * S * I/N
	deriv[7] <- progRt * I4
	return(list(deriv))
})

#  SI MODEL WITH INTERVENTION ---------------------------------------------
SI4control <- function(t, y, parms) {
	with(c(as.list(y), parms), {
		N <- sum(y[1:5])
		I <- I1 + I2 + I3 + I4
		lambdaHat <- Beta * exp(-alpha * I/N)
		g <- 4 * progRt
		control_effect <- min(1, 1 - cMax / (1 + exp(-cRate * (t - cHalf))))
		birth <- birthRt * N
		infection <- control_effect * lambdaHat * S * I/N
		death <- deathRt * y[1:5]
		progression <- g * c(I1, I2, I3, I4)
		dSdt <- birth - infection - death[1]
		dIdt <- c(infection, progression[1:3]) - progression - death[2:5]
		dCIdt <- infection
		dCDdt <- progression[4]
		return(list(c(dSdt, dIdt, dCIdt, dCDdt)))
	})
}


# Function to run the deterministic model simulation ----------------------
## Based on the ODE system defined in SImod().
simEpidemic <- function(init, tseq = tseqMonth, modFunction=SImod, parms = disease_params()) {
	simDat <- as.data.frame(lsoda(init, tseq, modFunction, parms=parms))
	simDat$I <- rowSums(simDat[, Is])
	simDat$N <- rowSums(simDat[, c('S',Is)])
	simDat$P <- with(simDat, I/N)
	return(simDat)
}


# ===============================================
# Estimation functions
# ===============================================


# Parameter substitute function -------------------------------------------
# Description: Substitute estimated parameters into the model parameter list
# Purpose: Converts log-transformed estimated parameters back to normal scale
# Inputs:
#   - fit.params: Named vector of estimated parameters (some in log scale)
#   - fixed.params: List of all model parameters (including fixed ones)
# Outputs:
#   - A complete parameter list with estimated parameters substituted in, ready for simulation
subsParms <- function(fit.params, fixed.params = disease_params()) {
	
	within(fixed.params, {
		# Identify which parameters were estimated on log scale
		loggedParms <- names(fit.params)[grepl('log_', names(fit.params))]
		# Identify which were estimated on natural scale
		unloggedParms <- names(fit.params)[!grepl('log_', names(fit.params))]
		
		# Update parameters estimated on natural scale
		for (nm in unloggedParms) assign(nm, as.numeric(fit.params[nm]))
		
		# Update parameters estimated on log scale (convert back by exponentiating)
		for (nm in loggedParms) assign(gsub('log_', '', nm), exp(as.numeric(fit.params[nm])))
		
		rm(nm, loggedParms, unloggedParms) # Clean up workspace variables
	})
}



#  Negative log-likelihood (MLE) function ---------------------------------
# Description: Calculate negative log-likelihood for prevalence data
# Purpose: Quantifies how well the model (with given parameters) fits the observed data
# Inputs:
#   - parms: List of parameters for the simulation
#   - obsDat: Observed (sampled) prevalence data
#   - simEpidemicFunction: Simulation function to generate model outputs
# Outputs:
#   - The negative log-likelihood value (to be minimized)
nllikelihood <- function(parms, obsDat, simEpidemicFunction) {
	
	# Simulate the epidemic with the provided parameters
	simDat <- simEpidemicFunction(init, tseqMonth, parms)
	
	# Match times in simulation with observation times
	matchedTimes <- simDat$time %in% obsDat$time
	
	# Calculate log-likelihood for each observation (Binomial distribution)
	nlls <- -dbinom(obsDat$numPos, obsDat$numSamp, prob = simDat$P[matchedTimes], log = TRUE)
	
	# Return total negative log-likelihood
	return(sum(nlls))
}


# Maximum Likelihood Estimation function ----------------------------------
# Description: Estimate intervention parameters using Maximum Likelihood Estimation (MLE)
# Inputs:
#   - myDat: Observed (sampled) prevalence data
#   - simEpidemicFunction: Epidemic simulation function
#   - init_guess: Starting guess for the parameters to be estimated (on log scale)
#   - fixed_params: List of other fixed model parameters
# Outputs:
#   - A list containing:
#       - params: Estimated parameter values (on real scale)
#       - se: Standard errors of the estimated parameters
#       - loglik: Final negative log-likelihood value
#       - fit: Full optim() output for diagnostics
estimate_MLE <- function(myDat, simEpidemicFunction,
												 init_guess = c(log_cMax = log(0.5), log_cRate = log(1)),
												 fixed_params = disease_params()) {
	
	# Create an objective function that returns the NLL for given parameter values
	objFXN_MLE <- function(fit.params, fixed.params = fixed_params, obsDat = myDat) {
		parms <- subsParms(fit.params, fixed.params) # Substitute parameters into full parameter list
		nllikelihood(parms, obsDat, simEpidemicFunction) # Calculate NLL
	}
	
	# Run optimization using Nelder-Mead method to find parameters that minimize NLL
	mle_fit <- optim(par = init_guess, objFXN_MLE, fixed.params = fixed_params, obsDat = myDat,
									 method = "Nelder-Mead", control = list(trace = 0, maxit = 1000, reltol = 1e-7), hessian = TRUE)
	
	# Extract estimated parameters (convert from log scale)
	est_params <- exp(mle_fit$par)
	
	# Estimate covariance matrix from the Hessian
	cov_matrix <- solve(mle_fit$hessian)
	
	# Calculate standard errors
	se_params <- sqrt(diag(cov_matrix))
	
	return(list(params = est_params, se = se_params, loglik = mle_fit$value, fit = mle_fit))
}

# Least Squares Estimation function ---------------------------------------
## Purpose: This function calculates the sum of squared differences 
## between the observed prevalence (numPos/numSamp) from sampled data and the prevalence 
## predicted by the epidemic model using a given set of parameters.
## We return the sum of squared residual 
# Inputs:
#   - parms: List of model parameters to be used for simulation.
#            Typically includes transmission rate, progression rate, etc.
#   - obsDat: Data frame containing observed prevalence data with the following columns:
#             * time: time points of observation
#             * numPos: number of positive cases at each time point
#             * numSamp: total number of individuals sampled at each time point
#
# Outputs:
#   - A single numeric value: the sum of squared residuals, which quantifies the overall 
#     discrepancy between the observed data and model predictions. This is the value minimized 
#     during least squares parameter estimation.

sum_squares <- function(parms = disease_params(), obsDat=myDat) {
	simDat <- simEpidemic(init, parms=parms)
	
	## What are the rows from our simulation at which we have observed data?
	matchedTimes <- simDat$time %in% obsDat$time
	
	## What is bserved prevalence (proportion positive)
	observed_prev <- obsDat$numPos / obsDat$numSamp
	
	## What are the model predicted prevalences from the simulation at matched time points
	predicted_prev <- simDat$P[matchedTimes]
	
	## Then we calculate squared residuals
	squared_residuals <- (observed_prev - predicted_prev)^2
	
	return(sum(squared_residuals))
}
