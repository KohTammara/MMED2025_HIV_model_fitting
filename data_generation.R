# HIV Model fitting for evaluating interventions: A simulation exercise 
#SIMULATING HIV PREVALENCE DATA FROM 1985 TO 2024

#clear the environment
rm(list = ls())

#load required packages

require(dplyr); require(boot); require(deSolve); require(ellipse); require(ggplot2)

#simulate HIV prevalence data from 1985 to 2024 with "function"
#Using Lab6 tutorials

## Function that makes a list of disease parameters with default values
disease_params <- function(Beta = 0.9 ## transmission coefficient when prevalence is 0 
													 , alpha = 8 ## for transmission coefficient: decline with prevalence
													 , progRt = (1/10)*4 ## rate of of progression through each of the I classes, for 10 years total
													 , birthRt = .03 ## birth rate, 3% of people give birth per year
													 , deathRt = 1/60 ## 60 year natural life expectancy
){
	return(as.list(environment()))
}

disease_params()

initPrev <- exp(-7) ## infected at start
tseqMonth <- seq(1985, 2024, by = 1/12)
init <- c(S=1-initPrev, I1=initPrev, I2=0, I3=0, I4=0, CI = 0, CD = 0) ## modeling proportion of population
Is <- paste0('I',1:4) ## for easy indexing

## Define the SI ODE model. This model is equivalent to the third model in the HIV in Harare tutorial
## 	(though some parameters may have different names)
SImod <- function(tt, yy, parms) with(c(parms,as.list(yy)), {
	## State variables are: S, I1, I2, I3, I4
	## derived quantitties
	I <- I1+I2+I3+I4           ## total infected
	N <- I + S                 ## total population
	transmissionCoef <- Beta * exp(-alpha * I/N) ## Infectious contact rate
	## state variable derivatives (ODE system)
	deriv <- rep(NA,7)
	deriv[1] <-	birthRt*N - deathRt*S - transmissionCoef*S*I/N ## Instantaneous rate of change: Susceptibles
	deriv[2] <-	transmissionCoef*S*I/N - progRt*I1 - deathRt*I1 ## Instantaneous rate of change: Infection class I1
	deriv[3] <-	progRt*I1 - progRt*I2 - deathRt*I2 ## Instantaneous rate of change:  Infection class I2
	deriv[4] <-	progRt*I2 - progRt*I3 - deathRt*I3 ## Instantaneous rate of change: Infection class I3 
	deriv[5] <-	progRt*I3 - progRt*I4 - deathRt*I4 ## Instantaneous rate of change: Infection class I4
	deriv[6] <-	transmissionCoef*S*I/N ## Instantaneous rate of change: Cumulative incidence
	deriv[7] <-	progRt*I4 ## Instantaneous rate of change: Cumulative mortality
	return(list(deriv))
})

## Function to run the deterministic model simulation, based on the ODE system defined in SImod().
simEpidemic <- function(init, tseq = tseqMonth, modFunction=SImod, parms = disease_params()) {
	simDat <- as.data.frame(lsoda(init, tseq, modFunction, parms=parms))
	simDat$I <- rowSums(simDat[, Is])
	simDat$N <- rowSums(simDat[, c('S',Is)])
	simDat$P <- with(simDat, I/N)
	return(simDat)
}

## Function to 'sample' the population:
## From a simulated epidemic, measure prevalence at several time points by drawing
## cross-sectional samples of individuals at each time, testing them, and then calculating sample
## prevalence and associated binomial confidence intervals

sampleEpidemic <- function(simDat # Simulated "data" which we treat as real 
													 , sampleDates = seq(1985, 2024, by = 3) # Sample every 3 years 
													 , numSamp = rep(100, length(sampleDates)) # Number of individuals sampled at each time point
){
	prev_at_sample_times <- simDat[simDat$time %in% sampleDates, 'P']
	numPos <- rbinom(length(numSamp), numSamp, prev_at_sample_times)
	lci <- mapply(function(x,n) binom.test(x,n)$conf.int[1], x = numPos, n = numSamp)
	uci <- mapply(function(x,n) binom.test(x,n)$conf.int[2], x = numPos, n = numSamp)    
	return(data.frame(time = sampleDates, numPos, numSamp, sampPrev =  numPos/numSamp,
										lci = lci, uci = uci))
}

## Run system of ODEs for "true" parameter values
trueParms <- disease_params() # Default model parameters are defined in lines 20-26
simDat <- simEpidemic(init, parms = trueParms) # Simulated epidemic (underlying process)
View(simDat)

par(bty='n', lwd = 2)
# Plot simulated prevalence through time:
with(simDat, plot(time, P, xlab = '', ylab = 'prevalence', type = 'l', ylim = c(0,.4), col='red', las = 1))
## Take cross-sectional sample of individuals to get prevalence estimates at multiple time points:
set.seed(1) # Initiate the random number generator (so everyone's simulation looks the same)
myDat <- sampleEpidemic(simDat) # Simulate data from the sampling process (function defined above)


