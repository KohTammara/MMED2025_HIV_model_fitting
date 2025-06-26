
#SIMULATING 2 SCENARIOS OF 20 DATASETS EACH
#Armand Mutwadi, Evidence, Folk

# 1. Clear the workspace =======================================================

rm(list = ls())

# 2. Load required libraries ==================================================

require(dplyr)
require(boot)
require(deSolve)
require(ellipse)
require(ggplot2)

# 3. Global Parameters and Initial Conditions ===================================

initial_prevalence <- exp(-9.5)
time_sequence <- seq(1976, 2025, by = 1)
initial_conditions <- c(S = 1 - initial_prevalence, I1 = initial_prevalence, I2 = 0, I3 = 0, I4 = 0, CI = 0, CD = 0)
infection_states <- paste0('I', 1:4)

# 4. Disease Parameters (C(t)) =========================================================

create_disease_parameters <- function(Beta = 0.6 ## transmission coefficient when prevalence is 0 
																			, alpha = 3.5 ## for transmission coefficient: decline with prevalence
																			, progRt = (1/15) ## rate of of progression through each of the I classes, for 10 years total
																			, birthRt = .03 ## birth rate, 3% of people give birth per year
																			, deathRt = 1/60 ## 60 year natural life expectancy
																			, cMax = 0.7 ## maximum control effect, 80% reduction in transmission
																			, cRate = 0.5 ## rate of control effect increase, 2.8 years to reach cMax
																			, cStart = 1998 ## year when control effect started
)
{
	return(as.list(environment())) 
}
true_params <- create_disease_parameters()

# 5. SI simulation function (SI4control) =====================================================

SI4control <- function(time, state, parameters) {
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
# 6. Sampling Functions (sample_epidemic) ==========================================

simulate_epidemic <- function(initial_conditions, time_sequence, model_function, parameters) {
	sim_data <- as.data.frame(lsoda(initial_conditions, time_sequence, model_function, parms = parameters))
	sim_data$I <- rowSums(sim_data[, infection_states])
	sim_data$N <- rowSums(sim_data[, c('S', infection_states)])
	sim_data$P <- sim_data$I / sim_data$N
	return(sim_data)
}

simulated_data <- simulate_epidemic(initial_conditions,time_sequence,SI4control,true_params)

sample_epidemic1 <- function(simulated_data, sample_times = seq(1985, 2024, by = 2), sample_sizes = rep(250, length(sample_times))) {
	prevalence <- simulated_data$P[round(simulated_data$time, 4) %in% round(sample_times, 4)]
	num_positive <- rbinom(length(sample_sizes), sample_sizes, prevalence)
	lci <- mapply(function(x, n) binom.test(x, n)$conf.int[1], x = num_positive, n = sample_sizes)
	uci <- mapply(function(x, n) binom.test(x, n)$conf.int[2], x = num_positive, n = sample_sizes)
	
	data.frame(time = sample_times, numPos = num_positive, numSamp = sample_sizes, sampPrev = num_positive / sample_sizes, lci = lci, uci = uci)
}

sample_epidemic2 <- function(simulated_data,sample_times = seq(1985, 2024, by = 5),sample_sizes = c(200, 200, 100, 100, 50, 50, 50, 50)) {
	prevalence <- simulated_data$P[round(simulated_data$time, 4) %in% round(sample_times, 4)]
	num_positive <- rbinom(length(sample_sizes), sample_sizes, prevalence)
	lci <- mapply(function(x, n) binom.test(x, n)$conf.int[1], x = num_positive, n = sample_sizes)
	uci <- mapply(function(x, n) binom.test(x, n)$conf.int[2], x = num_positive, n = sample_sizes)
	
	data.frame(time = sample_times,numPos = num_positive,numSamp = sample_sizes,sampPrev = num_positive / sample_sizes,lci = lci,uci = uci)
}

# 7. Data Simulation: 1st SCENARIO =============================================================== 

# Number of simulated datasets
n_simPData1 <- 20 #number of datasets to simulate

# Create an empty list to store each dataset
simPdata_list1 <- vector("list", n_simPData1) #list to store simulated datasets

# Loop to generate and store each dataset
for (i in 1:n_simPData1) {
	simPdata_list1[[i]] <- sample_epidemic1(simulated_data)
}

# 8. Data Simulation: 2nd SCENARIO ================================================================ 

# Number of simulated datasets
n_simPData2 <- 20 #number of datasets to simulate
simPdata_list2 <- replicate(n_simPData2, sample_epidemic2(simulated_data), simplify = FALSE)

# 9. Print simulated_Data ================================================================ 

print(simPdata_list1)
print(simPdata_list2)

saveRDS(simPdata_list1, file = "simPdata_list1.rds")
saveRDS(simPdata_list2, file = "simPdata_list2.rds")

# 10.Plotting the simulated datasets =================================================

# 1st Scenario
# Combine all into one data.frame with a simulation ID
simPdata_df1 <- bind_rows(simPdata_list1, .id = "sim_id")

ggplot(simPdata_df1, aes(x = time, y = sampPrev, group = sim_id)) +
	geom_line(alpha = 0.3, color = "darkblue") +
	labs(title = "Simulated HIV Prevalence Over Time (Multiple Samples 1st Scenario)",
			 x = "Year", y = "Sampled Prevalence") +
	theme_minimal(base_size = 12)

# 2nd Scenario
# Combine all into one data.frame with a simulation ID
simPdata_df2 <- bind_rows(simPdata_list2, .id = "sim_id")

ggplot(simPdata_df2, aes(x = time, y = sampPrev, group = sim_id)) +
	geom_line(alpha = 0.3, color = "darkgreen") +
	labs(title = "Simulated HIV Prevalence Over Time (Multiple Samples 2nd Scenario)",
			 x = "Year", y = "Sampled Prevalence") +
	theme_minimal(base_size = 12)



