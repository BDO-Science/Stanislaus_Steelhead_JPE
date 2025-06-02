library(tidyverse)
library(cmdstanr)
set_cmdstan_path("C:/CMDSTAN/cmdstan-2.34.1")
library(posterior)

#######import and establish necessary data
age <- read.csv('age_structure.csv') #data frame of count of fish by age class in the Stanislaus for 3 years
n_sim <- 500
alpha <- 0.213 #for fecundity function from hodge et al. 2016
beta <- 2.43 #for fecundity function from hodge et al. 2016


#######calculate egg production
sizeProp <- data.frame(Age = c(1,2,3,4,5,6), #setting up size-age structure using FL data from Cramer study
                       minFL = c(104,181,228,324,720,720), 
                       maxFL = c(284,310,400,690,720,720)) %>%
  left_join(age, by = 'Age') %>% #join age data to this
  group_by(Year) %>%
  mutate(prop = prop.table(Count)) #calculate proportion of each age by year

# Function to simulate fecundity per fish
simulate_egg_production <- function(sizeProp, minFL, maxFL, alpha, beta) {
  year <- sample(min(sizeProp$Year):max(sizeProp$Year), 1) # Random year sample
  females <- sample(2075:2450, 1) * 0.56  # Spawner estimate * 0.56 for females
  
  sizeProp %>%
    filter(Year == year) %>%
    mutate(females = round(prop * females, 0)) %>%
    uncount(females) %>%  # Expand data for individual fish
    mutate(
      FL = if_else(Age %in% c(5, 6), 720, sample(minFL:maxFL, n(), replace = TRUE)),  # Sample FL
      Fecundity = round((alpha * FL)^beta, 0)  # Compute fecundity
    )
}

# Run simulations and store individual fecundities
egg_data <- map_dfr(1:n_sim, ~simulate_egg_production(sizeProp, minFL, maxFL, alpha, beta))

# Calculate mean and SD for fecundity per fish
fecundity_stats <- egg_data %>%
  summarize(mean_fecundity = mean(Fecundity), sd_fecundity = sd(Fecundity))

# Print results
print(fecundity_stats)

fecundity_stats |> summarize(Mean = mean(mean_fecundity), SD = sd(sd_fecundity))

###############
##STAN MODEL
##############

# Define the data list based on priors
stan_data <- list(
  N = 1,  # Example: Number of years (adjust as needed)
  
  # Number of redds
  mean_Rd = 1267,  
  sd_Rd = 105,  
  
  # Fecundity
  mean_Fe = 11675,  
  sd_Fe = 11702,
  
  # Egg-to-fry survival probability
  mean_SEg = 0.34,
  sd_SEg = 0.17,
  
  # Fry-to-age 0+ parr survival probability
  mean_SFr = 0.26,
  sd_SFr = 0.11,
  
  # 0+ parr to age 1+ survival and remain
  mean_SR_01 = 0.33,
  sd_SR_01 = 0.2,
  
  # Age 1+ to age 2+ survival and remain
  mean_SR_12 = 0.34,
  sd_SR_12 = 0.09,
  
  # Probability of smolting and surviving age 2+
  mean_SmS = 0.16,
  sd_SmS = 0.04
)

# Compile the model
mod <- cmdstan_model("bayesian_jpe/adult_steelhead_jpe_abbrev.stan")

# Run the model
fit <- mod$sample(
  data = stan_data,
  seed = 123,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 2000,
  iter_sampling = 4000
)

# Print summary of posterior samples
print(fit$summary())

# Extract and plot results
posterior_samples <- fit$draws(format = "df")

ggplot(posterior_samples, aes(x = F_e)) +
  geom_histogram(bins = 30, fill = "blue", alpha = 0.6) +
  labs(title = "Posterior Distribution of Fecundity", x = "Fecundity", y = "Frequency") +
  theme_minimal()
