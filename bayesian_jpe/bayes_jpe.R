##############################
#Juvenile JPE Simulation
##############################

# Load libraries
library(cmdstanr)
library(tidyverse)

check_cmdstan_toolchain()
set_cmdstan_path("C:/CMDSTAN/cmdstan-2.34.1")

# Define simulation parameters
set.seed(123)
N <- 5
mean_P_r1 <- rnorm(N, mean = 1000, sd = 200)  # Simulated means for P_r1
sd_P_r1 <- runif(N, min = 50, max = 100)     # Simulated standard deviations for P_r1
mean_SP_r1 <- 0.6
sd_SP_r1 <- 0.1
mean_PS_m2 <- 0.5
sd_PS_m2 <- 0.1
mean_SO_2 <- 0.7
sd_SO_2 <- 0.1

# Compile the model
model <- cmdstan_model("juvenile_steelhead_jpe.stan")

# Create the data list
data_list <- list(
  N = N,
  mean_P_r1 = mean_P_r1,
  sd_P_r1 = sd_P_r1,
  mean_SP_r1 = mean_SP_r1,
  sd_SP_r1 = sd_SP_r1,
  mean_PS_m2 = mean_PS_m2,
  sd_PS_m2 = sd_PS_m2,
  mean_SO_2 = mean_SO_2,
  sd_SO_2 = sd_SO_2
)

# Fit the model
fit <- model$sample(
  data = data_list,
  seed = 123,
  chains = 4,
  parallel_chains = 4,
  iter_sampling = 2000,
  iter_warmup = 1000
)

# Summarize results
fit$summary() |> print()

##############################
#Adult JPE Simulation
##############################

# Define simulation parameters
set.seed(123)
N <- 5
mean_R_d <- rnorm(N, mean = 5000, sd = 1000)  # Simulated means for R_d
sd_R_d <- runif(N, min = 100, max = 500)     # Simulated standard deviations for R_d
mean_F_e <- 5000
sd_F_e <- 500
mean_SE_g <- 0.4
sd_SE_g <- 0.1
mean_SF_r <- 0.3
sd_SF_r <- 0.1
mean_SP_r0 <- 0.2
sd_SP_r0 <- 0.1
mean_PS_m1 <- 0.4
sd_PS_m1 <- 0.1
mean_SP_r1 <- 0.6
sd_SP_r1 <- 0.1
mean_PS_m2 <- 0.6
sd_PS_m2 <- 0.1
mean_SO_2 <- 0.7
sd_SO_2 <- 0.1

# Compile the model
model <- cmdstan_model("adult_steelhead_jpe.stan")

# Prepare the updated data list
data_list <- list(
  N = N,
  mean_R_d = mean_R_d,
  sd_R_d = sd_R_d,
  mean_F_e = mean_F_e,
  sd_F_e = sd_F_e,
  mean_SE_g = mean_SE_g,
  sd_SE_g = sd_SE_g,
  mean_SF_r = mean_SF_r,
  sd_SF_r = sd_SF_r,
  mean_SP_r0 = mean_SP_r0,
  sd_SP_r0 = sd_SP_r0,
  mean_PS_m1 = mean_PS_m1,
  sd_PS_m1 = sd_PS_m1,
  mean_SP_r1 = mean_SP_r1,
  sd_SP_r1 = sd_SP_r1,
  mean_PS_m2 = mean_PS_m2,
  sd_PS_m2 = sd_PS_m2,
  mean_SO_2 = mean_SO_2,
  sd_SO_2 = sd_SO_2
)

fit <- model$sample(
  data = data_list,
  seed = 123,
  chains = 4,
  parallel_chains = 4,
  iter_sampling = 2000,
  iter_warmup = 1000,
  init = function() list(
    F_e = 5000,
    SE_g = 0.3,
    SF_r = 0.2,
    SP_r0 = 0.1,
    PS_m1 = 0.5,
    SP_r1 = 0.4,
    PS_m2 = 0.5,
    SO_2 = 0.6
  )
)

# Summarize results
fit$summary() |> print()
