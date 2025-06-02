data {
  int<lower=0> N; // number of years
  vector[N] mean_Rd; // mean estimated number of redds for each year
  vector[N] sd_Rd;  // standard deviation for Rd
  real<lower=0> mean_Fe; // mean fecundity from literature
  real<lower=0> sd_Fe;   // standard deviation for fecundity
  real<lower=0, upper=1> mean_SEg; // mean egg-to-fry survival probability
  real<lower=0, upper=1> sd_SEg;  // standard deviation for SEg
  real<lower=0, upper=1> mean_SFr; // mean fry-to-age 0+ parr survival probability
  real<lower=0, upper=1> sd_SFr;  // standard deviation for SFr
  real<lower=0, upper=1> mean_SR_01; // mean probability of surviving and remaining from 0+ parr to age 1+
  real<lower=0, upper=1> sd_SR_01;  // standard deviation for SR_12
  real<lower=0, upper=1> mean_SR_12; // mean probability of surviving and remaining from age 1+ to age 2+
  real<lower=0, upper=1> sd_SR_12;  // standard deviation for SR_12
  real<lower=0, upper=1> mean_SmS; // mean probability of smolting and surviving age 2+
  real<lower=0, upper=1> sd_SmS;  // standard deviation for Sm
}

parameters {
  real<lower=0> Fe;    // fecundity
  vector<lower=0>[N] Rd; //Number of redds
  real<lower=0, upper=1> SEg;   // egg-to-fry survival probability
  real<lower=0, upper=1> SFr;   // fry-to-age 0+ parr survival probability
  real<lower=0, upper=1> SR_01;  // mean probability of surviving and remaining from 0+ parr to age 1+
  real<lower=0, upper=1> SR_12;  // mean probability of surviving and remaining from age 1+ to age 2+
  real<lower=0, upper=1> SmS;   // probability of smolting and surviving age 2+
  real<lower=0> sigma;  // standard deviation for the likelihood
}

model {
  // Priors for parameters based on provided estimates

  // Fecundity (Hodge et al. 2016)
  Fe ~ normal(11675, 11702)

  // Female Spawner Abundance (Zeug et al. 2024)
  //Rd ~ normal(1267, 105); // Mean = (1162+1372)/2, SD estimated from range

  // Egg to Fry Survival (Zeug et al. 2024)
  //SEg ~ beta(6.3, 12.2); // Mean = 0.34, SD = 0.17 (approximated)

  // Fry to Parr Survival (Baxter 1997)
  SFr ~ beta(4.4, 12.6); // Mean = 0.26, SD = 0.11 (approximated)

  // Parr to 1+ Survival and Remain (Cuniack et al. 1998)
  SR_01 ~ beta(6.6, 13.4); // Mean = 0.33, SD = 0.2 (approximated)

  // 1+ to 2+ Survival and Remain (Cuniack et al. 1998)
  SR_12 ~ beta(6.8, 13.2); // Mean = 0.34, SD = 0.09 (approximated)

  // 2+ smolting and surviving (Zeug et al. 2024)
  //SmS ~ beta(4.8, 25.2); // Mean = 0.16, SD = 0.04 (approximated)

  // Process error
  //sigma ~ normal(0, 1);

  for (i in 1:N_obs) {
   // Likelihood: Data from submodels are treated as observations with uncertainty
 mean_Fe[i] ~ normal(Fe, sd_Fe); 
 mean_Rd[i] ~ normal(Rd, sd_Rd); 
  mean_SEg[i] ~ normal(SEg, sd_SEg); 
  mean_SFr[i] ~ normal(SFr, sd_SFr); 
  mean_SR_01[i] ~ normal(SR_01, sd_SR_01); 
  mean_SR_12[i] ~ normal(SR_12, sd_SR_12); 
  mean_SmS[i] ~ normal(SmS, sd_SmS); 

  // Likelihood
    real JPE_i = Rd[i] * Fe * SEg * SFr * SR_01 * SR_12 * SmS;
    target += normal_lpdf(JPE_i | Rd[i] * Fe * SEg * SFr * SR_01 * SR_12 * SmS, sigma);
  }
}

generated quantities {
  vector[N] JPE;
  for (t in 1:N_pred) {
    JPE[i] = Rd[i] * Fe * SEg * SFr * SR_01 * SR_12 * SmS;
  }
}




