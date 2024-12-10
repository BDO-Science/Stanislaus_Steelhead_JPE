data {
  int<lower=0> N; // number of years
  vector[N] mean_R_d; // mean estimated number of redds for each year
  vector[N] sd_R_d;  // standard deviation for R_d
  //real<lower=0> mean_F_e; // mean fecundity from literature
  //real<lower=0> sd_F_e;   // standard deviation for fecundity
  real<lower=0> mean_SE_g; // mean egg-to-fry survival probability
  real<lower=0> sd_SE_g;  // standard deviation for SE_g
  real<lower=0> mean_SF_r; // mean fry-to-age 0+ parr survival probability
  real<lower=0> sd_SF_r;  // standard deviation for SF_r
  real<lower=0> mean_SP_r0; // mean survival probability from egg to age 0+ parr
  real<lower=0> sd_SP_r0;  // standard deviation for SP_r0
  real<lower=0> mean_PS_m1; // mean probability of smolting at age 1+
  real<lower=0> sd_PS_m1;  // standard deviation for PS_m1
  real<lower=0> mean_SP_r1; // mean survival probability from age 0+ to age 1+
  real<lower=0> sd_SP_r1;  // standard deviation for SP_r1
  real<lower=0> mean_PS_m2; // mean probability of smolting at age 2+
  real<lower=0> sd_PS_m2;  // standard deviation for PS_m2
  real<lower=0> mean_SO_2;  // mean survival probability of outmigration at age 2+
  real<lower=0> sd_SO_2;   // standard deviation for SO_2
}

parameters {
  real<lower=0> F_e;    // fecundity
  vector<lower=0>[N] R_d; //Number of redds
  real<lower=0, upper=1> SE_g;   // egg-to-fry survival probability
  real<lower=0, upper=1> SF_r;   // fry-to-age 0+ parr survival probability
  real<lower=0, upper=1> SP_r0;  // survival probability from egg to age 0+ parr
  real<lower=0, upper=1> PS_m1;  // probability of smolting at age 1+
  real<lower=0, upper=1> SP_r1;  // survival probability from age 0+ to age 1+
  real<lower=0, upper=1> PS_m2;  // probability of smolting at age 2+
  real<lower=0, upper=1> SO_2;   // survival probability of outmigration at age 2+
  real<lower=0> sigma;  // standard deviation for the likelihood
}

model {
  // Priors for parameters with beta and normal distributions
  F_e ~ normal(5000, 500);
  R_d ~ normal(300, 100);
  SE_g ~ beta(3, 7);   // Mean = 0.3, SD ~ 0.13
  SF_r ~ beta(2, 8);   // Mean = 0.2, SD ~ 0.12
  SP_r0 ~ beta(2, 18); // Mean = 0.1, SD ~ 0.06
  PS_m1 ~ beta(2.5, 2.5);  // Mean = 0.5, SD ~ 0.22
  SP_r1 ~ beta(4, 6);   // Mean = 0.4, SD ~ 0.15
  PS_m2 ~ beta(2.5, 2.5);  // Mean = 0.5, SD ~ 0.22
  SO_2 ~ beta(6, 4);   // Mean = 0.6, SD ~ 0.15
  sigma ~ normal(0, 1); // Keep sigma as normal

  // Data-informed priors for each parameter
 // mean_F_e ~ normal(F_e, sd_F_e); 
 mean_R_d ~ normal(R_d, sd_R_d); 
  mean_SE_g ~ normal(SE_g, sd_SE_g); 
  mean_SF_r ~ normal(SF_r, sd_SF_r); 
  mean_SP_r0 ~ normal(SP_r0, sd_SP_r0); 
  mean_PS_m1 ~ normal(PS_m1, sd_PS_m1); 
  mean_SP_r1 ~ normal(SP_r1, sd_SP_r1); 
  mean_PS_m2 ~ normal(PS_m2, sd_PS_m2); 
  mean_SO_2 ~ normal(SO_2, sd_SO_2); 

  // Likelihood
  for (i in 1:N) {
    real JPE_i = R_d[i] * F_e * SE_g * SF_r * SP_r0 * (1 - PS_m1) * SP_r1 * PS_m2 * SO_2;
    target += normal_lpdf(JPE_i | R_d[i] * F_e * SE_g * SF_r * SP_r0 * (1 - PS_m1) * SP_r1 * PS_m2 * SO_2, sigma);
  }
}

generated quantities {
  vector[N] JPE;
  for (i in 1:N) {
    JPE[i] = R_d[i] * F_e * SE_g * SF_r * SP_r0 * (1 - PS_m1) * SP_r1 * PS_m2 * SO_2;
  }
}




