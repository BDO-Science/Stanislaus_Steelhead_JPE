data {
  int<lower=0> N; // number of years
  vector[N] mean_P_r1; // mean estimated number of age 1+ parr for each year
  vector[N] sd_P_r1;  // standard deviation for P_r1
  real<lower=0> mean_SP_r1; // mean survival probability to reach age 2+
  real<lower=0> sd_SP_r1; // standard deviation for SP_r1
  real<lower=0> mean_PS_m2; // mean probability of smolting at age 2+
  real<lower=0> sd_PS_m2; // standard deviation for PS_m2
  real<lower=0> mean_SO_2;  // mean survival probability of outmigration at age 2+
  real<lower=0> sd_SO_2;  // standard deviation for SO_2
}

parameters {
  vector<lower=0>[N] P_r1; // estimated number of age 1+ parr for each year
  real<lower=0, upper=1> SP_r1; // survival probability to reach age 2+
  real<lower=0, upper=1> PS_m2; // probability of smolting at age 2+
  real<lower=0, upper=1> SO_2;  // survival probability of outmigration at age 2+
  real<lower=0> sigma; // standard deviation for the likelihood
}

model {
  // Priors
  P_r1 ~ normal(6000, 1000);
  SP_r1 ~ beta(3, 7);   // Mean = 0.3, SD ~ 0.13
  PS_m2 ~ beta(2.5, 2.5);  // Mean = 0.5, SD ~ 0.22
  SO_2 ~ beta(2, 8);   // Mean = 0.2, SD ~ 0.12
  sigma ~ normal(0, 1);

  // Data-informed priors
  mean_P_r1 ~ normal(P_r1, sd_P_r1); 
  mean_SP_r1 ~ normal(SP_r1, sd_SP_r1);
  mean_PS_m2 ~ normal(PS_m2, sd_PS_m2);
  mean_SO_2 ~ normal(SO_2, sd_SO_2);

  // Likelihood
  for (i in 1:N) {
    real JPE_i = P_r1[i] * SP_r1 * PS_m2 * SO_2;
    target += normal_lpdf(JPE_i | P_r1[i] * SP_r1 * PS_m2 * SO_2, sigma);
  }
}


generated quantities {
  vector[N] JPE;
  for (i in 1:N) {
    JPE[i] = P_r1[i] * SP_r1 * PS_m2 * SO_2;
  }
}


