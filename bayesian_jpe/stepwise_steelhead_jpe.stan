data {
  int<lower=0> N; // number of years
  vector[N] R_d; // estimated number of redds for each year
  real F_e; // average fecundity
  real SE_g; // egg-to-fry survival probability
  real SF_r; // fry-to-age 0+ parr survival probability
  vector[N] P_r1; // estimated number of age 1+ parr for each year
  real SP_r0; // survival probability from egg to age 0+ parr
  real PS_m1; // probability of smolting at age 1+
  real SP_r1; // survival probability from age 0+ to age 1+
  real PS_m2; // probability of smolting at age 2+
  real SO_2; // survival probability of outmigration at age 2+
}

parameters {
  real<lower=0> sigma_adult; // standard deviation for the adult likelihood
  real<lower=0> sigma_juvenile; // standard deviation for the juvenile likelihood
}

model {
  // Priors for the sigma parameters
  sigma_adult ~ normal(0, 1);
  sigma_juvenile ~ normal(0, 1);

  // Likelihood for the Adult Approach
  for (i in 1:N) {
    real log_JPE_adult = log(R_d[i] * F_e * SE_g * SF_r * SP_r0 * (1 - PS_m1) * SP_r1 * PS_m2 * SO_2);
    target += normal_lpdf(log_JPE_adult | 0, sigma_adult);
  }

  // Likelihood for the Juvenile Approach
  for (i in 1:N) {
    real log_JPE_juvenile = log(P_r1[i] * SP_r1 * PS_m2 * SO_2);
    target += normal_lpdf(log_JPE_juvenile | log_JPE_adult, sigma_juvenile);
  }
}

generated quantities {
  vector[N] JPE_stepwise;
  for (i in 1:N) {
    real JPE_adult = R_d[i] * F_e * SE_g * SF_r * SP_r0 * (1 - PS_m1) * SP_r1 * PS_m2 * SO_2;
    JPE_stepwise[i] = P_r1[i] * SP_r1 * PS_m2 * SO_2; // Updated with Juvenile Approach
  }
}


