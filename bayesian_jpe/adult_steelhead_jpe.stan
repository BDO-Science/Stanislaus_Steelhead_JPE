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
  real<lower=0, upper=1> mean_SPr_01; // mean survival probability from egg to age 0+ parr
  real<lower=0, upper=1> sd_SPr_01;  // standard deviation for SPr_12
  real<lower=0, upper=1> mean_PSm; // mean probability of smolting at age 1+
  real<lower=0, upper=1> sd_PSm;  // standard deviation for PSm
  real<lower=0, upper=1> mean_SPr_12; // mean survival probability from age 0+ to age 1+
  real<lower=0, upper=1> sd_SPr_12;  // standard deviation for SPr_12
  real<lower=0, upper=1> mean_PSm; // mean probability of smolting at age 2+
  real<lower=0, upper=1> sd_PSm;  // standard deviation for PSm
  real<lower=0, upper=1> mean_SO;  // mean survival probability of outmigration at age 2+
  real<lower=0, upper=1> sd_SO;   // standard deviation for SO
}

parameters {
  real<lower=0> Fe;    // fecundity
  vector<lower=0>[N] Rd; //Number of redds
  real<lower=0, upper=1> SEg;   // egg-to-fry survival probability
  real<lower=0, upper=1> SFr;   // fry-to-age 0+ parr survival probability
  real<lower=0, upper=1> SPr_01;  // survival probability from egg to age 0+ parr
  real<lower=0, upper=1> PSm;  // probability of smolting at age 1+
  real<lower=0, upper=1> SPr_12;  // survival probability from age 0+ to age 1+
  real<lower=0, upper=1> PSm;  // probability of smolting at age 2+
  real<lower=0, upper=1> SO;   // survival probability of outmigration at age 2+
  real<lower=0> sigma;  // standard deviation for the likelihood
}

model {
  // Priors for parameters based on provided estimates

  // Fecundity (Hodge et al. 2016)
  //F_e ~ normal(

  // Female Spawner Abundance (Zeug et al. 2024)
  //Rd ~ normal(1267, 105); // Mean = (1162+1372)/2, SD estimated from range

  // Egg to Fry Survival (Zeug et al. 2024)
  //SEg ~ beta(6.3, 12.2); // Mean = 0.34, SD = 0.17 (approximated)

  // Fry to Parr Survival (Baxter 1997)
  //SFr ~ beta(4.4, 12.6); // Mean = 0.26, SD = 0.11 (approximated)

  // Parr to 1+ Survival and Remain (Cuniack et al. 1998)
  //SPr_01 ~ beta(6.6, 13.4); // Mean = 0.33, SD = 0.2 (approximated)

  // 1+ to 2+ Survival and Remain (Cuniack et al. 1998)
  //SPr_12 ~ beta(6.8, 13.2); // Mean = 0.34, SD = 0.09 (approximated)

  // 2+ smolting and surviving (Zeug et al. 2024)
  //SO ~ beta(4.8, 25.2); // Mean = 0.16, SD = 0.04 (approximated)

  // Process error
  //sigma ~ normal(0, 1);

  // Data-informed priors for each parameter
 // mean_F_e ~ normal(F_e, sd_F_e); 
 mean_Rd ~ normal(Rd, sd_Rd); 
  mean_SEg ~ normal(SEg, sd_SEg); 
  mean_SFr ~ normal(SFr, sd_SFr); 
  mean_SPr_01 ~ normal(SPr_12, sd_SPr_12); 
  mean_PSm ~ normal(PSm, sd_PSm); 
  mean_SPr_12 ~ normal(SPr_12, sd_SPr_12); 
  mean_PSm ~ normal(PSm, sd_PSm); 
  mean_SO ~ normal(SO, sd_SO); 

  // Likelihood
  for (i in 1:N) {
    real JPE_i = Rd[i] * F_e * SEg * SFr * SPr_01 * (1 - PSm) * SPr_12 * PSm * SO;
    target += normal_lpdf(JPE_i | Rd[i] * F_e * SEg * SFr * SPr_01 * (1 - PSm) * SPr_12 * PSm * SO, sigma);
  }
}

generated quantities {
  vector[N] JPE;
  for (i in 1:N) {
    JPE[i] = Rd[i] * F_e * SEg * SFr * SPr_01 * (1 - PSm) * SPr_12 * PSm * SO;
  }
}




