//BT_estimation.stan
data {
  int<lower=0> N;                    		// Cells Number
  int<lower=0> nr;                   		// Number of Raters
  int<lower=0,upper=1>  BCM[N*nr];   		// Boolean Classification Matrices (BCMs)
  real<lower=0,upper=1> PUS_BMM[N];      	// Mean Raters' BMM prior (Unique Starting Belonging Measure Matrix)
  real<lower=0,upper=1> sigma[N*nr]; 		// Standard Deviation Raters' BMM prior
}

parameters {
  real<lower=0,upper=1> R_BMM[N*nr]; 		// Raters' Belonging Measure Matrix (R_BMM)
  real<lower=0,upper=1> US_BMM[N];    		// Unique Starting Belonging Measure Matrix (US_BMM)
  real<lower=0.1,upper=1> BT[nr];   		// Raters' Belonging Thresholds (BTs)
}


model {
  real lambda;
  int cntr;

  for (n in 1:N) {
      for (r in 1:nr) {

      cntr = ((r-1)*N)+n;
      US_BMM[n]   ~ normal(PUS_BMM[n],.0005);
      R_BMM[cntr] ~ normal( US_BMM[n], sigma[cntr])T[0,1];  // Raters' BMM prior

      if (R_BMM[cntr]>BT[r]) {
         lambda = .9999;
         } else {
         lambda = .0001;
      }
      BCM[cntr] ~ bernoulli(lambda);
    }
  }
}
