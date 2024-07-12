// R0_est_original_betabinom.stan
// Purpose: Stan Estimation of R0 using serological survey data. 
// Version: Estimate R0 simultaneously with estimation of serological profiles (fitting curve). 
//   -- This is an updated version of the original analysis using a fit of the sero curve and average age of infection.
// Author: Shaun Truelove, satruelove@gmail.com
// Date updated: 13 Jan 2020

functions {
  // Functional form for cumulative seroprevalence. Fits monotonically increasing curve.
  real x_a_funct(real a, real b0, real b1, real b2){
    return 1 - exp((b0/b1)*a*exp(-b1*a) + (1/b1)*((b0/b1)-b2)*(exp(-b1*a)-1) - b2*a);
  }
}

data {
  int<lower=1> C;                           // Number of countries
  int<lower=1> S;                           // Number of separate samples/serosurveys
  int<lower=1> SA;                          // Number of separate samples/serosurveys x age groups (n rows in data)
  int<lower=1> N[SA];                       // Number of indivs in each age grp in each serosurvey/sample
  int<lower=0> Y[SA];                       // Number of seropositive in each age in each serosurvey/sample
  vector<lower=0,upper=100>[SA] age;        // Age value of each age group
  vector<lower=0>[S] mu;                    // Birth rate for that country during the year of the survey
  int<lower=1> country[SA];                 // Country id number for each survey/age
  int<lower=1> survey[SA];                  // Survey id number for each survey/age
  int<lower=1> country_of_survey[S];        // Country id for each survey
  int<lower=1> age_pred_length;
  vector<lower=0>[age_pred_length] age_pred;
  matrix[S, age_pred_length] N_age;         // Matrix of population for 0-99 years for each survey 
  real<lower=0> M;                          // Duration of maternal antibodies
  //vector<lower=0>[S] L;                   // Birth rate for that country during the year of the survey
  vector<lower=0>[S] plus1;                 // Whether the average age equation needs a plus 1
  vector<lower=0>[S] r_ratio;               // Ratio of R to 1980 expected based on the contact matrix and age composition
}

parameters {
  // Serology Model Parameters
  vector<lower=0,upper=1>[S] b0;             // seroprev function parameter
  vector<lower=.00001,upper=1>[S] b1;        // seroprev function parameter
  vector<lower=0,upper=0.05>[S] b2;          // seroprev function parameter
  
  // R0 Parameters
  real<lower=1,upper=3> logR0_g;             // Global R0, log-tansformed
  real<lower=0,upper=2> sigma_g;            // Standard deviation of global R0 distribution
  real<lower=0,upper=1> sigma_c;            // Standard deviation of country R0 distribution
  vector<lower=0,upper=4>[C] logR0_c;        // log R0 for each country
  
  // Dispersion
  real<lower=0.0000001, upper=.3> gamma;
}      

transformed parameters {
  vector<lower=0,upper=1>[SA] p;                                // binom prob seropositive from each age group in each sample
  vector[S] R0_s;                                               // R0 for each survey/sample
  vector[S] Ave_age;                                            // estimated average age of infection
  vector[S] logR0_s;                                            // log R0 for each survey/sample
  vector[S] R0_s_adj;                                           // Adjusted R0 for each survey/sample
  
  
  // Fit functional form to seroprev data for each survey
  for (sa in 1:SA) {
    p[sa] = x_a_funct(age[sa], b0[survey[sa]], b1[survey[sa]], b2[survey[sa]]);
  } 
  
  // Loop through surveys and get estimates for each, assuming the samples values of b0:b2
  // -- Use fitted curves/betas to return probability of infection given age
  // -- calculate average age of infection
  // -- estimate R0 from average age of infection and birth rate
  for (s in 1:S) {
    vector[age_pred_length] p_inf_given_age = rep_vector(0.5, age_pred_length);
    vector[age_pred_length] n_new_inf_age = rep_vector(0, age_pred_length);
    real tot_new_inf = 0;
    
    for (a in 1:age_pred_length){
      p_inf_given_age[a] = x_a_funct(age_pred[a]+1, b0[s], b1[s], b2[s]) - x_a_funct(age_pred[a], b0[s], b1[s], b2[s]); 
      n_new_inf_age[a] = p_inf_given_age[a] * N_age[s, a];
    }
    tot_new_inf = sum(n_new_inf_age);
    Ave_age[s] = sum((n_new_inf_age / tot_new_inf) .* age_pred); // Average age of infection
  }
  
  // Vectorized calculations
  R0_s = plus1 + 1 ./ ((Ave_age - M) .* mu);                          // R0 
  //R0_s_2 = plus1 + L ./ Ave_age;
  
  // Adjust by the r_ratio to standardize to standard year
  R0_s_adj = R0_s .* r_ratio;
  logR0_s = log(R0_s_adj);                                            // log-trans R0s
} 
  
model { 
  // Weakly informative priors (uniform priors are implicitly implied)
  // b0 ~ uniform(0,1);
  // b1 ~ uniform(0.00001,1);
  // b2 ~ uniform(0,0.05);
  // gamma ~ uniform(0.0000001,1);  
  // moderately informative prior on global logR0
  logR0_g ~ normal(1.75, 0.5);                  // informative prior for mean of log(global R0)  -- R0=5
  logR0_c ~ normal(1.75, 1);                  // informative prior for mean of log(country R0)  -- R0=5
  sigma_g ~ normal(0, 1);                  // informative prior for mean of log(global R0)  -- R0=5
  sigma_c ~ normal(0, 1);                  // informative prior for mean of log(global R0)  -- R0=5

  target += normal_lpdf(logR0_s | logR0_c[country_of_survey], sigma_c); // Centered version
  target += normal_lpdf(logR0_c | logR0_g, sigma_g); // Centered version
  
  // likelihood of number seropositive per survey given number per survey (N) and binom prob seropositive(p) 
  target += beta_binomial_lpmf(Y | N, (1/gamma-1)*p, (1/gamma - 1)*(1-p)); // to reduce the number of parameters, putting gamma and p in this
  //target += binomial_lpmf(Y | N, p);   
}

generated quantities {
  real R0c[C];
  // Use the current param estimates to generate samples of R0c and R0s
  R0c = exp(normal_rng(logR0_c, sigma_c));
}

