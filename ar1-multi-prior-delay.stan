// No mean-reversion; non-stationary
// Refactored projections for T + 1
// Prior on total abundance 2023
data{
  // Max obserbed time-steps
  int<lower=0> T;
  // Number of stocks modeled
  int<lower=0> S;
  // Vector of indices
  matrix[T, S] x;
  // Mean of stock for z-score
  vector[S] stock_mean;
  // SD of stock for z-score
  vector<lower=0>[S] stock_sd;
  // Number of future projections
  int<lower=0> n_proj;
  // Upriver Brights stock ID
  int<lower=0> increase_stock_id;
  // Duration of increase in years
  int<lower=0> increase_years;
  // Delay duration for mean increase
  int<lower=0> increase_delay;
  // Pct increase in mean
  real<lower=0> mean_increase;
  // 2023 total abundance prior mean
  real prior_mean;
  // 2023 total abundance prior sd
  real<lower=0> prior_sd;
}

parameters{
  // Initial value
  vector[S] eta_init;
  // Intercept
  vector[S] mu;
  // Logit-scaled autoregression parameter
  real phi_logit;
  // Global sigma
  real<lower=0> sigma_global;
  // Global sigma SD
  real<lower=0> sigma_sigma;
  // Standard deviation
  vector<lower=0>[S] sigma;
}

transformed parameters{
  // Vector of predicted log index
  matrix[(T + 1), S] eta;
  matrix<lower=0>[(T + 1), S] abs_eta;
  // Original-scale autoregression parameter
  real phi = inv_logit(phi_logit) * 4 - 2;
  for(s in 1:S){
    // AR(1) initial condition & process
    eta[1, s] = eta_init[s];
    abs_eta[1, s] = exp(eta_init[s] * stock_sd[s] + stock_mean[s]);
    for(t in 2:(T + 1)) {
      eta[t, s] = mu[s] + phi * x[(t-1), s];
      abs_eta[t, s] = exp(eta[t, s] * stock_sd[s] + stock_mean[s]);
    }
  }
  // New mean, after increase
  real new_mu = (log(exp(mu[increase_stock_id] * stock_sd[increase_stock_id] + stock_mean[increase_stock_id]) * (1.0 + mean_increase)) - stock_mean[increase_stock_id]) / stock_sd[increase_stock_id];
  // Total abundances during data phase
  vector<lower=0>[T + 1] total_abs_eta;
  for(t in 1:(T + 1)) total_abs_eta[t] = sum(abs_eta[t, ]);
}

model{
  // Priors (weakly informative)
  phi_logit ~ std_normal();
  for(s in 1:S){
    eta_init[s] ~ std_normal();
    mu[s] ~ std_normal();
    sigma[s] ~ normal(sigma_global, sigma_sigma);
  }
  sigma_global ~ lognormal(0, 1);
  sigma_sigma ~ std_normal();
  // Prior on 2023 total abundance
  total_abs_eta[T + 1] ~ normal(prior_mean, prior_sd);
  // Likelihood statement
  for(s in 1:S){
    for (t in 1:T) x[t, s] ~ normal(eta[t, s], sigma[s]);
  }
}

generated quantities{
  // Log likelihood calcs
  matrix[T, S] log_lik;
  for(s in 1:S){
    for(t in 1:T) log_lik[t, s] = normal_lpdf(x[t, s] | eta[t, s], sigma[s]);
  }
  /* Posterior predictive distribution */
  // Index prediction
  matrix[(T + 1), S] x_rep;
  // Log mean abundance index prediction
  matrix[(T + 1), S] eta_rep;
  // AR(1) posterior prediction process
  for(s in 1:S){
    eta_rep[1, s] = eta_init[s];
    x_rep[1, s] = normal_rng(eta_rep[1, s], sigma[s]);
    for(t in 2:(T + 1)){
      eta_rep[t, s] = mu[s] + phi * x_rep[t - 1, s];
      x_rep[t, s] = normal_rng(eta_rep[t, s], sigma[s]);
    }
  }
  // Future predictions
  matrix[n_proj, S] eta_rep_proj;
  matrix[n_proj, S] x_rep_proj;
  matrix[n_proj, S] raw_indices;
  vector[n_proj] mean_value;
  // Fill the vector of mean values for each time-step
  for(t in 1:increase_delay) mean_value[t] = mu[increase_stock_id];
  for(t in (increase_delay + 1):(increase_delay + increase_years)){
    mean_value[t] = mu[increase_stock_id] + (new_mu - mu[increase_stock_id]) * ((t - increase_delay) * 1.0 / increase_years);
  }
  for(t in (increase_delay + increase_years + 1):n_proj) mean_value[t] = new_mu;
  for(s in 1:S){
    eta_rep_proj[1, s] = mu[s] + phi * x_rep[(T + 1), s];
    x_rep_proj[1, s] = normal_rng(eta_rep_proj[1, s], sigma[s]);
    raw_indices[1, s] = x_rep_proj[1, s] * stock_sd[s] + stock_mean[s];
    for(t in 2:n_proj){
      if(s != increase_stock_id){eta_rep_proj[t, s] = mu[s] + phi * x_rep_proj[t - 1, s];} 
      else{eta_rep_proj[t, s] = mean_value[t] + phi * x_rep_proj[t - 1, s];}
      x_rep_proj[t, s] = normal_rng(eta_rep_proj[t, s], sigma[s]);
      raw_indices[t, s] = x_rep_proj[t, s] * stock_sd[s] + stock_mean[s];
    }
  }
  // Total abundance time-series calcs
  vector[n_proj] total_indices;
  for(t in 1:n_proj){
    total_indices[t] = sum(exp(raw_indices[t, 1:S]));
  }
}
