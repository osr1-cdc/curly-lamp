// Multinomial logistic regression nowcast model for SARS-CoV-2 variant proportions
//
// This model fits a time-series structure to variant proportions using
// multinomial regression with random walk priors on time dynamics.
//
// It produces:
// 1. Smoothed proportions for all historical time periods (retrospective)
// 2. Predicted proportions for future time periods (nowcast)

data {
  int<lower=1> N_weeks;                    // Total number of time periods
  int<lower=1> N_variants;                 // Number of variants
  int<lower=1> N_fit;                      // Number of weeks with data (< N_weeks for prediction)

  matrix<lower=0>[N_weeks, N_variants] prop;  // Observed proportions
  matrix<lower=0>[N_weeks, N_variants] se;    // Standard errors of proportions

  real<lower=0, upper=1> credible_interval;   // For posterior intervals (e.g., 0.95)
}

parameters {
  // True log-odds (latent proportions)
  matrix[N_weeks, N_variants - 1] theta;     // (N_variants-1) for identifiability

  // Trend smoothness priors
  real<lower=0> tau_variant[N_variants];     // Per-variant time-series SD
}

transformed parameters {
  // Convert log-odds back to probability scale
  matrix<lower=0, upper=1>[N_weeks, N_variants] p;

  for (w in 1:N_weeks) {
    real denominator = 1;
    for (v in 1:(N_variants - 1)) {
      denominator += exp(theta[w, v]);
    }

    for (v in 1:(N_variants - 1)) {
      p[w, v] = exp(theta[w, v]) / denominator;
    }
    p[w, N_variants] = 1 / denominator;
  }
}

model {
  // Priors
  for (v in 1:(N_variants - 1)) {
    // Random walk prior on log-odds
    for (w in 2:N_fit) {
      theta[w, v] ~ normal(theta[w - 1, v], tau_variant[v]);
    }
    theta[1, v] ~ normal(0, 2);
  }

  // Hyperpriors on time-series SD
  for (v in 1:(N_variants - 1)) {
    tau_variant[v] ~ exponential(5);
  }

  // Likelihood
  for (w in 1:N_fit) {
    for (v in 1:N_variants) {
      // Observed proportions with measurement error
      prop[w, v] ~ normal(p[w, v], se[w, v]);
    }
  }
}

generated quantities {
  // Posterior predictions for future periods
  matrix[N_weeks, N_variants - 1] theta_pred;
  matrix<lower=0, upper=1>[N_weeks, N_variants] p_pred;

  // Copy fitted periods
  theta_pred[1:N_fit, ] = theta[1:N_fit, ];

  // Predict future periods via random walk continuation
  for (w in (N_fit + 1):N_weeks) {
    for (v in 1:(N_variants - 1)) {
      theta_pred[w, v] = normal_rng(theta_pred[w - 1, v], tau_variant[v]);
    }
  }

  // Convert predicted log-odds to probabilities
  for (w in 1:N_weeks) {
    real denom = 1;
    for (v in 1:(N_variants - 1)) {
      denom += exp(theta_pred[w, v]);
    }
    for (v in 1:(N_variants - 1)) {
      p_pred[w, v] = exp(theta_pred[w, v]) / denom;
    }
    p_pred[w, N_variants] = 1 / denom;
  }
}
