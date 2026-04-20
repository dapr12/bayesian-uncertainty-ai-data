// =============================================================================
// selection_model.stan
// Bayesian Partial Identification Framework for AI-Augmented Social Science Data
//
// Pérez Ruiz, D.A. (2026). Bayesian Uncertainty Quantification for
// AI-Augmented Social Science Data. RSS: Data Science and Artificial
// Intelligence. Special Issue: Uncertainty in the Era of AI.
//
// Model description:
//   Combines a probability sample (unbiased, known design) with a
//   non-probability sample (unknown AI-mediated selection mechanism).
//   The selection bias parameter delta ~ N(0, sigma_delta^2) encodes
//   uncertainty about the AI selection mechanism.
//
// Computational approximation:
//   The conditional expectation E[Z_i | R_i = 1] is approximated as
//   delta * phi(0) = delta * 0.3989 (first-order Taylor approximation,
//   exact at delta = 0). The latent covariate Z_i is integrated out
//   analytically, avoiding the funnel geometry that causes 100% HMC
//   divergences when Z_i is sampled explicitly.
//
// For full documentation see the paper and README.md.
// =============================================================================

data {
  int<lower=1> n_p;              // probability sample size
  int<lower=1> n_np;             // non-probability sample size
  vector[n_p]  Y_p;              // outcomes, probability sample
  vector[n_p]  X_p;              // covariates, probability sample
  vector[n_np] Y_np;             // outcomes, non-probability sample
  vector[n_np] X_np;             // covariates, non-probability sample
  real<lower=0> sigma_delta;     // prior SD on selection bias parameter
}

parameters {
  real beta_0;                   // intercept (= population mean mu)
  real beta_1;                   // effect of observed covariate X
  real<lower=0> sigma_y;         // outcome noise SD
  real delta;                    // selection bias parameter
                                 // delta = 1: 1 SD increase in Z
                                 // doubles log-odds of selection
}

transformed parameters {
  // Marginal SD of Y integrating out Z ~ N(0,1) with loading gamma = 1
  // Var(Y | X) = beta_1^2 * Var(X) + gamma^2 * Var(Z) + sigma_y^2
  //            = beta_1^2 + 1 + sigma_y^2
  real sigma_marg = sqrt(square(beta_1) + 1.0 + square(sigma_y));

  // First-order approximation: E[Z | R=1] ~ delta * phi(0) = delta * 0.3989
  // This is exact at delta = 0 and accurate in the neighbourhood of delta = 0,
  // where the model is locally identifiable from the probability sample.
  real bias_np = delta * 0.3989;
}

model {
  // ------- Priors -------
  beta_0  ~ normal(0, 5);
  beta_1  ~ normal(0, 2);
  sigma_y ~ exponential(1);

  // KEY prior: encodes researcher uncertainty about AI selection bias.
  // Small sigma_delta: confident the AI source is approximately representative.
  // Large sigma_delta: agnostic or sceptical about AI selection properties.
  delta ~ normal(0, sigma_delta);

  // ------- Probability sample likelihood -------
  // Unbiased marginal likelihood (Z integrated out analytically).
  // No selection correction needed: inclusion probabilities are known by design.
  Y_p ~ normal(beta_0 + beta_1 * X_p, sigma_marg);

  // ------- Non-probability sample likelihood -------
  // Bias-corrected marginal likelihood.
  // The bias_np term shifts the conditional mean to account for
  // selection on the unobserved latent covariate Z_i.
  Y_np ~ normal(beta_0 + beta_1 * X_np + bias_np, sigma_marg);
}

generated quantities {
  // Population mean: E[Y] = beta_0 since E[X] = E[Z] = 0 in the population.
  real mu = beta_0;
}
