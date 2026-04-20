# =============================================================================
# 02_sensitivity.R
# Sensitivity analysis: varying the prior on delta
#
# Pérez Ruiz, D.A. (2026). Bayesian Uncertainty Quantification for
# AI-Augmented Social Science Data. RSS: Data Science and Artificial
# Intelligence. Special Issue: Uncertainty in the Era of AI.
#
# Reproduces: Section 4.3 (sensitivity analysis paragraph)
#
# Requires: 01_simulation.R to have been run first (loads pop and
#           stan_model_compiled from the environment).
# If running standalone, set run_standalone <- TRUE below.
# =============================================================================

library(cmdstanr)
library(posterior)
library(tidyverse)

run_standalone <- FALSE   # set TRUE to run without 01_simulation.R

if (run_standalone) {
  set.seed(42)
  N <- 10000; beta_0 <- 1.0; beta_1 <- 0.5; gamma <- 0.8; sigma <- 1.0
  pop <- data.frame(X = rnorm(N), Z = rnorm(N)) |>
    transform(Y = beta_0 + beta_1 * X + gamma * Z + rnorm(N, 0, sigma))
  mu_true <- mean(pop$Y)
  stan_model_compiled <- cmdstan_model("stan/selection_model.stan")

  draw_samples <- function(pop, delta, n_p = 200, n_np = 2000,
                           alpha_0 = -2.0, alpha_1 = 0.3) {
    s_p      <- pop[sample(nrow(pop), n_p), ]
    prob_np  <- plogis(alpha_0 + alpha_1 * pop$X + delta * pop$Z)
    selected <- rbinom(nrow(pop), 1, prob_np)
    s_np_all <- pop[selected == 1, ]
    s_np     <- s_np_all[seq_len(min(n_np, nrow(s_np_all))), ]
    list(s_p = s_p, s_np = s_np)
  }
}

cat("\n--- Running sensitivity analysis (moderate scenario, delta = 1.5) ---\n")

sigma_delta_grid <- c(0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5)
sensitivity_list <- list()
samples_mod      <- draw_samples(pop, delta = 1.5)

for (sd_val in sigma_delta_grid) {

  cat("  sigma_delta =", sd_val, "\n")

  stan_data <- list(
    n_p         = nrow(samples_mod$s_p),
    n_np        = nrow(samples_mod$s_np),
    Y_p         = samples_mod$s_p$Y,
    X_p         = samples_mod$s_p$X,
    Y_np        = samples_mod$s_np$Y,
    X_np        = samples_mod$s_np$X,
    sigma_delta = sd_val
  )

  fit <- stan_model_compiled$sample(
    data            = stan_data,
    iter_warmup     = 1000,
    iter_sampling   = 1000,
    chains          = 4,
    parallel_chains = 4,
    seed            = 42,
    adapt_delta     = 0.95,
    refresh         = 0,
    show_messages   = FALSE
  )

  mu_draws <- as_draws_df(fit$draws())$mu
  lower    <- quantile(mu_draws, 0.025)
  upper    <- quantile(mu_draws, 0.975)

  sensitivity_list[[as.character(sd_val)]] <- data.frame(
    sigma_delta = sd_val,
    mean        = mean(mu_draws),
    lower       = as.numeric(lower),
    upper       = as.numeric(upper),
    width       = as.numeric(upper - lower)
  )
}

sensitivity_df <- do.call(rbind, sensitivity_list)
rownames(sensitivity_df) <- NULL

cat("\n--- Sensitivity analysis results ---\n")
print(sensitivity_df %>% mutate(across(where(is.numeric), ~round(.x, 3))))

# Save
dir.create("data/simulation_results", recursive = TRUE, showWarnings = FALSE)
saveRDS(sensitivity_df, "data/simulation_results/results_sensitivity.rds")
write.csv(sensitivity_df, "data/simulation_results/results_sensitivity.csv",
          row.names = FALSE)
cat("Saved to data/simulation_results/results_sensitivity.rds\n")
cat("\n=== Sensitivity analysis complete ===\n")
cat("Next: run R/03_visualisation.R\n")
