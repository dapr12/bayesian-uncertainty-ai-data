# =============================================================================
# 01_simulation.R
# Main simulation study
#
# Pérez Ruiz, D.A. (2026). Bayesian Uncertainty Quantification for
# AI-Augmented Social Science Data. RSS: Data Science and Artificial
# Intelligence. Special Issue: Uncertainty in the Era of AI.
#
# Reproduces: Table 1, coverage summary
#
# Usage:
#   setwd("path/to/repo")
#   source("R/01_simulation.R")
#
# For a quick test run, set n_sims <- 20 below.
# Full run (n_sims = 200) takes approximately 1-2 hours.
# =============================================================================

library(cmdstanr)
library(posterior)
library(tidyverse)

set.seed(42)

# =============================================================================
# USER SETTINGS
# =============================================================================
n_sims <- 200          # number of simulation replicates (set to 20 for quick test)
n_chains <- 4          # HMC chains per model fit
n_iter   <- 1000       # post-warmup iterations per chain
n_warmup <- 1000       # warmup iterations per chain

# =============================================================================
# SECTION 1: HELPER FUNCTIONS
# =============================================================================

#' Draw probability and non-probability samples from population
draw_samples <- function(pop, delta, n_p = 200, n_np = 2000,
                         alpha_0 = -2.0, alpha_1 = 0.3) {
  s_p      <- pop[sample(nrow(pop), n_p), ]
  prob_np  <- plogis(alpha_0 + alpha_1 * pop$X + delta * pop$Z)
  selected <- rbinom(nrow(pop), 1, prob_np)
  s_np_all <- pop[selected == 1, ]
  s_np     <- s_np_all[seq_len(min(n_np, nrow(s_np_all))), ]
  list(s_p = s_p, s_np = s_np)
}

#' Naive pooling estimator
naive_estimate <- function(s_p, s_np, mu_true) {
  dat    <- rbind(s_p, s_np)
  mu_hat <- mean(dat$Y)
  se_hat <- sd(dat$Y) / sqrt(nrow(dat))
  list(mean     = mu_hat,
       lower    = mu_hat - 1.96 * se_hat,
       upper    = mu_hat + 1.96 * se_hat,
       width    = 2 * 1.96 * se_hat,
       coverage = as.numeric(mu_hat - 1.96 * se_hat <= mu_true &
                             mu_hat + 1.96 * se_hat >= mu_true))
}

#' Horvitz-Thompson estimator (probability sample only)
ht_estimate <- function(s_p, mu_true) {
  mu_hat <- mean(s_p$Y)
  se_hat <- sd(s_p$Y) / sqrt(nrow(s_p))
  list(mean     = mu_hat,
       lower    = mu_hat - 1.96 * se_hat,
       upper    = mu_hat + 1.96 * se_hat,
       width    = 2 * 1.96 * se_hat,
       coverage = as.numeric(mu_hat - 1.96 * se_hat <= mu_true &
                             mu_hat + 1.96 * se_hat >= mu_true))
}

#' Fit Stan model and extract posterior summary for mu
fit_bayes <- function(stan_model, s_p, s_np, sigma_delta,
                      seed = 42, mu_true) {
  x_shift   <- mean(s_np$X) - mean(s_p$X)
  stan_data <- list(
    n_p         = nrow(s_p),
    n_np        = nrow(s_np),
    Y_p         = s_p$Y,
    X_p         = s_p$X,
    Y_np        = s_np$Y,
    X_np        = s_np$X,
    sigma_delta = sigma_delta
  )
  fit <- stan_model$sample(
    data            = stan_data,
    iter_warmup     = n_warmup,
    iter_sampling   = n_iter,
    chains          = n_chains,
    parallel_chains = n_chains,
    seed            = seed,
    adapt_delta     = 0.95,
    refresh         = 0,
    show_messages   = FALSE
  )
  mu_draws <- as_draws_df(fit$draws())$mu
  lower    <- quantile(mu_draws, 0.025)
  upper    <- quantile(mu_draws, 0.975)
  list(mean     = mean(mu_draws),
       lower    = as.numeric(lower),
       upper    = as.numeric(upper),
       width    = as.numeric(upper - lower),
       coverage = as.numeric(lower <= mu_true & upper >= mu_true))
}

# =============================================================================
# SECTION 2: POPULATION DATA GENERATION
# =============================================================================
cat("\n--- Generating population ---\n")

N      <- 10000
beta_0 <- 1.0
beta_1 <- 0.5
gamma  <- 0.8
sigma  <- 1.0

pop <- data.frame(X = rnorm(N), Z = rnorm(N)) |>
  transform(Y = beta_0 + beta_1 * X + gamma * Z + rnorm(N, 0, sigma))

mu_true <- mean(pop$Y)
cat("True population mean mu:", round(mu_true, 4), "\n")

# =============================================================================
# SECTION 3: COMPILE STAN MODEL
# =============================================================================
cat("\n--- Compiling Stan model ---\n")
stan_model_compiled <- cmdstan_model("stan/selection_model.stan")
cat("Stan model compiled successfully.\n")

# =============================================================================
# SECTION 4: MAIN SIMULATION — ONE DRAW PER SCENARIO
# =============================================================================
cat("\n--- Running main simulation fits ---\n")

delta_scenarios  <- c(mild = 0.5, moderate = 1.5, severe = 3.0)
sigma_delta_vals <- c(0.5, 3.0)
main_results     <- list()

for (scenario_name in names(delta_scenarios)) {

  delta_true <- delta_scenarios[scenario_name]
  cat("\nScenario:", scenario_name, "| delta =", delta_true, "\n")

  samples <- draw_samples(pop, delta_true)
  s_p     <- samples$s_p
  s_np    <- samples$s_np
  cat("  Prob sample n =", nrow(s_p),
      "| Non-prob sample n =", nrow(s_np), "\n")

  # Naive pooling
  naive_res <- naive_estimate(s_p, s_np, mu_true)
  main_results[[paste0(scenario_name, "_naive")]] <- data.frame(
    scenario = scenario_name, method = "Naive pooling",
    sigma_delta = NA_real_, mean = naive_res$mean,
    lower = naive_res$lower, upper = naive_res$upper,
    width = naive_res$width, coverage = naive_res$coverage)

  # Horvitz-Thompson
  ht_res <- ht_estimate(s_p, mu_true)
  main_results[[paste0(scenario_name, "_ht")]] <- data.frame(
    scenario = scenario_name, method = "Horvitz-Thompson",
    sigma_delta = NA_real_, mean = ht_res$mean,
    lower = ht_res$lower, upper = ht_res$upper,
    width = ht_res$width, coverage = ht_res$coverage)

  # Bayesian framework
  for (sd_val in sigma_delta_vals) {
    cat("  Fitting Bayesian model | sigma_delta =", sd_val, "...\n")
    bayes_res <- fit_bayes(stan_model_compiled, s_p, s_np,
                           sigma_delta = sd_val, seed = 42,
                           mu_true = mu_true)
    main_results[[paste0(scenario_name, "_bayes_", sd_val)]] <- data.frame(
      scenario = scenario_name,
      method   = paste0("Bayesian (sigma_delta=", sd_val, ")"),
      sigma_delta = sd_val, mean = bayes_res$mean,
      lower = bayes_res$lower, upper = bayes_res$upper,
      width = bayes_res$width, coverage = bayes_res$coverage)
  }
}

main_df <- do.call(rbind, main_results)
rownames(main_df) <- NULL
cat("\n--- Main simulation results ---\n")
print(main_df %>% mutate(across(where(is.numeric), ~round(.x, 3))))

# =============================================================================
# SECTION 5: COVERAGE SIMULATION — n_sims REPLICATES
# =============================================================================
cat("\n--- Running coverage simulation (", n_sims, "replicates) ---\n")
cat("    Tip: set n_sims <- 20 at the top for a quick test run.\n\n")

coverage_records <- list()

for (i in seq_len(n_sims)) {

  if (i %% 10 == 0) cat("  Replicate", i, "of", n_sims, "\n")

  for (scenario_name in names(delta_scenarios)) {

    samples   <- draw_samples(pop, delta_scenarios[scenario_name])
    s_p       <- samples$s_p
    s_np      <- samples$s_np

    naive_res <- naive_estimate(s_p, s_np, mu_true)
    ht_res    <- ht_estimate(s_p, mu_true)
    bayes_05  <- fit_bayes(stan_model_compiled, s_p, s_np,
                           sigma_delta = 0.5, seed = i, mu_true = mu_true)
    bayes_30  <- fit_bayes(stan_model_compiled, s_p, s_np,
                           sigma_delta = 3.0, seed = i, mu_true = mu_true)

    coverage_records[[length(coverage_records) + 1]] <- data.frame(
      replicate        = i,
      scenario         = scenario_name,
      covered_naive    = naive_res$coverage,
      covered_ht       = ht_res$coverage,
      covered_bayes_05 = bayes_05$coverage,
      covered_bayes_30 = bayes_30$coverage,
      width_naive      = naive_res$width,
      width_ht         = ht_res$width,
      width_bayes_05   = bayes_05$width,
      width_bayes_30   = bayes_30$width
    )
  }
}

coverage_df <- do.call(rbind, coverage_records)

coverage_summary <- coverage_df %>%
  group_by(scenario) %>%
  summarise(
    coverage_naive    = mean(covered_naive),
    coverage_ht       = mean(covered_ht),
    coverage_bayes_05 = mean(covered_bayes_05),
    coverage_bayes_30 = mean(covered_bayes_30),
    se_naive          = sqrt(mean(covered_naive) *
                             (1 - mean(covered_naive)) / n()),
    se_ht             = sqrt(mean(covered_ht) *
                             (1 - mean(covered_ht)) / n()),
    se_bayes_05       = sqrt(mean(covered_bayes_05) *
                             (1 - mean(covered_bayes_05)) / n()),
    se_bayes_30       = sqrt(mean(covered_bayes_30) *
                             (1 - mean(covered_bayes_30)) / n()),
    width_naive       = mean(width_naive),
    width_ht          = mean(width_ht),
    width_bayes_05    = mean(width_bayes_05),
    width_bayes_30    = mean(width_bayes_30),
    .groups = "drop"
  )

cat("\n--- Coverage summary ---\n")
print(coverage_summary %>% mutate(across(where(is.numeric), ~round(.x, 3))))

# =============================================================================
# SECTION 6: SAVE RESULTS
# =============================================================================
cat("\n--- Saving results ---\n")
dir.create("data/simulation_results", recursive = TRUE, showWarnings = FALSE)

saveRDS(main_df,          "data/simulation_results/results_main.rds")
saveRDS(coverage_df,      "data/simulation_results/results_coverage_full.rds")
saveRDS(coverage_summary, "data/simulation_results/results_coverage_summary.rds")

write.csv(main_df,          "data/simulation_results/results_main.csv",
          row.names = FALSE)
write.csv(coverage_summary, "data/simulation_results/results_coverage_summary.csv",
          row.names = FALSE)

cat("Saved to data/simulation_results/\n")
cat("\n=== Simulation complete ===\n")
cat("Next: run R/02_sensitivity.R\n")
