# =============================================================================
# 03_visualisation.R
# Uncertainty dashboard plots
#
# Pérez Ruiz, D.A. (2026). Bayesian Uncertainty Quantification for
# AI-Augmented Social Science Data. RSS: Data Science and Artificial
# Intelligence. Special Issue: Uncertainty in the Era of AI.
#
# Reproduces: Figure components described in Section 2.5
# Requires:   data/simulation_results/ to exist (run 01 and 02 first)
# =============================================================================

library(tidyverse)
library(patchwork)

# Load pre-computed results
main_df        <- readRDS("data/simulation_results/results_main.rds")
sensitivity_df <- readRDS("data/simulation_results/results_sensitivity.rds")
mu_true        <- 1.003

# Colour palette
method_colours <- c(
  "Naive pooling"                  = "#e74c3c",
  "Horvitz-Thompson"               = "#3498db",
  "Bayesian (sigma_delta=0.5)"     = "#27ae60",
  "Bayesian (sigma_delta=3)"       = "#1a5c35"
)

# =============================================================================
# PANEL 1: Layered interval display
# Inner band = HT interval (sampling uncertainty)
# Outer band = Bayesian interval (sampling + selection uncertainty)
# =============================================================================

# Construct layered data: pair each Bayesian estimate with its HT inner band
ht_bands <- main_df %>%
  filter(method == "Horvitz-Thompson") %>%
  select(scenario, ht_lower = lower, ht_upper = upper)

layered_df <- main_df %>%
  left_join(ht_bands, by = "scenario") %>%
  mutate(
    inner_lower = ifelse(grepl("Bayesian", method), ht_lower, lower),
    inner_upper = ifelse(grepl("Bayesian", method), ht_upper, upper),
    outer_lower = lower,
    outer_upper = upper,
    scenario    = factor(scenario, levels = c("mild", "moderate", "severe"))
  )

p1 <- layered_df %>%
  ggplot(aes(y = method, colour = method)) +
  geom_linerange(aes(xmin = outer_lower, xmax = outer_upper),
                 linewidth = 3, alpha = 0.20) +
  geom_linerange(aes(xmin = inner_lower, xmax = inner_upper),
                 linewidth = 3, alpha = 0.70) +
  geom_point(aes(x = mean), size = 3) +
  geom_vline(xintercept = mu_true, linetype = "dashed",
             colour = "black", linewidth = 0.7) +
  facet_wrap(~scenario, ncol = 1, labeller = label_both) +
  scale_colour_manual(values = method_colours, guide = "none") +
  scale_x_continuous(limits = c(0.5, 2.1)) +
  labs(title = "Panel 1: Layered interval display",
       subtitle = paste("Dark band = sampling uncertainty (inner, HT).",
                        "Light band = total uncertainty (outer, Bayesian).",
                        "Dashed = true μ."),
       x = "Estimate of population mean μ", y = NULL) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor  = element_blank(),
        plot.subtitle      = element_text(size = 9, colour = "grey40"),
        strip.text         = element_text(face = "bold"))

# =============================================================================
# PANEL 2: Sensitivity curve
# =============================================================================
sigma_delta_grid <- c(0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5)

p2 <- sensitivity_df %>%
  ggplot(aes(x = sigma_delta)) +
  geom_ribbon(aes(ymin = lower, ymax = upper),
              fill = "#27ae60", alpha = 0.20) +
  geom_line(aes(y = mean), colour = "#27ae60", linewidth = 1.2) +
  geom_point(aes(y = mean), colour = "#27ae60", size = 3) +
  geom_hline(yintercept = mu_true, linetype = "dashed", linewidth = 0.7) +
  annotate("text", x = 3.6, y = mu_true + 0.02,
           label = "True μ", size = 3.5, hjust = 0) +
  scale_x_continuous(breaks = sigma_delta_grid) +
  labs(title = "Panel 2: Sensitivity to prior on δ",
       subtitle = "Moderate bias scenario (δ = 1.5). Posterior mean and 95% CI.",
       x = "Prior SD on δ (σ_δ)",
       y = "Posterior mean and 95% CI for μ") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        plot.subtitle     = element_text(size = 9, colour = "grey40"))

# =============================================================================
# PANEL 3: Width decomposition
# =============================================================================
width_df <- main_df %>%
  filter(grepl("Bayesian", method)) %>%
  left_join(ht_bands, by = "scenario") %>%
  mutate(
    sampling_width  = ht_upper - ht_lower,
    selection_width = pmax(0, (upper - lower) - sampling_width),
    scenario        = factor(scenario, levels = c("mild", "moderate", "severe"))
  ) %>%
  select(scenario, method, sampling_width, selection_width) %>%
  pivot_longer(cols = c(sampling_width, selection_width),
               names_to = "component", values_to = "width") %>%
  mutate(component = recode(component,
    "sampling_width"  = "Sampling uncertainty",
    "selection_width" = "Selection uncertainty (AI)"
  ))

p3 <- width_df %>%
  ggplot(aes(x = scenario, y = width, fill = component)) +
  geom_col(width = 0.5, alpha = 0.88) +
  facet_wrap(~method, ncol = 1) +
  scale_fill_manual(values = c(
    "Sampling uncertainty"       = "#3498db",
    "Selection uncertainty (AI)" = "#e67e22"
  )) +
  labs(title = "Panel 3: Posterior width decomposition",
       subtitle = "Blue = sampling; orange = selection uncertainty.",
       x = "Bias scenario", y = "95% CI width", fill = NULL) +
  theme_minimal(base_size = 12) +
  theme(legend.position  = "bottom",
        panel.grid.minor  = element_blank(),
        plot.subtitle     = element_text(size = 9, colour = "grey40"),
        strip.text        = element_text(size = 9))

# =============================================================================
# COMBINED FIGURE
# =============================================================================
combined <- (p1 | (p2 / p3)) +
  plot_annotation(
    title   = "Uncertainty Dashboard",
    caption = paste(
      "Left: layered intervals separate sampling (dark) from selection (light) uncertainty.",
      "Top right: sensitivity of posterior to prior on δ (moderate scenario).",
      "Bottom right: decomposition of interval width by uncertainty source."
    ),
    theme = theme(
      plot.title   = element_text(size = 14, face = "bold"),
      plot.caption = element_text(size = 8, colour = "grey40")
    )
  )

dir.create("data/simulation_results", recursive = TRUE, showWarnings = FALSE)
ggsave("data/simulation_results/figure_dashboard.pdf",
       combined, width = 14, height = 10, dpi = 300)
ggsave("data/simulation_results/figure_dashboard.png",
       combined, width = 14, height = 10, dpi = 300)

cat("Figures saved to data/simulation_results/\n")
cat("\n=== Visualisation complete ===\n")
