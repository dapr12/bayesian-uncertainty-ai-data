---

# Bayesian Uncertainty for AI-Augmented Data

This repository accompanies the paper:

**Pérez Ruiz, D.A. (2026)**
*Bayesian Uncertainty Quantification for AI-Augmented Social Science Data: A Framework for Combining Probability and Non-Probability Sources.*

---

## Abstract

AI-mediated data sources—such as synthetic survey respondents, web-scraped behavioural traces, and algorithmically mediated panels—are increasingly combined with traditional probability samples. These sources introduce **selection uncertainty**: bias arising from unknown and unobservable data-generating mechanisms.

This repository implements a **Bayesian partial identification framework** that treats the selection mechanism as an explicit object of inference. A prior is placed on a selection bias parameter, and posterior inference marginalises over this uncertainty, producing credible intervals that reflect both sampling variability and structural uncertainty in the data source.

---

## Methodological Summary

We consider two samples from a population:

* A probability sample with known inclusion mechanism
* A non-probability (AI-mediated) sample with unknown selection mechanism

The framework introduces:

* A latent variable ( Z_i ) capturing unobserved selection drivers
* A selection model:
  [
  P(R_i = 1 \mid X_i, Z_i) = \mathrm{logistic}(X_i^\top \alpha + \delta Z_i)
  ]
* A prior on the selection bias parameter:
  [
  \delta \sim \mathcal{N}(0, \sigma_\delta^2)
  ]

The parameter ( \delta ) is not identified from the observed data. Instead, inference proceeds via **partial identification**, with posterior uncertainty reflecting both data and assumptions.

---

## Computational Approach

Direct inference is computationally challenging due to the latent structure. To ensure tractability, the implementation uses a **first-order mean-shift approximation**:

[
\mathbb{E}[Z_i \mid R_i = 1] \approx \delta \cdot \phi(0),
]

where ( \phi(0) ) is the standard normal density at zero.

This approximation:

* Avoids sampling latent variables (preventing funnel pathologies in HMC)
* Preserves the dependence of posterior inference on the selection parameter
* Is locally accurate around ( \delta = 0 ), where the model is anchored by the probability sample

---

## Repository Structure

```
.
├── R/                # Simulation and analysis scripts
├── stan/             # Stan model implementation
├── data/             # Simulated datasets (if included)
└── README.md
```

---

## Reproducibility

All results in the paper are reproducible from this repository.

### Requirements

* R (≥ 4.0)
* Packages:

  * `rstan`
  * `tidyverse`
  * `posterior`


---

## Simulation Design

The simulation study evaluates performance under varying levels of selection bias:

* Mild bias: ( \delta = 0.5 )
* Moderate bias: ( \delta = 1.5 )
* Severe bias: ( \delta = 3.0 )

Estimators compared:

* Naive pooling
* Horvitz–Thompson estimator (probability sample only)
* Bayesian framework with varying prior assumptions on ( \delta )

Performance is assessed via:

* Bias
* Coverage of 95% intervals
* Interval width

---

## Uncertainty Decomposition

The framework distinguishes between two sources of uncertainty:

* **Sampling uncertainty**: due to finite probability sample size
* **Selection uncertainty**: due to unknown AI selection mechanism

This distinction is operationalised in an interactive dashboard.

---

## Interactive Dashboard

An interactive Shiny application is available at:

[https://dapr12.shinyapps.io/Bayesian_Uncertainty/](https://dapr12.shinyapps.io/Bayesian_Uncertainty/)

The dashboard provides:

* Layered interval visualisation
* Sensitivity analysis over ( \sigma_\delta )
* Decomposition of posterior uncertainty

---

## Interpretation

The key principle of the framework is:

> Increasing the size of an AI-mediated dataset does not reduce uncertainty if the selection mechanism is unknown.

Naive pooling implicitly assumes ( \delta = 0 ). The proposed approach instead makes this assumption explicit and propagates its uncertainty through inference.

---

## Limitations

* The selection mechanism is modelled using a low-dimensional parametric form
* The mean-shift approximation simplifies the full likelihood
* Prior specification for ( \delta ) requires substantive judgement

Extensions include:

* Higher-dimensional latent structures
* Non-parametric selection models
* Integration of AI audit information into prior elicitation

---

## Citation

If you use this code, please cite:

```
Pérez Ruiz, D.A. (2026).
Bayesian Uncertainty Quantification for AI-Augmented Social Science Data.
```

---

## License

Specify license (e.g. MIT)

---

## Contact

Diego Andrés Pérez Ruiz
Department of Social Statistics
University of Manchester
