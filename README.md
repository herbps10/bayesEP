# bayesEP

<!-- badges: start -->
[![R-CMD-check](https://github.com/herbps10/bayesEP/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/herbps10/bayesEP/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

## Overview
`bayesEP` implements a flexible framework for distributed Bayesian inference using Expectation Propagation (EP) based on the algorithm described by [Vehtari et al. 2020](https://jmlr.org/papers/v21/18-817.html).

Fitting large Bayesian models can be computationally expensive and time consuming. EP is an alternative method that approximates the full joint posterior by divide and conquer. This package focuses on Bayesian hierarchical models with a natural distinction between group-specific and global shared parameters. EP partitions the data by group to multiple sites, fits each site independently against a _cavity distribution_ that carries information from all other sites, and interatively refines a global approximation to the posterior over the shared (hierarchical) parameters. 

This package has a modular design: you supply a model fitting function that can use any MCM backend (e.g., [Stan](https://mc-stan.org) via [cmdstanr](https://mc-stan.github.io/cmdstanr/)), and `bayesEP` handles the EP algorithm.

## Installation

You can install the development version of `bayesEP` from [Github](https://github.com/herbps10/bayesEP) with:

```r
# install.packages("pak")
pak::pak("herbps10/bayesEP")
```

For MCMC sampling, we recommend using [cmdstanr](https://mc-stan.github.io/cmdstanr), an R interface to [Stan](https://mc-stan.org/). See the [cmdstanr documentation](https://mc-stan.org/cmdstanr/articles/cmdstanr.html) for installation instructions.

## Getting started
See the [Getting started with BayesEP](http://herbsusmann.com/bayesEP/getting-started.html) vignette for a complete example.

## Quick example
```{r}
library(bayesEP

# Model fitting function
# Must accept: data, use_cavity, cavity_mu, cavity_Sigma
# Must return: NULL on failure, or a list with a `phi` element
#              containing a matrix of shared parameter draws.
fit_model <- function(data, use_cavity = FALSE,
                      cavity_mu = NULL, cavity_Sigma = NULL) {

  # Run MCMC, for example using cmdstanr
  samples <- model$samples(
    # ...
  )

  # Return matrix of posterior draws for the
  # vector of shared parameters
  list(
    phi = posterior::as_draws_matrix(samples$draws("phi"))
  )
}

# Run EP
result <- fit_ep(
  data = my_data,
  group_column = "group",
  K = 10 # Number of sites
  d = 3  # Number of shared parameters
  fit_model = fit_model,
  max_iter = 20,
  conv_tol = 10
)

# EP multivariate normal approximation to the 
# posterior of the shared parameters
result$mu    # Approximated posterior mean
result$Sigma # Approximated posterior covariance
```

## References
- Vehtari, A., Gelman, A., Sivula, T., Jylänki, P., Tran, D., Sahai, S., Blomstedt, P., Cunningham, J. P., Schiminovich, D., and Robert, C. P. (2020).
  Expectation Propagation as a Way of Life: A Framework for Bayesian Inference
  on Partitioned Data. Journal of Machine Learning Research, 21(17), 1–53. [http://jmlr.org/papers/v21/18-817.html](http://jmlr.org/papers/v21/18-817.html)
