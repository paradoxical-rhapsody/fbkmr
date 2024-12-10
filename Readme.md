# FBKMR: Fast Bayesian Kernel Machine Regression

FBKMR is an R package that implements a scalable version of Bayesian Kernel Machine Regression (BKMR) for analyzing the health effects of complex environmental mixtures.

## Installation

```R
# Install from GitHub
devtools::install_github("username/fbkmr")
```

Note: This package requires Gurobi optimization software. Please visit https://www.gurobi.com/ for installation instructions.

## Overview

BKMR is a statistical method for studying the joint effects of multiple exposures on a health outcome, while allowing for nonlinear relationships and interactions. The main function `kmbayes()` fits the BKMR model using MCMC methods.

## Basic Usage

```R
library(fbkmr)

# Fit BKMR model
fit <- kmbayes(
  y = outcome,              # Outcome vector 
  Z = exposures,           # Exposure matrix
  X = covariates,          # Optional covariate matrix
  iter = 1000,             # Number of MCMC iterations
  family = "gaussian",     # "gaussian" or "binomial"
  varsel = FALSE          # Whether to do variable selection
)
```

### Key Parameters

- `y`: Vector of outcome data (n x 1)
- `Z`: Matrix of exposure variables (n x M) to be included in the h() function
- `X`: Optional matrix of covariates (n x K)
- `iter`: Number of MCMC iterations
- `family`: "gaussian" for continuous outcomes or "binomial" for binary outcomes
- `varsel`: Whether to perform variable selection on the exposures
- `groups`: Optional vector for grouped variable selection
- `knots`: Optional matrix for implementing Gaussian predictive process approximation
- `est.h`: Whether to sample the subject-specific effects h within MCMC

### Additional Options

- `control.params`: List of parameters for priors and MCMC tuning
- `starting.values`: List of starting values for parameters
- `verbose`: Whether to print progress during model fitting

## Example

```R
# Generate example data
n <- 100  # Sample size
M <- 5    # Number of exposures

set.seed(123)
Z <- matrix(rnorm(n*M), n, M)
X <- matrix(rnorm(n*2), n, 2) 
y <- rnorm(n)

# Fit model
fit <- kmbayes(
  y = y,
  Z = Z, 
  X = X,
  iter = 1000,
  verbose = TRUE
)

# Print summary
summary(fit)
```

## Scalable Version

For large datasets, use the `skmbayes()` function which implements a scalable version using data sketching:

```R
fit <- skmbayes(
  Z = Z,
  X = X, 
  y = y,
  n_subset = 5,    # Number of data subsets
  n_samp = 200,    # Posterior samples per subset
  iter = 1000
)
```

## References

- Bobb, JF, et al. (2015). Bayesian Kernel Machine Regression for Estimating the Health Effects of Multi-Pollutant Mixtures. Biostatistics 16(3): 493-508.

- For guided examples, visit: https://jenfb.github.io/bkmr/overview.html

## License

[License details]

## Contributing

[Contributing guidelines]

## Support

[Support information]

Let me know if you would like me to expand on any section or add additional details!
