# eiv-fin

R functions and datasets for errors-in-variables models

## Overview

`eivfin` implements classical errors-in-variables (measurement error)
regression estimators, developed as a companion to the dissertation
*Errors-in-Variables Models and their Applications in Finance*:

- `tls()` — total least squares (orthogonal) regression via the SVD of the
  centered data matrix, with pairs-bootstrap covariance estimation.
- `deming()` — Deming regression with known error variance ratio
  (λ = σ²_y-error / σ²_x-error) and delete-one jackknife standard errors.

Both return `eivfit` objects with `coef()`, `vcov()`, `confint()` and
`summary()` methods.

Planned estimators: regression calibration, SIMEX, instrumental-variable
and grouped-data corrections, and maximum likelihood for the linear EIV
model with replication.

## Installation

```r
# install.packages("remotes")
remotes::install_github("erkind/eiv-fin")
```

## Example

OLS understates the slope when the regressor is measured with error
(attenuation bias); the EIV estimators correct it:

```r
library(eivfin)

set.seed(2024)
n <- 240
xi <- rnorm(n, sd = 0.04)                  # true factor
r  <- 0.005 + 1.2 * xi + rnorm(n, sd = 0.01)   # excess return
m  <- xi + rnorm(n, sd = 0.01)             # noisy factor proxy

coef(tls(m, r, seed = 2024))               # recovers ~1.2
coef(deming(m, r, lambda = 1))             # orthogonal regression, ~1.2
cov(m, r) / var(m)                         # OLS, attenuated to ~1.13
```

See `vignette("eiv-fin")` for a full treatment of attenuation bias in
financial beta estimation.
