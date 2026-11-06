# Calibrate the efficacy threshold to a one-sided type I error

Simulates `n_sim` trials under the null hypothesis (ignoring futility,
i.e. non-binding) and returns the (1 - alpha)-quantile of the maximum
posterior probability over all analyses.

## Usage

``` r
calibrate_bayes(..., n_sim, alpha = 0.025, looks = c(124, 185, 247), a0, b0)
```

## Arguments

- ...:

  Arguments passed to
  [`simulate_trial()`](https://adrianwurm.github.io/adaptsurv/reference/simulate_trial.md)
  (scenario under H0).

- n_sim:

  Number of simulated trials.

- alpha:

  Target one-sided type I error.

- looks:

  Increasing event numbers of the analyses (last = final analysis).

- a0, b0:

  Shape and rate of the Gamma prior (same prior in both arms).

## Value

Calibrated threshold `c_eff`.
