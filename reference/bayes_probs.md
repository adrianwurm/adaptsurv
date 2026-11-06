# Posterior probabilities P(HR \< 1 \| data) at each analysis

Posterior probabilities P(HR \< 1 \| data) at each analysis

## Usage

``` r
bayes_probs(trial, looks, a0, b0)
```

## Arguments

- trial:

  Output of
  [`simulate_trial()`](https://adrianwurm.github.io/adaptsurv/reference/simulate_trial.md).

- looks:

  Increasing event numbers of the analyses (last = final analysis).

- a0, b0:

  Shape and rate of the Gamma prior (same prior in both arms).

## Value

Numeric vector, one probability per element of `looks`.
