# Bayesian analysis of a data cut

Bayesian analysis of a data cut

## Usage

``` r
analyse_bayes(data, a0, b0)
```

## Arguments

- data:

  Data frame with columns `time`, `status`, `arm`.

- a0, b0:

  Shape and rate of the Gamma prior (same prior in both arms).

## Value

Named numeric vector: `prob` = P(HR \< 1 \| data), `hr_median` =
posterior median of the hazard ratio.
