# Gamma posterior parameters of the exponential hazard rate per arm

Gamma posterior parameters of the exponential hazard rate per arm

## Usage

``` r
posterior_gamma(data, a0, b0)
```

## Arguments

- data:

  Data frame with columns `time`, `status`, `arm`.

- a0, b0:

  Shape and rate of the Gamma prior (same prior in both arms).

## Value

List with `shape` and `rate`, each of length 2 (control, treatment).
