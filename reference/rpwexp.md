# Simulate piecewise exponential survival times

Simulate piecewise exponential survival times

## Usage

``` r
rpwexp(n, rates, breaks = numeric(0))
```

## Arguments

- n:

  Number of observations.

- rates:

  Positive hazard rates, one per interval.

- breaks:

  Interior change points (increasing, \> 0); `length(breaks)` must equal
  `length(rates) - 1`.

## Value

Numeric vector of length `n`.

## Examples

``` r
rpwexp(5, rates = c(0.1, 0.07), breaks = 6)  # effect starts after 6 months
#> [1] 31.306339  6.980944  2.015218 20.913741 71.361566
```
