# Simulate latent data for a two-arm survival trial

Patients are recruited uniformly over `[0, accrual_time]` and randomised
1:1. Event times follow piecewise exponential distributions (see
[`rpwexp()`](https://adrianwurm.github.io/adaptsurv/reference/rpwexp.md)),
drop-out times an exponential distribution. All times except `entry` are
measured from the patient's entry.

## Usage

``` r
simulate_trial(
  n,
  accrual_time,
  rates_c,
  rates_t,
  breaks_c = numeric(0),
  breaks_t = numeric(0),
  dropout_rate = 0
)
```

## Arguments

- n:

  Total number of patients (even).

- accrual_time:

  Length of the recruitment period.

- rates_c, breaks_c:

  Hazard rates and change points, control arm.

- rates_t, breaks_t:

  Hazard rates and change points, treatment arm.

- dropout_rate:

  Exponential drop-out rate (0 = no drop-out).

## Value

A data frame with columns `id`, `arm` (0 = control, 1 = treatment),
`entry`, `event_time`, `dropout_time`.
