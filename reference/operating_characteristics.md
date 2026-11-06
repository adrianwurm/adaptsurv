# Operating characteristics of the fixed and the Bayesian adaptive design

Simulates `n_sim` trials per scenario and applies both designs to the
same simulated trials.

## Usage

``` r
operating_characteristics(
  scenarios,
  n_sim,
  c_eff,
  a0,
  b0,
  looks = c(124, 185, 247),
  f_fut = 0.2,
  n = 400,
  accrual_time = 24,
  dropout_rate = 0.005
)
```

## Arguments

- scenarios:

  Named list of scenarios; each a list with `rates_c`, `breaks_c`,
  `rates_t`, `breaks_t` (see
  [`simulate_trial()`](https://adrianwurm.github.io/adaptsurv/reference/simulate_trial.md)).

- n_sim:

  Number of simulated trials per scenario.

- c_eff:

  Efficacy threshold: stop and reject if P(HR \< 1 \| data) \> `c_eff`.

- a0, b0:

  Shape and rate of the Gamma prior (same prior in both arms).

- looks:

  Increasing event numbers of the analyses (last = final analysis).

- f_fut:

  Futility threshold at interim analyses: stop if P(HR \< 1 \| data) \<
  `f_fut`.

- n, accrual_time, dropout_rate:

  Trial parameters, see
  [`simulate_trial()`](https://adrianwurm.github.io/adaptsurv/reference/simulate_trial.md).

## Value

Data frame with one row per scenario and design.
