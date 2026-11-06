# Run the Bayesian adaptive design on one simulated trial

Run the Bayesian adaptive design on one simulated trial

## Usage

``` r
run_bayes_adaptive(
  trial,
  looks = c(124, 185, 247),
  c_eff,
  f_fut = 0.2,
  futility = TRUE,
  a0,
  b0
)
```

## Arguments

- trial:

  Output of
  [`simulate_trial()`](https://adrianwurm.github.io/adaptsurv/reference/simulate_trial.md).

- looks:

  Increasing event numbers of the analyses (last = final analysis).

- c_eff:

  Efficacy threshold: stop and reject if P(HR \< 1 \| data) \> `c_eff`.

- f_fut:

  Futility threshold at interim analyses: stop if P(HR \< 1 \| data) \<
  `f_fut`.

- futility:

  Apply the futility rule? (`FALSE` for non-binding calibration.)

- a0, b0:

  Shape and rate of the Gamma prior (same prior in both arms).

## Value

Named numeric vector: `reject` (0/1), `stop_look` (index of the analysis
at which the trial stopped), `stop_futility` (0/1), `hr` (Cox HR at the
stopping analysis), `n_events`, `t_cut`, `n_recruited`.
