# Run the fixed design on one simulated trial

Run the fixed design on one simulated trial

## Usage

``` r
run_fixed(trial, n_events, alpha = 0.025)
```

## Arguments

- trial:

  Output of
  [`simulate_trial()`](https://adrianwurm.github.io/adaptsurv/reference/simulate_trial.md).

- n_events:

  Planned number of events for the single analysis.

- alpha:

  One-sided significance level.

## Value

Named numeric vector: `reject` (0/1), `hr`, `t_cut`, `n_recruited`
,`n_events`
