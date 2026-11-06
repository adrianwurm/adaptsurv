# Cut trial data at the calendar time of the k-th observed event

Cut trial data at the calendar time of the k-th observed event

## Usage

``` r
cut_at_events(trial, n_events)
```

## Arguments

- trial:

  Output of
  [`simulate_trial()`](https://adrianwurm.github.io/adaptsurv/reference/simulate_trial.md).

- n_events:

  Number of events triggering the analysis.

## Value

A list with `data` (columns `id`, `arm`, `time`, `status` for patients
recruited before the cut) and `t_cut` (calendar time of the cut).
