# Operating characteristics of both designs (precomputed)

Result of
[`operating_characteristics()`](https://adrianwurm.github.io/adaptsurv/reference/operating_characteristics.md)
for six scenarios and two calibrations of the Bayesian efficacy
threshold, 10,000 simulated trials each (see `dev/run_oc.R` for the
exact calls).

## Usage

``` r
oc_results
```

## Format

A data frame with 24 rows, one per combination of calibration, scenario
and design. The first three columns identify the row; all other columns
are summaries over the simulated trials:

- calibration:

  `"constant"` (c = 0.9873, calibrated under constant hazards) or
  `"robust"` (c = 0.9918, also calibrated under a hazard with a change
  point).

- scenario:

  Scenario name.

- design:

  `"fixed"` or `"bayes"`.

- reject:

  Rejection rate (type I error under H0, power otherwise).

- reject_se:

  Monte Carlo standard error of `reject`.

- hr:

  Geometric mean of the Cox hazard ratio estimates (treatment vs.
  control) at the last analysis performed.

- n_events:

  Mean number of events at the last analysis performed.

- duration:

  Mean calendar time (months) from the start of recruitment to the last
  analysis performed.

- n_recruited:

  Mean number of patients recruited by the last analysis performed.

- early_stop:

  Proportion of trials stopped at an interim analysis (efficacy or
  futility); always 0 for the fixed design.

- stop_futility:

  Proportion of trials stopped for futility; `NA` for the fixed design.
