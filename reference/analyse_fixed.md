# Analyse a data cut with the log-rank test and a Cox model

Analyse a data cut with the log-rank test and a Cox model

## Usage

``` r
analyse_fixed(data)
```

## Arguments

- data:

  Data frame with columns `time`, `status`, `arm` (e.g.
  `cut_at_events()$data`).

## Value

Named numeric vector: `z` (signed log-rank statistic, negative favours
treatment), `p_one_sided`, `hr` (Cox hazard ratio, treatment vs
control).
