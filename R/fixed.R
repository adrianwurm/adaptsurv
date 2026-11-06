#' Analyse a data cut with the log-rank test and a Cox model
#'
#' @param data Data frame with columns `time`, `status`, `arm`
#'   (e.g. `cut_at_events()$data`).
#' @return Named numeric vector: `z` (signed log-rank statistic, negative
#'   favours treatment), `p_one_sided`, `hr` (Cox hazard ratio, treatment vs control).
#' @export
analyse_fixed <- function(data) {
  lr  <- survival::survdiff(survival::Surv(time, status) ~ arm, data = data)
  z <- (lr$obs[2]-lr$exp[2])/sqrt(lr$var[2,2])
  cox <- survival::coxph(survival::Surv(time, status) ~ arm, data = data)
  c(z = z, p_one_sided = unname(stats::pnorm(z)), hr = unname(exp(stats::coef(cox))))
}

#' Run the fixed design on one simulated trial
#'
#' @param trial Output of [simulate_trial()].
#' @param n_events Planned number of events for the single analysis.
#' @param alpha One-sided significance level.
#' @return Named numeric vector: `reject` (0/1), `hr`, `t_cut`, `n_recruited` ,`n_events`
#' @export
run_fixed <- function(trial, n_events, alpha = 0.025) {
  cut <- cut_at_events(trial, n_events)
  res <- analyse_fixed(cut$data)
  c(reject      = res[["z"]] < stats::qnorm(alpha),                 # reject if z < qnorm(alpha)
    hr          = unname(res["hr"]),
    t_cut       = cut$t_cut,
    n_recruited = nrow(cut$data),
    n_events = n_events)
}
