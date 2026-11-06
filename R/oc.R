#' Operating characteristics of the fixed and the Bayesian adaptive design
#'
#' Simulates `n_sim` trials per scenario and applies both designs to the
#' same simulated trials.
#'
#' @param scenarios Named list of scenarios; each a list with `rates_c`,
#'   `breaks_c`, `rates_t`, `breaks_t` (see [simulate_trial()]).
#' @param n_sim Number of simulated trials per scenario.
#' @param n,accrual_time,dropout_rate Trial parameters, see [simulate_trial()].
#' @inheritParams run_bayes_adaptive
#' @return Data frame with one row per scenario and design.
#' @export
operating_characteristics <- function(scenarios, n_sim, c_eff, a0, b0,
                                      looks = c(124, 185, 247), f_fut = 0.2,
                                      n = 400, accrual_time = 24, dropout_rate = 0.005) {
  D <- looks[length(looks)]
  rows <- lapply(names(scenarios), function(nm) {
    sc <- scenarios[[nm]]
    sims <- replicate(n_sim, {
      tr <- simulate_trial(n = n, accrual_time = accrual_time,
                           rates_c = sc$rates_c, breaks_c = sc$breaks_c,
                           rates_t = sc$rates_t, breaks_t = sc$breaks_t,
                           dropout_rate = dropout_rate)
      c(fixed = run_fixed(tr, n_events = D),
        bayes = run_bayes_adaptive(tr,looks = looks, f_fut = f_fut, c_eff = c_eff, futility = TRUE, a0 = a0, b0 = b0))
    })
    rbind(summarise_sims(sims, "fixed", nm, D),
          summarise_sims(sims, "bayes", nm, D))
  })
  do.call(rbind, rows)
}

#' Summarise simulated trials of one design (internal)
#' @noRd
summarise_sims <- function(sims, design, scenario, D) {
  get <- function(x) sims[paste0(design, ".", x), ]     # one row of the result matrix
  p <- mean(get("reject"))
  data.frame(
    scenario      = scenario,
    design        = design,
    reject        = p,
    reject_se     = sqrt(p*(1-p)/ncol(sims)),       # Monte Carlo standard error
    hr            = exp(mean(log(get("hr")))),        # geometric mean (averaged on the log scale)
    n_events      = mean(get("n_events")),
    duration      = mean(get("t_cut")),
    n_recruited   = mean(get("n_recruited")),
    early_stop    = mean(get("n_events") < D),   # proportion of trials stopped at an interim
    stop_futility = if (design == "bayes") mean(get("stop_futility")) else NA
  )
}
