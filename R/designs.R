#' Run the Bayesian adaptive design on one simulated trial
#'
#' @param trial Output of [simulate_trial()].
#' @param looks Increasing event numbers of the analyses (last = final analysis).
#' @param c_eff Efficacy threshold: stop and reject if P(HR < 1 | data) > `c_eff`.
#' @param f_fut Futility threshold at interim analyses: stop if P(HR < 1 | data) < `f_fut`.
#' @param futility Apply the futility rule? (`FALSE` for non-binding calibration.)
#' @inheritParams posterior_gamma
#' @return Named numeric vector: `reject` (0/1), `stop_look` (index of the
#'   analysis at which the trial stopped), `stop_futility` (0/1), `hr` (Cox HR
#'   at the stopping analysis), `n_events`, `t_cut`, `n_recruited`.
#' @export
run_bayes_adaptive <- function(trial, looks = c(124, 185, 247), c_eff,
                               f_fut = 0.2, futility = TRUE, a0, b0) {
  K <- length(looks)
  for (k in seq_len(K)) {
    cut  <- cut_at_events(trial, looks[k])
    prob <- analyse_bayes(cut$data,a0,b0)[["prob"]]

    eff  <- prob > c_eff                                    # efficacy
    fut  <- futility && k < K && prob < f_fut                # futility (interim analyses only)

    if (eff || fut || k == K) {
      return(c(reject        = as.numeric(eff) ,
               stop_look     = k,
               stop_futility = as.numeric(fut && !eff),
               hr            = analyse_fixed(cut$data)[["hr"]],
               n_events      = looks[k],
               t_cut         = cut$t_cut,
               n_recruited   = nrow(cut$data)))
    }
  }
}
