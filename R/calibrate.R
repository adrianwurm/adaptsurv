#' Posterior probabilities P(HR < 1 | data) at each analysis
#'
#' @inheritParams run_bayes_adaptive
#' @return Numeric vector, one probability per element of `looks`.
#' @export
bayes_probs <- function(trial, looks, a0, b0) {
  vapply(looks, function(d) analyse_bayes(cut_at_events(trial, d)$data, a0, b0)[["prob"]], numeric(1))
}

#' Calibrate the efficacy threshold to a one-sided type I error
#'
#' Simulates `n_sim` trials under the null hypothesis (ignoring futility,
#' i.e. non-binding) and returns the (1 - alpha)-quantile of the maximum
#' posterior probability over all analyses.
#'
#' @param n_sim Number of simulated trials.
#' @param alpha Target one-sided type I error.
#' @param ... Arguments passed to [simulate_trial()] (scenario under H0).
#' @inheritParams run_bayes_adaptive
#' @return Calibrated threshold `c_eff`.
#' @export
calibrate_bayes <- function(..., n_sim, alpha = 0.025, looks = c(124, 185, 247), a0, b0) {
  max_prob <- vapply(seq_len(n_sim), function(i) {
    tr <- simulate_trial(...)
    max(bayes_probs(tr, looks, a0, b0))
  }, numeric(1))
  unname(stats::quantile(max_prob, 1-alpha))
}

