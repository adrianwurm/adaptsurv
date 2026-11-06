#' Gamma posterior parameters of the exponential hazard rate per arm
#'
#' @param data Data frame with columns `time`, `status`, `arm`.
#' @param a0,b0 Shape and rate of the Gamma prior (same prior in both arms).
#' @return List with `shape` and `rate`, each of length 2 (control, treatment).
#' @export
posterior_gamma <- function(data, a0, b0) {
  d  <- tapply(X=data$status, INDEX=data$arm, sum)     # events per arm
  ex <- tapply(X=data$time, INDEX=data$arm, sum)      # total exposure time per arm
  list(shape = unname(a0 + d), rate = unname(b0 + ex))
}

#' Bayesian analysis of a data cut
#'
#' @inheritParams posterior_gamma
#' @return Named numeric vector: `prob` = P(HR < 1 | data),
#'   `hr_median` = posterior median of the hazard ratio.
#' @export
analyse_bayes <- function(data, a0, b0) {
  post <- posterior_gamma(data, a0, b0)
  a_c <- post$shape[1]; b_c <- post$rate[1]
  a_t <- post$shape[2]; b_t <- post$rate[2]

  prob <- stats::pbeta(b_t / (b_t + b_c), a_t, a_c)
  # HR = lambda_t / lambda_c follows a scaled F distribution:
  # (lambda_t / lambda_c) * (b_t * a_c) / (a_t * b_c) ~ F(2 a_t, 2 a_c)
  hr_median <- stats::qf(0.5, 2*a_t,2*a_c)*(a_t * b_c)/(b_t * a_c)  # posterior median
  c(prob = prob, hr_median = hr_median)
}
