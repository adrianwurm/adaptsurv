#' Simulate latent data for a two-arm survival trial
#'
#' Patients are recruited uniformly over `[0, accrual_time]` and randomised 1:1.
#' Event times follow piecewise exponential distributions (see [rpwexp()]),
#' drop-out times an exponential distribution. All times except `entry`
#' are measured from the patient's entry.
#'
#' @param n Total number of patients (even).
#' @param accrual_time Length of the recruitment period.
#' @param rates_c,breaks_c Hazard rates and change points, control arm.
#' @param rates_t,breaks_t Hazard rates and change points, treatment arm.
#' @param dropout_rate Exponential drop-out rate (0 = no drop-out).
#' @return A data frame with columns `id`, `arm` (0 = control, 1 = treatment),
#'   `entry`, `event_time`, `dropout_time`.
#' @export
simulate_trial <- function(n, accrual_time, rates_c, rates_t,
                           breaks_c = numeric(0), breaks_t = numeric(0),
                           dropout_rate = 0) {
  stopifnot(n %% 2 == 0, accrual_time > 0, dropout_rate >= 0)

  arm   <- sample(rep(0:1, each = n / 2))         # exactly 1:1, random order
  entry <- sort(stats::runif(n, 0, accrual_time))

  event_time <- numeric(n)
  event_time[arm == 0] <- rpwexp(n/2, rates_c, breaks_c)
  event_time[arm == 1] <- rpwexp(n/2, rates_t, breaks_t)

  dropout_time <- if (dropout_rate > 0) stats::rexp(n, rate = dropout_rate) else rep(Inf, n)

  data.frame(id = seq_len(n), arm, entry, event_time, dropout_time)
}
