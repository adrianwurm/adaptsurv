#' Simulate piecewise exponential survival times
#'
#' @param n Number of observations.
#' @param rates Positive hazard rates, one per interval.
#' @param breaks Interior change points (increasing, > 0);
#'   `length(breaks)` must equal `length(rates) - 1`.
#' @return Numeric vector of length `n`.
#' @export
#' @examples
#' rpwexp(5, rates = c(0.1, 0.07), breaks = 6)  # effect starts after 6 months
rpwexp <- function(n, rates, breaks = numeric(0)) {
  stopifnot(length(rates) == length(breaks) + 1, all(rates > 0))

  starts <- c(0, breaks)                          # s_1, ..., s_k
  # cumulative hazard at the interval starts: H(s_1) = 0, H(s_2) = r_1 (s_2 - s_1), ...
  H_starts <- c(0, cumsum(rates[-length(rates)] * (starts[-1] - starts[-length(starts)])))

  E <- stats::rexp(n)                             # Exp(1)
  j <- findInterval(E, H_starts)                       # interval in which E falls
  starts[j]+(E-H_starts[j])/rates[j]          # inverse of the cumulative hazard
}
