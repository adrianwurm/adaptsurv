#' Cut trial data at the calendar time of the k-th observed event
#'
#' @param trial Output of [simulate_trial()].
#' @param n_events Number of events triggering the analysis.
#' @return A list with `data` (columns `id`, `arm`, `time`, `status` for
#'   patients recruited before the cut) and `t_cut` (calendar time of the cut).
#' @export
cut_at_events <- function(trial, n_events) {
  # an event is only observable if it occurs before drop-out
  observed <- (trial$event_time<trial$dropout_time)
  calendar_event <- ifelse(observed, trial$entry + trial$event_time, Inf)

  if (sum(observed) < n_events) {
    stop("Only ", sum(observed), " events can ever be observed.")
  }

  t_cut <- sort(calendar_event)[n_events]

  in_trial <- trial$entry < t_cut                 # recruited before the cut?
  d <- trial[in_trial, ]
  cal <- calendar_event[in_trial]

  d$time   <- pmin(d$entry+d$event_time,d$entry+d$dropout_time, t_cut)-d$entry  # event, drop-out or end of follow-up at t_cut
  d$status <- as.integer(cal<=t_cut)

  list(data = d[, c("id", "arm", "time", "status")], t_cut = t_cut)
}
