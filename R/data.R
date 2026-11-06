#' Operating characteristics of both designs (precomputed)
#'
#' Result of [operating_characteristics()] for six scenarios and two
#' calibrations of the Bayesian efficacy threshold, 10,000 simulated trials
#' each (see `dev/run_oc.R` for the exact calls).
#'
#' @format A data frame with 24 rows, one per combination of calibration,
#'   scenario and design. The first three columns identify the row; all
#'   other columns are summaries over the simulated trials:
#' \describe{
#'   \item{calibration}{`"constant"` (c = 0.9873, calibrated under constant
#'     hazards) or `"robust"` (c = 0.9918, also calibrated under a hazard with
#'     a change point).}
#'   \item{scenario}{Scenario name.}
#'   \item{design}{`"fixed"` or `"bayes"`.}
#'   \item{reject}{Rejection rate (type I error under H0, power otherwise).}
#'   \item{reject_se}{Monte Carlo standard error of `reject`.}
#'   \item{hr}{Geometric mean of the Cox hazard ratio estimates
#'     (treatment vs. control) at the last analysis performed.}
#'   \item{n_events}{Mean number of events at the last analysis performed.}
#'   \item{duration}{Mean calendar time (months) from the start of
#'     recruitment to the last analysis performed.}
#'   \item{n_recruited}{Mean number of patients recruited by the last
#'     analysis performed.}
#'   \item{early_stop}{Proportion of trials stopped at an interim analysis
#'     (efficacy or futility); always 0 for the fixed design.}
#'   \item{stop_futility}{Proportion of trials stopped for futility;
#'     `NA` for the fixed design.}
#' }
"oc_results"
