# Exploratory comparison of both designs under non-proportional hazards
# (constant calibration only, seed 2030). Superseded by dev/run_oc.R and the
# dataset oc_results; kept for reference. Result stored in dev/tab_nph.rds.

devtools::load_all()
lambda_c <- log(2) / 12
a0 <- 1; b0 <- 12 / log(2)
c_eff <- 0.9873009
D <- 247

# each scenario: hazard rates and change points for both arms
scenarios <- list(
  h0       = list(rates_c = lambda_c, breaks_c = numeric(0),
                  rates_t = lambda_c, breaks_t = numeric(0)),
  h0_break = list(rates_c = c(1, 0.6) * lambda_c , breaks_c = 6,    # both arms identical, with change point
                  rates_t = c(1, 0.6) * lambda_c, breaks_t = 6),
  ph_07    = list(rates_c = lambda_c, breaks_c = numeric(0),
                  rates_t = 0.7 * lambda_c, breaks_t = numeric(0)),
  delayed  = list(rates_c = lambda_c , breaks_c = numeric(0),
                   rates_t = c(1, 0.6) * lambda_c, breaks_t = 6),
  diminish = list(rates_c = lambda_c , breaks_c = numeric(0),
                  rates_t = c(0.5, 1) * lambda_c, breaks_t = 6),
  crossing = list(rates_c = lambda_c , breaks_c = numeric(0),
                  rates_t = c(1.3, 0.5) * lambda_c, breaks_t = 6)
)

one_trial <- function(sc) {
  tr <- simulate_trial(n = 400, accrual_time = 24,
                       rates_c = sc$rates_c, breaks_c = sc$breaks_c,
                       rates_t = sc$rates_t, breaks_t = sc$breaks_t,
                       dropout_rate = 0.005)
  c(fixed = run_fixed(tr, n_events = D),
    bayes = run_bayes_adaptive(tr, c_eff = c_eff, futility = TRUE, a0 = a0, b0 = b0))
}

set.seed(2030)
res <- lapply(scenarios, function(sc) rowMeans(replicate(10000, one_trial(sc))))
tab <- do.call(rbind, res)
round(tab[, c("fixed.reject", "bayes.reject", "bayes.stop_futility",
              "bayes.n_events", "fixed.t_cut", "bayes.t_cut")], 3)

