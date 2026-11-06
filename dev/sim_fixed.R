# Validation of the fixed design (1,000 trials each, seed 2026).
# Results: type I error 0.032 (within Monte Carlo error, SE about 0.005),
#          power under HR 0.7: 0.796 (Schoenfeld: 0.80).

devtools::load_all()
set.seed(2026)
lambda_c <- log(2) / 12          # control median 12 months
D        <- ceiling(4*(qnorm(0.975) + qnorm(0.8))^2 / (log (7/10))^2)  # Schoenfeld: 247 events

one_trial <- function(hr) {
  tr <- simulate_trial(n = 400, accrual_time = 24,
                       rates_c = lambda_c, rates_t = hr * lambda_c,
                       dropout_rate = 0.005)
  run_fixed(tr, n_events = D)
}


sims_h0 <- replicate(1000, one_trial(hr = 1))     # matrix: rows = summaries, columns = trials
sims_h1 <- replicate(1000, one_trial(hr = 0.7))
rowMeans(sims_h0)   # reject = simulated type I error
rowMeans(sims_h1)   # reject = simulated power

