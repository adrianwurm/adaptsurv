# Uncalibrated adaptive design: c = 0.975 at all three looks, no futility
# (2,000 trials each, seed 2026). Result: type I error 0.048, i.e. the
# repeated looks inflate the type I error; hence the calibration in
# dev/calibrate.R.

devtools::load_all()
set.seed(2026)
lambda_c <- log(2) / 12          # control median 12 months

a0 <- 1
b0 <- 12 / log(2)
one_trial <- function(hr) {
  tr <- simulate_trial(n = 400, accrual_time = 24,
                       rates_c = lambda_c, rates_t = hr * lambda_c,
                       dropout_rate = 0.005)
  run_bayes_adaptive(tr, c_eff = 0.975, futility = FALSE, a0 = a0, b0 = b0)
}


sims_h0 <- replicate(2000, one_trial(hr = 1))     # matrix: rows = summaries, columns = trials
sims_h1 <- replicate(2000, one_trial(hr = 0.7))
rowMeans(sims_h0)   # reject = simulated type I error
rowMeans(sims_h1)   # reject = simulated power

