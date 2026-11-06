# Calibration of the efficacy threshold c of the Bayesian adaptive design.
# Each step simulates 10,000 trials and takes several minutes.
#
# Results:
#   Step 1  c_eff    = 0.9873009     (calibrated under constant hazards)
#   Step 2  type I error with c_eff on independent trials: 0.0237
#           (identical with and without futility for this seed)
#   Step 3  c_break  = 0.9918179602  (calibrated under h0_break, i.e. both
#                                     arms share a hazard with a change point)
#           c_robust = max(c_eff, c_break) = 0.9918179602

devtools::load_all()
lambda_c <- log(2) / 12          # control median 12 months
a0 <- 1
b0 <- 12 / log(2)

# Step 1: calibration under constant hazards
set.seed(2026)
c_eff <- calibrate_bayes(n_sim = 10000, a0 = a0, b0 = b0,
                         n = 400, accrual_time = 24,
                         rates_c = lambda_c, rates_t = lambda_c,
                         dropout_rate = 0.005)

# Step 2: validation of c_eff on independent trials under H0 (no futility)
one_trial <- function(hr, futility) {
  tr <- simulate_trial(n = 400, accrual_time = 24,
                       rates_c = lambda_c, rates_t = hr * lambda_c,
                       dropout_rate = 0.005)
  run_bayes_adaptive(tr, c_eff = c_eff, futility = futility, a0 = a0, b0 = b0)
}
set.seed(2027)
val_h0 <- replicate(10000, one_trial(hr = 1, futility = FALSE))
mean(val_h0["reject", ])         # type I error, target 0.025

# Step 3: calibration under h0_break and robust threshold
set.seed(2032)
c_break <- calibrate_bayes(n_sim = 10000, a0 = a0, b0 = b0,
                           n = 400, accrual_time = 24,
                           rates_c = c(1, 0.6) * lambda_c, breaks_c = 6,
                           rates_t = c(1, 0.6) * lambda_c, breaks_t = 6,
                           dropout_rate = 0.005)
c_robust <- max(c_eff, c_break)
print(c(c_eff = c_eff, c_break = c_break, c_robust = c_robust), digits = 10)
