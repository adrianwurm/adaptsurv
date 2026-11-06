test_that("prob is 0.5 for identical arms", {
  set.seed(42)
  a0 <- 1
  b0 <- 12/log(2)
  n_arm <- 10

  time   <- round(stats::runif(n_arm, 1, 24), 1)
  status <- rep(c(1, 1, 0), length.out = n_arm)

  data <- data.frame(
    time   = rep(time,   times = 2),
    status = rep(status, times = 2),
    arm    = rep(0:1, each = n_arm)
  )
  x <- analyse_bayes(data, a0, b0)
  expect_equal(x[["prob"]], 0.5)
})

test_that("prob matches Monte Carlo and F distribution", {
  set.seed(42)
  a0 <- 1
  b0 <- 12/log(2)

  data <- cut_at_events(simulate_trial(400, 30, 0.4, 0.5),100)$data



  post <- posterior_gamma(data, a0, b0)
  a_c <- post$shape[1]; b_c <- post$rate[1]
  a_t <- post$shape[2]; b_t <- post$rate[2]
  x <- analyse_bayes(data, a0, b0)

  mc <- mean(stats::rgamma(1e5, a_t, b_t) < stats::rgamma(1e5, a_c, b_c))
  expect_lt(abs(mc - x[["prob"]]), 0.01)

  expect_equal(stats::pf((b_t * a_c) / (a_t * b_c), 2 * a_t, 2 * a_c), x[["prob"]])
})

test_that("strong effect gives prob close to 1", {
  set.seed(42)

  rc <- c(0.10, 0.06)
  trial <- simulate_trial(n = 400, accrual_time = 12,
                          rates_c = rc,       breaks_c = 6,
                          rates_t = 0.5 * rc, breaks_t = 6)
  x <- cut_at_events(trial,240)

  y <- analyse_bayes(x$data, 1, 12/log(2))[["prob"]]
  expect_gt(y, 0.99)
})



