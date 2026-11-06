
test_that("Extremal efficacy", {
  set.seed(42)
  a0 <- 1
  b0 <- 12/log(2)

  rc <- c(0.10, 0.06)
  trial <- simulate_trial(n = 400, accrual_time = 12,
                          rates_c = rc,       breaks_c = 6,
                          rates_t = 1*rc, breaks_t = 6)

  x <- run_bayes_adaptive(trial, looks = c(124, 185, 247), c_eff=1,
                                 f_fut = 0.2, futility = FALSE, a0 = a0, b0 = b0)

  expect_equal(x[["stop_look"]], 3)
  expect_equal(x[["reject"]], 0)

  x <- run_bayes_adaptive(trial, looks = c(124, 185, 247), c_eff=0,
                          f_fut = 0.2, futility = FALSE, a0 = a0, b0 = b0)
  expect_equal(x[["stop_look"]], 1)
  expect_equal(x[["reject"]], 1)
})

test_that("Futility works", {
  set.seed(42)
  a0 <- 1
  b0 <- 12/log(2)

  rc <- c(0.10, 0.06)
  trial <- simulate_trial(n = 400, accrual_time = 12,
                          rates_c = rc,       breaks_c = 6,
                          rates_t = 1.5*rc, breaks_t = 6)
  looks = c(124, 185, 247)
  x <- run_bayes_adaptive(trial, looks = looks, c_eff=0.99,
                          f_fut = 0.2, futility = TRUE, a0 = a0, b0 = b0)

  expect_equal(x[["stop_futility"]], 1)
  expect_equal(x[["reject"]], 0)
  expect_equal(x[["n_events"]], looks[x[["stop_look"]]])

})

