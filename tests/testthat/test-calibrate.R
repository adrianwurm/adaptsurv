test_that("length correct", {
  set.seed(2026)
  lambda_c <- log(2) / 12
  looks    <- c(124, 185, 247)
  a0 <- 1
  b0 <- 12 / log(2)
  trial <- simulate_trial(n = 400, accrual_time = 24,
                          rates_c = lambda_c, rates_t = lambda_c,
                          dropout_rate = 0.005)
  pr <- bayes_probs(trial, looks, a0, b0)
  expect_equal(length(pr),3)
  expect_equal(pr[3], analyse_bayes(cut_at_events(trial, 247)$data, a0, b0)[["prob"]])
})
