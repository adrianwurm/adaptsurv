test_that("sizes match", {
  set.seed(2)
  lambda_c <- log(2) / 12
  a0 <- 1; b0 <- 12 / log(2)
  c_eff <- 0.9873009
  sc <- list(
    h0      = list(rates_c = lambda_c, breaks_c = numeric(0),
                   rates_t = lambda_c, breaks_t = numeric(0)),
    delayed = list(rates_c = lambda_c, breaks_c = numeric(0),
                   rates_t = c(1, 0.6) * lambda_c, breaks_t = 6)
  )


  x <- operating_characteristics(scenarios=c(sc), n_sim=20,
        c_eff = c_eff, a0 = a0, b0 = b0, looks = c(124, 185, 247),
        f_fut = 0.2, n = 400, accrual_time = 24, dropout_rate = 0.005)


  expect_equal(nrow(x), 4)
  expect_true(all(x$reject<=1) & all(x$reject>=0))
  fx <- x[x$design == "fixed", ]
  expect_true(all(fx$n_events == 247))
  expect_true(all(fx$early_stop == 0))
  expect_true(all(is.na(fx$stop_futility)))
})


