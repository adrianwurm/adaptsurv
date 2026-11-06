test_that("sum(status) is correct", {
  set.seed(42)

  tr <- simulate_trial(40, accrual_time=20, rates_c=3, rates_t=5,
                       breaks_c = numeric(0), breaks_t = numeric(0),
                       dropout_rate = 2)
  x <- cut_at_events (tr,8)


  expect_equal(sum(x$data$status), 8)


  expect_true(all(x$data$time>=0))

  expect_true(all(x$t_cut-tr$entry[x$data$id]>=x$data$time))
  expect_gt(nrow(x$data), 0)
  expect_error(cut_at_events(tr, 1000))

})

