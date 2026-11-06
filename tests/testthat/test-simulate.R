test_that("data frame sizes fit", {
  set.seed(42)
  x <- simulate_trial(18, accrual_time=20, rates_c=3, rates_t=5,
                             breaks_c = numeric(0), breaks_t = numeric(0),
                             dropout_rate = 2)

  expect_equal(nrow(x),18)
  expect_equal(sum(x$arm==0),9)
  expect_equal(sum(x$arm==1),9)

  expect_true(all(x$entry<20 & x$entry>0))
})

test_that("dropout never happens at rate 0", {
  set.seed(42)
  x <- simulate_trial(40, accrual_time=20, rates_c=3, rates_t=5,
                      breaks_c = numeric(0), breaks_t = numeric(0),
                      dropout_rate = 0)

  expect_true(all(is.infinite(x$dropout_time)))
})
