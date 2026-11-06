test_that("z matches survdiff chi-square", {
  set.seed(42)

  tr <- simulate_trial(400, accrual_time=20, rates_c=8, rates_t=2,
                       breaks_c = numeric(0), breaks_t = numeric(0),
                       dropout_rate = 0.5)
  y <- cut_at_events (tr,200)
  x <- analyse_fixed(y$data)

  lr <- survival::survdiff(survival::Surv(time, status) ~ arm, data = y$data)
  expect_equal(x[["z"]]^2, lr$chisq)

})
test_that("strong effect gives z < 0 and hr < 1", {

  set.seed(42)

  tr <- simulate_trial(400, accrual_time=20, rates_c=8, rates_t=2,
                       breaks_c = numeric(0), breaks_t = numeric(0),
                       dropout_rate = 0.5)
  y <- cut_at_events (tr,200)
  x <- analyse_fixed(y$data)

  a <- run_fixed(tr, 200, alpha = 0.025)


  expect_lt(x[["z"]], 0)
  expect_lt(a[["hr"]], 1)
  expect_true(a[["reject"]] %in% c(0, 1))
  expect_gt(a[["t_cut"]], 0)


})
