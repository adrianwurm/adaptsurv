test_that("single rate reduces to exponential", {
  set.seed(42)
  x <- rpwexp(1e5, rates = 0.1)
  expect_equal(mean(x), 1 / 0.1, tolerance = 0.02)   # relative tolerance
})

test_that("survival function matches exp(-H(t))", {
  set.seed(42)
  t <- 5
  t2 <- 26
  x <- rpwexp(1e5, rates = c(0.1, 0.05), breaks = 6)
  expect_equal(mean(x>t), exp(-1/2) , tolerance = 0.02)
  expect_equal((mean(x>t2)/mean(x>6)), exp(-1) , tolerance = 0.02)
})

test_that("no missing values and correct length", {
  x <- rpwexp(1000, rates = c(0.1, 0.05, 0.02), breaks = c(6, 12))
  expect_length(x, 1000)
  expect_false(anyNA(x))
})
