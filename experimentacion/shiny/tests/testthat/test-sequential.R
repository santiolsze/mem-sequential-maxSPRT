test_that("Poisson MaxSPRT uses the one-sided generalized likelihood ratio", {
  source(file.path("..", "..", "R", "sequential.R"))

  expect_equal(
    poisson_maxsprt_llr(c(5, 3), c(2, 3)),
    c(5 * log(5 / 2) + 2 - 5, 0)
  )
})

test_that("classical SPRT exposes Wald boundaries", {
  source(file.path("..", "..", "R", "sequential.R"))

  result <- classical_poisson_sprt(
    cumulative_observed = c(5, 8),
    cumulative_expected = c(2, 4),
    rr = 2, alpha = 0.05, beta = 0.2
  )

  expect_equal(result$upper[1], log(0.8 / 0.05))
  expect_equal(result$lower[1], log(0.2 / 0.95))
  expect_equal(result$llr[2], 8 * log(2) + 4 * (1 - 2))
})

test_that("monthly Poisson calibration uses a conservative discrete boundary", {
  source(file.path("..", "..", "R", "sequential.R"))

  boundary <- calibrate_monthly_poisson_boundary(
    expected = rep(1, 12), alpha = 0.001, reps = 1000, seed = 42
  )

  expect_lte(mean(boundary$null_maxima >= boundary$critical_value), 0.001)
  expect_equal(calibration_reps(0.001, 1000), 100000L)
})

test_that("binomial MaxSPRT returns zero without excess exposure", {
  source(file.path("..", "..", "R", "sequential.R"))

  expect_equal(binomial_maxsprt_llr(2, 5, matching_ratio = 1), 0)
  expect_gt(binomial_maxsprt_llr(4, 5, matching_ratio = 1), 0)
  expect_equal(binomial_critical_value(100, 1), 3.46574)
})
