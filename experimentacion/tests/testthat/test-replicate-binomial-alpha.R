script_path <- file.path("..", "..", "scripts", "replicate_binomial_alpha.R")
source(script_path, local = TRUE)

test_that("binomial MaxSPRT LLR is one-sided", {
  expect_equal(binomial_maxsprt_llr_exact(2, 4, z = 1), 0)
  expect_equal(binomial_maxsprt_llr_exact(4, 4, z = 1), 4 * log(2))
  expect_gt(binomial_maxsprt_llr_exact(4, 5, z = 1), 0)
})

test_that("Markov recursion conserves null probability", {
  result <- binomial_alpha_markov(V = 2.77259, N = 10, z = 1)
  totals <- with(result$history, alive_probability + absorbed_probability)

  expect_equal(totals, rep(1, 10), tolerance = 1e-14)
  expect_equal(result$alpha, tail(result$history$absorbed_probability, 1))
})

test_that("Markov alpha agrees with direct enumeration of all event paths", {
  V <- 2.77259
  markov <- binomial_alpha_markov(V = V, N = 10, z = 1)$alpha
  direct <- binomial_alpha_enumeration(V = V, N = 10, z = 1)

  expect_equal(markov, direct, tolerance = 1e-14)
})

test_that("calibration reproduces Table 4 for matching 1 to 1 and N 10", {
  calibration <- calibrate_binomial_boundary(alpha = 0.05, N = 10, z = 1)

  expect_equal(calibration$critical_value_display, 2.77259)
  expect_lte(calibration$actual_alpha, 0.05)
  expect_gt(calibration$alpha_if_equality_included, 0.05)
})
