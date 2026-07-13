script_path <- file.path("..", "..", "scripts", "replicate_paper_tables_1_3.R")

test_that("lambert_w solves both real branches", {
  source(script_path, local = environment())

  principal <- lambert_w(1)
  lower <- lambert_w(-0.1, branch = -1)

  expect_equal(principal * exp(principal), 1, tolerance = 1e-11)
  expect_lt(lower, -1)
  expect_equal(lower * exp(lower), -0.1, tolerance = 1e-11)
  expect_equal(lambert_w(-1 / exp(1), branch = 0), -1)
  expect_equal(lambert_w(-1 / exp(1), branch = -1), -1)
})

test_that("lambert_w rejects values outside a requested real branch", {
  source(script_path, local = environment())

  expect_error(lambert_w(-0.5), "domain")
  expect_error(lambert_w(0.1, branch = -1), "domain")
  expect_error(lambert_w(-0.1, branch = 2), "branch")
})

test_that("MaxSPRT crossing times satisfy the likelihood-ratio equation", {
  source(script_path, local = environment())

  n <- c(1, 2, 10, 100)
  critical_value <- 3.952321
  crossing <- maxsprt_crossing_time(n, critical_value)

  expect_true(all(crossing > 0 & crossing < n))
  expect_equal(
    poisson_maxsprt_llr(n, crossing),
    rep(critical_value, length(n)),
    tolerance = 1e-10
  )
})

test_that("exact recursion conserves probability and stopping-time mass", {
  source(script_path, local = environment())

  result <- maxsprt_exact(critical_value = 3.5, horizon = 10, rr = 1.5)

  expect_equal(
    result$signal_probability + result$no_signal_probability,
    1,
    tolerance = 1e-11
  )
  expect_equal(
    result$expected_surveillance_time,
    result$signal_time_numerator + 10 * result$no_signal_probability,
    tolerance = 1e-11
  )
  expect_equal(
    result$conditional_signal_time,
    result$signal_time_numerator / result$signal_probability,
    tolerance = 1e-11
  )
})

test_that("exact rejection probability has the expected monotonicity", {
  source(script_path, local = environment())

  low_boundary <- maxsprt_exact(3, horizon = 10, rr = 1)$signal_probability
  high_boundary <- maxsprt_exact(4, horizon = 10, rr = 1)$signal_probability
  null <- maxsprt_exact(3.5, horizon = 10, rr = 1)$signal_probability
  alternative <- maxsprt_exact(3.5, horizon = 10, rr = 2)$signal_probability

  expect_gt(low_boundary, high_boundary)
  expect_gt(alternative, null)
})

test_that("exact calibration reproduces selected Table 1 boundaries", {
  source(script_path, local = environment())

  cases <- data.frame(
    horizon = c(0.1, 10, 1000),
    alpha = c(0.05, 0.01, 0.001),
    published = c(2.044069, 5.260513, 8.577253)
  )
  calibrated <- Map(
    calibrate_maxsprt_boundary,
    horizon = cases$horizon,
    alpha = cases$alpha
  )

  expect_equal(
    vapply(calibrated, `[[`, numeric(1), "critical_value"),
    cases$published,
    tolerance = 5e-7
  )
  expect_lt(
    max(abs(vapply(calibrated, `[[`, numeric(1), "achieved_alpha") -
      cases$alpha)),
    1e-8
  )
})

test_that("exact recursion reproduces representative Table 2 power cells", {
  source(script_path, local = environment())

  cases <- data.frame(
    horizon = c(0.1, 10, 1000),
    rr = c(1.2, 2, 1.2),
    boundary = c(2.044069, 3.467952, 4.324917),
    published = c(0.060, 0.685, 0.9996),
    tolerance = c(5e-4, 5e-4, 5e-5)
  )
  reproduced <- mapply(
    function(horizon, rr, boundary) {
      maxsprt_exact(boundary, horizon, rr)$signal_probability
    },
    cases$horizon, cases$rr, cases$boundary
  )

  expect_true(all(abs(reproduced - cases$published) < cases$tolerance))
})

test_that("exact recursion reproduces representative Table 3 time cells", {
  source(script_path, local = environment())

  short <- maxsprt_exact(2.853937, horizon = 1, rr = 1)
  middle <- maxsprt_exact(3.628123, horizon = 20, rr = 1.5)
  long <- maxsprt_exact(4.324917, horizon = 1000, rr = 1.2)

  expect_equal(short$conditional_signal_time, 0.22, tolerance = 0.005)
  expect_equal(short$expected_surveillance_time, 0.96, tolerance = 0.005)
  expect_equal(middle$conditional_signal_time, 8.68, tolerance = 0.005)
  expect_equal(middle$expected_surveillance_time, 14.91, tolerance = 0.005)
  expect_equal(long$conditional_signal_time, 196.27, tolerance = 0.005)
  expect_equal(long$expected_surveillance_time, 196.59, tolerance = 0.005)
})

test_that("table replication functions expose all published dimensions", {
  source(script_path, local = environment())

  table_2 <- replicate_table_2(paper_table_1)
  table_3 <- replicate_table_3(paper_table_1)

  expect_equal(dim(table_2), c(33, 7))
  expect_equal(dim(table_3$signal_time), c(10, 8))
  expect_equal(dim(table_3$surveillance_time), c(10, 8))
})

test_that("display-precision comparison identifies matches and failures", {
  source(script_path, local = environment())

  reference <- data.frame(horizon = 1, value = 0.060)
  display <- data.frame(horizon = "1", value = ".060")
  attr(reference, "display_values") <- display

  matching <- compare_paper_table(reference, reference, "example")
  failing <- reference
  failing$value <- 0.061
  mismatch <- compare_paper_table(failing, reference, "example")

  expect_true(all(matching$pass))
  expect_false(mismatch$pass[mismatch$column == "value"])
})

test_that("sourcing the replication script has no output side effects", {
  temporary_output <- tempfile("maxsprt-output-")
  old_option <- options(maxsprt.output_dir = temporary_output)
  on.exit(options(old_option), add = TRUE)

  source(script_path, local = environment())

  expect_false(dir.exists(temporary_output))
  expect_true(is.function(run_paper_table_replication))
})
