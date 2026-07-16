module_path <- file.path("..", "..", "R", "discretization.R")
sequential_path <- file.path("..", "..", "R", "sequential.R")

test_that("look grid is equally spaced in cumulative expected events", {
  source(module_path, local = environment())

  expect_equal(maxsprt_look_grid(12, 4), c(3, 6, 9, 12))
})

test_that("Wilson intervals remain valid at binomial extremes", {
  source(module_path, local = environment())

  none <- wilson_interval(0, 100)
  all <- wilson_interval(100, 100)

  expect_gte(none[1], 0)
  expect_lte(none[2], 1)
  expect_gte(all[1], 0)
  expect_lte(all[2], 1)
  expect_equal(unname(none[1]), 0, tolerance = 1e-12)
  expect_equal(unname(all[2]), 1, tolerance = 1e-12)
})

test_that("very fine grids cap Monte Carlo repetitions according to workload", {
  source(module_path, local = environment())

  expect_equal(discretization_effective_reps(c(1, 1000), 20000), 20000L)
  expect_equal(discretization_effective_reps(c(1, 10000), 50000), 50000L)
  expect_equal(discretization_effective_reps(c(1, 10000), 200000), 100000L)
  expect_equal(discretization_effective_reps(c(1, 25000), 20000), 5000L)
  expect_equal(discretization_effective_reps(c(1, 50000), 50000), 5000L)
  expect_equal(discretization_effective_reps(c(1, 50000), 1000), 1000L)
})

test_that("discrete crossing simulation is reproducible", {
  source(sequential_path, local = environment())
  source(module_path, local = environment())

  first <- simulate_discrete_crossing(
    horizon = 10, looks = 12, rr = 1, boundary = 3.5,
    reps = 1000, seed = 42
  )
  second <- simulate_discrete_crossing(
    horizon = 10, looks = 12, rr = 1, boundary = 3.5,
    reps = 1000, seed = 42
  )

  expect_equal(first, second)
  expect_named(
    first,
    c("rate", "standard_error", "lower", "upper", "signals", "reps")
  )
})

test_that("experiment uses one exact boundary and exposes comparison columns", {
  source(module_path, local = environment())
  calls <- new.env(parent = emptyenv())
  calls$calibration <- 0L
  calls$power <- 0L
  exact_engine <- list(
    calibrate_maxsprt_boundary = function(horizon, alpha) {
      calls$calibration <- calls$calibration + 1L
      list(critical_value = 3.5, achieved_alpha = alpha)
    },
    maxsprt_exact = function(critical_value, horizon, rr) {
      calls$power <- calls$power + 1L
      list(signal_probability = 0.4)
    }
  )

  result <- run_discretization_experiment(
    horizon = 10, alpha = 0.05, rr = 1.5,
    looks = c(12, 2, 4), reps = 500, seed = 7,
    exact_engine = exact_engine
  )

  expect_equal(calls$calibration, 1L)
  expect_equal(calls$power, 1L)
  expect_equal(result$boundary, 3.5)
  expect_equal(result$continuous_power, 0.4)
  expect_equal(result$results$looks, c(2L, 4L, 12L))
  expect_named(result$results, c(
    "looks", "alpha_experimental", "alpha_se", "alpha_lower", "alpha_upper",
    "alpha_difference", "power_experimental", "power_se", "power_lower",
    "power_upper", "power_difference"
  ))
})

test_that("null stream is independent of the power RR", {
  source(module_path, local = environment())
  exact_engine <- list(
    calibrate_maxsprt_boundary = function(horizon, alpha) {
      list(critical_value = 3.5, achieved_alpha = alpha)
    },
    maxsprt_exact = function(critical_value, horizon, rr) {
      list(signal_probability = 0.4)
    }
  )

  low_rr <- run_discretization_experiment(
    10, 0.05, 1.2, c(2, 12), 1000, 91, exact_engine
  )
  high_rr <- run_discretization_experiment(
    10, 0.05, 3, c(2, 12), 1000, 91, exact_engine
  )

  expect_equal(
    low_rr$results$alpha_experimental,
    high_rr$results$alpha_experimental
  )
  expect_false(identical(
    low_rr$results$power_experimental,
    high_rr$results$power_experimental
  ))
})

test_that("experiment rejects invalid parameter combinations", {
  source(module_path, local = environment())
  engine <- list(
    calibrate_maxsprt_boundary = function(...) list(critical_value = 3),
    maxsprt_exact = function(...) list(signal_probability = 0.5)
  )
  run <- function(horizon = 10, alpha = 0.05, rr = 1.5,
                  looks = c(2, 4), reps = 100, seed = 1) {
    run_discretization_experiment(
      horizon, alpha, rr, looks, reps, seed, engine
    )
  }

  expect_error(run(horizon = 0), "horizon")
  expect_error(run(alpha = 1), "alpha")
  expect_error(run(rr = 1), "rr")
  expect_error(run(looks = 2), "at least two")
  expect_error(run(looks = c(2, 2)), "distinct")
  expect_error(run(looks = c(2, 3.5)), "integers")
  expect_error(run(reps = 0), "reps")
  expect_error(run(seed = 1.2), "seed")
})

test_that("exact experiment accepts alpha 0.10", {
  source(module_path, local = environment())
  exact_engine <- list(
    calibrate_maxsprt_boundary = function(horizon, alpha) {
      expect_equal(alpha, 0.10)
      list(critical_value = 2.5, achieved_alpha = alpha)
    },
    maxsprt_exact = function(...) list(signal_probability = 0.6)
  )

  result <- run_discretization_experiment(
    10, 0.10, 1.5, c(2, 4), 100, 1, exact_engine
  )

  expect_equal(result$alpha, 0.10)
})

test_that("finer looks approach continuous alpha and power within Monte Carlo error", {
  source(module_path, local = environment())
  exact_engine <- new.env(parent = globalenv())
  sys.source(
    file.path("..", "..", "..", "scripts", "replicate_paper_tables_1_3.R"),
    envir = exact_engine
  )

  experiment <- run_discretization_experiment(
    horizon = 10, alpha = 0.05, rr = 1.5,
    looks = c(1, 100), reps = 10000, seed = 20260715,
    exact_engine = exact_engine
  )
  coarse <- experiment$results[1, ]
  fine <- experiment$results[2, ]

  expect_lte(fine$alpha_experimental, 0.05 + 4 * fine$alpha_se)
  expect_lte(
    abs(fine$alpha_difference),
    abs(coarse$alpha_difference) + 4 * (fine$alpha_se + coarse$alpha_se)
  )
  expect_lte(
    abs(fine$power_difference),
    abs(coarse$power_difference) + 4 * (fine$power_se + coarse$power_se)
  )
})
