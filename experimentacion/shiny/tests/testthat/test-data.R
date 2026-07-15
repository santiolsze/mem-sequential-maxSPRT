test_that("real-data presets are complete monthly series", {
  source(file.path("..", "..", "R", "sequential.R"))
  source(file.path("..", "..", "R", "data.R"))

  processed_root <- file.path("..", "..", "..", "data", "processed")
  fda <- load_real_series("eliquis", processed_root)
  hpv <- load_real_series("hpv9_syncope", processed_root)
  menb <- load_real_series("menb_pyrexia", processed_root)

  expect_equal(nrow(fda), 120L)
  expect_equal(nrow(hpv), 120L)
  expect_equal(nrow(menb), 120L)
  expect_true(all(fda$expected > 0))
  expect_true(all(hpv$expected > 0))
  expect_true(all(menb$expected > 0))
  expect_equal(fda$observed, fda$eliquis_haemorrhage_reports)
  expect_equal(
    fda$expected,
    fda$eliquis_reports * fda$non_eliquis_haemorrhage_reports /
      fda$non_eliquis_reports,
    tolerance = 1e-12
  )
  expect_equal(hpv$observed, hpv$events_target)
  expect_equal(
    hpv$expected,
    hpv$denominator * (hpv$events_comparator + 0.5) / (hpv$denominator_comparator + 1),
    tolerance = 1e-12
  )
  expect_equal(menb$observed, menb$events_target)
  expect_equal(
    menb$expected,
    menb$denominator * (menb$events_comparator + 0.5) / (menb$denominator_comparator + 1),
    tolerance = 1e-12
  )
  # MNQ is age-matched to both (same adolescent-visit vaccines): HPV9/Syncope
  # shows no excess against it, while MENB/Pyrexia still shows a clear excess.
  expect_lt(max(hpv$maxsprt_llr), 1)
  expect_true(max(menb$maxsprt_llr) > 100)
})

test_that("restarting surveillance rebuilds cumulative evidence", {
  source(file.path("..", "..", "R", "sequential.R"))
  source(file.path("..", "..", "R", "data.R"))

  series <- append_sequential_columns(data.frame(
    month_date = as.Date(c("2020-01-01", "2020-02-01", "2020-03-01")),
    observed = c(5, 2, 4),
    expected = c(2, 3, 2)
  ))
  restarted <- restart_real_series(series, as.Date("2020-02-01"))

  expect_equal(restarted$month_date, as.Date(c("2020-02-01", "2020-03-01")))
  expect_equal(restarted$cumulative_observed, c(2, 6))
  expect_equal(restarted$cumulative_expected, c(3, 5))
  expect_equal(restarted$maxsprt_llr[1], 0)
})

test_that("one-month Poisson diagnostic has no lag correlation", {
  source(file.path("..", "..", "R", "sequential.R"))
  source(file.path("..", "..", "R", "data.R"))

  gof <- poisson_gof(data.frame(observed = 2, expected = 1.5))

  expect_true(is.na(gof$lag1_residual_correlation))
})
