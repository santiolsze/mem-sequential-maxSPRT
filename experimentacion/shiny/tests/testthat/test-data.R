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
    hpv$denominator * hpv$non_target_events / hpv$non_target_denominator,
    tolerance = 1e-12
  )
  expect_equal(menb$observed, menb$events_target)
  expect_equal(
    menb$expected,
    menb$denominator * menb$non_target_events / menb$non_target_denominator,
    tolerance = 1e-12
  )
  # Contemporaneous non-target-vaccine comparator: both VAERS presets show excess reporting.
  expect_true(max(hpv$maxsprt_llr) > 100)
  expect_true(max(menb$maxsprt_llr) > 100)
})
