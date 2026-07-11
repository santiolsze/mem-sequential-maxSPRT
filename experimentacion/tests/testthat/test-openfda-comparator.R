library(testthat)
library(arrow)

test_that("the openFDA ELIQUIS comparator is a complete monthly Parquet dataset", {
  root <- normalizePath(getwd(), mustWork = TRUE)
  while (!dir.exists(file.path(root, "experimentacion"))) {
    parent <- dirname(root)
    if (identical(parent, root)) {
      stop("Could not locate the repository root.")
    }
    root <- parent
  }
  path <- file.path(
    root, "experimentacion", "data", "processed", "fda",
    "openfda_eliquis_haemorrhage_monthly"
  )

  expect_true(dir.exists(path), info = path)
  if (!dir.exists(path)) {
    return(invisible())
  }

  monthly <- open_dataset(path, format = "parquet") |>
    dplyr::collect()

  expect_equal(nrow(monthly), 120L)
  expect_equal(
    sort(monthly$month),
    format(seq(as.Date("2016-01-01"), by = "month", length.out = 120), "%Y-%m")
  )
  expect_true(all(monthly$non_eliquis_reports > 0))
  expect_true(all(monthly$expected_under_no_eliquis_association > 0))
  expect_equal(
    monthly$expected_under_no_eliquis_association,
    monthly$eliquis_reports * monthly$non_eliquis_haemorrhage_reports /
      monthly$non_eliquis_reports,
    tolerance = 1e-12
  )
})
