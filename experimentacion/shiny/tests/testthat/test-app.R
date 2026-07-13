test_that("the Shiny server exposes the selected local series", {
  source(file.path("..", "..", "app.R"), local = TRUE)

  shiny::testServer(app_server, {
    session$setInputs(preset = "menb_pyrexia")

    expect_equal(nrow(selected_series()), 120L)
    expect_gt(max(selected_series()$maxsprt_llr), 100)
  })
})
