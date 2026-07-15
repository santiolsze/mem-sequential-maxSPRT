test_that("the Shiny server exposes the selected local series", {
  source(file.path("..", "..", "app.R"), local = TRUE)

  shiny::testServer(app_server, {
    session$setInputs(preset = "menb_pyrexia")

    expect_equal(nrow(selected_series()), 120L)
    expect_gt(max(selected_series()$maxsprt_llr), 100)
  })
})

test_that("start month restarts every cumulative real-data calculation", {
  source(file.path("..", "..", "app.R"), local = TRUE)

  ui_text <- as.character(ui)
  expect_match(ui_text, 'id="start_month"', fixed = TRUE)
  expect_match(ui_text, "Mes inicial de vigilancia", fixed = TRUE)

  shiny::testServer(app_server, {
    session$setInputs(preset = "menb_pyrexia")
    session$setInputs(start_month = "2018-01-01")

    series <- selected_series()
    expect_equal(series$month_date[1], as.Date("2018-01-01"))
    expect_equal(series$cumulative_observed[1], series$observed[1])
    expect_equal(series$cumulative_expected[1], series$expected[1])
    expect_equal(nrow(series), 96L)

    session$setInputs(start_month = "2025-12-01")
    expect_equal(nrow(selected_series()), 1L)
    expect_true(all(nzchar(output$diagnostic_interpretation)))
  })
})

test_that("the app exposes an independent discretization experiment", {
  source(file.path("..", "..", "app.R"), local = TRUE)

  ui_text <- as.character(ui)
  expect_match(ui_text, "Discretizacion", fixed = TRUE)
  expect_match(ui_text, "disc_looks", fixed = TRUE)
  expect_match(ui_text, "run_discretization", fixed = TRUE)
  expect_match(ui_text, 'value="1000" checked="checked"', fixed = TRUE)
  expect_match(ui_text, 'value="2500"', fixed = TRUE)
  expect_match(ui_text, 'value="5000"', fixed = TRUE)
  expect_match(ui_text, 'value="10000"', fixed = TRUE)
  expect_match(ui_text, 'value="50000"', fixed = TRUE)
  expect_false(grepl('value="10000" checked="checked"', ui_text, fixed = TRUE))

  shiny::testServer(app_server, {
    session$setInputs(
      disc_t = "1",
      disc_alpha = "0.05",
      disc_rr = 1.5,
      disc_looks = c("1", "2"),
      disc_reps = "100",
      disc_seed = 20260715
    )
    session$setInputs(run_discretization = 1)

    result <- discretization_experiment()
    expect_equal(result$results$looks, c(1L, 2L))
    expect_equal(result$horizon, 1)
    expect_equal(result$alpha, 0.05)
    expect_equal(result$reps, 100L)
  })
})

test_that("simulated power labels rejection time in months", {
  app_source <- paste(
    readLines(file.path("..", "..", "app.R"), warn = FALSE),
    collapse = "\n"
  )

  expect_match(
    app_source, "Tiempo medio hasta el rechazo (meses)", fixed = TRUE
  )
})

test_that("real data exposes labelled boundaries and rejection comparison", {
  app_source <- paste(
    readLines(file.path("..", "..", "app.R"), warn = FALSE),
    collapse = "\n"
  )

  expect_match(app_source, "V MaxSPRT =", fixed = TRUE)
  expect_match(app_source, "Frontera SPRT =", fixed = TRUE)
  expect_match(app_source, "crossing_summary", fixed = TRUE)
  expect_match(app_source, "truncate_after_decision", fixed = TRUE)
  expect_match(app_source, "Resultado secuencial", fixed = TRUE)
  expect_match(app_source, "Rechaza en %s", fixed = TRUE)
  expect_match(
    app_source,
    "Finaliza sin rechazo al alcanzar T = %.2f en %s",
    fixed = TRUE
  )
})

test_that("simulated power distinguishes real and assumed risk ratios", {
  source(file.path("..", "..", "app.R"), local = TRUE)

  ui_text <- as.character(ui)
  expect_match(
    ui_text, "RR real usado para generar las simulaciones", fixed = TRUE
  )
  expect_match(ui_text, "RR supuesto por SPRT clasico 1", fixed = TRUE)
  expect_match(ui_text, "RR supuesto por SPRT clasico 2", fixed = TRUE)
})

test_that("the binomial simulation is absent from the app interface", {
  source(file.path("..", "..", "app.R"), local = TRUE)

  ui_text <- as.character(ui)
  expect_false(grepl("Binomial simulado", ui_text, fixed = TRUE))
  expect_false(grepl('id="binomial_n"', ui_text, fixed = TRUE))
  expect_false(grepl('id="matching_ratio"', ui_text, fixed = TRUE))
  expect_false(grepl('id="binomial_rr"', ui_text, fixed = TRUE))
})

test_that("Poisson diagnostics explains its objective and interpretation", {
  source(file.path("..", "..", "app.R"), local = TRUE)

  ui_text <- as.character(ui)
  expect_match(ui_text, "Objetivo del diagnostico", fixed = TRUE)
  expect_match(
    ui_text,
    "variabilidad mensual deberia ser aproximadamente igual a la media",
    fixed = TRUE
  )
  expect_match(
    ui_text, "puede no controlar correctamente el alfa nominal", fixed = TRUE
  )
  expect_match(
    ui_text, "Conclusion sobre el supuesto Poisson", fixed = TRUE
  )
})

test_that("the method tab is absent from the app interface", {
  source(file.path("..", "..", "app.R"), local = TRUE)

  ui_text <- as.character(ui)
  expect_false(grepl('data-value="Metodo"', ui_text, fixed = TRUE))
  expect_false(grepl("Metodo: ", ui_text, fixed = TRUE))
})

test_that("alpha 0.10 is available in general and discretization controls", {
  source(file.path("..", "..", "app.R"), local = TRUE)

  ui_text <- as.character(ui)
  expect_equal(
    lengths(regmatches(ui_text, gregexpr('value="0.10"', ui_text, fixed = TRUE))),
    2L
  )
})

test_that("Poisson diagnostics is the final navigation tab", {
  source(file.path("..", "..", "app.R"), local = TRUE)

  ui_text <- as.character(ui)
  method_position <- regexpr('data-value="Metodo"', ui_text, fixed = TRUE)[1]
  diagnostic_position <- regexpr(
    'data-value="Diagnostico Poisson"', ui_text, fixed = TRUE
  )[1]

  expect_gt(diagnostic_position, method_position)
})
