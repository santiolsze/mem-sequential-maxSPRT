library(shiny)
library(bslib)
library(ggplot2)
library(plotly)
library(DT)

app_file <- tryCatch(sys.frame(1)$ofile, error = function(error) NULL)
candidate_app_dirs <- unique(c(
  if (!is.null(app_file)) dirname(normalizePath(app_file, mustWork = FALSE)),
  normalizePath(getwd(), mustWork = FALSE),
  normalizePath(file.path(getwd(), "..", ".."), mustWork = FALSE),
  normalizePath("experimentacion/shiny", mustWork = FALSE)
))
app_dir <- candidate_app_dirs[
  dir.exists(file.path(candidate_app_dirs, "R"))
][1]

if (is.na(app_dir)) {
  stop("Could not locate the Shiny application directory.")
}

or_else <- function(value, fallback) {
  if (is.null(value) || length(value) == 0) fallback else value
}

source(file.path(app_dir, "R", "sequential.R"))
source(file.path(app_dir, "R", "data.R"))

processed_root <- normalizePath(
  file.path(app_dir, "..", "data", "processed"),
  mustWork = TRUE
)
presets <- real_series_presets()

metric <- function(label, value, detail = NULL) {
  tags$div(
    class = "metric",
    tags$div(class = "metric-label", label),
    tags$div(class = "metric-value", value),
    if (!is.null(detail)) tags$div(class = "metric-detail", detail)
  )
}

ui <- page_sidebar(
  title = tags$div(
    class = "app-title",
    tags$span("MaxSPRT Lab"),
    tags$small("Vigilancia secuencial con series reales y simuladas")
  ),
  sidebar = sidebar(
    width = 320,
    tags$div(class = "sidebar-section", "Serie real"),
    selectInput("preset", "Preset", choices = stats::setNames(presets$id, presets$label)),
    radioButtons("alpha", "Alpha", choices = c("0.05", "0.01", "0.001"), selected = "0.05", inline = TRUE),
    sliderInput("beta", "Beta (SPRT clasico)", min = 0.05, max = 0.40, value = 0.20, step = 0.05),
    numericInput("rr_low", "RR alternativo bajo", value = 1.2, min = 1.01, step = 0.05),
    numericInput("rr_high", "RR alternativo alto", value = 3, min = 1.05, step = 0.1),
    sliderInput("true_rr", "RR verdadero simulado", min = 1.05, max = 3, value = 1.5, step = 0.05),
    selectInput("simulation_t", "T simulado", choices = c(50, 100, 500), selected = 100),
    selectInput("simulation_reps", "Repeticiones de potencia", choices = c("1,000" = 1000, "5,000" = 5000), selected = 1000),
    actionButton("run_power", "Simular potencia", class = "btn-primary"),
    tags$hr(),
    tags$div(class = "sidebar-section", "Binomial simulado"),
    selectInput("binomial_n", "Eventos observados N", choices = c(20, 50, 100, 200, 500, 1000), selected = 100),
    selectInput("matching_ratio", "Ratio expuesto:no expuesto", choices = c("1:1" = 1, "1:2" = 2, "1:3" = 3), selected = 1),
    sliderInput("binomial_rr", "RR verdadero", min = 1, max = 4, value = 1.8, step = 0.1)
  ),
  tags$head(includeCSS(file.path(app_dir, "www", "app.css"))),
  tags$div(
    class = "caveat-strip",
    "FAERS y VAERS son sistemas de reportes espontaneos. Los paneles reales describen reportes; no estiman incidencia ni causalidad."
  ),
  navset_tab(
    nav_panel(
      "Guia",
      tags$p(class = "tab-intro", "Como usar la app y que representa cada preset de datos reales."),
      layout_columns(
        card(
          card_header("Uso basico"),
          tags$ol(
            tags$li("Elegi un preset de serie real en la barra lateral (Eliquis, HPV9 o MENB)."),
            tags$li("Ajusta Alpha (tasa de falsos positivos nominal del MaxSPRT) y Beta (para el SPRT clasico)."),
            tags$li("Defini RR alternativo bajo y alto: son las hipotesis alternativas fijas del SPRT clasico."),
            tags$li("La pestana 'Datos reales' se recalcula sola al cambiar cualquiera de esos parametros."),
            tags$li("En 'Potencia simulada' elegi el RR verdadero simulado, T simulado y repeticiones, y apreta 'Simular potencia'."),
            tags$li("En 'Binomial simulado' usa los controles propios de esa seccion (N, ratio expuesto:no expuesto, RR verdadero).")
          )
        ),
        card(
          card_header("Que mirar en cada pestana"),
          tags$ul(
            tags$li(tags$strong("Datos reales: "), "observados vs esperados, trayectorias LLR (MaxSPRT y SPRT clasico) y la tabla mensual."),
            tags$li(tags$strong("Diagnostico Poisson: "), "chequea si el nulo Poisson simple es razonable para la serie elegida."),
            tags$li(tags$strong("Potencia simulada: "), "ilustra formalmente alpha, beta, V y el stopping time con datos simulados."),
            tags$li(tags$strong("Binomial simulado: "), "ejemplo didactico expuesto/no expuesto, independiente de los presets reales."),
            tags$li(tags$strong("Metodo: "), "formulas y supuestos de MaxSPRT y del SPRT clasico.")
          )
        ),
        col_widths = c(6, 6)
      ),
      card(
        card_header("Presets de serie real"),
        tags$table(
          class = "table table-sm",
          tags$thead(tags$tr(
            tags$th("Preset"), tags$th("Fuente"), tags$th("Observado"), tags$th("Esperado")
          )),
          tags$tbody(
            tags$tr(
              tags$td("ELIQUIS / Haemorrhage"),
              tags$td("FAERS/openFDA"),
              tags$td("Reportes mensuales que mencionan ELIQUIS y Hemorragia"),
              tags$td("Fraccion de Hemorragia en reportes no-ELIQUIS del mismo mes, aplicada a los reportes de ELIQUIS")
            ),
            tags$tr(
              tags$td("HPV9 (GARDASIL 9) / Syncope"),
              tags$td("VAERS"),
              tags$td("Reportes de Sincope tras HPV9 por mes"),
              tags$td("Fraccion de Sincope en reportes de otras vacunas del mismo mes, aplicada a los reportes de HPV9")
            ),
            tags$tr(
              tags$td("MENB / Pyrexia"),
              tags$td("VAERS"),
              tags$td("Reportes de Pyrexia tras MENB por mes"),
              tags$td("Fraccion de Pyrexia en reportes de otras vacunas del mismo mes, aplicada a los reportes de MENB")
            )
          )
        ),
        tags$p(
          class = "diagnostic-warning",
          "FAERS y VAERS son sistemas de reportes espontaneos: 'esperado' es una comparacion contemporanea contra la propia tasa de reportes de otra vacuna o droga, no una incidencia clinica ni un denominador poblacional real. Sincope y Pyrexia son eventos frecuentes tras cualquier vacunacion, no especificos de HPV9 o MENB."
        )
      )
    ),
    nav_panel(
      "Datos reales",
      uiOutput("series_summary"),
      layout_columns(
        card(card_header("Reportes observados y esperados"), plotlyOutput("observed_plot", height = "330px")),
        card(card_header("Trayectorias secuenciales"), plotlyOutput("llr_plot", height = "330px")),
        col_widths = c(6, 6)
      ),
      card(card_header("Serie mensual"), DTOutput("monthly_table"))
    ),
    nav_panel(
      "Diagnostico Poisson",
      tags$p(class = "tab-intro", "Chequeo retrospectivo del nulo Poisson simple usado por la serie. Si el ajuste falla, el alpha nominal del MaxSPRT no queda garantizado."),
      uiOutput("gof_summary"),
      layout_columns(
        card(card_header("Residuos de Pearson por mes"), plotOutput("residual_plot", height = "340px")),
        card(card_header("Lectura"), uiOutput("diagnostic_interpretation")),
        col_widths = c(7, 5)
      )
    ),
    nav_panel(
      "Potencia simulada",
      tags$p(class = "tab-intro", "Simulacion Poisson con esperado conocido. Este es el panel donde se ilustran formalmente V, alpha, beta y el efecto de elegir alternativas clasicas altas o bajas."),
      uiOutput("power_results")
    ),
    nav_panel(
      "Binomial simulado",
      tags$p(class = "tab-intro", "Ejemplo didactico de casos expuestos/no expuestos. No se deriva de VAERS."),
      uiOutput("binomial_summary"),
      card(card_header("MaxSPRT binomial"), plotOutput("binomial_plot", height = "380px"))
    ),
    nav_panel(
      "Metodo",
      layout_columns(
        card(
          card_header("Poisson"),
          tags$p("MaxSPRT maximiza el LLR sobre RR > 1."),
          tags$pre("LLR_max = C log(C / Mu) + Mu - C, si C > Mu"),
          tags$p("V se calibra bajo H0 para controlar la probabilidad de cruzar al menos una vez. En esta app, para observaciones mensuales, se obtiene por Monte Carlo con semilla fija.")
        ),
        card(
          card_header("SPRT clasico"),
          tags$pre("LLR_r = C log(r) + Mu (1 - r)"),
          tags$p("Frontera superior: log((1-beta)/alpha). Frontera inferior: log(beta/(1-alpha)). A diferencia de MaxSPRT, puede aceptar H0 temprano para una alternativa fija.")
        ),
        card(
          card_header("Amenazas a la validez"),
          tags$ul(
            tags$li("Los esperados locales se estiman desde reportes y no son tasas clinicas conocidas."),
            tags$li("Sobredispersion o dependencia temporal invalidan la calibracion Poisson simple."),
            tags$li("Los presets se eligieron de forma exploratoria; alpha no controla esa seleccion multiple.")
          )
        ),
        col_widths = c(4, 4, 4)
      )
    )
  )
)

app_server <- function(input, output, session) {
  selected_series <- reactive({
    load_real_series(or_else(input$preset, "eliquis"), processed_root)
  })

  boundary <- reactive({
    series <- selected_series()
    alpha <- as.numeric(or_else(input$alpha, "0.05"))
    calibrate_monthly_poisson_boundary(
      series$expected,
      alpha = alpha,
      reps = calibration_reps(
        alpha,
        as.integer(or_else(input$simulation_reps, "1000"))
      ),
      seed = 20260711L
    )
  })

  classical_paths <- reactive({
    rr_low <- as.numeric(or_else(input$rr_low, 1.2))
    rr_high <- as.numeric(or_else(input$rr_high, 3))
    alpha <- as.numeric(or_else(input$alpha, "0.05"))
    beta <- as.numeric(or_else(input$beta, 0.2))
    validate(need(rr_high > rr_low, "El RR alto debe ser mayor que el RR bajo."))
    series <- selected_series()
    data.frame(
      low = classical_poisson_sprt(
        series$cumulative_observed, series$cumulative_expected,
        rr_low, alpha, beta
      )$llr,
      high = classical_poisson_sprt(
        series$cumulative_observed, series$cumulative_expected,
        rr_high, alpha, beta
    )$llr
    )
  })

  output$series_summary <- renderUI({
    series <- selected_series()
    bound <- boundary()
    signal <- first_crossing(series$maxsprt_llr, bound$critical_value)
    tagList(
      tags$div(
        class = "metric-grid",
        metric("Fuente", series$source[1]),
        metric("Observados", format(tail(series$cumulative_observed, 1), big.mark = ",")),
        metric("Esperados", format(round(tail(series$cumulative_expected, 1), 1), big.mark = ",")),
        metric("LLR maximo", round(max(series$maxsprt_llr), 2)),
        metric("V mensual", round(bound$critical_value, 2), paste(bound$reps, "simulaciones nulas")),
        metric(
          "Resultado mecanico",
          if (signal$decision == "reject") "Cruza V*" else "No cruza V",
          "* sin garantia formal de alpha"
        )
      ),
      tags$p(class = "metric-detail", series$interpretation[1])
    )
  })

  output$observed_plot <- renderPlotly({
    series <- selected_series()
    plot <- ggplot(series, aes(x = month_date)) +
      geom_line(aes(y = observed, color = "Observados"), linewidth = 0.8) +
      geom_line(aes(y = expected, color = "Esperados"), linewidth = 0.8, linetype = "dashed") +
      scale_color_manual(values = c("Observados" = "#0e7490", "Esperados" = "#b45309")) +
      labs(x = NULL, y = "Reportes por mes", color = NULL) +
      theme_minimal(base_size = 12) + theme(legend.position = "top")
    ggplotly(plot, tooltip = c("x", "y", "colour"))
  })

  output$llr_plot <- renderPlotly({
    series <- selected_series()
    paths <- classical_paths()
    bound <- boundary()
    plot_data <- data.frame(
      month_date = rep(series$month_date, 3),
      llr = c(series$maxsprt_llr, paths$low, paths$high),
      method = rep(c("MaxSPRT", "SPRT RR bajo", "SPRT RR alto"), each = nrow(series))
    )
    plot <- ggplot(plot_data, aes(month_date, llr, color = method)) +
      geom_hline(yintercept = bound$critical_value, linetype = "dashed", color = "#b91c1c") +
      geom_hline(
        yintercept = log((1 - as.numeric(or_else(input$beta, 0.2))) /
          as.numeric(or_else(input$alpha, "0.05"))),
        linetype = "dotted", color = "#475569"
      ) +
      geom_line(linewidth = 0.8) +
      scale_color_manual(values = c("MaxSPRT" = "#0e7490", "SPRT RR bajo" = "#6d28d9", "SPRT RR alto" = "#b45309")) +
      labs(x = NULL, y = "LLR acumulado", color = NULL) +
      theme_minimal(base_size = 12) + theme(legend.position = "top")
    ggplotly(plot, tooltip = c("x", "y", "colour"))
  })

  output$monthly_table <- renderDT({
    series <- selected_series()
    datatable(
      transform(
        series[, c("month_date", "observed", "expected", "cumulative_rr", "maxsprt_llr")],
        expected = round(expected, 2), cumulative_rr = round(cumulative_rr, 3), maxsprt_llr = round(maxsprt_llr, 3)
      ),
      rownames = FALSE,
      options = list(pageLength = 10, dom = "tip")
    )
  })

  output$gof_summary <- renderUI({
    gof <- poisson_gof(selected_series())
    tags$div(
      class = "metric-grid",
      metric("Deviance", round(gof$deviance, 1)),
      metric("Pearson", round(gof$pearson, 1)),
      metric("Dispersion", round(gof$dispersion, 2), "Poisson ideal: cerca de 1"),
      metric("Correlacion lag-1", round(gof$lag1_residual_correlation, 2))
    )
  })

  output$residual_plot <- renderPlot({
    series <- selected_series()
    residuals <- (series$observed - series$expected) / sqrt(series$expected)
    plot_data <- data.frame(month_date = series$month_date, residual = residuals)
    ggplot(plot_data, aes(month_date, residual)) +
      geom_hline(yintercept = c(-2, 0, 2), linetype = c("dashed", "solid", "dashed"), color = c("#94a3b8", "#334155", "#94a3b8")) +
      geom_line(color = "#0e7490", linewidth = 0.7) +
      geom_point(color = "#0e7490", size = 1.4) +
      labs(x = NULL, y = "Residuo de Pearson") +
      theme_minimal(base_size = 12)
  })

  output$diagnostic_interpretation <- renderUI({
    gof <- poisson_gof(selected_series())
    if (gof$dispersion > 1.5 || abs(gof$lag1_residual_correlation) > 0.25) {
      tags$p(class = "diagnostic-warning", "El nulo Poisson simple no describe bien esta serie. El cruce de un umbral se debe interpretar como demostracion mecanica, no como una garantia formal de alpha.")
    } else {
      tags$p("El ajuste Poisson simple no muestra una desviacion grande en este resumen. Sigue siendo un analisis de reportes, no de incidencia clinica.")
    }
  })

  simulations <- eventReactive(input$run_power, {
    series <- selected_series()
    simulation_expected <- series$expected / sum(series$expected) *
      as.numeric(or_else(input$simulation_t, "100"))
    simulation_boundary <- calibrate_monthly_poisson_boundary(
      simulation_expected,
      alpha = as.numeric(or_else(input$alpha, "0.05")),
      reps = calibration_reps(
        as.numeric(or_else(input$alpha, "0.05")),
        as.integer(or_else(input$simulation_reps, "1000"))
      ),
      seed = 20260711L
    )
    simulation <- simulate_poisson_comparison(
      expected = simulation_expected,
      true_rr = as.numeric(or_else(input$true_rr, 1.5)),
      alpha = as.numeric(or_else(input$alpha, "0.05")),
      beta = as.numeric(or_else(input$beta, 0.2)),
      rr_low = as.numeric(or_else(input$rr_low, 1.2)),
      rr_high = as.numeric(or_else(input$rr_high, 3)),
      critical_value = simulation_boundary$critical_value,
      reps = as.integer(or_else(input$simulation_reps, "1000")),
      seed = 20260711L
    )
    summary <- rbind(
      summarise_simulation(simulation, "maxsprt"),
      summarise_simulation(simulation, "low"),
      summarise_simulation(simulation, "high")
    )
    attr(summary, "critical_value") <- simulation_boundary$critical_value
    attr(summary, "simulation_t") <- sum(simulation_expected)
    summary
  }, ignoreInit = TRUE)

  output$power_results <- renderUI({
    validate(need(!is.null(simulations()), "Elegí los parametros y presioná Simular potencia."))
    result <- simulations()
    tagList(
      tags$p(
        class = "tab-intro",
        sprintf(
          "Perfil temporal del preset reescalado a T = %.0f; V mensual = %.2f.",
          attr(result, "simulation_t"), attr(result, "critical_value")
        )
      ),
      tags$div(
        class = "power-table",
        tags$table(
          class = "table table-sm",
          tags$thead(tags$tr(tags$th("Metodo"), tags$th("Prob. rechazo"), tags$th("Aceptacion temprana"), tags$th("Look medio de detencion"))),
          tags$tbody(lapply(seq_len(nrow(result)), function(index) {
            row <- result[index, ]
            tags$tr(
              tags$td(row$method),
              tags$td(sprintf("%.1f%%", 100 * row$rejection_rate)),
              tags$td(sprintf("%.1f%%", 100 * row$early_acceptance_rate)),
              tags$td(ifelse(is.nan(row$mean_stopping_look), "-", round(row$mean_stopping_look, 1)))
            )
          }))
        )
      )
    )
  })

  binomial_trajectory <- reactive({
    simulate_binomial_trajectory(
      total_cases = as.integer(or_else(input$binomial_n, "100")),
      matching_ratio = as.integer(or_else(input$matching_ratio, "1")),
      true_rr = as.numeric(or_else(input$binomial_rr, 1.8)),
      seed = 20260711L
    )
  })

  output$binomial_summary <- renderUI({
    trajectory <- binomial_trajectory()
    total_cases <- as.integer(or_else(input$binomial_n, "100"))
    matching_ratio <- as.integer(or_else(input$matching_ratio, "1"))
    critical <- binomial_critical_value(total_cases, matching_ratio)
    tags$div(
      class = "metric-grid",
      metric("N maximo", total_cases),
      metric("RR verdadero", as.numeric(or_else(input$binomial_rr, 1.8))),
      metric("LLR maximo", round(max(trajectory$llr), 2)),
      metric("V Tabla 4", critical, "alpha = 0.05")
    )
  })

  output$binomial_plot <- renderPlot({
    trajectory <- binomial_trajectory()
    critical <- binomial_critical_value(
      as.integer(or_else(input$binomial_n, "100")),
      as.integer(or_else(input$matching_ratio, "1"))
    )
    ggplot(trajectory, aes(case, llr)) +
      geom_hline(yintercept = critical, linetype = "dashed", color = "#b91c1c") +
      geom_line(color = "#6d28d9", linewidth = 0.8) +
      labs(x = "Eventos observados", y = "LLR binomial acumulado") +
      theme_minimal(base_size = 12)
  })
}

shinyApp(ui, app_server)
