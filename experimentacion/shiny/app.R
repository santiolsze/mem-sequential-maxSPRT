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
source(file.path(app_dir, "R", "discretization.R"))

exact_maxsprt <- new.env(parent = globalenv())
sys.source(
  file.path(app_dir, "..", "scripts", "replicate_paper_tables_1_3.R"),
  envir = exact_maxsprt
)

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
    selectInput("preset", "Serie", choices = stats::setNames(presets$id, presets$label)),
    selectInput(
      "start_month", "Mes inicial de vigilancia", choices = character()
    ),
    radioButtons("alpha", "Alpha", choices = c("0.10", "0.05", "0.01", "0.001"), selected = "0.05", inline = TRUE),
    sliderInput("beta", "Beta (SPRT clasico)", min = 0.05, max = 0.40, value = 0.20, step = 0.05),
    numericInput("rr_low", "RR supuesto por SPRT clasico 1", value = 1.2, min = 1.01, step = 0.05),
    numericInput("rr_high", "RR supuesto por SPRT clasico 2", value = 3, min = 1.05, step = 0.1),
    sliderInput("true_rr", "RR real usado para generar las simulaciones", min = 1.05, max = 3, value = 1.5, step = 0.05),
    selectInput("simulation_t", "T simulado", choices = c(50, 100, 500), selected = 100),
    selectInput("simulation_reps", "Repeticiones de potencia", choices = c("1,000" = 1000, "5,000" = 5000), selected = 1000),
    actionButton("run_power", "Simular potencia", class = "btn-primary")
  ),
  tags$head(includeCSS(file.path(app_dir, "www", "app.css"))),
  navset_tab(
    nav_panel(
      "Guia",
      tags$p(class = "tab-intro", "Como usar la app y que representa cada serie de datos reales."),
      layout_columns(
        card(
          card_header("Uso basico"),
          tags$ol(
            tags$li("Elegi una serie real en la barra lateral (Eliquis, HPV9 o MENB)."),
            tags$li("Elegi el mes inicial: los acumulados y las fronteras se recalculan desde ese mes."),
            tags$li("Ajusta Alpha (tasa de falsos positivos nominal del MaxSPRT) y Beta (para el SPRT clasico)."),
            tags$li("Defini los RR supuestos por los dos SPRT clasicos: son sus hipotesis alternativas fijas."),
            tags$li("La pestana 'Datos reales' se recalcula sola al cambiar cualquiera de esos parametros."),
            tags$li("En 'Potencia simulada' elegi el RR real usado para generar los datos, T simulado y repeticiones, y apreta 'Simular potencia'.")
          )
        ),
        card(
          card_header("Que mirar en cada pestana"),
          tags$ul(
            tags$li(tags$strong("Datos reales: "), "observados vs esperados, trayectorias LLR (MaxSPRT y SPRT clasico) y la tabla mensual."),
            tags$li(tags$strong("Binomial MENB vs MNQ: "), "compara la proporcion de reportes con Pyrexia entre ambas vacunas, condicionando por mes."),
            tags$li(tags$strong("Potencia simulada: "), "ilustra formalmente alpha, beta, V y el stopping time con datos simulados."),
            tags$li(tags$strong("Discretizacion: "), "compara alpha y potencia continuos exactos contra observaciones en una cantidad finita de looks, sin usar datos reales.")
          )
        ),
        col_widths = c(6, 6)
      ),
      card(
        card_header("Series reales"),
        tags$table(
          class = "table table-sm",
          tags$thead(tags$tr(
            tags$th("Serie"), tags$th("Fuente"), tags$th("Observado"), tags$th("Esperado")
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
              tags$td("Fraccion de Sincope en reportes de MNQ del mismo mes, aplicada a los reportes de HPV9")
            ),
            tags$tr(
              tags$td("MENB / Pyrexia"),
              tags$td("VAERS"),
              tags$td("Reportes de Pyrexia tras MENB por mes"),
              tags$td("Fraccion de Pyrexia en reportes de MNQ del mismo mes, aplicada a los reportes de MENB")
            )
          )
        ),
        tags$p(
          class = "diagnostic-warning",
          "FAERS y VAERS son sistemas de reportes espontaneos: 'esperado' es una comparacion contemporanea contra la propia tasa de reportes de otra vacuna, no una incidencia clinica ni un denominador poblacional real. MNQ (meningococica ACWY) se eligio como comparador de HPV9 y MENB por tener una distribucion de edad similar (misma visita de vacunacion adolescente); comparar contra un pool de todas las vacunas mezclaria el efecto de la vacuna con diferencias de edad entre programas de vacunacion."
        )
      )
    ),
    nav_panel(
      "Datos reales",
      uiOutput("series_summary"),
      layout_columns(
        card(card_header("Reportes observados y esperados"), plotlyOutput("observed_plot", height = "330px")),
        card(
          card_header("Trayectorias secuenciales"),
          plotlyOutput("llr_plot", height = "330px"),
          uiOutput("crossing_summary")
        ),
        col_widths = c(6, 6)
      ),
      card(card_header("Serie mensual"), DTOutput("monthly_table"))
    ),
    nav_panel(
      "Binomial MENB vs MNQ",
      tags$div(
        class = "tab-intro",
        tags$h4("Comparacion fija de Pyrexia en VAERS"),
        tags$p(
          "Contrasta si la proporcion de reportes con Pyrexia es mayor para MENB que para MNQ. En cada mes condiciona por el total de reportes de Pyrexia y ajusta la probabilidad nula por el volumen de reportes de cada vacuna."
        ),
        tags$p(
          "La frontera Monte Carlo controla alpha condicionalmente a esos margenes mensuales. Este analisis describe desproporcionalidad de reporte; no estima incidencia, riesgo clinico ni causalidad."
        )
      ),
      uiOutput("binomial_summary"),
      layout_columns(
        card(
          card_header("Proporcion mensual de reportes con Pyrexia"),
          plotlyOutput("binomial_proportion_plot", height = "340px")
        ),
        card(
          card_header("MaxSPRT binomial condicional"),
          plotlyOutput("binomial_llr_plot", height = "340px"),
          uiOutput("binomial_decision")
        ),
        col_widths = c(6, 6)
      ),
      card(
        card_header("Comparacion mensual MENB vs MNQ"),
        DTOutput("binomial_monthly_table")
      )
    ),
    nav_panel(
      "Potencia simulada",
      tags$p(class = "tab-intro", "Simulacion Poisson con esperado conocido. Este es el panel donde se ilustran formalmente V, alpha, beta y el efecto de elegir alternativas clasicas altas o bajas."),
      uiOutput("power_results")
    ),
    nav_panel(
      "Discretizacion",
      tags$p(
        class = "tab-intro",
        "Experimento Poisson independiente de FAERS y VAERS. Los K looks estan igualmente espaciados en eventos esperados acumulados: mu_k = kT/K. La misma frontera continua se usa para todo K, de modo que se vea el efecto conservador de observar con menor frecuencia."
      ),
      layout_sidebar(
        sidebar = sidebar(
          width = 300,
          selectInput("disc_t", "Horizonte T", choices = c(10, 50, 100), selected = 50),
          radioButtons(
            "disc_alpha", "Alpha nominal",
            choices = c("0.10", "0.05", "0.01", "0.001"),
            selected = "0.05", inline = TRUE
          ),
          sliderInput(
            "disc_rr", "RR verdadero para potencia",
            min = 1.05, max = 3, value = 1.5, step = 0.05
          ),
          checkboxGroupInput(
            "disc_looks", "Cantidad de looks K",
            choices = c(1, 2, 4, 12, 24, 52, 100, 365, 1000, 2500, 5000, 10000, 25000, 50000),
            selected = c(1, 2, 4, 12, 24, 52, 100, 365, 1000)
          ),
          selectInput(
            "disc_reps", "Replicas Monte Carlo",
            choices = c(
              "1,000" = 1000,
              "5,000" = 5000,
              "10,000" = 10000,
              "20,000" = 20000,
              "50,000" = 50000
            ),
            selected = 5000
          ),
          numericInput(
            "disc_seed", "Semilla", value = 20260715,
            min = 0, step = 1
          ),
          actionButton(
            "run_discretization", "Comparar discretizaciones",
            class = "btn-primary"
          )
        ),
        uiOutput("discretization_workload"),
        uiOutput("discretization_metrics"),
        layout_columns(
          card(
            card_header("Alpha: continuo vs. looks discretos"),
            plotlyOutput("discretization_alpha_plot", height = "350px")
          ),
          card(
            card_header("Potencia: continuo vs. looks discretos"),
            plotlyOutput("discretization_power_plot", height = "350px")
          ),
          col_widths = c(6, 6)
        ),
        card(
          card_header("Resultados por cantidad de looks"),
          DTOutput("discretization_table")
        )
      )
    )
  )
)

app_server <- function(input, output, session) {
  binomial_series <- reactive({
    load_menb_mnq_pyrexia_binomial(processed_root)
  })

  binomial_analysis <- reactive({
    series <- binomial_series()
    stratified_binomial_maxsprt(
      series$target_events,
      series$comparator_events,
      series$target_reports,
      series$comparator_reports
    )
  })

  binomial_boundary <- reactive({
    series <- binomial_series()
    alpha <- as.numeric(or_else(input$alpha, "0.05"))
    calibrate_stratified_binomial_boundary(
      total_events = series$target_events + series$comparator_events,
      target_reports = series$target_reports,
      comparator_reports = series$comparator_reports,
      alpha = alpha,
      reps = calibration_reps(
        alpha,
        as.integer(or_else(input$simulation_reps, "1000"))
      ),
      seed = 20260715L
    )
  })

  raw_series <- reactive({
    load_real_series(or_else(input$preset, "eliquis"), processed_root)
  })

  observeEvent(raw_series(), {
    series <- raw_series()
    month_values <- format(series$month_date, "%Y-%m-%d")
    month_labels <- format(series$month_date, "%Y-%m")
    current <- isolate(input$start_month)
    selected <- if (!is.null(current) && current %in% month_values) {
      current
    } else {
      month_values[1]
    }
    updateSelectInput(
      session, "start_month",
      choices = stats::setNames(month_values, month_labels),
      selected = selected
    )
  }, ignoreInit = FALSE)

  selected_series <- reactive({
    series <- raw_series()
    month_values <- format(series$month_date, "%Y-%m-%d")
    start_month <- or_else(input$start_month, month_values[1])
    if (!start_month %in% month_values) {
      start_month <- month_values[1]
    }
    restart_real_series(series, as.Date(start_month))
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

  real_crossing_summary <- reactive({
    series <- selected_series()
    paths <- classical_paths()
    maxsprt_boundary <- boundary()$critical_value
    alpha <- as.numeric(or_else(input$alpha, "0.05"))
    beta <- as.numeric(or_else(input$beta, 0.2))
    sprt_upper <- log((1 - beta) / alpha)
    sprt_lower <- log(beta / (1 - alpha))
    rr_low <- as.numeric(or_else(input$rr_low, 1.2))
    rr_high <- as.numeric(or_else(input$rr_high, 3))

    definitions <- list(
      list(
        method = "MaxSPRT", values = series$maxsprt_llr,
        upper = maxsprt_boundary, lower = NULL
      ),
      list(
        method = sprintf("SPRT supuesto RR %g", rr_low), values = paths$low,
        upper = sprt_upper, lower = sprt_lower
      ),
      list(
        method = sprintf("SPRT supuesto RR %g", rr_high), values = paths$high,
        upper = sprt_upper, lower = sprt_lower
      )
    )

    rows <- lapply(definitions, function(definition) {
      result <- sequential_decision_summary(
        definition$values, series$month_date,
        definition$upper, definition$lower
      )
      outcome <- if (result$decision == "reject") {
        sprintf("Rechaza en %s", format(result$first_rejection, "%Y-%m"))
      } else if (definition$method == "MaxSPRT") {
        sprintf(
          "Finaliza sin rechazo al alcanzar T = %.2f en %s",
          tail(series$cumulative_expected, 1),
          format(tail(series$month_date, 1), "%Y-%m")
        )
      } else if (result$decision == "accept") {
        sprintf("No rechaza (acepta H0 en %s)", format(result$decision_date, "%Y-%m"))
      } else {
        sprintf("Finaliza sin decision en %s", format(tail(series$month_date, 1), "%Y-%m"))
      }
      data.frame(
        method = definition$method,
        boundary = definition$upper,
        outcome = outcome,
        crossing_index = result$crossing_index,
        crossing_value = result$crossing_value,
        stringsAsFactors = FALSE
      )
    })
    do.call(rbind, rows)
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
    rr_low <- as.numeric(or_else(input$rr_low, 1.2))
    rr_high <- as.numeric(or_else(input$rr_high, 3))
    method_names <- c(
      "MaxSPRT",
      sprintf("SPRT supuesto RR %g", rr_low),
      sprintf("SPRT supuesto RR %g", rr_high)
    )
    sprt_boundary <- log((1 - as.numeric(or_else(input$beta, 0.2))) /
      as.numeric(or_else(input$alpha, "0.05")))
    sprt_lower <- log(as.numeric(or_else(input$beta, 0.2)) /
      (1 - as.numeric(or_else(input$alpha, "0.05"))))
    plot_data <- data.frame(
      month_date = rep(series$month_date, 3),
      llr = c(
        truncate_after_decision(
          series$maxsprt_llr, upper = bound$critical_value
        ),
        truncate_after_decision(
          paths$low, upper = sprt_boundary, lower = sprt_lower
        ),
        truncate_after_decision(
          paths$high, upper = sprt_boundary, lower = sprt_lower
        )
      ),
      method = rep(method_names, each = nrow(series))
    )
    crossings <- real_crossing_summary()
    crossing_points <- crossings[!is.na(crossings$crossing_index), ]
    crossing_points$month_date <- series$month_date[crossing_points$crossing_index]
    boundary_labels <- data.frame(
      month_date = max(series$month_date),
      llr = c(bound$critical_value, sprt_boundary),
      label = c(
        sprintf("V MaxSPRT = %.2f", bound$critical_value),
        sprintf("Frontera SPRT = %.2f", sprt_boundary)
      )
    )
    method_colors <- stats::setNames(
      c("#0e7490", "#6d28d9", "#b45309"), method_names
    )
    plot <- ggplot(plot_data, aes(month_date, llr, color = method)) +
      geom_hline(yintercept = bound$critical_value, linetype = "dashed", color = "#b91c1c") +
      geom_hline(
        yintercept = sprt_boundary,
        linetype = "dotted", color = "#475569"
      ) +
      geom_line(linewidth = 0.8) +
      geom_point(
        data = crossing_points,
        aes(month_date, crossing_value, color = method),
        inherit.aes = FALSE, size = 3.2, shape = 21, fill = "white", stroke = 1.2
      ) +
      geom_text(
        data = boundary_labels,
        aes(month_date, llr, label = label),
        inherit.aes = FALSE, hjust = 1, vjust = -0.4,
        color = "#334155", size = 3.2
      ) +
      scale_color_manual(values = method_colors) +
      labs(x = NULL, y = "LLR acumulado", color = NULL) +
      theme_minimal(base_size = 12) + theme(legend.position = "top")
    ggplotly(plot, tooltip = c("x", "y", "colour"))
  })

  output$crossing_summary <- renderUI({
    crossings <- real_crossing_summary()
    tags$div(
      class = "power-table",
      tags$table(
        class = "table table-sm",
        tags$thead(tags$tr(
          tags$th("Metodo"), tags$th("Frontera"), tags$th("Resultado secuencial")
        )),
        tags$tbody(lapply(seq_len(nrow(crossings)), function(index) {
          row <- crossings[index, ]
          tags$tr(
            tags$td(row$method),
            tags$td(sprintf("%.2f", row$boundary)),
            tags$td(row$outcome)
          )
        }))
      )
    )
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

  output$binomial_summary <- renderUI({
    series <- binomial_series()
    analysis <- binomial_analysis()
    bound <- binomial_boundary()
    target_events <- sum(series$target_events)
    comparator_events <- sum(series$comparator_events)
    target_reports <- sum(series$target_reports)
    comparator_reports <- sum(series$comparator_reports)
    tagList(
      tags$div(
        class = "metric-grid",
        metric("MENB: reportes", format(target_reports, big.mark = ",")),
        metric("MNQ: reportes", format(comparator_reports, big.mark = ",")),
        metric(
          "MENB: Pyrexia",
          sprintf("%s (%.1f%%)", format(target_events, big.mark = ","), 100 * target_events / target_reports)
        ),
        metric(
          "MNQ: Pyrexia",
          sprintf("%s (%.1f%%)", format(comparator_events, big.mark = ","), 100 * comparator_events / comparator_reports)
        ),
        metric("RR ajustado", sprintf("%.2f", tail(analysis$rr_hat, 1)), "razon relativa de reporte"),
        metric("Frontera V", sprintf("%.2f", bound$critical_value), paste(bound$reps, "simulaciones nulas"))
      )
    )
  })

  output$binomial_proportion_plot <- renderPlotly({
    series <- binomial_series()
    plot_data <- rbind(
      data.frame(month_date = series$month_date, proportion = series$target_proportion, vaccine = "MENB"),
      data.frame(month_date = series$month_date, proportion = series$comparator_proportion, vaccine = "MNQ")
    )
    plot <- ggplot(plot_data, aes(month_date, proportion, color = vaccine)) +
      geom_line(linewidth = 0.75) +
      scale_color_manual(values = c(MENB = "#0e7490", MNQ = "#b45309")) +
      scale_y_continuous(labels = scales::label_percent(accuracy = 1)) +
      labs(x = NULL, y = "Pyrexia / reportes", color = NULL) +
      theme_minimal(base_size = 12) + theme(legend.position = "top")
    ggplotly(plot, tooltip = c("x", "y", "colour"))
  })

  output$binomial_llr_plot <- renderPlotly({
    series <- binomial_series()
    analysis <- binomial_analysis()
    bound <- binomial_boundary()
    crossing <- first_crossing(analysis$llr, bound$critical_value)
    plot_data <- data.frame(
      month_date = series$month_date,
      llr = truncate_after_decision(analysis$llr, bound$critical_value)
    )
    plot <- ggplot(plot_data, aes(month_date, llr)) +
      geom_hline(yintercept = bound$critical_value, linetype = "dashed", color = "#b91c1c") +
      geom_line(color = "#0e7490", linewidth = 0.85) +
      annotate(
        "text", x = max(series$month_date), y = bound$critical_value,
        label = sprintf("V = %.2f", bound$critical_value),
        hjust = 1, vjust = -0.4, color = "#334155", size = 3.2
      ) +
      labs(x = NULL, y = "LLR acumulado") +
      theme_minimal(base_size = 12)
    if (!is.na(crossing$look)) {
      plot <- plot + geom_point(
        data = plot_data[crossing$look, , drop = FALSE],
        color = "#b91c1c", size = 3
      )
    }
    ggplotly(plot, tooltip = c("x", "y"))
  })

  output$binomial_decision <- renderUI({
    series <- binomial_series()
    analysis <- binomial_analysis()
    bound <- binomial_boundary()
    crossing <- first_crossing(analysis$llr, bound$critical_value)
    outcome <- if (crossing$decision == "reject") {
      sprintf("Rechaza H0 por primera vez en %s.", format(series$month_date[crossing$look], "%Y-%m"))
    } else {
      sprintf("Finaliza sin rechazo al alcanzar N = %s eventos de Pyrexia.", format(sum(series$target_events + series$comparator_events), big.mark = ","))
    }
    tags$p(
      class = "diagnostic-warning",
      outcome,
      " El alfa es condicional a los margenes mensuales observados; la conclusion se refiere a proporciones de reportes VAERS."
    )
  })

  output$binomial_monthly_table <- renderDT({
    series <- binomial_series()
    analysis <- binomial_analysis()
    display <- data.frame(
      Mes = format(series$month_date, "%Y-%m"),
      `Reportes MENB` = series$target_reports,
      `Pyrexia MENB` = series$target_events,
      `% MENB` = sprintf("%.1f%%", 100 * series$target_proportion),
      `Reportes MNQ` = series$comparator_reports,
      `Pyrexia MNQ` = series$comparator_events,
      `% MNQ` = sprintf("%.1f%%", 100 * series$comparator_proportion),
      `RR ajustado acumulado` = round(analysis$rr_hat, 3),
      `LLR acumulado` = round(analysis$llr, 3),
      check.names = FALSE
    )
    datatable(display, rownames = FALSE, options = list(pageLength = 10, dom = "tip"))
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
          "Perfil temporal de la serie reescalada a T = %.0f; V mensual = %.2f.",
          attr(result, "simulation_t"), attr(result, "critical_value")
        )
      ),
      tags$div(
        class = "power-table",
        tags$table(
          class = "table table-sm",
          tags$thead(tags$tr(tags$th("Metodo"), tags$th("Prob. rechazo"), tags$th("Aceptacion temprana"), tags$th("Tiempo medio hasta el rechazo (meses)"))),
          tags$tbody(lapply(seq_len(nrow(result)), function(index) {
            row <- result[index, ]
            tags$tr(
              tags$td(row$method),
              tags$td(sprintf("%.1f%%", 100 * row$rejection_rate)),
              tags$td(sprintf("%.1f%%", 100 * row$early_acceptance_rate)),
              tags$td(ifelse(is.nan(row$mean_rejection_look), "-", round(row$mean_rejection_look, 1)))
            )
          }))
        )
      )
    )
  })

  output$discretization_workload <- renderUI({
    looks <- as.integer(or_else(input$disc_looks, character()))
    requested_reps <- as.integer(or_else(input$disc_reps, "5000"))
    reps <- if (length(looks) > 0 && all(is.finite(looks))) {
      discretization_effective_reps(looks, requested_reps)
    } else {
      requested_reps
    }
    workload <- reps * sum(looks, na.rm = TRUE) * 2
    messages <- list()
    if (reps < requested_reps) {
      messages <- c(messages, list(tags$p(
        class = "diagnostic-warning",
        if (max(looks) >= 25000) {
          "Con K >= 25000, las replicas se limitan automaticamente a 5,000 para mantener un tiempo interactivo razonable."
        } else {
          "Con K >= 10000, las replicas se limitan automaticamente a 100,000 para mantener un tiempo interactivo razonable."
        }
      )))
    }
    if (is.finite(workload) && workload > 5e6) {
      messages <- c(messages, list(tags$p(
          class = "diagnostic-warning",
          sprintf(
            "Calculo grande: aproximadamente %s actualizaciones Poisson entre alpha y potencia.",
            format(workload, big.mark = ",", scientific = FALSE)
          )
      )))
    }
    tagList(messages)
  })

  discretization_experiment <- eventReactive(input$run_discretization, {
    looks <- as.integer(or_else(input$disc_looks, character()))
    validate(
      need(length(looks) >= 2, "Elegi al menos dos cantidades de looks."),
      need(all(is.finite(looks)) && all(looks > 0), "Los looks deben ser positivos.")
    )
    requested_reps <- as.integer(or_else(input$disc_reps, "5000"))
    effective_reps <- discretization_effective_reps(looks, requested_reps)
    run_discretization_experiment(
      horizon = as.numeric(or_else(input$disc_t, "50")),
      alpha = as.numeric(or_else(input$disc_alpha, "0.05")),
      rr = as.numeric(or_else(input$disc_rr, 1.5)),
      looks = looks,
      reps = effective_reps,
      seed = as.integer(or_else(input$disc_seed, 20260715)),
      exact_engine = exact_maxsprt
    )
  }, ignoreInit = TRUE)

  output$discretization_metrics <- renderUI({
    validate(need(
      !is.null(discretization_experiment()),
      "Elegi los parametros y presiona Comparar discretizaciones."
    ))
    result <- discretization_experiment()
    tags$div(
      class = "metric-grid",
      metric("Frontera continua V", sprintf("%.4f", result$boundary)),
      metric("Alpha continuo", sprintf("%.3f", result$alpha)),
      metric("Potencia continua", sprintf("%.1f%%", 100 * result$continuous_power)),
      metric("Horizonte T", result$horizon, "eventos esperados bajo H0"),
      metric("RR para potencia", result$rr),
      metric("Replicas por escenario", format(result$reps, big.mark = ","))
    )
  })

  output$discretization_alpha_plot <- renderPlotly({
    result <- discretization_experiment()
    data <- result$results
    plot <- ggplot(data, aes(looks, alpha_experimental)) +
      geom_hline(
        yintercept = result$alpha, linetype = "dashed", color = "#b91c1c"
      ) +
      geom_errorbar(
        aes(ymin = alpha_lower, ymax = alpha_upper),
        width = 0.08, color = "#0e7490"
      ) +
      geom_line(color = "#0e7490", linewidth = 0.8) +
      geom_point(color = "#0e7490", size = 2.2) +
      scale_x_log10(breaks = data$looks) +
      scale_y_continuous(labels = scales::label_percent(accuracy = 0.1)) +
      labs(x = "Cantidad de looks K (escala log)", y = "Alpha experimental") +
      theme_minimal(base_size = 12)
    ggplotly(plot, tooltip = c("x", "y"))
  })

  output$discretization_power_plot <- renderPlotly({
    result <- discretization_experiment()
    data <- result$results
    plot <- ggplot(data, aes(looks, power_experimental)) +
      geom_hline(
        yintercept = result$continuous_power,
        linetype = "dashed", color = "#6d28d9"
      ) +
      geom_errorbar(
        aes(ymin = power_lower, ymax = power_upper),
        width = 0.08, color = "#b45309"
      ) +
      geom_line(color = "#b45309", linewidth = 0.8) +
      geom_point(color = "#b45309", size = 2.2) +
      scale_x_log10(breaks = data$looks) +
      scale_y_continuous(labels = scales::label_percent(accuracy = 1)) +
      labs(x = "Cantidad de looks K (escala log)", y = "Potencia experimental") +
      theme_minimal(base_size = 12)
    ggplotly(plot, tooltip = c("x", "y"))
  })

  output$discretization_table <- renderDT({
    result <- discretization_experiment()
    table <- result$results
    display <- data.frame(
      K = table$looks,
      `Alpha experimental` = sprintf("%.3f%%", 100 * table$alpha_experimental),
      `IC 95% alpha` = sprintf(
        "[%.3f%%, %.3f%%]", 100 * table$alpha_lower, 100 * table$alpha_upper
      ),
      `Diferencia alpha` = sprintf("%+.3f pp", 100 * table$alpha_difference),
      `Potencia experimental` = sprintf("%.2f%%", 100 * table$power_experimental),
      `IC 95% potencia` = sprintf(
        "[%.2f%%, %.2f%%]", 100 * table$power_lower, 100 * table$power_upper
      ),
      `Diferencia potencia` = sprintf("%+.2f pp", 100 * table$power_difference),
      check.names = FALSE
    )
    datatable(display, rownames = FALSE, options = list(dom = "t", pageLength = 10))
  })

}

shinyApp(ui, app_server)
