# Exact numerical replication of Poisson MaxSPRT Tables 1-3 in paper.pdf.
#
# This file intentionally uses base R only. Sourcing it defines reusable
# functions; direct execution additionally builds and writes the tables.

lambert_w_scalar <- function(x, branch = 0L, tol = 1e-12) {
  branch <- as.integer(branch)
  branch_point <- -1 / exp(1)

  if (length(x) != 1L || is.na(x) || !is.finite(x)) {
    stop("x must be one finite numeric value", call. = FALSE)
  }
  if (!branch %in% c(-1L, 0L)) {
    stop("branch must be 0 or -1", call. = FALSE)
  }
  if (x < branch_point - tol) {
    stop("x is outside the real Lambert W domain", call. = FALSE)
  }
  if (abs(x - branch_point) <= tol) {
    return(-1)
  }

  equation <- function(w) w * exp(w) - x

  if (branch == 0L) {
    if (x < 0) {
      return(uniroot(equation, c(-1, 0), tol = tol)$root)
    }
    if (x == 0) {
      return(0)
    }
    upper <- max(1, log(x) + 1)
    while (equation(upper) < 0) {
      upper <- upper * 2
    }
    return(uniroot(equation, c(0, upper), tol = tol)$root)
  }

  if (x > 0) {
    stop("x is outside the real domain of Lambert W branch -1", call. = FALSE)
  }
  if (x == 0) {
    return(-Inf)
  }

  lower <- -2
  while (equation(lower) < 0) {
    lower <- lower * 2
  }
  uniroot(equation, c(lower, -1), tol = tol)$root
}

lambert_w <- function(x, branch = 0L, tol = 1e-12) {
  vapply(
    x,
    lambert_w_scalar,
    numeric(1),
    branch = branch,
    tol = tol,
    USE.NAMES = FALSE
  )
}

poisson_maxsprt_llr <- function(events, expected) {
  if (length(events) != length(expected) &&
      length(events) != 1L && length(expected) != 1L) {
    stop("events and expected must have compatible lengths", call. = FALSE)
  }
  if (any(events < 0) || any(expected < 0)) {
    stop("events and expected must be nonnegative", call. = FALSE)
  }

  events <- rep_len(events, max(length(events), length(expected)))
  expected <- rep_len(expected, length(events))
  result <- numeric(length(events))
  excess <- events > expected & expected > 0
  result[excess] <- expected[excess] - events[excess] +
    events[excess] * log(events[excess] / expected[excess])
  result[events > 0 & expected == 0] <- Inf
  result
}

maxsprt_crossing_time <- function(n, critical_value, tol = 1e-12) {
  if (any(n <= 0) || any(n != floor(n))) {
    stop("n must contain positive integers", call. = FALSE)
  }
  if (length(critical_value) != 1L || critical_value <= 0) {
    stop("critical_value must be one positive number", call. = FALSE)
  }

  argument <- -exp(-1 - critical_value / n)
  -n * lambert_w(argument, branch = 0L, tol = tol)
}

advance_nonabsorbed_states <- function(state, delta, threshold, rr) {
  increment_probability <- dpois(0:(threshold - 1L), lambda = rr * delta)
  propagated <- convolve(
    state,
    rev(increment_probability),
    type = "open"
  )[seq_len(threshold)]
  propagated[abs(propagated) < .Machine$double.eps * 100] <- 0
  pmax(propagated, 0)
}

signal_time_contribution <- function(state, start_time, delta, threshold, rr) {
  previous_counts <- seq_along(state) - 1L
  arrivals_needed <- threshold - previous_counts
  mean_increment <- rr * delta
  signal_probability <- ppois(
    arrivals_needed - 1L,
    lambda = mean_increment,
    lower.tail = FALSE
  )
  truncated_arrival_time <- arrivals_needed / rr * ppois(
    arrivals_needed,
    lambda = mean_increment,
    lower.tail = FALSE
  )
  sum(state * (start_time * signal_probability + truncated_arrival_time))
}

maxsprt_exact <- function(critical_value, horizon, rr = 1, tol = 1e-12) {
  if (length(horizon) != 1L || !is.finite(horizon) || horizon <= 0) {
    stop("horizon must be one positive finite number", call. = FALSE)
  }
  if (length(rr) != 1L || !is.finite(rr) || rr <= 0) {
    stop("rr must be one positive finite number", call. = FALSE)
  }

  state <- 1
  previous_time <- 0
  signal_probability <- 0
  signal_time_numerator <- 0
  crossing_rows <- list()
  n <- 1L

  repeat {
    theoretical_crossing <- maxsprt_crossing_time(n, critical_value, tol)
    boundary_time <- min(theoretical_crossing, horizon)
    delta <- boundary_time - previous_time
    if (delta < -tol) {
      stop("crossing times are not ordered", call. = FALSE)
    }
    delta <- max(delta, 0)

    mass_before <- sum(state)
    time_contribution <- signal_time_contribution(
      state, previous_time, delta, n, rr
    )
    next_state <- advance_nonabsorbed_states(state, delta, n, rr)
    absorbed <- mass_before - sum(next_state)
    if (absorbed < -tol) {
      stop("probability mass increased during absorption", call. = FALSE)
    }
    absorbed <- max(absorbed, 0)

    signal_probability <- signal_probability + absorbed
    signal_time_numerator <- signal_time_numerator + time_contribution
    crossing_rows[[length(crossing_rows) + 1L]] <- data.frame(
      event = n,
      time = boundary_time,
      theoretical_time = theoretical_crossing,
      signal_mass = absorbed
    )
    state <- next_state
    previous_time <- boundary_time

    if (theoretical_crossing >= horizon - tol) {
      break
    }
    n <- n + 1L
  }

  no_signal_probability <- sum(state)
  total_probability <- signal_probability + no_signal_probability
  if (abs(total_probability - 1) > max(tol * 10, 1e-10)) {
    stop(sprintf("probability conservation failed: %.16g", total_probability),
         call. = FALSE)
  }

  list(
    signal_probability = signal_probability,
    signal_time_numerator = signal_time_numerator,
    conditional_signal_time = signal_time_numerator / signal_probability,
    expected_surveillance_time = signal_time_numerator +
      horizon * no_signal_probability,
    no_signal_probability = no_signal_probability,
    crossings = do.call(rbind, crossing_rows)
  )
}

calibrate_maxsprt_boundary <- function(horizon, alpha, tol = 1e-10) {
  if (length(alpha) != 1L || alpha <= 0 || alpha >= 1) {
    stop("alpha must be between zero and one", call. = FALSE)
  }

  objective <- function(critical_value) {
    maxsprt_exact(
      critical_value,
      horizon = horizon,
      rr = 1,
      tol = min(tol, 1e-12)
    )$signal_probability - alpha
  }
  lower <- 0.25
  upper <- 10
  while (objective(lower) < 0) {
    lower <- lower / 2
  }
  while (objective(upper) > 0) {
    upper <- upper * 1.5
  }

  critical_value <- uniroot(
    objective,
    interval = c(lower, upper),
    tol = tol
  )$root
  achieved_alpha <- maxsprt_exact(
    critical_value,
    horizon = horizon,
    rr = 1,
    tol = min(tol, 1e-12)
  )$signal_probability
  list(
    critical_value = critical_value,
    achieved_alpha = achieved_alpha
  )
}

paper_table_1 <- data.frame(
  horizon = c(
    0.1, 0.2, 0.5, 1, 1.5, 2, 2.5, 3, 4, 5, 6, 8, 10, 12, 15,
    20, 25, 30, 40, 50, 60, 80, 100, 120, 150, 200, 250, 300, 400,
    500, 600, 800, 1000
  ),
  alpha_0.05 = c(
    2.044069, 2.266893, 2.637928, 2.853937, 2.964971, 3.046977,
    3.110419, 3.162106, 3.245004, 3.297183, 3.342729, 3.413782,
    3.467952, 3.511749, 3.562591, 3.628123, 3.676320, 3.715764,
    3.774663, 3.819903, 3.855755, 3.910853, 3.952321, 3.985577,
    4.025338, 4.074828, 4.112234, 4.142134, 4.188031, 4.222632,
    4.250310, 4.292829, 4.324917
  ),
  alpha_0.01 = c(
    4.119293, 4.179630, 4.483740, 4.670428, 4.778944, 4.862223,
    4.924475, 4.971792, 5.040311, 5.091907, 5.136461, 5.206326,
    5.260513, 5.302914, 5.351279, 5.414770, 5.463382, 5.502563,
    5.561620, 5.605972, 5.642209, 5.697631, 5.738974, 5.772435,
    5.812121, 5.862113, 5.899824, 5.929897, 5.976241, 6.011088,
    6.039013, 6.081871, 6.114225
  ),
  alpha_0.001 = c(
    6.579669, 6.754862, 7.034472, 7.172614, 7.278202, 7.341453,
    7.397851, 7.445736, 7.518319, 7.569312, 7.608607, 7.673013,
    7.724863, 7.767520, 7.814719, 7.877573, 7.924478, 7.962688,
    8.022182, 8.067072, 8.102340, 8.157530, 8.199403, 8.232827,
    8.272692, 8.322983, 8.360938, 8.391288, 8.438008, 8.473183,
    8.501314, 8.544590, 8.577253
  ),
  check.names = FALSE
)
attr(paper_table_1, "display_values") <- data.frame(
  horizon = format(paper_table_1$horizon, scientific = FALSE, trim = TRUE),
  alpha_0.05 = sprintf("%.6f", paper_table_1$alpha_0.05),
  alpha_0.01 = sprintf("%.6f", paper_table_1$alpha_0.01),
  alpha_0.001 = sprintf("%.6f", paper_table_1$alpha_0.001),
  check.names = FALSE
)

replicate_table_1 <- function(tol = 1e-10) {
  alpha_values <- c(0.05, 0.01, 0.001)
  result <- data.frame(horizon = paper_table_1$horizon)
  for (alpha in alpha_values) {
    column <- paste0("alpha_", alpha)
    result[[column]] <- vapply(
      paper_table_1$horizon,
      function(horizon) calibrate_maxsprt_boundary(
        horizon, alpha, tol
      )$critical_value,
      numeric(1)
    )
  }
  result
}

paper_table_2_text <- "horizon 1.2 1.5 2 3 5 10
0.1 .060 .075 .100 .148 .242 .449
0.2 .062 .081 .115 .187 .337 .648
0.5 .066 .093 .147 .273 .532 .899
1 .070 .107 .185 .379 .729 .987
1.5 .073 .118 .221 .475 .852 .9987
2 .076 .130 .255 .561 .924 .9999
2.5 .078 .140 .289 .637 .962 .999990
3 .081 .151 .323 .703 .981 .999999
4 .086 .172 .390 .809 .996 1
5 .089 .190 .447 .876 .9992 1
6 .093 .208 .500 .920 .9998 1
8 .100 .244 .600 .970 .9999943 1
10 .107 .280 .685 .989 .9999998 1
12 .114 .315 .756 .996 1 1
15 .123 .367 .836 .9993 1 1
20 .138 .450 .921 .99997 1 1
25 .153 .526 .963 .999999 1 1
30 .167 .596 .984 1 1 1
40 .196 .713 .997 1 1 1
50 .225 .803 .9996 1 1 1
60 .254 .868 .99994 1 1 1
80 .311 .944 .999999 1 1 1
100 .368 .978 1 1 1 1
120 .424 .992 1 1 1 1
150 .505 .998 1 1 1 1
200 .623 .99990 1 1 1 1
250 .722 .9999952 1 1 1 1
300 .800 .9999998 1 1 1 1
400 .903 1 1 1 1 1
500 .956 1 1 1 1 1
600 .981 1 1 1 1 1
800 .997 1 1 1 1 1
1000 .9996 1 1 1 1 1"

paper_table_3_signal_text <- "horizon 1 1.2 1.5 2 3 5 10
1 0.22 0.25 0.30 0.35 0.39 0.37 0.22
2 0.40 0.51 0.63 0.75 0.79 0.62 0.24
5 0.96 1.38 1.82 2.09 1.78 0.83 0.26
10 1.83 2.99 4.02 4.13 2.45 0.87 0.27
20 3.48 6.70 8.68 6.96 2.67 0.91 0.28
50 8.06 19.76 20.45 8.94 2.82 0.96 0.30
100 15.13 43.78 29.93 9.30 2.92 0.99 0.31
200 28.37 89.73 33.00 9.62 3.01 1.02 0.32
500 64.95 171.79 34.40 10.01 3.12 1.06 0.33
1000 121.49 196.27 35.37 10.27 3.20 1.08 0.33"

paper_table_3_surveillance_text <- "horizon 1 1.2 1.5 2 3 5 10
1 0.96 0.95 0.93 0.88 0.77 0.54 0.23
2 1.92 1.89 1.82 1.68 1.32 0.72 0.24
5 4.80 4.68 4.40 3.70 2.18 0.83 0.26
10 9.59 9.25 8.33 5.98 2.53 0.87 0.27
20 19.17 18.17 14.91 8.00 2.67 0.91 0.28
50 47.90 43.20 26.28 8.96 2.82 0.96 0.30
100 95.76 79.28 31.46 9.30 2.92 0.99 0.31
200 191.42 131.25 33.02 9.62 3.01 1.02 0.32
500 478.25 186.15 34.40 10.01 3.12 1.06 0.33
1000 956.07 196.59 35.37 10.27 3.20 1.08 0.33"

read_paper_table <- function(text) {
  character_table <- read.table(
    text = text,
    header = TRUE,
    colClasses = "character",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  numeric_table <- as.data.frame(lapply(character_table, as.numeric),
                                 check.names = FALSE)
  attr(numeric_table, "display_values") <- character_table
  numeric_table
}

paper_table_2 <- read_paper_table(paper_table_2_text)
paper_table_3_signal <- read_paper_table(paper_table_3_signal_text)
paper_table_3_surveillance <- read_paper_table(paper_table_3_surveillance_text)

replicate_table_2 <- function(boundaries) {
  rr_values <- as.numeric(names(paper_table_2)[-1])
  result <- data.frame(horizon = paper_table_2$horizon, check.names = FALSE)
  for (rr in rr_values) {
    result[[as.character(rr)]] <- mapply(
      function(horizon, boundary) {
        maxsprt_exact(boundary, horizon, rr)$signal_probability
      },
      paper_table_2$horizon,
      boundaries$alpha_0.05
    )
  }
  result
}

replicate_table_3 <- function(boundaries) {
  horizons <- paper_table_3_signal$horizon
  boundary_index <- match(horizons, boundaries$horizon)
  if (anyNA(boundary_index)) {
    stop("boundaries do not contain every Table 3 horizon", call. = FALSE)
  }
  critical_values <- boundaries$alpha_0.05[boundary_index]
  rr_values <- as.numeric(names(paper_table_3_signal)[-1])
  signal_time <- data.frame(horizon = horizons, check.names = FALSE)
  surveillance_time <- data.frame(horizon = horizons, check.names = FALSE)

  for (rr in rr_values) {
    exact <- Map(
      maxsprt_exact,
      critical_value = critical_values,
      horizon = horizons,
      MoreArgs = list(rr = rr)
    )
    signal_time[[as.character(rr)]] <- vapply(
      exact, `[[`, numeric(1), "conditional_signal_time"
    )
    surveillance_time[[as.character(rr)]] <- vapply(
      exact, `[[`, numeric(1), "expected_surveillance_time"
    )
  }
  list(signal_time = signal_time, surveillance_time = surveillance_time)
}

display_decimal_places <- function(value) {
  if (!grepl("\\.", value)) {
    return(0L)
  }
  nchar(sub("^[^.]*\\.", "", value))
}

compare_paper_table <- function(reproduced, reference, table_name) {
  if (!identical(dim(reproduced), dim(reference)) ||
      !identical(names(reproduced), names(reference))) {
    stop("reproduced and reference tables must have matching dimensions and columns",
         call. = FALSE)
  }
  display <- attr(reference, "display_values")
  if (is.null(display)) {
    display <- as.data.frame(lapply(reference, as.character),
                             check.names = FALSE)
  }

  rows <- lapply(seq_along(reference), function(column_index) {
    column <- names(reference)[column_index]
    decimals <- vapply(
      display[[column_index]], display_decimal_places, integer(1)
    )
    reproduced_value <- reproduced[[column_index]]
    published_value <- reference[[column_index]]
    data.frame(
      table = table_name,
      horizon = reference$horizon,
      column = column,
      published_display = display[[column_index]],
      published_value = published_value,
      reproduced_value = reproduced_value,
      difference = reproduced_value - published_value,
      absolute_difference = abs(reproduced_value - published_value),
      decimals = decimals,
      pass = mapply(
        function(value, published, digits) {
          abs(value - published) < 10^(-digits)
        },
        reproduced_value, published_value, decimals
      ),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

run_paper_table_replication <- function(
  output_dir = getOption(
    "maxsprt.output_dir",
    file.path("experimentacion", "results", "paper_tables")
  )
) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  message("Calibrating all Table 1 critical values...")
  table_1 <- replicate_table_1()
  message("Calculating Table 2 power...")
  table_2 <- replicate_table_2(table_1)
  message("Calculating Table 3 stopping times...")
  table_3 <- replicate_table_3(table_1)

  comparisons <- list(
    table_1 = compare_paper_table(table_1, paper_table_1, "table_1"),
    table_2 = compare_paper_table(table_2, paper_table_2, "table_2"),
    table_3_signal_time = compare_paper_table(
      table_3$signal_time, paper_table_3_signal, "table_3_signal_time"
    ),
    table_3_surveillance_time = compare_paper_table(
      table_3$surveillance_time,
      paper_table_3_surveillance,
      "table_3_surveillance_time"
    )
  )

  outputs <- list(
    table_1_critical_values = table_1,
    table_2_power = table_2,
    table_3_signal_time = table_3$signal_time,
    table_3_surveillance_time = table_3$surveillance_time
  )
  for (name in names(outputs)) {
    write.csv(
      outputs[[name]],
      file.path(output_dir, paste0(name, ".csv")),
      row.names = FALSE
    )
  }
  for (name in names(comparisons)) {
    write.csv(
      comparisons[[name]],
      file.path(output_dir, paste0(name, "_comparison.csv")),
      row.names = FALSE
    )
  }

  comparison <- do.call(rbind, comparisons)
  failures <- comparison[!comparison$pass, , drop = FALSE]
  message(sprintf(
    "Compared %d cells: %d pass, %d fail at displayed precision.",
    nrow(comparison), sum(comparison$pass), nrow(failures)
  ))

  list(
    table_1 = table_1,
    table_2 = table_2,
    table_3 = table_3,
    comparisons = comparisons,
    failures = failures,
    output_dir = normalizePath(output_dir, mustWork = TRUE)
  )
}

if (sys.nframe() == 0L) {
  replication <- run_paper_table_replication()
  if (nrow(replication$failures) > 0L) {
    print(replication$failures)
    stop("Paper-table replication has displayed-precision failures",
         call. = FALSE)
  }
}
