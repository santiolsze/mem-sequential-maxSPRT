#!/usr/bin/env Rscript

# Exact replication of binomial MaxSPRT alpha calibration in Kulldorff et al.
# Base R only. The default walkthrough reproduces Table 4 for z=1, N=10,
# nominal alpha=0.05.

binomial_maxsprt_llr_exact <- function(c, n, z = 1) {
  if (length(c) != length(n) && length(c) != 1L && length(n) != 1L) {
    stop("c and n must have compatible lengths", call. = FALSE)
  }
  if (any(!is.finite(c)) || any(!is.finite(n)) || !is.finite(z) ||
      z <= 0 || any(n < 0) || any(c < 0) || any(c > n)) {
    stop("invalid binomial counts or matching ratio", call. = FALSE)
  }

  p0 <- 1 / (z + 1)
  phat <- ifelse(n == 0, 0, c / n)
  result <- numeric(max(length(c), length(n)))
  excess <- n > 0 & phat > p0

  x_log_ratio <- function(x, denominator) {
    ifelse(x == 0, 0, x * log(x / denominator))
  }
  result[excess] <- n[excess] * (
    x_log_ratio(phat[excess], p0) +
      x_log_ratio(1 - phat[excess], 1 - p0)
  )
  result
}

binomial_alpha_markov <- function(V, N, z = 1, keep_states = TRUE) {
  if (!is.finite(V) || V < 0 || length(N) != 1L || N < 1 || N != as.integer(N) ||
      !is.finite(z) || z <= 0) {
    stop("V, N, and z must define a valid calibration", call. = FALSE)
  }

  N <- as.integer(N)
  p0 <- 1 / (z + 1)
  alive <- 1
  absorbed <- 0
  history <- vector("list", N)
  state_history <- vector("list", N)

  for (n in seq_len(N)) {
    incoming <- numeric(n + 1L)
    incoming[seq_len(n)] <- incoming[seq_len(n)] + alive * (1 - p0)
    incoming[seq_len(n) + 1L] <- incoming[seq_len(n) + 1L] + alive * p0

    c_values <- 0:n
    llr <- binomial_maxsprt_llr_exact(c_values, rep(n, n + 1L), z)
    reject <- llr >= V & incoming > 0
    newly_absorbed <- sum(incoming[reject])
    absorbed <- absorbed + newly_absorbed
    alive <- incoming
    alive[reject] <- 0

    history[[n]] <- data.frame(
      n = n,
      newly_absorbed = newly_absorbed,
      absorbed_probability = absorbed,
      alive_probability = sum(alive)
    )
    if (keep_states) {
      state_history[[n]] <- data.frame(
        n = n,
        c = c_values,
        llr = llr,
        incoming_probability = incoming,
        rejects = reject,
        alive_probability = alive
      )
    }
  }

  list(
    alpha = absorbed,
    history = do.call(rbind, history),
    states = if (keep_states) do.call(rbind, state_history) else NULL,
    p0 = p0,
    V = V,
    N = N,
    z = z
  )
}

binomial_alpha_enumeration <- function(V, N, z = 1) {
  if (N > 20) {
    stop("direct enumeration is limited to N <= 20", call. = FALSE)
  }
  p0 <- 1 / (z + 1)
  paths <- as.matrix(expand.grid(rep(list(0:1), N)))
  cumulative_exposed <- t(apply(paths, 1, cumsum))
  crossed <- logical(nrow(paths))

  for (n in seq_len(N)) {
    crossed <- crossed | binomial_maxsprt_llr_exact(
      cumulative_exposed[, n], rep(n, nrow(paths)), z
    ) >= V
  }
  path_probability <- apply(
    paths,
    1,
    function(path) prod(ifelse(path == 1, p0, 1 - p0))
  )
  sum(path_probability[crossed])
}

calibrate_binomial_boundary <- function(alpha, N, z = 1, digits = 5L) {
  if (!is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("alpha must lie strictly between zero and one", call. = FALSE)
  }

  reachable <- unique(unlist(lapply(seq_len(N), function(n) {
    binomial_maxsprt_llr_exact(0:n, rep(n, n + 1L), z)
  })))
  reachable <- sort(reachable[reachable > 0])

  candidates <- do.call(rbind, lapply(reachable, function(boundary) {
    epsilon <- 8 * .Machine$double.eps * max(1, abs(boundary))
    internal_V <- boundary + epsilon
    data.frame(
      excluded_llr = boundary,
      internal_V = internal_V,
      actual_alpha = binomial_alpha_markov(
        internal_V, N, z, keep_states = FALSE
      )$alpha,
      alpha_if_equality_included = binomial_alpha_markov(
        boundary, N, z, keep_states = FALSE
      )$alpha
    )
  }))

  safe <- which(candidates$actual_alpha <= alpha)
  if (length(safe) == 0L) {
    stop("no conservative binomial boundary exists", call. = FALSE)
  }
  selected <- candidates[safe[1], , drop = FALSE]
  scale <- 10^as.integer(digits)
  display_value <- ceiling(selected$excluded_llr * scale) / scale

  list(
    critical_value = selected$internal_V,
    critical_value_display = display_value,
    excluded_llr = selected$excluded_llr,
    actual_alpha = selected$actual_alpha,
    alpha_if_equality_included = selected$alpha_if_equality_included,
    nominal_alpha = alpha,
    N = N,
    z = z,
    p0 = 1 / (z + 1),
    candidates = candidates
  )
}

run_binomial_alpha_walkthrough <- function(
  alpha = 0.05,
  N = 10L,
  z = 1,
  output_dir = file.path("experimentacion", "results", "binomial_alpha")
) {
  calibration <- calibrate_binomial_boundary(alpha, N, z)
  markov <- binomial_alpha_markov(
    calibration$critical_value_display, N, z, keep_states = TRUE
  )
  direct_alpha <- binomial_alpha_enumeration(
    calibration$critical_value_display, N, z
  )

  cat("REPLICA EXACTA: ALPHA DEL MAXSPRT BINOMIAL\n")
  cat("==========================================\n\n")
  cat(sprintf("1. Diseno: z = %g, N = %d, alpha nominal = %.3f\n", z, N, alpha))
  cat(sprintf("   Bajo H0: p0 = 1/(z+1) = %.6f\n\n", calibration$p0))
  cat("2. Estado (n,c): n eventos totales y c en el lado expuesto.\n")
  cat("   El siguiente evento lleva a (n+1,c+1) con probabilidad p0\n")
  cat("   o a (n+1,c) con probabilidad 1-p0.\n\n")
  cat("3. En cada estado se calcula el LLR binomial unilateral.\n")
  cat("   Si LLR >= V, el estado absorbe su probabilidad: ya hubo senal.\n\n")
  cat(sprintf(
    "4. 4*log(2) = %.9f daria alpha = %.6f al incluir la igualdad.\n",
    calibration$excluded_llr,
    calibration$alpha_if_equality_included
  ))
  cat(sprintf(
    "   La tabla publica V = %.5f, apenas mayor, y el alpha efectivo es %.6f.\n\n",
    calibration$critical_value_display,
    markov$alpha
  ))
  cat("5. Propagacion exacta por numero de evento:\n")
  print(markov$history, row.names = FALSE, digits = 8)
  cat(sprintf("\n6. Comprobacion enumerando los 2^%d caminos: %.8f\n", N, direct_alpha))
  cat(sprintf("   Diferencia Markov - enumeracion: %.3g\n", markov$alpha - direct_alpha))

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(markov$history, file.path(output_dir, "step_history.csv"), row.names = FALSE)
  write.csv(markov$states, file.path(output_dir, "state_history.csv"), row.names = FALSE)
  write.csv(calibration$candidates, file.path(output_dir, "candidate_boundaries.csv"), row.names = FALSE)

  invisible(list(calibration = calibration, markov = markov, direct_alpha = direct_alpha))
}

if (sys.nframe() == 0L) {
  run_binomial_alpha_walkthrough()
}
