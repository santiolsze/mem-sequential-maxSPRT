poisson_maxsprt_llr <- function(cumulative_observed, cumulative_expected) {
  ifelse(
    cumulative_observed > cumulative_expected & cumulative_expected > 0,
    cumulative_observed * log(cumulative_observed / cumulative_expected) +
      cumulative_expected - cumulative_observed,
    0
  )
}

classical_poisson_sprt <- function(
  cumulative_observed,
  cumulative_expected,
  rr,
  alpha,
  beta
) {
  stopifnot(rr > 1, alpha > 0, alpha < 1, beta > 0, beta < 1)

  upper <- log((1 - beta) / alpha)
  lower <- log(beta / (1 - alpha))
  llr <- cumulative_observed * log(rr) + cumulative_expected * (1 - rr)
  decision <- ifelse(llr >= upper, "reject", ifelse(llr <= lower, "accept", "continue"))

  data.frame(
    llr = llr,
    upper = rep(upper, length(llr)),
    lower = rep(lower, length(llr)),
    decision = decision
  )
}

with_seed <- function(seed, code) {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) {
    previous_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  }
  on.exit({
    if (had_seed) {
      assign(".Random.seed", previous_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed)
  force(code)
}

first_crossing <- function(values, upper, lower = NULL) {
  if (is.null(lower)) {
    index <- which(values >= upper)[1]
    return(if (is.na(index)) list(decision = "no_signal", look = NA_integer_) else list(decision = "reject", look = index))
  }

  index <- which(values >= upper | values <= lower)[1]
  if (is.na(index)) {
    return(list(decision = "no_decision", look = NA_integer_))
  }
  list(decision = if (values[index] >= upper) "reject" else "accept", look = index)
}

sequential_decision_summary <- function(values, dates, upper, lower = NULL) {
  if (length(values) != length(dates) || length(values) == 0L) {
    stop("values and dates must have the same positive length", call. = FALSE)
  }

  crossing <- first_crossing(values, upper, lower)
  has_decision <- !is.na(crossing$look)
  rejected <- identical(crossing$decision, "reject")

  list(
    decision = crossing$decision,
    decision_date = if (has_decision) dates[crossing$look] else as.Date(NA),
    first_rejection = if (rejected) dates[crossing$look] else as.Date(NA),
    crossing_index = if (rejected) as.integer(crossing$look) else NA_integer_,
    crossing_value = if (rejected) values[crossing$look] else NA_real_
  )
}

truncate_after_decision <- function(values, upper, lower = NULL) {
  crossing <- first_crossing(values, upper, lower)
  if (!is.na(crossing$look) && crossing$look < length(values)) {
    values[seq.int(crossing$look + 1L, length(values))] <- NA_real_
  }
  values
}

calibration_reps <- function(alpha, requested_reps, expected_tail_events = 100L) {
  stopifnot(alpha > 0, alpha < 1, requested_reps >= 1, expected_tail_events >= 1)
  as.integer(max(requested_reps, ceiling(expected_tail_events / alpha)))
}

conservative_critical_value <- function(null_maxima, alpha) {
  candidates <- sort(unique(null_maxima))
  valid <- candidates[vapply(
    candidates,
    function(candidate) mean(null_maxima >= candidate) <= alpha,
    logical(1)
  )]
  if (length(valid) == 0L) {
    return(Inf)
  }
  valid[1]
}

calibrate_monthly_poisson_boundary <- function(expected, alpha, reps = 5000L, seed = 20260711L) {
  stopifnot(all(expected > 0), alpha > 0, alpha < 1, reps >= 100)

  null_maxima <- with_seed(seed, {
    looks <- length(expected)
    counts <- matrix(
      rpois(reps * looks, lambda = rep(expected, each = reps)),
      nrow = reps,
      ncol = looks
    )
    cumulative_counts <- t(apply(counts, 1, cumsum))
    cumulative_expected <- matrix(
      rep(cumsum(expected), each = reps),
      nrow = reps,
      ncol = looks
    )
    llr <- ifelse(
      cumulative_counts > cumulative_expected,
      cumulative_counts * log(cumulative_counts / cumulative_expected) +
        cumulative_expected - cumulative_counts,
      0
    )
    apply(llr, 1, max)
  })

  list(
    critical_value = conservative_critical_value(null_maxima, alpha),
    null_maxima = null_maxima,
    reps = reps,
    alpha = alpha
  )
}

simulate_poisson_comparison <- function(
  expected,
  true_rr,
  alpha,
  beta,
  rr_low,
  rr_high,
  critical_value,
  reps = 2000L,
  seed = 20260711L
) {
  stopifnot(all(expected > 0), true_rr > 0, rr_low > 1, rr_high > rr_low)

  with_seed(seed, {
    simulations <- lapply(seq_len(reps), function(index) {
      observed <- rpois(length(expected), expected * true_rr)
      cumulative_observed <- cumsum(observed)
      cumulative_expected <- cumsum(expected)
      maxsprt <- poisson_maxsprt_llr(cumulative_observed, cumulative_expected)
      low <- classical_poisson_sprt(cumulative_observed, cumulative_expected, rr_low, alpha, beta)
      high <- classical_poisson_sprt(cumulative_observed, cumulative_expected, rr_high, alpha, beta)
      max_result <- first_crossing(maxsprt, critical_value)
      low_result <- first_crossing(low$llr, low$upper[1], low$lower[1])
      high_result <- first_crossing(high$llr, high$upper[1], high$lower[1])

      data.frame(
        simulation = index,
        maxsprt_decision = max_result$decision,
        maxsprt_look = max_result$look,
        low_decision = low_result$decision,
        low_look = low_result$look,
        high_decision = high_result$decision,
        high_look = high_result$look
      )
    })
    do.call(rbind, simulations)
  })
}

summarise_simulation <- function(simulations, method) {
  decision_col <- paste0(method, "_decision")
  look_col <- paste0(method, "_look")
  decisions <- simulations[[decision_col]]
  looks <- simulations[[look_col]]
  data.frame(
    method = method,
    rejection_rate = mean(decisions == "reject"),
    early_acceptance_rate = mean(decisions == "accept"),
    mean_rejection_look = mean(looks[decisions == "reject"], na.rm = TRUE)
  )
}

binomial_maxsprt_llr <- function(exposed_cases, total_cases, matching_ratio = 1) {
  unexposed_cases <- total_cases - exposed_cases
  excess <- unexposed_cases > 0 &
    matching_ratio * exposed_cases / unexposed_cases > 1

  result <- rep(0, length(total_cases))
  result[excess] <- with(
    list(c = exposed_cases[excess], n = total_cases[excess], u = unexposed_cases[excess]),
    c * log(c / n) + u * log(u / n) -
      c * log(1 / (matching_ratio + 1)) -
      u * log(matching_ratio / (matching_ratio + 1))
  )
  result
}

stratified_binomial_profile <- function(total_events, p0, target_event_totals) {
  cumulative_total <- sum(total_events)
  cumulative_target <- as.integer(target_event_totals)
  expected_target <- sum(total_events * p0)
  null_log_probability <- sum(total_events * log(p0))

  llr <- numeric(length(cumulative_target))
  rr_hat <- rep(1, length(cumulative_target))
  above_null <- cumulative_target > expected_target

  for (index in which(above_null)) {
    target_total <- cumulative_target[index]
    if (target_total == cumulative_total) {
      llr[index] <- -null_log_probability
      rr_hat[index] <- Inf
      next
    }

    score <- function(theta) {
      target_total - sum(total_events * plogis(qlogis(p0) + theta))
    }
    upper <- 1
    while (score(upper) > 0 && upper < 32) {
      upper <- upper * 2
    }
    theta_hat <- uniroot(score, interval = c(0, upper), tol = 1e-10)$root
    log_normalizer <- log1p(p0 * expm1(theta_hat))
    llr[index] <- theta_hat * target_total - sum(total_events * log_normalizer)
    rr_hat[index] <- exp(theta_hat)
  }

  list(llr = llr, rr_hat = rr_hat)
}

stratified_binomial_maxsprt <- function(
  target_events,
  comparator_events,
  target_reports,
  comparator_reports
) {
  lengths <- lengths(list(
    target_events, comparator_events, target_reports, comparator_reports
  ))
  if (length(unique(lengths)) != 1L || lengths[1] == 0L) {
    stop("all inputs must have the same positive length", call. = FALSE)
  }
  if (any(!is.finite(c(
    target_events, comparator_events, target_reports, comparator_reports
  ))) || any(target_events < 0) || any(comparator_events < 0) ||
      any(target_reports <= 0) || any(comparator_reports <= 0) ||
      any(target_events > target_reports) ||
      any(comparator_events > comparator_reports)) {
    stop("counts and report margins must be valid", call. = FALSE)
  }

  total_events <- target_events + comparator_events
  p0 <- target_reports / (target_reports + comparator_reports)
  llr <- numeric(length(total_events))
  rr_hat <- rep(1, length(total_events))
  cumulative_target <- cumsum(target_events)

  for (look in seq_along(total_events)) {
    profile <- stratified_binomial_profile(
      total_events[seq_len(look)], p0[seq_len(look)], cumulative_target[look]
    )
    llr[look] <- profile$llr
    rr_hat[look] <- profile$rr_hat
  }

  list(llr = llr, rr_hat = rr_hat, p0 = p0)
}

calibrate_stratified_binomial_boundary <- function(
  total_events,
  target_reports,
  comparator_reports,
  alpha,
  reps = 5000L,
  seed = 20260715L
) {
  if (length(total_events) == 0L ||
      length(target_reports) != length(total_events) ||
      length(comparator_reports) != length(total_events) ||
      any(total_events < 0) || any(target_reports <= 0) ||
      any(comparator_reports <= 0) || alpha <= 0 || alpha >= 1 || reps < 100) {
    stop("invalid conditional binomial calibration inputs", call. = FALSE)
  }

  p0 <- target_reports / (target_reports + comparator_reports)
  null_maxima <- with_seed(seed, {
    maxima <- numeric(reps)
    cumulative_target <- integer(reps)
    for (look in seq_along(total_events)) {
      cumulative_target <- cumulative_target + rbinom(
        reps, size = total_events[look], prob = p0[look]
      )
      possible_targets <- 0:sum(total_events[seq_len(look)])
      profile <- stratified_binomial_profile(
        total_events[seq_len(look)], p0[seq_len(look)], possible_targets
      )
      maxima <- pmax(maxima, profile$llr[cumulative_target + 1L])
    }
    maxima
  })

  list(
    critical_value = conservative_critical_value(null_maxima, alpha),
    null_maxima = null_maxima,
    reps = reps,
    alpha = alpha
  )
}

simulate_binomial_trajectory <- function(total_cases, matching_ratio, true_rr, seed = 20260711L) {
  stopifnot(total_cases > 1, matching_ratio %in% c(1, 2, 3), true_rr > 0)
  exposed_probability <- true_rr / (true_rr + matching_ratio)
  exposed <- with_seed(seed, rbinom(total_cases, size = 1, prob = exposed_probability))
  cumulative_exposed <- cumsum(exposed)
  cases <- seq_len(total_cases)
  data.frame(
    case = cases,
    exposed_cases = cumulative_exposed,
    total_cases = cases,
    llr = binomial_maxsprt_llr(cumulative_exposed, cases, matching_ratio)
  )
}

binomial_critical_value <- function(total_cases, matching_ratio) {
  table <- data.frame(
    total_cases = c(20, 50, 100, 200, 500, 1000),
    `1` = c(3.09884, 3.46574, 3.46574, 3.68065, 3.95630, 4.12966),
    `2` = c(3.29584, 3.29584, 3.43691, 3.59470, 3.87392, 4.09409),
    `3` = c(2.77259, 3.31895, 3.45219, 3.64449, 3.92779, 4.11634),
    check.names = FALSE
  )
  table[table$total_cases == total_cases, as.character(matching_ratio)][[1]]
}
