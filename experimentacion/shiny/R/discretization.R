maxsprt_look_grid <- function(horizon, looks) {
  if (length(horizon) != 1L || !is.finite(horizon) || horizon <= 0) {
    stop("horizon must be one positive finite number", call. = FALSE)
  }
  if (length(looks) != 1L || !is.finite(looks) ||
      looks < 1 || looks != floor(looks)) {
    stop("looks must be one positive integer", call. = FALSE)
  }

  seq_len(as.integer(looks)) * horizon / looks
}

wilson_interval <- function(successes, trials, level = 0.95) {
  if (length(successes) != 1L || length(trials) != 1L ||
      !is.finite(successes) || !is.finite(trials) || trials <= 0 ||
      successes < 0 || successes > trials ||
      successes != floor(successes) || trials != floor(trials)) {
    stop("successes and trials must define valid binomial counts", call. = FALSE)
  }
  if (length(level) != 1L || !is.finite(level) || level <= 0 || level >= 1) {
    stop("level must be between zero and one", call. = FALSE)
  }

  proportion <- successes / trials
  z <- qnorm(1 - (1 - level) / 2)
  denominator <- 1 + z^2 / trials
  center <- (proportion + z^2 / (2 * trials)) / denominator
  half_width <- z / denominator * sqrt(
    proportion * (1 - proportion) / trials + z^2 / (4 * trials^2)
  )
  c(lower = max(0, center - half_width), upper = min(1, center + half_width))
}

discretization_effective_reps <- function(looks, requested_reps) {
  if (length(looks) == 0L || any(!is.finite(looks)) || any(looks < 1)) {
    stop("looks must contain positive finite values", call. = FALSE)
  }
  if (length(requested_reps) != 1L || !is.finite(requested_reps) ||
      requested_reps < 1 || requested_reps != floor(requested_reps)) {
    stop("requested_reps must be one positive integer", call. = FALSE)
  }

  as.integer(if (max(looks) >= 10000) min(requested_reps, 100000) else requested_reps)
}

with_discretization_seed <- function(seed, code) {
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

simulate_discrete_crossing <- function(
  horizon,
  looks,
  rr,
  boundary,
  reps,
  seed
) {
  expected_grid <- maxsprt_look_grid(horizon, looks)
  if (length(rr) != 1L || !is.finite(rr) || rr <= 0) {
    stop("rr must be one positive finite number", call. = FALSE)
  }
  if (length(boundary) != 1L || !is.finite(boundary) || boundary <= 0) {
    stop("boundary must be one positive finite number", call. = FALSE)
  }
  if (length(reps) != 1L || !is.finite(reps) ||
      reps < 1 || reps != floor(reps)) {
    stop("reps must be one positive integer", call. = FALSE)
  }
  if (length(seed) != 1L || !is.finite(seed) || seed != floor(seed)) {
    stop("seed must be one finite integer", call. = FALSE)
  }

  signals <- with_discretization_seed(seed, {
    cumulative_observed <- integer(reps)
    maximum_llr <- numeric(reps)
    increment_mean <- rr * horizon / looks
    for (look_index in seq_len(looks)) {
      cumulative_observed <- cumulative_observed + rpois(reps, increment_mean)
      expected <- expected_grid[look_index]
      llr <- ifelse(
        cumulative_observed > expected,
        cumulative_observed * log(cumulative_observed / expected) +
          expected - cumulative_observed,
        0
      )
      maximum_llr <- pmax(maximum_llr, llr)
    }
    sum(maximum_llr >= boundary)
  })

  rate <- signals / reps
  interval <- wilson_interval(signals, reps)
  list(
    rate = rate,
    standard_error = sqrt(rate * (1 - rate) / reps),
    lower = unname(interval["lower"]),
    upper = unname(interval["upper"]),
    signals = signals,
    reps = reps
  )
}

run_discretization_experiment <- function(
  horizon,
  alpha,
  rr,
  looks,
  reps,
  seed,
  exact_engine
) {
  if (length(horizon) != 1L || !is.finite(horizon) || horizon <= 0) {
    stop("horizon must be one positive finite number", call. = FALSE)
  }
  if (length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("alpha must be between zero and one", call. = FALSE)
  }
  if (length(rr) != 1L || !is.finite(rr) || rr <= 1) {
    stop("rr must be greater than one for the power comparison", call. = FALSE)
  }
  if (length(looks) < 2L) {
    stop("looks must contain at least two values", call. = FALSE)
  }
  if (any(!is.finite(looks)) || any(looks < 1) || any(looks != floor(looks))) {
    stop("looks must contain positive integers", call. = FALSE)
  }
  if (anyDuplicated(looks)) {
    stop("looks values must be distinct", call. = FALSE)
  }
  if (length(reps) != 1L || !is.finite(reps) ||
      reps < 1 || reps != floor(reps)) {
    stop("reps must be one positive integer", call. = FALSE)
  }
  if (length(seed) != 1L || !is.finite(seed) || seed != floor(seed) ||
      seed < 0 || seed > .Machine$integer.max - 2L * length(looks) - 1L) {
    stop("seed must be a nonnegative integer with room for stream offsets",
         call. = FALSE)
  }
  if (is.null(exact_engine) ||
      !is.function(exact_engine$calibrate_maxsprt_boundary) ||
      !is.function(exact_engine$maxsprt_exact)) {
    stop("exact_engine must provide calibration and exact power functions",
         call. = FALSE)
  }

  looks <- sort(as.integer(looks))
  calibration <- exact_engine$calibrate_maxsprt_boundary(horizon, alpha)
  boundary <- calibration$critical_value
  continuous_power <- exact_engine$maxsprt_exact(
    boundary, horizon, rr
  )$signal_probability

  rows <- lapply(seq_along(looks), function(index) {
    null_result <- simulate_discrete_crossing(
      horizon = horizon,
      looks = looks[index],
      rr = 1,
      boundary = boundary,
      reps = reps,
      seed = seed + 2L * index
    )
    power_result <- simulate_discrete_crossing(
      horizon = horizon,
      looks = looks[index],
      rr = rr,
      boundary = boundary,
      reps = reps,
      seed = seed + 2L * index + 1L
    )

    data.frame(
      looks = looks[index],
      alpha_experimental = null_result$rate,
      alpha_se = null_result$standard_error,
      alpha_lower = null_result$lower,
      alpha_upper = null_result$upper,
      alpha_difference = null_result$rate - alpha,
      power_experimental = power_result$rate,
      power_se = power_result$standard_error,
      power_lower = power_result$lower,
      power_upper = power_result$upper,
      power_difference = power_result$rate - continuous_power
    )
  })

  list(
    boundary = boundary,
    continuous_power = continuous_power,
    alpha = alpha,
    rr = rr,
    horizon = horizon,
    reps = as.integer(reps),
    results = do.call(rbind, rows)
  )
}
