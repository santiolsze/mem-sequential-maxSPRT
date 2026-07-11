library(arrow)
library(dplyr)

real_series_presets <- function() {
  data.frame(
    id = c("eliquis", "hpv9_syncope", "menb_pyrexia"),
    label = c(
      "ELIQUIS / Haemorrhage",
      "HPV9 (GARDASIL 9) / Syncope",
      "MENB / Pyrexia"
    ),
    source = c("FAERS/openFDA", "VAERS", "VAERS"),
    stringsAsFactors = FALSE
  )
}

append_sequential_columns <- function(series) {
  series <- series[order(series$month_date), ]
  series$cumulative_observed <- cumsum(series$observed)
  series$cumulative_expected <- cumsum(series$expected)
  series$cumulative_rr <- series$cumulative_observed / series$cumulative_expected
  series$maxsprt_llr <- poisson_maxsprt_llr(
    series$cumulative_observed,
    series$cumulative_expected
  )
  series
}

load_openfda_series <- function(processed_root) {
  dataset <- open_dataset(
    file.path(processed_root, "fda", "openfda_eliquis_haemorrhage_monthly"),
    format = "parquet"
  )
  frame <- dataset |>
    select(
      month, all_reports, eliquis_reports, haemorrhage_reports,
      eliquis_haemorrhage_reports, non_eliquis_reports,
      non_eliquis_haemorrhage_reports,
      expected_under_no_eliquis_association
    ) |>
    collect() |>
    arrange(month)

  append_sequential_columns(data.frame(
    month_date = as.Date(paste0(frame$month, "-01")),
    observed = frame$eliquis_haemorrhage_reports,
    expected = frame$expected_under_no_eliquis_association,
    denominator = frame$eliquis_reports,
    all_reports = frame$all_reports,
    eliquis_reports = frame$eliquis_reports,
    haemorrhage_reports = frame$haemorrhage_reports,
    eliquis_haemorrhage_reports = frame$eliquis_haemorrhage_reports,
    non_eliquis_reports = frame$non_eliquis_reports,
    non_eliquis_haemorrhage_reports = frame$non_eliquis_haemorrhage_reports,
    label = "ELIQUIS / Haemorrhage",
    source = "FAERS/openFDA",
    interpretation = "Contemporaneous non-ELIQUIS report fraction comparator; report-level only, not an untreated patient comparator."
  ))
}

load_vaers_series <- function(target_vax_type, target_symptom, label, processed_root) {
  vaccines <- open_dataset(
    file.path(processed_root, "vaers", "vaccines"),
    format = "parquet"
  )
  vaccine_symptoms <- open_dataset(
    file.path(processed_root, "vaers", "vaccine_symptoms"),
    format = "parquet"
  )

  denominator_target <- vaccines |>
    filter(month > 0, vax_type == target_vax_type) |>
    distinct(source_year, vaers_id, vax_record_order, year, month) |>
    count(year, month, name = "denominator_target") |>
    collect()
  denominator_all <- vaccines |>
    filter(month > 0) |>
    distinct(source_year, vaers_id, vax_record_order, year, month) |>
    count(year, month, name = "denominator_all") |>
    collect()
  events_target <- vaccine_symptoms |>
    filter(month > 0, vax_type == target_vax_type, symptom == target_symptom) |>
    distinct(source_year, vaers_id, vax_record_order, year, month) |>
    count(year, month, name = "events_target") |>
    collect()
  events_all <- vaccine_symptoms |>
    filter(month > 0, symptom == target_symptom) |>
    distinct(source_year, vaers_id, vax_record_order, year, month) |>
    count(year, month, name = "events_all") |>
    collect()

  months <- expand.grid(year = 2016:2025, month = 1:12)
  frame <- merge(months, denominator_target, by = c("year", "month"), all.x = TRUE)
  frame <- merge(frame, denominator_all, by = c("year", "month"), all.x = TRUE)
  frame <- merge(frame, events_target, by = c("year", "month"), all.x = TRUE)
  frame <- merge(frame, events_all, by = c("year", "month"), all.x = TRUE)
  frame$denominator_target[is.na(frame$denominator_target)] <- 0
  frame$denominator_all[is.na(frame$denominator_all)] <- 0
  frame$events_target[is.na(frame$events_target)] <- 0
  frame$events_all[is.na(frame$events_all)] <- 0

  frame$non_target_denominator <- frame$denominator_all - frame$denominator_target
  frame$non_target_events <- frame$events_all - frame$events_target
  stopifnot(all(frame$non_target_denominator > 0))
  frame$expected <- frame$denominator_target *
    frame$non_target_events / frame$non_target_denominator

  append_sequential_columns(data.frame(
    month_date = as.Date(sprintf("%d-%02d-01", frame$year, frame$month)),
    observed = frame$events_target,
    expected = frame$expected,
    denominator = frame$denominator_target,
    denominator_all = frame$denominator_all,
    events_target = frame$events_target,
    events_all = frame$events_all,
    non_target_denominator = frame$non_target_denominator,
    non_target_events = frame$non_target_events,
    label = label,
    source = "VAERS",
    interpretation = sprintf(
      "Contemporaneous non-%s report fraction comparator; report-level only, not an unexposed patient comparator.",
      target_vax_type
    )
  ))
}

load_real_series <- function(id, processed_root) {
  switch(
    id,
    eliquis = load_openfda_series(processed_root),
    hpv9_syncope = load_vaers_series(
      "HPV9", "Syncope", "HPV9 (GARDASIL 9) / Syncope", processed_root
    ),
    menb_pyrexia = load_vaers_series(
      "MENB", "Pyrexia", "MENB / Pyrexia", processed_root
    ),
    stop("Unknown real-data preset: ", id)
  )
}

poisson_gof <- function(series) {
  observed <- series$observed
  expected <- series$expected
  pearson_residual <- (observed - expected) / sqrt(expected)
  deviance <- 2 * sum(ifelse(
    observed == 0,
    expected,
    observed * log(observed / expected) - (observed - expected)
  ))
  pearson <- sum(pearson_residual^2)
  data.frame(
    deviance = deviance,
    pearson = pearson,
    dispersion = pearson / nrow(series),
    lag1_residual_correlation = cor(
      pearson_residual[-1], pearson_residual[-length(pearson_residual)]
    )
  )
}
