#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(arrow)
  library(jsonlite)
})

script_path <- sub(
  "^--file=", "", commandArgs(trailingOnly = FALSE)[
    grep("^--file=", commandArgs(trailingOnly = FALSE))
  ][1]
)
root <- normalizePath(file.path(dirname(script_path), "..", ".."))

api_url <- "https://api.fda.gov/drug/event.json"
start_date <- as.Date("2016-01-01")
end_date <- as.Date("2025-12-31")
date_query <- "receivedate:[20160101 TO 20251231]"

queries <- c(
  all_reports = date_query,
  eliquis_reports = paste(
    "patient.drug.medicinalproduct:\"ELIQUIS\" AND", date_query
  ),
  haemorrhage_reports = paste(
    "patient.reaction.reactionmeddrapt:\"Haemorrhage\" AND", date_query
  ),
  eliquis_haemorrhage_reports = paste(
    "patient.drug.medicinalproduct:\"ELIQUIS\" AND",
    "patient.reaction.reactionmeddrapt:\"Haemorrhage\" AND", date_query
  )
)

data_dir <- file.path(root, "experimentacion", "data")
raw_path <- file.path(
  data_dir, "raw", "openfda", "eliquis_haemorrhage_daily"
)
processed_path <- file.path(
  data_dir, "processed", "fda", "openfda_eliquis_haemorrhage_monthly"
)
metadata_path <- file.path(
  data_dir, "raw", "openfda", "eliquis_haemorrhage_metadata.json"
)
manifest_path <- file.path(data_dir, "processed", "parquet_manifest.json")
legacy_paths <- c(
  file.path(data_dir, "raw", "openfda", "eliquis_haemorrhage_comparator_daily"),
  file.path(data_dir, "processed", "fda", "openfda_eliquis_haemorrhage_comparator_monthly"),
  file.path(data_dir, "raw", "openfda", "eliquis_haemorrhage_comparator_metadata.json")
)

fetch_daily_counts <- function(query) {
  request_url <- paste0(
    api_url,
    "?search=", URLencode(query, reserved = TRUE),
    "&count=receivedate"
  )
  response <- fromJSON(request_url)

  if (is.null(response$results)) {
    stop("openFDA returned no daily count results for query: ", query)
  }

  list(
    data = data.frame(
      received_date = as.Date(response$results$time, format = "%Y%m%d"),
      count = as.numeric(response$results$count)
    ),
    meta = response$meta
  )
}

responses <- lapply(queries, fetch_daily_counts)
calendar <- data.frame(received_date = seq(start_date, end_date, by = "day"))
daily <- calendar

for (name in names(responses)) {
  counts <- responses[[name]]$data
  names(counts)[2] <- name
  daily <- merge(daily, counts, by = "received_date", all.x = TRUE)
}

count_columns <- names(queries)
daily[count_columns] <- lapply(daily[count_columns], function(column) {
  column[is.na(column)] <- 0
  as.numeric(column)
})
daily$year <- as.integer(format(daily$received_date, "%Y"))
daily$month_num <- as.integer(format(daily$received_date, "%m"))
daily$month <- format(daily$received_date, "%Y-%m")

monthly <- aggregate(daily[count_columns], by = list(month = daily$month), FUN = sum)
monthly$non_eliquis_reports <-
  monthly$all_reports - monthly$eliquis_reports
monthly$non_eliquis_haemorrhage_reports <-
  monthly$haemorrhage_reports - monthly$eliquis_haemorrhage_reports
monthly$expected_under_no_eliquis_association <- with(
  monthly,
  eliquis_reports * non_eliquis_haemorrhage_reports / non_eliquis_reports
)
monthly$year <- as.integer(substr(monthly$month, 1, 4))
monthly$month_num <- as.integer(substr(monthly$month, 6, 7))
monthly <- monthly[, c(
  "month", "all_reports", "eliquis_reports", "haemorrhage_reports",
  "eliquis_haemorrhage_reports", "non_eliquis_reports",
  "non_eliquis_haemorrhage_reports",
  "expected_under_no_eliquis_association", "year", "month_num"
)]

if (nrow(monthly) != 120L || any(monthly$non_eliquis_reports <= 0) ||
    any(monthly$expected_under_no_eliquis_association <= 0)) {
  stop("The downloaded comparator aggregate failed its monthly integrity checks.")
}

write.csv(
  monthly,
  file.path(data_dir, "processed", "openfda_eliquis_haemorrhage_monthly.csv"),
  row.names = FALSE
)

write_parquet_dataset <- function(frame, path) {
  if (dir.exists(path)) {
    unlink(path, recursive = TRUE, force = TRUE)
  }
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  write_dataset(
    Table$create(frame),
    path = path,
    format = "parquet",
    partitioning = c("year", "month_num"),
    compression = "zstd"
  )
}

write_parquet_dataset(daily, raw_path)
write_parquet_dataset(monthly, processed_path)

api_last_updated <- unique(vapply(
  responses,
  function(response) response$meta$last_updated %||% NA_character_,
  character(1)
))
metadata <- list(
  source = "openFDA drug adverse-event API, sourced from FAERS",
  api_url = api_url,
  fetched_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  api_last_updated = unname(api_last_updated),
  period = "2016-01-01 through 2025-12-31",
  aggregation = "Daily report-level API counts rolled up to calendar month",
  queries = as.list(queries),
  raw_parquet = normalizePath(raw_path, winslash = "/", mustWork = FALSE),
  processed_parquet = normalizePath(processed_path, winslash = "/", mustWork = FALSE),
  limitation = paste(
    "The expected count compares ELIQUIS-tagged FAERS reports with",
    "non-ELIQUIS-tagged reports in the same month. It is not a patient-level",
    "exposed-versus-unexposed comparison and cannot establish incidence or causality."
  )
)
dir.create(dirname(metadata_path), recursive = TRUE, showWarnings = FALSE)
write_json(metadata, metadata_path, pretty = TRUE, auto_unbox = TRUE)

if (file.exists(manifest_path)) {
  manifest <- read_json(manifest_path, simplifyVector = FALSE)
  dataset_entry <- list(
    name = "openfda_eliquis_haemorrhage_monthly",
    path = "experimentacion/data/processed/fda/openfda_eliquis_haemorrhage_monthly",
    rows = nrow(monthly),
    columns = names(monthly),
    partitions = c("year", "month_num"),
    compression = "zstd"
  )
  manifest$datasets <- Filter(
    function(item) item$name != "openfda_eliquis_haemorrhage_comparator_monthly",
    manifest$datasets
  )
  dataset_names <- vapply(manifest$datasets, function(item) item$name, character(1))
  existing <- which(dataset_names == dataset_entry$name)
  if (length(existing) == 0L) {
    manifest$datasets[[length(manifest$datasets) + 1L]] <- dataset_entry
  } else {
    manifest$datasets[[existing[1]]] <- dataset_entry
  }
  write_json(manifest, manifest_path, pretty = TRUE, auto_unbox = TRUE)
}

for (legacy_path in legacy_paths) {
  if (dir.exists(legacy_path)) {
    unlink(legacy_path, recursive = TRUE, force = TRUE)
  } else if (file.exists(legacy_path)) {
    unlink(legacy_path, force = TRUE)
  }
}

message(
  "Wrote ", nrow(monthly), " monthly rows to ", processed_path,
  " (", min(monthly$month), " through ", max(monthly$month), ")."
)
