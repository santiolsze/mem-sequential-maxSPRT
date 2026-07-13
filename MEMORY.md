# Project Memory

- 2026-07-11: Converted local openFDA and VAERS data to compressed Parquet with
  `zstd` using `experimentacion/scripts/build_parquet_datasets.py`.
- Processed openFDA control series is at
  `experimentacion/data/processed/fda/openfda_eliquis_haemorrhage_monthly/`,
  partitioned by `product_slug`, `event_slug`, `year`, `month_num`.
- Processed VAERS datasets are under `experimentacion/data/processed/vaers/`:
  `reports`, `vaccines`, `symptoms_long`, and `vaccine_symptoms`.
- VAERS joins must use the composite record key. The Parquet conversion attaches
  symptoms to vaccines using `source_year`, `vaers_id`, and `ORDER`
  (`vax_record_order` / `symptoms_record_order`) to avoid crossing follow-up
  records.
- For R analysis, prefer `vaers/vaccine_symptoms` when filtering by vaccine
  product and symptom term. It is partitioned by `year`, `month`,
  `vax_type_slug`, and `medicine_slug`; filter by `medicine_slug` for partition
  pruning, then use `medicine_name` for display.
- Manifest with row counts, columns, partitions, and compression is at
  `experimentacion/data/processed/parquet_manifest.json`.
- 2026-07-11: The retained data horizon is 2016--2025. VAERS archives and
  Parquet partitions for 2013--2015 were removed. The openFDA series was
  rebased to 2016--2017 for expected counts, so the earlier 2013--2014 positive
  control result is historical context only and cannot be reproduced locally.
- 2026-07-11: Approved scope for an interactive Shiny demonstration:
  - real-data Poisson MaxSPRT with local Parquet datasets;
  - openFDA ELIQUIS/Haemorrhage as the signal example. This must use a
    contemporaneous FAERS report-level comparator, not ELIQUIS's own historical
    haemorrhage fraction: expected monthly reports are
    `ELIQUIS reports * (non-ELIQUIS haemorrhage reports / non-ELIQUIS reports)`;
  - VAERS HPV9 (GARDASIL 9)/Syncope as the signal example, using the same
    baseline and surveillance split (1,614 observed vs 1,229 expected,
    cumulative RR 1.313, maximum LLR 54.82, paper-table crossing in 2018-08);
  - optional VAERS MENB/Pyrexia no-signal preset (910 observed vs 1,597
    expected; maximum LLR 0);
  - classic Poisson SPRT comparisons with low and high fixed alternatives;
  - a separate simulated binomial MaxSPRT panel. VAERS does not provide the
    exposure-time or population denominator needed for a valid binomial
    surveillance interpretation; report counts must not be presented as all
    vaccine recipients, so that panel remains explicitly didactic.
- The Shiny interface must explain that alpha is chosen before surveillance,
  the MaxSPRT critical value V is numerically calibrated from alpha and the
  surveillance horizon T (Monte Carlo for the app's monthly looks), beta sets
  the classic-SPRT boundary, and all FAERS/VAERS findings are report signals,
  not causal or incidence claims.
- 2026-07-11: Verified the openFDA comparison aggregate via the official API
  for 2016--2025. It gives 16,908 observed ELIQUIS/Haemorrhage reports versus
  5,073.77 expected from the non-ELIQUIS report fraction (report RR 3.33,
  maximum LLR 8,517.97, first crossing in 2016-01). This is a contrast among
  spontaneous report records, not evidence of a clinical effect in people
  taking ELIQUIS versus people not taking it.
- 2026-07-11: Built the local Shiny application in `experimentacion/shiny/`.
  It reads only existing local Parquet datasets at runtime and has four panels:
  real report series, Poisson goodness-of-fit diagnostics, formal Poisson
  simulations comparing MaxSPRT with low/high classical-SPRT alternatives, and
  a fully simulated binomial MaxSPRT panel.
- The real-data presets are ELIQUIS/Haemorrhage, HPV9/Syncope, and
  MENB/Pyrexia. Their displayed MaxSPRT result is explicitly labelled
  mechanical because their local expected counts are estimated from reporting
  patterns and the simple Poisson model can fail. The diagnostics panel exposes
  deviance, Pearson dispersion, residual autocorrelation, and residuals.
- `testthat` 3.3.2 was installed locally. Run tests with
  `Rscript experimentacion/shiny/tests/run_tests.R`.
  Run the app with
  `Rscript -e 'shiny::runApp("experimentacion/shiny", host = "127.0.0.1", port = 3838)'`.
- 2026-07-11: ELIQUIS/Haemorrhage has one canonical local aggregate only:
  `experimentacion/data/processed/fda/openfda_eliquis_haemorrhage_monthly/`.
  Refresh it with `Rscript experimentacion/scripts/fetch_openfda_eliquis_comparator.R`.
  The app reads this Parquet directly; no alternate comparator Parquet remains.
- The monthly Poisson calibration chooses a discrete critical value with
  simulated tail probability at most alpha and uses at least `100 / alpha`
  null simulations (100,000 for alpha = 0.001).
