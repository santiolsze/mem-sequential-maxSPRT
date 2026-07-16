# Project Memory

- 2026-07-16: Removed the active `Diagnostico Poisson` Shiny tab, its guide entry, and its server outputs. The reusable `poisson_gof()` implementation and its unit tests remain available. Also removed the global FAERS/VAERS caveat strip requested for the live demo. Extended the discretization selector with optional `K = 25000` and `K = 50000`; selecting any `K >= 25000` caps effective Monte Carlo repetitions at 5,000 and reports the adjustment in the UI.

- 2026-07-16: Reorganized the presentation's experimental handoff into three slides placed before the Shiny separator: (1) two public spontaneous-report sources and three 2016--2025 target series, (2) separate contemporaneous comparator constructions for ELIQUIS and the VAERS vaccine series plus the distinct binomial MENB--MNQ panel, and (3) what report-level signals can and cannot establish. Accumulated observed/expected totals were removed from the deck because they are shown interactively in Shiny. The Shiny separator now follows these slides, and the closing separator follows the live-demo handoff.
- 2026-07-16: Reorganized the Typst deck into seven explicit narrative sections with consistent divider slides and `01 / 07`--`07 / 07` counters: public-health context, classical SPRT, Poisson MaxSPRT, binomial variant, Poisson threshold calculation, periodic/discrete looks, and the real-data demonstration. Short technical labels remain in the top bar while the centered divider titles state the audience-facing question or takeaway. Existing technical slides remain in their prior order.
- 2026-07-16: Reworked the Typst deck introduction into two slides. The first establishes the CDC as a major U.S. federal public-health institution, gives the FY2025 requested budget and VSD scale, and accurately states that MaxSPRT arose from a concrete need in the CDC-sponsored VSD rather than calling the paper a CDC publication. The second motivates post-market pharmacovigilance through rare, delayed, and subgroup-specific adverse events and closes with the repeated-monitoring/false-alarm question. CDC sources are listed only in the references slide, not in the audience-facing narrative.

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
- 2026-07-13: Added an exact, deterministic base-R replication of Poisson
  MaxSPRT Tables 1-3 from `paper.pdf` at
  `experimentacion/scripts/replicate_paper_tables_1_3.R`. It implements both
  real Lambert W branches with `uniroot()`, uses the principal branch for the
  rejection crossing time, and propagates nonabsorbed Poisson probabilities at
  event times. Run it from the repository root with
  `Rscript experimentacion/scripts/replicate_paper_tables_1_3.R`.
- Replicated table and comparison CSVs are written under
  `experimentacion/results/paper_tables/`. All 523 compared cells pass the
  project tolerance of less than one unit in the paper's last displayed
  decimal. Three cells differ by about 0.5-0.6 units in their last displayed
  decimal; tighter roots and an independent direct convolution confirmed these
  are paper-level numerical/rounding differences rather than implementation
  errors.
- 2026-07-15: Added a dataset-independent `Discretizacion` Shiny tab design and
  implementation. For fixed `T`, nominal alpha, and RR, it calibrates one exact
  continuous Poisson MaxSPRT boundary and compares continuous alpha/power with
  Monte Carlo estimates observed only at `K` equally spaced expected-count
  looks `mu_k = kT/K`. The continuous boundary is deliberately not recalibrated
  for each `K`, so the conservative effect of sparse looks remains visible.
- Discretization simulation helpers live in
  `experimentacion/shiny/R/discretization.R`; focused tests are in
  `experimentacion/shiny/tests/testthat/test-discretization.R`. The experiment
  has its own Shiny inputs and random seeds and does not access FAERS, VAERS,
  Parquet data, presets, or real-series reactives.
- The discretization selector includes K values through 10,000. K=1,000 is
  selected by default; 2,500, 5,000, and 10,000 are optional. Monte Carlo
  repetition choices include 1,000, 5,000, 10,000, 20,000, and 50,000. When
  K=10,000 is selected, requested repetitions are capped at 100,000 and the UI
  explains any adjustment, since runtime is proportional to repetitions times
  total looks.
- In the `Potencia simulada` table, `Tiempo medio hasta el rechazo (meses)` is
  conditional on rejection for every method. Early acceptances and simulations
  without a decision are excluded from that average; each look is one month in
  the retained monthly profile.
- 2026-07-15: Removed the simulated binomial tab, its sidebar controls, and its
  Shiny server outputs from the active app. The reusable binomial functions and
  mathematical tests remain in `R/sequential.R` and the test suite.
- 2026-07-15: The real-data sequential plot labels the monthly Monte Carlo
  MaxSPRT boundary and the common theoretical Wald upper boundary for the two
  classical SPRTs. It marks each valid first rejection and includes a table of
  method, boundary, and first rejection month. A later upper crossing is not
  counted when a classical SPRT already accepted H0 at its lower boundary.
- 2026-07-15: The Poisson diagnostics tab now explains that it assesses whether
  monthly variance is compatible with the expected Poisson mean and whether
  residuals show temporal dependence. It explicitly connects overdispersion or
  autocorrelation with loss of the nominal MaxSPRT alpha guarantee.
- 2026-07-15: Removed the `Metodo` tab and its guide entry from the active Shiny
  interface. Statistical implementations and documentation elsewhere remain.
- 2026-07-15: Added alpha 0.10 to both the general Shiny alpha control and the
  independent discretization experiment control.
- 2026-07-15: Real-data LLR trajectories are truncated immediately after their
  first valid sequential decision: rejection for MaxSPRT, and rejection or
  early H0 acceptance for each classical SPRT.
- The real-data decision table reports `Rechaza en YYYY-MM` when MaxSPRT
  signals. Otherwise it reports the formal surveillance completion as
  `Finaliza sin rechazo al alcanzar T = ... en YYYY-MM`, where T is the final
  cumulative expected count under H0 rather than an observed-event limit.
- 2026-07-15: Added `Mes inicial de vigilancia` to the Shiny sidebar. It filters
  the chosen real series and rebuilds cumulative observed/expected counts, RR,
  MaxSPRT LLR, Monte Carlo boundary, sequential decisions, remaining T, plots,
  monthly table, simulation profile, and Poisson diagnostics from that month.
  The selector uses YYYY-MM labels, defaults to the first available month, and
  refreshes when the real-data preset changes.
- 2026-07-15: Added a fixed `Binomial MENB vs MNQ` Shiny analysis for Pyrexia.
  It conditions on each month's combined MENB/MNQ Pyrexia reports, uses the
  monthly allocation probability `MENB reports / (MENB reports + MNQ reports)`,
  fits one common reporting RR above one, and calibrates the monthly GLR
  boundary by conditional Monte Carlo. Its alpha statement is conditional on
  the observed monthly event totals and report margins; it compares reporting
  proportions, not incidence or causality. The presentation now defines `z`
  as the control/exposed opportunity ratio, states that `P(exposed)=1/2` only
  for `z=1`, and explicitly labels the later `T` example as a return to the
  Poisson case.
- 2026-07-15: Expanded the Typst presentation's Table 1--3 construction section.
  It now derives the crossing point from the LLR equation through the
  substitution `x = mu(s_n)/n`, puts it into Lambert W canonical form, explains
  why the principal branch is used, shows the Poisson non-absorbed-state
  recursion, inverts `alpha(V,T)` numerically to obtain the critical value, and
  maps the same recursion to power and expected-time quantities in Tables 2--3.
- 2026-07-15: Added a standalone, base-R pedagogical replication of the exact
  binomial MaxSPRT alpha calibration at
  `experimentacion/scripts/replicate_binomial_alpha.R`, with focused tests at
  `experimentacion/tests/testthat/test-replicate-binomial-alpha.R`. It uses an
  absorbing Markov chain and an independent enumeration of all 1,024 paths for
  `z = 1`, `N = 10`. It reproduces Table 4's displayed boundary `V = 2.77259`
  and effective alpha `0.041015625`; including equality at `4 log(2)` would give
  alpha `0.0703125`. Runtime CSV walkthroughs are written under
  `experimentacion/results/binomial_alpha/`. This is a private teaching script
  and is deliberately not integrated into Shiny.
- 2026-07-15: Replaced the presentation's generic critical-value trajectory
  image with a native Typst/CeTZ conceptual diagram. The new slide shows that
  the Poisson MaxSPRT LLR decreases while cumulative expected events grow
  between observed events, jumps when an event arrives, rejects at the
  horizontal boundary `V`, and otherwise stops at the vertical horizon
  `T = mu(t_max)`. It explicitly states that when `lambda_0(t) > 0`, calendar
  time and cumulative expected events are in bijective correspondence.
