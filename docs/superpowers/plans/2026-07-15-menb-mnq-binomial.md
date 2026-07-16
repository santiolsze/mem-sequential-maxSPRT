# MENB vs MNQ Binomial MaxSPRT Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a fixed Shiny tab that performs a monthly-stratified conditional binomial MaxSPRT comparison of MENB versus MNQ Pyrexia reports, and correct the binomial explanation in the Typst deck.

**Architecture:** Pure statistical functions in `R/sequential.R` compute the stratified GLR path and a reproducible conditional Monte Carlo boundary. `R/data.R` supplies the fixed monthly MENB/MNQ margins; `app.R` only renders the analysis. The test suite covers numerical identities, calibration, Parquet integration, and Shiny outputs before production code is changed.

**Tech Stack:** R, Shiny, testthat, Arrow/Dplyr, ggplot2/plotly, DT, Typst.

## Global Constraints

- The estimand is relative Pyrexia reporting among VAERS reports, not clinical incidence or causality.
- Use 2016--2025 and a fixed MENB versus MNQ comparison independent of the Poisson preset/start month.
- Calibrate conditionally on observed monthly Pyrexia totals and MENB/MNQ report margins.
- Do not use `p0 = 1/2` or the paper's fixed-`z` critical-value table.

---

### Task 1: Stratified binomial engine

**Files:**
- Modify: `experimentacion/shiny/R/sequential.R`
- Test: `experimentacion/shiny/tests/testthat/test-sequential.R`

**Interfaces:**
- Produces: `stratified_binomial_maxsprt(target_events, comparator_events, target_reports, comparator_reports)` returning `llr`, `rr_hat`, and `p0`; `calibrate_stratified_binomial_boundary(..., alpha, reps, seed)` returning `critical_value`, `null_maxima`, `reps`, and `alpha`.

- [ ] **Step 1: Write failing numerical tests** for the constant-`p0` identity, zero evidence at the null allocation, positive evidence under MENB excess, zero-event months, deterministic seed, and conservative empirical crossing probability.
- [ ] **Step 2: Run** `Rscript -e 'testthat::test_file("experimentacion/shiny/tests/testthat/test-sequential.R")'` and verify failures report missing stratified functions.
- [ ] **Step 3: Implement the minimal engine.** For each prefix maximize
  `sum(x_i log(p_i(rr)) + (n_i-x_i) log1p(-p_i(rr))) - loglik(rr=1)`
  over `log(rr) >= 0` using `optimize`, with a finite data-derived upper bound and stable zero-count handling. Simulate monthly binomials under `p0_i` and reuse `conservative_critical_value` for calibration.
- [ ] **Step 4: Re-run the focused test** and verify it passes.

### Task 2: Fixed MENB/MNQ data series

**Files:**
- Modify: `experimentacion/shiny/R/data.R`
- Test: `experimentacion/shiny/tests/testthat/test-data.R`

**Interfaces:**
- Consumes: existing VAERS `vaccines` and `vaccine_symptoms` Parquet datasets.
- Produces: `load_menb_mnq_pyrexia_binomial(processed_root)` with `month_date`, `target_reports`, `comparator_reports`, `target_events`, `comparator_events`, and monthly proportions.

- [ ] **Step 1: Write a failing integration test** asserting 120 ordered months, nonnegative events, positive report margins, events no larger than reports, and totals identical to the existing MENB preset columns.
- [ ] **Step 2: Run** `Rscript -e 'testthat::test_file("experimentacion/shiny/tests/testthat/test-data.R")'` and verify the loader is missing.
- [ ] **Step 3: Extract the existing four-count VAERS aggregation into the fixed loader** without continuity correction, preserving composite-record deduplication.
- [ ] **Step 4: Re-run the focused data test** and verify it passes.

### Task 3: Shiny tab

**Files:**
- Modify: `experimentacion/shiny/app.R`
- Modify: `experimentacion/shiny/www/app.css` only if the existing responsive metric layout is insufficient.
- Test: `experimentacion/shiny/tests/testthat/test-app.R`

**Interfaces:**
- Consumes: fixed loader plus stratified path/calibration functions.
- Produces: tab `Binomial MENB vs MNQ`; reactives `binomial_series`, `binomial_analysis`, `binomial_boundary`; outputs `binomial_summary`, `binomial_proportion_plot`, `binomial_llr_plot`, `binomial_decision`, and `binomial_monthly_table`.

- [ ] **Step 1: Replace the obsolete absence test with failing UI/server tests** asserting the fixed tab copy, output IDs, 120-month series, finite path, calibrated boundary, and a coherent first-crossing result.
- [ ] **Step 2: Run** `Rscript -e 'testthat::test_file("experimentacion/shiny/tests/testthat/test-app.R")'` and verify failures are caused by the missing tab/reactives.
- [ ] **Step 3: Add the fixed tab and server reactives.** Reuse general `alpha`; use at least `calibration_reps(alpha, selected repetitions)` and seed `20260715L`; truncate the displayed LLR after first rejection; state conditional-alpha and VAERS interpretation caveats in the tab.
- [ ] **Step 4: Add metrics, plots, decision text, and a monthly DT** using the interfaces above.
- [ ] **Step 5: Re-run the focused app test** and verify it passes.

### Task 4: Slides and durable documentation

**Files:**
- Modify: `presentation/maxsprt.typ`
- Modify: `experimentacion/shiny/README.md`
- Modify: `MEMORY.md`
- Test: `experimentacion/shiny/tests/testthat/test-app.R`

**Interfaces:**
- Documents the same definitions implemented in Tasks 1--3.

- [ ] **Step 1: Add failing source assertions** that the slide identifies `z=1` as the special `1/2` case, uses `(alpha,N,z)` for the binomial boundary, and labels the later return to the Poisson `T` example.
- [ ] **Step 2: Run the focused app/source tests** and verify they fail on the old wording.
- [ ] **Step 3: Edit the Typst text** to define `z` as the control/exposed opportunity ratio, state the `z=1` special case, use `(alpha,N,z)` for the binomial boundary, and label the later `T` section as a return to Poisson.
- [ ] **Step 4: Document the new tab and its conditional/report-level interpretation** in README and MEMORY.
- [ ] **Step 5: Compile** `typst compile presentation/maxsprt.typ presentation/maxsprt.pdf` and verify exit code 0.

### Task 5: Full verification

**Files:**
- Verify all modified files; do not change unrelated user edits.

- [ ] **Step 1: Run** `Rscript experimentacion/shiny/tests/run_tests.R` and require zero failures.
- [ ] **Step 2: Run** `git diff --check` and require no whitespace errors.
- [ ] **Step 3: Inspect** `git diff --stat` and `git status --short`, separating pre-existing presentation/MEMORY changes from this feature in the handoff.
