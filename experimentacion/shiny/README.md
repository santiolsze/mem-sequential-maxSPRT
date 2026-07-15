# MaxSPRT Shiny Lab

## Run the app locally

From the repository root, start the Shiny server with:

```sh
Rscript -e 'shiny::runApp("experimentacion/shiny", host = "127.0.0.1", port = 3838)'
```

Then open <http://127.0.0.1:3838> in a web browser. Keep the terminal process
running while using the app; stop it with `Ctrl+C`.

Required R packages:

- `shiny`, `bslib`, `arrow`, `dplyr`, `ggplot2`, `plotly`, `DT`
- `testthat` for the test suite

Run the tests with:

```sh
Rscript experimentacion/shiny/tests/run_tests.R
```

Refresh the sole local ELIQUIS / Haemorrhage Parquet aggregate with:

```sh
Rscript experimentacion/scripts/fetch_openfda_eliquis_comparator.R
```

## Panels

- **Datos reales** reads only the local Parquet datasets and offers ELIQUIS /
  Haemorrhage, HPV9 / Syncope, and MENB / Pyrexia presets. The `Mes inicial de
  vigilancia` control restarts the analysis at the selected month: cumulative
  observed and expected counts, LLR paths, sequential decisions, the remaining
  horizon `T`, and the monthly Monte Carlo boundary are all recomputed from
  that point.
- **Diagnostico Poisson** reports deviance, Pearson dispersion, residual
  autocorrelation, and the monthly Pearson residuals for the selected simple
  null model over the same selected surveillance period.
- **Potencia simulada** uses generated Poisson trajectories with known expected
  counts. This is the formal illustration of MaxSPRT, classical SPRT, `alpha`,
  `beta`, the critical value, and stopping time.
- **Discretizacion** is a dataset-independent Poisson experiment. It keeps the
  exact continuous MaxSPRT boundary fixed and evaluates trajectories at
  `K` equally spaced cumulative-expected-count looks, `mu_k = kT/K`. It compares
  empirical alpha and power, with Wilson 95% intervals, against their exact
  continuous references as `K` increases. Available grids extend through
  `K = 10000`; values above 1000 are optional because runtime grows with the
  total number of simulated looks. Selecting `K = 10000` caps the effective
  Monte Carlo repetitions at 100,000. The repetition selector offers 1,000,
  5,000, 10,000, 20,000, and 50,000 simulations per scenario.
## Interpretation

FAERS/openFDA and VAERS contain spontaneous reports rather than all exposed
people or administered doses. The real-data panels are useful for method
mechanics and diagnostics, but do not estimate incidence or establish causality.
For the shipped real-data presets, expected counts are estimated from report
patterns; the nominal Poisson MaxSPRT error rate is therefore not a formal
guarantee. The diagnostics tab makes this limitation visible.

The discretization tab does not read FAERS, VAERS, Parquet files, or any real
series. Its continuous boundary and power come from the exact event-time
recursion used to reproduce Tables 1-3 of the paper. Because that same
continuous boundary is used at every finite `K`, sparse observation schedules
are intentionally conservative; recalibrating a separate boundary for each
`K` would hide the effect this experiment is designed to show.

The ELIQUIS expected count is the same-month non-ELIQUIS Haemorrhage reporting
fraction multiplied by ELIQUIS reports. The app reads only
`data/processed/fda/openfda_eliquis_haemorrhage_monthly/`; no alternate ELIQUIS
Parquet aggregate is used.

The HPV9/Syncope and MENB/Pyrexia expected counts use a contemporaneous,
age-matched single-vaccine comparator: MNQ (meningococcal ACWY conjugate).
VAERS report ages show HPV9 (median 14, IQR 11-17) and MNQ (median 15, IQR
11-17) are given at essentially the same adolescent visit, so MNQ's same-month
reporting fraction for the target symptom is a much closer match than a pool
of every other vaccine type, which spans all ages and would conflate the
target vaccine's effect with age or vaccination-setting differences (e.g.
Syncope is a generic post-injection event, more common in adolescents than in
infants or older adults, regardless of which vaccine is given). Using MNQ as
the comparator, HPV9's Syncope reporting rate is statistically indistinguishable
from MNQ's (9.6% vs 9.6%), while MENB's Pyrexia rate remains more than double
MNQ's (14.4% vs 6.3%) even after age-matching. A Haldane-Anscombe continuity
correction (`+0.5` events, `+1` denominator) is applied to the comparator rate
because a single comparator vaccine can have zero events in a given month.
