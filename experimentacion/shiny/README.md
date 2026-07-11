# MaxSPRT Shiny Lab

Launch from the repository root:

```r
shiny::runApp("experimentacion/shiny", host = "127.0.0.1", port = 3838)
```

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
  Haemorrhage, HPV9 / Syncope, and MENB / Pyrexia presets.
- **Diagnostico Poisson** reports deviance, Pearson dispersion, residual
  autocorrelation, and the monthly Pearson residuals for the selected simple
  null model.
- **Potencia simulada** uses generated Poisson trajectories with known expected
  counts. This is the formal illustration of MaxSPRT, classical SPRT, `alpha`,
  `beta`, the critical value, and stopping time.
- **Binomial simulado** reproduces the exposed/unexposed case design from the
  paper. It is not derived from VAERS.

## Interpretation

FAERS/openFDA and VAERS contain spontaneous reports rather than all exposed
people or administered doses. The real-data panels are useful for method
mechanics and diagnostics, but do not estimate incidence or establish causality.
For the shipped real-data presets, expected counts are estimated from report
patterns; the nominal Poisson MaxSPRT error rate is therefore not a formal
guarantee. The diagnostics tab makes this limitation visible.

The ELIQUIS expected count is the same-month non-ELIQUIS Haemorrhage reporting
fraction multiplied by ELIQUIS reports. The app reads only
`data/processed/fda/openfda_eliquis_haemorrhage_monthly/`; no alternate ELIQUIS
Parquet aggregate is used.

The HPV9/Syncope and MENB/Pyrexia expected counts use the same contemporaneous
design: the same-month reporting fraction of the target symptom among reports
for every *other* vaccine type, multiplied by the target vaccine's reports that
month. This avoids attributing a generic time trend in report volume or
coding practice to the target vaccine, but it does not control for age or
setting-specific confounders shared across vaccines (e.g. Syncope and Pyrexia
are common post-injection events regardless of vaccine).
