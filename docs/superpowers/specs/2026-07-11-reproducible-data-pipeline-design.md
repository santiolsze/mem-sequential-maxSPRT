# Reproducible Data Pipeline Design

## Goal

Version only the source data in `experimentacion/data/raw/` and regenerate every
file in `experimentacion/data/processed/` locally from those inputs.

## Version-control policy

- Keep all raw inputs in Git.
- Store the large VAERS yearly ZIP archives (`raw/vaers/*VAERSData.zip`) through
  Git LFS.
- Keep small raw files, including openFDA daily Parquet and its metadata JSON, as
  ordinary Git files.
- Ignore `experimentacion/data/processed/` completely. It is a local build output,
  including the Parquet datasets, manifest, and monthly openFDA CSV export.

This policy preserves the source material needed to reproduce the analysis while
avoiding ordinary-Git history growth from the large VAERS archives. It does not
use LFS for small files merely because they happen to be below `raw/`.

## Data flow

```text
openFDA API --fetch script--> raw/openfda/daily Parquet + metadata
                                     |
VAERS yearly ZIP archives ----------+--> build_parquet_datasets.py --> processed/
```

`fetch_openfda_eliquis_comparator.R` will only fetch and persist raw openFDA
daily counts and provenance metadata. `build_parquet_datasets.py` will read
those daily counts, calculate the monthly ELIQUIS comparator and expected count,
then write the openFDA CSV export and partitioned Parquet dataset under
`processed/`. The same builder will continue to derive the four VAERS Parquet
datasets from the raw ZIP files.

## Failure handling

- The builder fails with a direct error when the required raw openFDA dataset or
  VAERS ZIP archives are absent.
- The builder retains its current clean-build behavior: it removes stale
  processed outputs, then produces them again from raw data.
- The fetch script remains the explicit network-refresh step. Re-fetching may
  change raw data when the upstream API changes; rebuilding from existing raw
  data is local and deterministic.

## Tests and documentation

- Add tests that demonstrate openFDA monthly output is derived from the raw daily
  dataset, not a committed processed CSV.
- Retain coverage of the expected-count calculation and output integrity.
- Document prerequisites, the raw-to-processed rebuild command, the LFS policy,
  and the API-refresh command in `experimentacion/data/README.md`.

## Acceptance criteria

1. No `experimentacion/data/processed/**` path is tracked by Git.
2. The only raw files configured for LFS are the VAERS yearly ZIP archives.
3. With the raw files present, `python3 experimentacion/scripts/build_parquet_datasets.py`
   creates all required processed outputs, including the openFDA monthly CSV and
   Parquet dataset.
4. The app and its tests read only generated processed outputs.
5. Documentation describes the exact rebuild and refresh commands.
