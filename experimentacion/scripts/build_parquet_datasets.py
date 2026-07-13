#!/usr/bin/env python3
"""Build compressed Parquet datasets from local openFDA and VAERS files.

The script keeps VAERS normalized for reproducibility and also writes a
vaccine-symptom fact table for fast exploratory filters in R.
"""

from __future__ import annotations

import json
import re
import shutil
from pathlib import Path
from zipfile import ZipFile

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq


ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = ROOT / "experimentacion" / "data"
RAW_DIR = DATA_DIR / "raw"
PROCESSED_DIR = DATA_DIR / "processed"

COMPRESSION = "zstd"
ANALYSIS_START_YEAR = 2016
ANALYSIS_END_YEAR = 2025
FDA_CSV = PROCESSED_DIR / "openfda_eliquis_haemorrhage_monthly.csv"
FDA_OUT = PROCESSED_DIR / "fda" / "openfda_eliquis_haemorrhage_monthly"
FDA_COMPARATOR_COLUMNS = {
    "month",
    "all_reports",
    "eliquis_reports",
    "haemorrhage_reports",
    "eliquis_haemorrhage_reports",
    "non_eliquis_reports",
    "non_eliquis_haemorrhage_reports",
    "expected_under_no_eliquis_association",
}

VAERS_RAW = RAW_DIR / "vaers"
VAERS_OUT = PROCESSED_DIR / "vaers"

MANIFEST_OUT = PROCESSED_DIR / "parquet_manifest.json"


DATE_COLUMNS = {
    "received_date",
    "report_date",
    "date_died",
    "vax_date",
    "onset_date",
    "todays_date",
}

NUMERIC_COLUMNS = {
    "age_yrs",
    "cage_yr",
    "cage_mo",
    "hospdays",
    "numdays",
    "form_vers",
}

SERIOUS_FLAGS = [
    "died",
    "l_threat",
    "hospital",
    "x_stay",
    "disable",
    "birth_defect",
]


def slugify(value: object) -> str:
    text = "" if pd.isna(value) else str(value)
    text = text.lower().strip()
    text = re.sub(r"[^a-z0-9]+", "_", text)
    text = re.sub(r"_+", "_", text).strip("_")
    return text or "unknown"


def snake_case(name: str) -> str:
    return slugify(name)


def in_analysis_range(year: int) -> bool:
    return ANALYSIS_START_YEAR <= year <= ANALYSIS_END_YEAR


def read_zip_csv(zip_path: Path, member: str) -> pd.DataFrame:
    with ZipFile(zip_path) as zf:
        with zf.open(member) as handle:
            return pd.read_csv(
                handle,
                dtype=str,
                encoding="latin1",
                keep_default_na=False,
                low_memory=False,
            )


def parse_mmddyyyy(series: pd.Series) -> pd.Series:
    return pd.to_datetime(series.replace("", pd.NA), format="%m/%d/%Y", errors="coerce")


def write_dataset(
    frame: pd.DataFrame,
    root: Path,
    partition_cols: list[str],
) -> None:
    root.mkdir(parents=True, exist_ok=True)
    table = pa.Table.from_pandas(frame, preserve_index=False)
    pq.write_to_dataset(
        table,
        root_path=str(root),
        partition_cols=partition_cols,
        compression=COMPRESSION,
        use_dictionary=True,
    )


def reset_outputs() -> None:
    for path in [
        FDA_OUT,
        VAERS_OUT / "reports",
        VAERS_OUT / "vaccines",
        VAERS_OUT / "symptoms_long",
        VAERS_OUT / "vaccine_symptoms",
        MANIFEST_OUT,
    ]:
        if path.is_dir():
            shutil.rmtree(path)
        elif path.exists():
            path.unlink()


def add_received_partitions(frame: pd.DataFrame, fallback_year: int) -> pd.DataFrame:
    received = frame["received_date"]
    frame["source_year"] = fallback_year
    frame["year"] = received.dt.year.fillna(fallback_year).astype("int16")
    frame["month"] = received.dt.month.fillna(0).astype("int8")
    return frame


def build_fda() -> dict[str, object]:
    frame = pd.read_csv(FDA_CSV)
    missing_columns = FDA_COMPARATOR_COLUMNS.difference(frame.columns)
    if missing_columns:
        raise ValueError(
            "The canonical openFDA CSV is stale. Run "
            "experimentacion/scripts/fetch_openfda_eliquis_comparator.R first. "
            f"Missing columns: {sorted(missing_columns)}"
        )
    frame["month_date"] = pd.to_datetime(frame["month"] + "-01", errors="coerce")
    frame = frame[
        frame["month_date"].dt.year.between(ANALYSIS_START_YEAR, ANALYSIS_END_YEAR)
    ].copy()
    frame = frame.drop(columns="month_date")

    frame["month_date"] = pd.to_datetime(frame["month"] + "-01", errors="coerce")
    frame["year"] = frame["month_date"].dt.year.astype("int16")
    frame["month_num"] = frame["month_date"].dt.month.astype("int8")

    write_dataset(
        frame,
        FDA_OUT,
        ["year", "month_num"],
    )

    return {
        "name": "openfda_eliquis_haemorrhage_monthly",
        "path": str(FDA_OUT.relative_to(ROOT)),
        "rows": int(len(frame)),
        "columns": list(frame.columns),
        "partitions": ["year", "month_num"],
        "compression": COMPRESSION,
    }


def clean_reports(raw: pd.DataFrame, source_year: int) -> pd.DataFrame:
    raw.columns = [snake_case(c) for c in raw.columns]
    raw = raw.rename(
        columns={
            "recvdate": "received_date",
            "rpt_date": "report_date",
            "datedied": "date_died",
            "order": "record_order",
        }
    )

    for col in DATE_COLUMNS.intersection(raw.columns):
        raw[col] = parse_mmddyyyy(raw[col])

    for col in NUMERIC_COLUMNS.intersection(raw.columns):
        raw[col] = pd.to_numeric(raw[col].replace("", pd.NA), errors="coerce")

    for col in SERIOUS_FLAGS:
        if col not in raw.columns:
            raw[col] = ""

    raw["serious"] = raw[SERIOUS_FLAGS].eq("Y").any(axis=1)
    raw = add_received_partitions(raw, source_year)
    raw["received_month"] = raw["received_date"].dt.strftime("%Y-%m").fillna(
        f"{source_year}-00"
    )
    return raw


def clean_vaccines(
    raw: pd.DataFrame, reports_index: pd.DataFrame, source_year: int
) -> pd.DataFrame:
    raw.columns = [snake_case(c) for c in raw.columns]
    raw = raw.rename(columns={"order": "record_order"})
    raw["source_year"] = source_year
    raw = raw.merge(
        reports_index,
        on=["source_year", "vaers_id", "record_order"],
        how="left",
    )
    raw = raw.rename(columns={"record_order": "vax_record_order"})
    raw["medicine_name"] = raw["vax_name"].replace("", "UNKNOWN")
    raw["medicine_slug"] = raw["medicine_name"].map(slugify)
    raw["vax_type"] = raw["vax_type"].replace("", "UNKNOWN")
    raw["vax_type_slug"] = raw["vax_type"].map(slugify)
    raw["manufacturer_slug"] = raw["vax_manu"].map(slugify)
    raw["year"] = raw["year"].fillna(raw["source_year"]).astype("int16")
    raw["month"] = raw["month"].fillna(0).astype("int8")
    raw["serious"] = raw["serious"].fillna(False).astype(bool)
    return raw


def clean_symptoms(
    raw: pd.DataFrame, reports_index: pd.DataFrame, source_year: int
) -> pd.DataFrame:
    raw.columns = [snake_case(c) for c in raw.columns]
    raw = raw.rename(columns={"order": "record_order"})
    raw["source_year"] = source_year
    pieces = []
    for slot in range(1, 6):
        symptom_col = f"symptom{slot}"
        version_col = f"symptomversion{slot}"
        if symptom_col not in raw.columns:
            continue
        piece = raw[
            ["source_year", "vaers_id", "record_order", symptom_col, version_col]
        ].copy()
        piece = piece.rename(
            columns={
                symptom_col: "symptom",
                version_col: "symptom_version",
            }
        )
        piece["symptom_slot"] = slot
        piece = piece[piece["symptom"] != ""]
        pieces.append(piece)

    if pieces:
        long = pd.concat(pieces, ignore_index=True)
    else:
        long = pd.DataFrame(
            columns=[
                "vaers_id",
                "record_order",
                "symptom",
                "symptom_version",
                "symptom_slot",
            ]
        )

    long = long.merge(
        reports_index,
        on=["source_year", "vaers_id", "record_order"],
        how="left",
    )
    long = long.rename(columns={"record_order": "symptoms_record_order"})
    long["symptom_slug"] = long["symptom"].map(slugify)
    long["year"] = long["year"].fillna(long["source_year"]).astype("int16")
    long["month"] = long["month"].fillna(0).astype("int8")
    long["serious"] = long["serious"].fillna(False).astype(bool)
    return long


def build_vaers() -> list[dict[str, object]]:
    outputs = {
        "reports": {
            "path": VAERS_OUT / "reports",
            "partitions": ["year", "month", "serious"],
            "rows": 0,
            "columns": None,
        },
        "vaccines": {
            "path": VAERS_OUT / "vaccines",
            "partitions": ["year", "month", "vax_type_slug", "medicine_slug"],
            "rows": 0,
            "columns": None,
        },
        "symptoms_long": {
            "path": VAERS_OUT / "symptoms_long",
            "partitions": ["year", "month"],
            "rows": 0,
            "columns": None,
        },
        "vaccine_symptoms": {
            "path": VAERS_OUT / "vaccine_symptoms",
            "partitions": ["year", "month", "vax_type_slug", "medicine_slug"],
            "rows": 0,
            "columns": None,
        },
    }

    zip_files = [
        path
        for path in sorted(VAERS_RAW.glob("*VAERSData.zip"))
        if in_analysis_range(int(path.name[:4]))
    ]
    for zip_path in zip_files:
        source_year = int(zip_path.name[:4])
        report_member = f"{source_year}VAERSDATA.csv"
        symptom_member = f"{source_year}VAERSSYMPTOMS.csv"
        vaccine_member = f"{source_year}VAERSVAX.csv"

        reports = clean_reports(read_zip_csv(zip_path, report_member), source_year)
        reports_index = reports[
            [
                "vaers_id",
                "record_order",
                "source_year",
                "year",
                "month",
                "received_month",
                "received_date",
                "state",
                "sex",
                "serious",
            ]
        ].copy()

        vaccines = clean_vaccines(
            read_zip_csv(zip_path, vaccine_member), reports_index, source_year
        )
        symptoms = clean_symptoms(
            read_zip_csv(zip_path, symptom_member), reports_index, source_year
        )

        vaccine_symptoms = vaccines[
            [
                "vaers_id",
                "source_year",
                "vax_record_order",
                "year",
                "month",
                "received_month",
                "received_date",
                "state",
                "sex",
                "serious",
                "vax_type",
                "vax_type_slug",
                "vax_manu",
                "manufacturer_slug",
                "vax_name",
                "medicine_name",
                "medicine_slug",
                "vax_dose_series",
            ]
        ].merge(
            symptoms[
                [
                    "vaers_id",
                    "source_year",
                    "symptoms_record_order",
                    "symptom",
                    "symptom_slug",
                    "symptom_version",
                    "symptom_slot",
                ]
            ],
            left_on=["source_year", "vaers_id", "vax_record_order"],
            right_on=["source_year", "vaers_id", "symptoms_record_order"],
            how="inner",
        )

        write_dataset(reports, outputs["reports"]["path"], outputs["reports"]["partitions"])
        write_dataset(
            vaccines,
            outputs["vaccines"]["path"],
            outputs["vaccines"]["partitions"],
        )
        write_dataset(
            symptoms,
            outputs["symptoms_long"]["path"],
            outputs["symptoms_long"]["partitions"],
        )
        write_dataset(
            vaccine_symptoms,
            outputs["vaccine_symptoms"]["path"],
            outputs["vaccine_symptoms"]["partitions"],
        )

        for name, frame in [
            ("reports", reports),
            ("vaccines", vaccines),
            ("symptoms_long", symptoms),
            ("vaccine_symptoms", vaccine_symptoms),
        ]:
            outputs[name]["rows"] += int(len(frame))
            if outputs[name]["columns"] is None:
                outputs[name]["columns"] = list(frame.columns)

        print(
            f"{source_year}: reports={len(reports):,} vaccines={len(vaccines):,} "
            f"symptoms={len(symptoms):,} vaccine_symptoms={len(vaccine_symptoms):,}"
        )

    return [
        {
            "name": name,
            "path": str(info["path"].relative_to(ROOT)),
            "rows": info["rows"],
            "columns": info["columns"],
            "partitions": info["partitions"],
            "compression": COMPRESSION,
        }
        for name, info in outputs.items()
    ]


def main() -> None:
    reset_outputs()
    manifest = {
        "compression": COMPRESSION,
        "datasets": [build_fda(), *build_vaers()],
    }
    MANIFEST_OUT.write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"Wrote {MANIFEST_OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
