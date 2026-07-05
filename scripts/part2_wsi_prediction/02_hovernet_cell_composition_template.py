#!/usr/bin/env python3
"""
Part 2 skeleton: summarize HoVer-Net nuclei composition.

This script represents the downstream analysis used after HoVer-Net inference.
It does not run HoVer-Net itself. Instead, it reads per-slide HoVer-Net JSON
outputs, counts predicted nuclei classes, calculates cell proportions, and can
join the summaries to study labels such as low/high PC1.

Expected HoVer-Net classes:
    0 = unlabeled
    1 = epithelial
    2 = lymphocyte
    3 = macrophage
    4 = neutrophil

HoVer-Net output formats differ between versions. The parser below accepts the
common structure in which each nucleus is represented by a dictionary with a
"type" field.
"""

from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path
from typing import Any

import pandas as pd


CELL_TYPES = {
    0: "unlabeled",
    1: "epithelial",
    2: "lymphocyte",
    3: "macrophage",
    4: "neutrophil",
}


def extract_nucleus_records(payload: Any) -> list[dict[str, Any]]:
    """
    Extract nucleus records from common HoVer-Net JSON structures.

    Supported examples:
      {"nuc": {"0": {"type": 1}, "1": {"type": 2}}}
      {"0": {"type": 1}, "1": {"type": 2}}
      [{"type": 1}, {"type": 2}]
    """
    if isinstance(payload, dict) and "nuc" in payload:
        payload = payload["nuc"]

    if isinstance(payload, dict):
        records = list(payload.values())
    elif isinstance(payload, list):
        records = payload
    else:
        raise ValueError("Unsupported HoVer-Net JSON structure.")

    return [record for record in records if isinstance(record, dict)]


def summarize_hovernet_json(json_path: Path) -> dict[str, float | int | str]:
    """Count nuclei classes and calculate per-slide cell proportions."""
    payload = json.loads(json_path.read_text(encoding="utf-8"))
    records = extract_nucleus_records(payload)

    type_ids: list[int] = []

    for record in records:
        raw_type = record.get("type", 0)

        try:
            type_id = int(raw_type)
        except (TypeError, ValueError):
            type_id = 0

        type_ids.append(type_id)

    counts = Counter(type_ids)
    total_cells = sum(counts.values())

    summary: dict[str, float | int | str] = {
        "slide_id": json_path.stem,
        "total_cells": total_cells,
    }

    for type_id, cell_name in CELL_TYPES.items():
        count = counts.get(type_id, 0)
        proportion = count / total_cells if total_cells else 0.0

        summary[f"{cell_name}_count"] = count
        summary[f"{cell_name}_prop"] = proportion

    immune_count = (
        counts.get(2, 0)
        + counts.get(3, 0)
        + counts.get(4, 0)
    )
    epithelial_count = counts.get(1, 0)

    summary["immune_count"] = immune_count
    summary["immune_prop"] = immune_count / total_cells if total_cells else 0.0
    summary["immune_epithelial_ratio"] = (
        immune_count / epithelial_count if epithelial_count else float("nan")
    )

    return summary


def summarize_directory(json_dir: Path) -> pd.DataFrame:
    """Summarize all HoVer-Net JSON files in a directory."""
    json_files = sorted(json_dir.glob("*.json"))

    if not json_files:
        raise FileNotFoundError(f"No JSON files found in {json_dir}")

    rows = [summarize_hovernet_json(path) for path in json_files]
    return pd.DataFrame(rows)


def join_study_labels(
    composition: pd.DataFrame,
    metadata_csv: Path,
    slide_id_col: str,
    group_col: str,
) -> pd.DataFrame:
    """
    Join per-slide cell composition to a study grouping variable.

    In the thesis, this can be used to compare cell composition between groups
    derived from the TE-expression PC1 labels.
    """
    metadata = pd.read_csv(metadata_csv)

    required = [slide_id_col, group_col]
    missing = [column for column in required if column not in metadata.columns]

    if missing:
        raise ValueError(f"Metadata are missing required columns: {missing}")

    metadata = metadata[required].drop_duplicates()

    return composition.merge(
        metadata,
        left_on="slide_id",
        right_on=slide_id_col,
        how="left",
        validate="one_to_one",
    )


def summarize_by_group(
    joined: pd.DataFrame,
    group_col: str,
) -> pd.DataFrame:
    """Calculate median and mean cell proportions for each study group."""
    proportion_columns = [
        column
        for column in joined.columns
        if column.endswith("_prop") or column == "immune_epithelial_ratio"
    ]

    summary = (
        joined.groupby(group_col)[proportion_columns]
        .agg(["median", "mean"])
        .reset_index()
    )

    # Flatten the pandas multi-index for easier CSV use.
    summary.columns = [
        "_".join(str(part) for part in column if part)
        if isinstance(column, tuple)
        else column
        for column in summary.columns
    ]

    return summary


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Summarize HoVer-Net cell composition by slide."
    )
    parser.add_argument(
        "--json-dir",
        required=True,
        type=Path,
        help="Directory containing per-slide HoVer-Net JSON files.",
    )
    parser.add_argument(
        "--output-csv",
        required=True,
        type=Path,
        help="Output path for the per-slide composition table.",
    )
    parser.add_argument(
        "--metadata-csv",
        type=Path,
        help="Optional slide metadata containing a study group label.",
    )
    parser.add_argument(
        "--slide-id-col",
        default="slide_id",
        help="Slide identifier column in the metadata table.",
    )
    parser.add_argument(
        "--group-col",
        default="pc1_group",
        help="Study-group column used for group summaries.",
    )
    parser.add_argument(
        "--group-summary-csv",
        type=Path,
        help="Optional output path for mean/median summaries by group.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    composition = summarize_directory(args.json_dir)

    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    composition.to_csv(args.output_csv, index=False)

    if args.metadata_csv is not None:
        joined = join_study_labels(
            composition=composition,
            metadata_csv=args.metadata_csv,
            slide_id_col=args.slide_id_col,
            group_col=args.group_col,
        )
        joined.to_csv(args.output_csv, index=False)

        if args.group_summary_csv is not None:
            group_summary = summarize_by_group(
                joined=joined,
                group_col=args.group_col,
            )
            args.group_summary_csv.parent.mkdir(
                parents=True,
                exist_ok=True,
            )
            group_summary.to_csv(args.group_summary_csv, index=False)


if __name__ == "__main__":
    main()
