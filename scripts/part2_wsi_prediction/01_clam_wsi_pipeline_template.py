#!/usr/bin/env python3
"""
Part 2 skeleton: CLAM whole-slide image workflow.

This file is intentionally a readable workflow template rather than a drop-in
replacement for the external CLAM repository. It shows the sequence used in the
thesis:

1. prepare slide-level labels,
2. segment tissue and extract image patches,
3. extract patch embeddings with UNI2-h,
4. create patient-level cross-validation splits,
5. train/evaluate CLAM,
6. summarize fold-level predictions.

External CLAM and UNI2-h code, model weights, TCGA slides, and feature tensors
are not distributed with this repository.

Typical usage is to edit the configuration paths below and first run with
--dry-run to inspect the commands that would be executed.
"""

from __future__ import annotations

import argparse
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import pandas as pd


@dataclass
class WorkflowConfig:
    """Local paths and project-specific settings."""

    clam_root: Path
    slide_dir: Path
    patch_dir: Path
    feature_dir: Path
    metadata_csv: Path
    dataset_csv: Path
    split_dir: Path
    results_dir: Path

    slide_id_col: str = "slide_id"
    patient_id_col: str = "patient_id"
    label_col: str = "label"

    patch_size: int = 256
    slide_extension: str = ".svs"
    feature_model_name: str = "uni_v1"  # local CLAM alias used for UNI2-h
    feature_batch_size: int = 1536

    n_splits: int = 5
    random_seed: int = 42


def run_command(command: list[str], dry_run: bool = True) -> None:
    """Print a command and optionally execute it."""
    printable = " ".join(str(item) for item in command)
    print(f"\n$ {printable}")

    if not dry_run:
        subprocess.run(command, check=True)


def validate_metadata(
    metadata: pd.DataFrame,
    slide_id_col: str,
    patient_id_col: str,
    label_col: str,
) -> pd.DataFrame:
    """
    Keep the minimum columns needed by the WSI workflow and verify that labels
    are defined at slide level while splits can be grouped at patient level.
    """
    required = [slide_id_col, patient_id_col, label_col]
    missing = [column for column in required if column not in metadata.columns]

    if missing:
        raise ValueError(f"Metadata are missing required columns: {missing}")

    dataset = metadata[required].copy()
    dataset = dataset.dropna(subset=required)
    dataset[slide_id_col] = dataset[slide_id_col].astype(str)
    dataset[patient_id_col] = dataset[patient_id_col].astype(str)

    if dataset[slide_id_col].duplicated().any():
        duplicates = dataset.loc[
            dataset[slide_id_col].duplicated(), slide_id_col
        ].tolist()
        raise ValueError(f"Duplicate slide IDs found: {duplicates[:10]}")

    # Patient-level splitting prevents slides from the same patient from
    # appearing in both training and test sets.
    patient_label_counts = (
        dataset.groupby(patient_id_col)[label_col].nunique(dropna=True)
    )
    inconsistent = patient_label_counts[patient_label_counts > 1]

    if not inconsistent.empty:
        raise ValueError(
            "Some patients have conflicting slide labels: "
            f"{inconsistent.index.tolist()[:10]}"
        )

    return dataset


def create_patient_level_splits(
    dataset: pd.DataFrame,
    output_dir: Path,
    patient_id_col: str,
    label_col: str,
    n_splits: int,
    random_seed: int,
) -> None:
    """
    Create stratified patient-level folds.

    CLAM expects train/validation/test split CSVs. The exact file format varies
    across CLAM forks, so this function writes a transparent long-form table
    that can be converted to the format expected by a local CLAM installation.
    """
    try:
        from sklearn.model_selection import StratifiedKFold
    except ImportError as exc:
        raise ImportError(
            "scikit-learn is required to create patient-level folds."
        ) from exc

    output_dir.mkdir(parents=True, exist_ok=True)

    patient_table = (
        dataset[[patient_id_col, label_col]]
        .drop_duplicates()
        .reset_index(drop=True)
    )

    splitter = StratifiedKFold(
        n_splits=n_splits,
        shuffle=True,
        random_state=random_seed,
    )

    fold_rows: list[pd.DataFrame] = []

    for fold, (train_index, test_index) in enumerate(
        splitter.split(patient_table[patient_id_col], patient_table[label_col])
    ):
        train_patients = set(patient_table.iloc[train_index][patient_id_col])
        test_patients = set(patient_table.iloc[test_index][patient_id_col])

        fold_data = dataset.copy()
        fold_data["fold"] = fold
        fold_data["split"] = fold_data[patient_id_col].map(
            lambda patient: "train" if patient in train_patients else "test"
        )

        # A validation subset can be created from the training patients inside
        # the local CLAM implementation or by adding a second grouped split.
        fold_rows.append(fold_data)

    split_table = pd.concat(fold_rows, ignore_index=True)
    split_table.to_csv(output_dir / "patient_level_folds.csv", index=False)


def build_patch_command(config: WorkflowConfig) -> list[str]:
    """Build the CLAM tissue segmentation and patch-extraction command."""
    return [
        "python",
        str(config.clam_root / "create_patches_fp.py"),
        "--source",
        str(config.slide_dir),
        "--save_dir",
        str(config.patch_dir),
        "--patch_size",
        str(config.patch_size),
        "--seg",
        "--patch",
        "--stitch",
    ]


def build_feature_command(config: WorkflowConfig) -> list[str]:
    """Build the CLAM/UNI2-h feature-extraction command."""
    h5_list = config.patch_dir / "h5_list.csv"

    return [
        "python",
        str(config.clam_root / "extract_features_fp.py"),
        "--data_h5_dir",
        str(config.patch_dir),
        "--data_slide_dir",
        str(config.slide_dir),
        "--csv_path",
        str(h5_list),
        "--feat_dir",
        str(config.feature_dir),
        "--batch_size",
        str(config.feature_batch_size),
        "--slide_ext",
        config.slide_extension,
        "--model_name",
        config.feature_model_name,
    ]


def write_h5_list(patch_dir: Path) -> Path:
    """Create the list of patch files required by CLAM feature extraction."""
    patches_dir = patch_dir / "patches"
    h5_files = sorted(patches_dir.glob("*.h5"))

    if not h5_files:
        raise FileNotFoundError(f"No .h5 patch files found in {patches_dir}")

    output_path = patch_dir / "h5_list.csv"
    output_path.write_text(
        "\n".join(path.stem for path in h5_files) + "\n",
        encoding="utf-8",
    )
    return output_path


def build_clam_training_command(
    config: WorkflowConfig,
    task_name: str,
) -> list[str]:
    """
    Build a generic CLAM training command.

    Argument names differ slightly between CLAM forks. Adapt this command to the
    local version while keeping the same inputs, patient-level folds, and task
    labels documented here.
    """
    return [
        "python",
        str(config.clam_root / "main.py"),
        "--data_root_dir",
        str(config.feature_dir),
        "--csv_path",
        str(config.dataset_csv),
        "--split_dir",
        str(config.split_dir),
        "--results_dir",
        str(config.results_dir),
        "--task",
        task_name,
        "--model_type",
        "clam_sb",
        "--k",
        str(config.n_splits),
        "--seed",
        str(config.random_seed),
    ]


def summarize_prediction_files(
    prediction_files: Iterable[Path],
    output_csv: Path,
) -> pd.DataFrame:
    """
    Combine fold-level prediction tables.

    Each input table is expected to include at least:
      - slide_id
      - true_label
      - predicted_label

    Optional probability columns are retained.
    """
    frames: list[pd.DataFrame] = []

    for fold, prediction_file in enumerate(prediction_files):
        frame = pd.read_csv(prediction_file)
        frame["fold"] = fold
        frames.append(frame)

    if not frames:
        raise ValueError("No prediction files were provided.")

    combined = pd.concat(frames, ignore_index=True)
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    combined.to_csv(output_csv, index=False)
    return combined


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Readable CLAM workflow skeleton for the thesis repository."
    )
    parser.add_argument(
        "--run",
        action="store_true",
        help="Execute external commands. Default behavior is dry-run only.",
    )
    parser.add_argument(
        "--task-name",
        default="replace_with_task_name",
        help="CLAM task name used by the local installation.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    # Replace these placeholders with local paths before use.
    config = WorkflowConfig(
        clam_root=Path("<path_to_clam_repository>"),
        slide_dir=Path("<path_to_svs_slides>"),
        patch_dir=Path("<path_to_patch_outputs>"),
        feature_dir=Path("<path_to_uni2_features>"),
        metadata_csv=Path("<path_to_slide_metadata.csv>"),
        dataset_csv=Path("<path_to_clam_dataset.csv>"),
        split_dir=Path("<path_to_split_directory>"),
        results_dir=Path("<path_to_clam_results>"),
    )

    dry_run = not args.run

    # 1. Prepare the slide-level label table.
    metadata = pd.read_csv(config.metadata_csv)
    dataset = validate_metadata(
        metadata=metadata,
        slide_id_col=config.slide_id_col,
        patient_id_col=config.patient_id_col,
        label_col=config.label_col,
    )
    config.dataset_csv.parent.mkdir(parents=True, exist_ok=True)
    dataset.to_csv(config.dataset_csv, index=False)

    # 2. Create patient-level cross-validation folds.
    create_patient_level_splits(
        dataset=dataset,
        output_dir=config.split_dir,
        patient_id_col=config.patient_id_col,
        label_col=config.label_col,
        n_splits=config.n_splits,
        random_seed=config.random_seed,
    )

    # 3. Segment tissue and extract patches.
    run_command(build_patch_command(config), dry_run=dry_run)

    # 4. After patch extraction completes, create h5_list.csv:
    # write_h5_list(config.patch_dir)

    # 5. Extract UNI2-h patch embeddings through the local CLAM fork.
    run_command(build_feature_command(config), dry_run=dry_run)

    # 6. Train and evaluate CLAM.
    run_command(
        build_clam_training_command(config, task_name=args.task_name),
        dry_run=dry_run,
    )


if __name__ == "__main__":
    main()
