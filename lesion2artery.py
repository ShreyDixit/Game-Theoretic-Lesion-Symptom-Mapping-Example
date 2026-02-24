#!/usr/bin/env python3
"""Build artery-level lesion features from lesion masks and an arterial atlas."""

from pathlib import Path
import argparse
import glob

import nibabel as nib
import numpy as np
import pandas as pd


ROOT_DIR = Path(__file__).resolve().parent
DATA_DIR = ROOT_DIR / "data"
NIFTI_DIR = DATA_DIR / "NIfTI"


def load_atlas_labels(labels_path: Path) -> dict[int, str]:
    labels: dict[int, str] = {}
    with labels_path.open("r", encoding="utf-8") as f:
        for line in f:
            parts = [p.strip() for p in line.strip().split("|")]
            if len(parts) < 2:
                continue
            idx = int(parts[0])
            labels[idx] = f"{idx}_{parts[1]}"
    return labels


def compute_row_for_lesion(
    lesion_path: Path, atlas_data: np.ndarray, atlas_counts: dict[int, int], atlas_labels: dict[int, str]
) -> dict[str, float]:
    lesion_data = nib.load(str(lesion_path)).get_fdata()
    lesion_nonzero = lesion_data != 0

    subject = lesion_path.name.split("_")[0]
    if subject.startswith("bwsr"):
        subject = subject[4:]

    row: dict[str, float] = {
        "participant_id": subject,
        "lesion_volume": float(np.count_nonzero(lesion_nonzero)),
    }
    for region in range(1, int(atlas_data.max()) + 1):
        denom = atlas_counts[region] if atlas_counts[region] > 0 else 1
        numer = int(np.count_nonzero((atlas_data == region) & lesion_nonzero))
        row[atlas_labels.get(region, str(region))] = numer / denom
    return row


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Create data/artery.tsv from lesion masks and arterial atlas.")
    parser.add_argument("--nifti-dir", type=Path, default=NIFTI_DIR)
    parser.add_argument("--atlas-nifti", type=Path, default=NIFTI_DIR / "ArterialAtlas136.nii.gz")
    parser.add_argument("--atlas-labels", type=Path, default=NIFTI_DIR / "ArterialAtlas136.txt")
    parser.add_argument("--lesion-pattern", default="bwsr*_lesion.nii.gz")
    parser.add_argument("--output", type=Path, default=DATA_DIR / "artery.tsv")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)

    atlas_data = nib.load(str(args.atlas_nifti)).get_fdata().astype(int)
    atlas_labels = load_atlas_labels(args.atlas_labels)
    atlas_counts = {r: int(np.count_nonzero(atlas_data == r)) for r in range(1, int(atlas_data.max()) + 1)}

    lesion_files = sorted(glob.glob(str(args.nifti_dir / args.lesion_pattern)))
    print(f"{len(lesion_files)} files match {args.lesion_pattern} in {args.nifti_dir}")

    rows = [
        compute_row_for_lesion(Path(lesion_path), atlas_data, atlas_counts, atlas_labels)
        for lesion_path in lesion_files
    ]
    out_df = pd.DataFrame(rows)
    out_df.to_csv(args.output, sep="\t", index=False)
    print(f"Wrote {len(out_df)} rows to {args.output}")


if __name__ == "__main__":
    main()
