#!/usr/bin/env python3
"""Filter artery features by prevalence and write a cleaned TSV in data/."""

from pathlib import Path
import argparse

import pandas as pd


ROOT_DIR = Path(__file__).resolve().parent
DATA_DIR = ROOT_DIR / "data"


def filter_columns(input_file: Path, output_file: Path, threshold: float = 0.05, epsilon: float = 0.001) -> None:
    df = pd.read_csv(input_file, sep="\t", index_col=0)

    # Keep lesion_volume. Filter only artery columns by prevalence:
    # prevalence = fraction of rows where lesion percentage exceeds epsilon.
    artery_cols = [c for c in df.columns if c != "lesion_volume"]
    non_zero_fraction = (df[artery_cols] > epsilon).mean()
    kept_artery_cols = non_zero_fraction[non_zero_fraction >= threshold].index.tolist()

    selected_cols = []
    if "lesion_volume" in df.columns:
        selected_cols.append("lesion_volume")
    selected_cols.extend(kept_artery_cols)

    out_df = df[selected_cols]
    output_file.parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(output_file, sep="\t")
    print(f"Wrote {out_df.shape[0]} rows, {out_df.shape[1]} columns to {output_file}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Create data/artery_cleaned.tsv from data/artery.tsv.")
    parser.add_argument("--input-file", type=Path, default=DATA_DIR / "artery.tsv")
    parser.add_argument("--output-file", type=Path, default=DATA_DIR / "artery_cleaned.tsv")
    parser.add_argument("--threshold", type=float, default=0.05)
    parser.add_argument("--epsilon", type=float, default=0.001)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    filter_columns(
        input_file=args.input_file,
        output_file=args.output_file,
        threshold=args.threshold,
        epsilon=args.epsilon,
    )


if __name__ == "__main__":
    main()
