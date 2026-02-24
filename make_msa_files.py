#!/usr/bin/env python3
import argparse
from pathlib import Path
import glob

import nibabel as nib
import numpy as np
import pandas as pd

ROOT_DIR = Path(__file__).resolve().parent
DATA_DIR = ROOT_DIR / "data"
NIFTI_DIR = DATA_DIR / "NIfTI"


def load_atlas_labels(atlas_labels_path: Path) -> dict[int, str]:
    labels = {}
    with atlas_labels_path.open("r", encoding="utf-8") as f:
        for line in f:
            parts = [p.strip() for p in line.strip().split("|")]
            if len(parts) < 2:
                continue
            idx = int(parts[0])
            short = parts[1]
            labels[idx] = f"{idx}_{short}"
    return labels


def compute_voxel_counts(atlas_nifti_path: Path, atlas_labels: dict[int, str]) -> pd.DataFrame:
    atlas = nib.load(str(atlas_nifti_path)).get_fdata().astype(int)
    max_label = int(atlas.max())
    row = {}
    for region in range(1, max_label + 1):
        row[atlas_labels.get(region, str(region))] = int(np.count_nonzero(atlas == region))
    return pd.DataFrame([row])


def lesion_data_from_artery_tsv(artery_tsv_path: Path) -> pd.DataFrame:
    df = pd.read_csv(artery_tsv_path, sep="\t")
    if "participant_id" not in df.columns:
        raise ValueError("artery.tsv must include a participant_id column")
    roi_cols = [c for c in df.columns if c not in {"participant_id", "lesion_volume"}]
    out = df[["participant_id"] + roi_cols].copy()
    out[roi_cols] = out[roi_cols].astype(float) * 100.0
    return out


def lesion_data_from_nifti(lesion_dir: Path, atlas_nifti_path: Path, atlas_labels: dict[int, str]) -> pd.DataFrame:
    atlas = nib.load(str(atlas_nifti_path)).get_fdata().astype(int)
    max_label = int(atlas.max())
    atlas_counts = {r: int(np.count_nonzero(atlas == r)) for r in range(1, max_label + 1)}
    lesion_files = sorted(glob.glob(str(lesion_dir / "bwsr*_lesion.nii*")))
    rows = []
    for lesion_path in lesion_files:
        lesion_img = nib.load(lesion_path).get_fdata()
        lesion_bin = lesion_img != 0
        subj = Path(lesion_path).name.split("_")[0]
        subj = subj[4:] if subj.startswith("bwsr") else subj
        row = {"participant_id": subj}
        for region in range(1, max_label + 1):
            mask = atlas == region
            numer = int(np.count_nonzero(mask & lesion_bin))
            denom = atlas_counts[region] if atlas_counts[region] > 0 else 1
            row[atlas_labels.get(region, str(region))] = (numer / denom) * 100.0
        rows.append(row)
    return pd.DataFrame(rows)


def make_score_df(participants_tsv_path: Path, participant_ids: set[str]) -> pd.DataFrame:
    df = pd.read_csv(participants_tsv_path, sep="\t")
    if "participant_id" not in df.columns or "nihss" not in df.columns:
        raise ValueError("participants.tsv must include participant_id and nihss columns")
    out = df[["participant_id", "nihss"]].copy()
    out["nihss"] = pd.to_numeric(out["nihss"], errors="coerce")
    out = out[out["participant_id"].isin(participant_ids)]
    out = out.dropna(subset=["nihss"]).reset_index(drop=True)
    return out


def write_df(df: pd.DataFrame, out_path: Path, output_format: str):
    if output_format == "xlsx":
        df.to_excel(out_path, index=False)
    else:
        df.to_csv(out_path, index=False)


def format_voxel_file(voxel_df: pd.DataFrame, roi_columns: list[str]) -> pd.DataFrame:
    """Return voxels in two-column long format aligned to lesion ROI column order."""
    rows = []
    for roi_col in roi_columns:
        rows.append([roi_col, int(round(float(voxel_df.iloc[0][roi_col])))])
    return pd.DataFrame(rows)


def _connected_components_from_correlation(corr: pd.DataFrame, threshold: float) -> list[list[str]]:
    cols = list(corr.columns)
    neighbors = {c: set() for c in cols}
    for i, c1 in enumerate(cols):
        for j in range(i + 1, len(cols)):
            c2 = cols[j]
            if float(corr.iloc[i, j]) >= threshold:
                neighbors[c1].add(c2)
                neighbors[c2].add(c1)
    visited = set()
    groups: list[list[str]] = []
    for c in cols:
        if c in visited:
            continue
        stack = [c]
        comp = []
        visited.add(c)
        while stack:
            cur = stack.pop()
            comp.append(cur)
            for nb in neighbors[cur]:
                if nb not in visited:
                    visited.add(nb)
                    stack.append(nb)
        groups.append(sorted(comp))
    return groups


def merge_correlated_rois(
    lesion_df: pd.DataFrame, voxel_df: pd.DataFrame, threshold: float
) -> tuple[pd.DataFrame, pd.DataFrame, list[dict[str, str]]]:
    corr = lesion_df.corr(method="pearson")
    groups = _connected_components_from_correlation(corr, threshold)

    new_lesion = pd.DataFrame(index=lesion_df.index)
    new_voxels = {}
    mapping_rows: list[dict[str, str]] = []
    merge_idx = 1

    for g in groups:
        if len(g) == 1:
            c = g[0]
            new_lesion[c] = lesion_df[c]
            new_voxels[c] = float(voxel_df.iloc[0][c])
            continue

        merged_name = f"MERGED_{merge_idx:03d}"
        merge_idx += 1
        vox = voxel_df.iloc[0][g].astype(float)
        denom = float(vox.sum())
        # Weighted recomputation of lesion percent for combined ROI.
        numer = (lesion_df[g].astype(float) / 100.0).mul(vox, axis=1).sum(axis=1)
        new_lesion[merged_name] = (numer / denom) * 100.0 if denom > 0 else 0.0
        new_voxels[merged_name] = denom
        mapping_rows.append(
            {
                "merged_roi": merged_name,
                "num_members": str(len(g)),
                "members": ",".join(g),
            }
        )

    return new_lesion, pd.DataFrame([new_voxels]), mapping_rows


def parse_args():
    p = argparse.ArgumentParser(description="Generate lesion_data, NIHSS score, and ROI voxel files for MSA-App.")
    p.add_argument("--mode", choices=["artery", "nifti"], default="artery")
    p.add_argument("--artery-tsv", default=str(DATA_DIR / "artery.tsv"))
    p.add_argument("--participants-tsv", default=str(DATA_DIR / "participants.tsv"))
    p.add_argument("--lesion-dir", default=str(NIFTI_DIR))
    p.add_argument("--atlas-nifti", default=str(NIFTI_DIR / "ArterialAtlas136.nii.gz"))
    p.add_argument("--atlas-labels", default=str(NIFTI_DIR / "ArterialAtlas136.txt"))
    p.add_argument("--output-dir", default=str(DATA_DIR))
    p.add_argument("--output-format", choices=["csv", "xlsx"], default="csv")
    p.add_argument("--lesion-out", default="roi_data")
    p.add_argument("--scores-out", default="nihss_scores")
    p.add_argument("--voxels-out", default="num_voxels")
    p.add_argument("--merge-corr-threshold", type=float, default=None)
    p.add_argument("--merge-map-out", default="roi_merge_map")
    return p.parse_args()


def main():
    args = parse_args()
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    suffix = ".xlsx" if args.output_format == "xlsx" else ".csv"

    atlas_labels = load_atlas_labels(Path(args.atlas_labels))
    if args.mode == "artery":
        lesion_df = lesion_data_from_artery_tsv(Path(args.artery_tsv))
    else:
        lesion_df = lesion_data_from_nifti(Path(args.lesion_dir), Path(args.atlas_nifti), atlas_labels)

    score_df = make_score_df(Path(args.participants_tsv), set(lesion_df["participant_id"].tolist()))
    # Keep only participants with NIHSS and enforce identical row order between
    # lesion and score files, then drop IDs because downstream expects same order.
    keep_ids = score_df["participant_id"].tolist()
    lesion_df = lesion_df.set_index("participant_id").loc[keep_ids].reset_index()
    score_df = score_df.set_index("participant_id").loc[keep_ids].reset_index()
    lesion_df = lesion_df.drop(columns=["participant_id"])
    score_df = score_df.drop(columns=["participant_id"])
    voxel_df = compute_voxel_counts(Path(args.atlas_nifti), atlas_labels)

    merge_rows: list[dict[str, str]] = []
    if args.merge_corr_threshold is not None:
        lesion_df, voxel_df, merge_rows = merge_correlated_rois(
            lesion_df, voxel_df, args.merge_corr_threshold
        )

    lesion_path = output_dir / f"{args.lesion_out}{suffix}"
    score_path = output_dir / f"{args.scores_out}{suffix}"
    voxel_path = output_dir / f"{args.voxels_out}{suffix}"
    merge_map_path = output_dir / f"{args.merge_map_out}{suffix}"

    write_df(lesion_df, lesion_path, args.output_format)
    write_df(score_df, score_path, args.output_format)
    voxel_out_df = format_voxel_file(voxel_df, list(lesion_df.columns))
    if args.output_format == "xlsx":
        voxel_out_df.to_excel(voxel_path, index=False, header=False)
    else:
        voxel_out_df.to_csv(voxel_path, index=False, header=False)
    if args.merge_corr_threshold is not None:
        write_df(
            pd.DataFrame(merge_rows, columns=["merged_roi", "num_members", "members"]),
            merge_map_path,
            args.output_format,
        )

    print(f"Wrote lesion data: {lesion_path}")
    print(f"Wrote score data:  {score_path}")
    print(f"Wrote voxel data:  {voxel_path}")
    if args.merge_corr_threshold is not None:
        print(f"Wrote merge map:   {merge_map_path}")
    print(f"Lesion rows: {len(lesion_df)}")
    print(f"Score rows:  {len(score_df)}")
    print(f"Voxel rows:  {len(voxel_out_df)}")
    print(f"ROI columns: {len(lesion_df.columns)}")


if __name__ == "__main__":
    main()
