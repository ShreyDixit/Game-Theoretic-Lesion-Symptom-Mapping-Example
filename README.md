# Stroke MSA Data Preparation (Example Repository)

This repository documents and reproduces a practical pipeline for preparing stroke lesion data for **Multi-perturbation Shapley-value Analysis (MSA)**.

In plain terms, this repo takes lesion masks + patient metadata and converts them into the three tables that MSA-App expects:
- `roi_data.csv`
- `nihss_scores.csv`
- `num_voxels.csv`

To run lesion-symptom mapping itself, download the latest MSA-App release from:
- https://github.com/ShreyDixit/MSA-App/tree/master

This repository is the data-preparation example that feeds that app.

## What This Repository Produces

All generated files are written to `data/`:
- `data/artery.tsv`
- `data/artery_cleaned.tsv`
- `data/roi_data.csv`
- `data/nihss_scores.csv`
- `data/num_voxels.csv`
- optional: `data/roi_merge_map.csv`

Important:
- These generated files are already present in this repository.
- You can still rerun scripts to regenerate them.

## Data and Folder Layout

Expected layout:

```text
.
├─ lesion2artery.py
├─ clean_artery_tsv.py
├─ make_msa_files.py
└─ data/
   ├─ participants.tsv
   ├─ artery.tsv
   ├─ artery_cleaned.tsv
   ├─ roi_data.csv
   ├─ nihss_scores.csv
   ├─ num_voxels.csv
   └─ NIfTI/
      ├─ ArterialAtlas136.nii.gz
      ├─ ArterialAtlas136.txt
      └─ bwsrsub-*_lesion.nii.gz
```

## End-to-End Pipeline

Run in this order:

1. Build arterial-territory lesion features:

```bash
python lesion2artery.py
```

2. Optionally filter sparse arterial features:

```bash
python clean_artery_tsv.py
```

3. Build MSA input tables:

```bash
python make_msa_files.py --mode artery --output-format csv --output-dir data
```

After step 3, pass these files to MSA-App:
- `data/roi_data.csv`
- `data/nihss_scores.csv`
- `data/num_voxels.csv`

## File Definitions

### `data/artery.tsv`
- One row per participant.
- Columns:
  - `participant_id`
  - `lesion_volume` (count of non-zero lesion voxels)
  - arterial-territory ROI columns (`0..1`, fraction of each territory affected)

What “arterial-territory ROI” means:
- Regions are based on vascular supply territories (ACA/MCA/PCA branches, perforator territories), not fine anatomical parcels.

### `data/artery_cleaned.tsv`
- Same rows as `artery.tsv`.
- Optional filtered version where rare artery columns can be removed using prevalence rules.

### `data/roi_data.csv`
- Final lesion feature matrix for MSA-App.
- Rows: participants with valid NIHSS.
- Columns: ROI lesion percentages in `0..100`.
- No participant ID column (row alignment is used across MSA input files).

### `data/nihss_scores.csv`
- Final target/output table for MSA-App.
- Single column: `nihss`.
- Row order is aligned with `data/roi_data.csv`.

### `data/num_voxels.csv`
- Two-column ROI voxel table:
  - column 1: ROI name (must match ROI columns in `roi_data.csv`)
  - column 2: voxel count for that ROI

### `data/roi_merge_map.csv` (optional)
- Created only when ROI-correlation merge is enabled.
- Documents which original ROIs were merged.

## Script Documentation

## `lesion2artery.py`

Purpose:
- Converts lesion masks to arterial-territory feature values.

Default arguments:
- `--nifti-dir data/NIfTI`
- `--atlas-nifti data/NIfTI/ArterialAtlas136.nii.gz`
- `--atlas-labels data/NIfTI/ArterialAtlas136.txt`
- `--lesion-pattern bwsr*_lesion.nii.gz`
- `--output data/artery.tsv`

Example:

```bash
python lesion2artery.py --output data/artery.tsv
```

## `clean_artery_tsv.py`

Purpose:
- Optional feature filtering for `artery.tsv`.
- Keeps `lesion_volume`.
- Filters artery columns by prevalence:
  - a value is considered present if it is greater than `epsilon`
  - a column is kept if at least `threshold` fraction of rows is present

Default arguments:
- `--input-file data/artery.tsv`
- `--output-file data/artery_cleaned.tsv`
- `--threshold 0.05`
- `--epsilon 0.001`

Example:

```bash
python clean_artery_tsv.py --threshold 0.05 --epsilon 0.001
```

## `make_msa_files.py`

Purpose:
- Produces final MSA-App input files.
- Aligns rows between lesion and NIHSS tables.

Core behavior:
- Keeps only participants with numeric NIHSS values.
- Enforces identical row order for lesion and score tables.
- Removes participant IDs from final MSA input files by design.

Default arguments:
- `--mode artery`
- `--artery-tsv data/artery.tsv`
- `--participants-tsv data/participants.tsv`
- `--lesion-dir data/NIfTI`
- `--atlas-nifti data/NIfTI/ArterialAtlas136.nii.gz`
- `--atlas-labels data/NIfTI/ArterialAtlas136.txt`
- `--output-dir data`
- `--output-format csv`
- `--lesion-out roi_data`
- `--scores-out nihss_scores`
- `--voxels-out num_voxels`

Optional ROI correlation merge:
- `--merge-corr-threshold <float>` (example: `0.95`)
- `--merge-map-out roi_merge_map`

Examples:

```bash
# Standard export
python make_msa_files.py --mode artery --output-format csv --output-dir data

# Export with correlated ROI merging
python make_msa_files.py --mode artery --merge-corr-threshold 0.95 --output-dir data
```

## Example MSA Interpretation (Included Result Files)

Included MSA outputs:
- `results_iterative_20260224-152628.json`
- `shapley_values_iterative_20260224-152628.csv`

Interpretation:
- In this run, MSA identifies a subset of arterial-territory regions with higher Shapley contribution values.
- This indicates these regions contribute more strongly to deficit variation measured by NIHSS under the chosen analysis setup.
- Treat this as ranked contribution evidence within the model configuration, not standalone causal proof.

## Using a Different Parcellation / Atlas

If you want a different atlas/parcellation:
1. Download SOOP data from OpenNeuro:  
   https://openneuro.org/datasets/ds004889/versions/1.1.2
2. Follow preprocessing from the original SOOP demo repository:  
   https://github.com/neurolabusc/StrokeOutcomeOptimizationProjectDemo/tree/main
3. Make sure lesion masks and the new atlas are in the same spatial template/grid before feature extraction.

## Attribution

- Some code and workflow ideas are adapted from:  
  https://github.com/neurolabusc/StrokeOutcomeOptimizationProjectDemo/tree/main

- Dataset/publication source:  
  https://www.nature.com/articles/s41597-023-02457-9  
  *A large public dataset of annotated clinical MRIs and metadata of patients with acute stroke*

## Support

If something is unclear, or if code/data behavior does not match what you expect, please email:  
`dixit@cbs.mpg.de`
