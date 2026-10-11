# Tumor neoantigen training feasibility v1

User-selected scope: tumor neoantigens. Work is isolated from the production app.
Read [the findings](REPORT.zh-TW.md), [model card](MODEL_CARD.md) and
[frozen protocol](PROTOCOL.md). No new animal experiments were performed.

## What was executed

1. Read original papers, supplementary tables and full IEDB `tcell_export` records.
2. Correctly distinguish mutant assays from WT controls, and vaccination from
   untreated/tumor-only or post-challenge conditions.
3. Freeze study and tumor-model-separated train/validation/test assignments with
   mutant/WT near-sequence and source-protein overlap purging.
4. Fit fixed logistic and gradient boosting models using training data only.
5. Evaluate validation and the fixed test, without tuning or selecting on test.

The split was locked at the timestamp in `locked/split_lock.json` before the fitted
checkpoint timestamp in `run/model_manifest.json`. External evaluation timestamps
are in `run/evaluation_started.json` and `run/evaluation_completed.json`.
Some records appeared in the earlier broad benchmark. This is a holdout from NEW
model fitting, not a previously unexamined or prospective cohort.

## Reproduction

Python dependencies: numpy, pandas, scikit-learn, scipy, joblib, openpyxl,
matplotlib; the existing local netMHCpan 4.2d installation is required. Versions
are recorded in `run/artifact_audit.json`. Run from repository root:

```bash
python -m unittest discover -s mouse_benchmark/tumor_feasibility -p 'test_*.py' -v
python mouse_benchmark/tumor_feasibility/audit.py
python mouse_benchmark/tumor_feasibility/plot.py
```

The original analysis commands were:

```bash
python mouse_benchmark/tumor_feasibility/prepare.py
python mouse_benchmark/tumor_feasibility/train.py
python mouse_benchmark/tumor_feasibility/assess.py
```

Those three commands intentionally refuse to overwrite a completed lock, fitted
checkpoint or opened test evaluation. For an exact rerun, copy this experiment
directory to a **new sibling directory under `mouse_benchmark/`**, excluding only
`locked/`, `run/` and `__pycache__/`; preserve sources and code. Execute the same
three scripts there. Do not adjust hyperparameters against the existing test.
The parent IEDB snapshot must remain accessible. Numerical results should match;
timestamps and absolute command paths will differ.

## Source provenance

`sources/details_manifest.json` records exact URLs and checksums for 313 detailed
assays across eight candidate studies. Only four fully scoped studies are used.
`sources/manifest.json` records paper retrievals (including failures);
`sources/supplement_archive_manifest.json` records successful official supplementary
archives and extracted member checksums. `sources/smith_supp_metadata.json` contains
the publisher's file metadata and citation for the Smith spreadsheet. Some initial
download attempts failed; the subsequent successful archives are the files used.
Full-text articles and supplementary files are research source material, not
models to redistribute as part of the production service.

## Outputs

| Artifact | Purpose |
|---|---|
|`locked/assay_curation.csv`|Every assay in the four scoped studies and its inclusion/exclusion reason|
|`locked/before_purging.csv`|Context-cleaned pairs before cross-split isolation|
|`locked/purged.csv`|Development pairs removed for overlap, with the matching holdout|
|`locked/train.csv`, `validation.csv`, `test.csv`|Frozen data and provenance|
|`locked/split_lock.json`|File hashes, counts and zero-overlap checks|
|`run/*.joblib`|Research-only serialized models|
|`run/model_manifest.json`|Training settings, code/weights/binary hashes|
|`run/features.csv`|The 13 descriptive features, no outcome/context IDs|
|`run/scored_cohort.csv`|Mutant/WT binding predictions and all model scores|
|`run/external_metrics.csv`|Validation/test, per allele and per study metrics|
|`run/development_metrics.csv`|Purged leave-one-tumor-model-out diagnostics|
|`run/uncertainty.json`|Exploratory within-study grouped bootstrap intervals|
|`run/test_comparison.png`, `.svg`|Exportable comparison figure|
|`run/artifact_audit.json`|Integrity, fitted-scaler and serialized-prediction checks|
|`study_review.csv`|Reviewed source decisions and context restrictions|

The test set does not authorize a production replacement. Results did not establish
an advantage for the primary new model; this run is retained as a negative
feasibility finding rather than tuned further against the same test.
