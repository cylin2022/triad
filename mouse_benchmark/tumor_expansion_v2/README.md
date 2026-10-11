# Mouse Tumor Neoantigen Benchmark Expansion (v2)

This directory contains the v2 expanded benchmark for mouse tumor neoantigen immunogenicity prediction, following up on `tumor_feasibility/`.

## Directory Structure

- `PROTOCOL.md`: Immutable research protocol declared before training.
- `prepare.py`: Curates source cohorts, applies leakage detection, and writes locked tables.
- `train.py`: Computes 13 predeclared features with NetMHCpan 4.2d, fits models, and logs manifests.
- `assess.py`: Evaluates metrics across all holdout test cohorts with bootstrap intervals.
- `plot.py`: Generates PR and ROC curves for holdout test sets.
- `test_expansion.py`: 6 regression and unit tests verifying lock integrity, zero leakage, and reproducibility.
- `MODEL_CARD.md`: Detailed research model card with data lineage and performance tables.
- `REPORT.zh-TW.md`: Comprehensive traditional Chinese report with empirical analysis and decision rationale.
- `locked/`: Immutable split CSVs and `split_lock.json`.
- `sources/`: Downloaded publication supplements, including `JEM_20190179_DataS1.xlsx`.
- `run/`: Serialized models (`.joblib`), evaluation metrics, and curves (`*.png`).

## Reproduction Commands

```bash
# 1. Run unit tests
python3 mouse_benchmark/tumor_expansion_v2/test_expansion.py

# 2. View external evaluation metrics
cat mouse_benchmark/tumor_expansion_v2/run/external_metrics.csv | head -n 30
```
