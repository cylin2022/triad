# Mouse immunogenicity retrospective benchmark

Independent development directory on `feature/mouse-immunogenicity`. Existing
`benchmark/`, production scoring code and UI are unchanged. No new model was
trained, no entropy/free-energy calculations were performed, and no claim of
prospective validation is made.

Read [the Chinese findings and training decision](REPORT.zh-TW.md) and the
[frozen evaluation protocol](PROTOCOL.md). Download provenance comes from the
[IEDB Query API](https://query-api.iedb.org/docs/swagger/) and the
[official API explanation](https://discuss.iedb.org/t/immune-epitope-database-query-api-iq-api/154).

## Reproduce using this snapshot

Run from the repository root. Requires Python 3.10+, numpy, pandas, scikit-learn,
matplotlib, and the locally installed licensed netMHCpan 4.2 package. See
`results/environment.json` for the exact versions used. No web server is required.
`NETMHC_HOME` can select another installation; the default is the existing
`/home/cylin/NETMHC/netMHCpan-4.2`. Changing models changes the benchmark.

```bash
python mouse_benchmark/curate.py
python mouse_benchmark/predict.py
python mouse_benchmark/evaluate.py
python mouse_benchmark/plot.py
python mouse_benchmark/audit.py
python -m unittest discover -s mouse_benchmark -p 'test_*.py' -v
```

These commands overwrite generated results but not the downloaded snapshot.
The report is a reviewed interpretation of this snapshot, not an automatically
updated conclusion for future datasets. `predict.py` runs both alleles, saves
the raw binary output, requires complete predictions and joins on peptide AND
allele. It imports the current production Calis scorer. A parity test checks the
benchmark composite against `BatchJob._finalize`, including clipping boundaries.

## Download a fresh snapshot

The original retrieval command was `python mouse_benchmark/fetch_iedb.py`.
It refuses to overwrite an existing snapshot. To inspect a new snapshot without
altering this run:

```bash
python mouse_benchmark/fetch_iedb.py --outdir /tmp/iedb_mouse_new_snapshot
```

The downloader uses mouse host ancestry and explicit `H2-Kb`/`H2-Db` annotations,
not the antigen's species. It includes all outcomes, pages by increasing assay ID,
checks total count and saves query URLs, UTC timestamps and a checksum. Retrieval
errors fail the run rather than silently producing an empty dataset. To benchmark
another snapshot, use an isolated copy of this directory with that `data/` folder.
IEDB is live; a later download is not expected to reproduce this snapshot exactly.

## Artifacts

| Path | Contents |
|---|---|
| `data/iedb_mouse_kb_db.jsonl` | Full selected-field snapshot, one original assay per line |
| `data/manifest.json` | Query filters, page counts/URLs, timestamps, SHA-256 |
| `results/assay_eligibility.csv` | Every assay and its first exclusion reason |
| `results/curated_pairs.csv` | Deduplicated primary cohort and evidence/source IDs |
| `results/conflicting_pairs.csv` | Excluded conflicting peptide:allele pairs |
| `results/studies.csv` | Study titles, PMID, dates, eligible assays and retained labels |
| `results/inventory.json` | Host, endpoint, restriction, source and leakage-group inventory |
| `results/all_predictions.csv` | All 6,260 peptides evaluated against both alleles |
| `results/scored_pairs.csv` | All 6,379 eligible experimental pairs with baseline scores |
| `results/netmhcpan_raw/` | Input peptide lists and unmodified predictor output |
| `results/prediction_manifest.json` | Binary/source/data hashes and exact commands |
| `results/metrics.csv` | All models × all prespecified analysis cohorts |
| `results/within_study_metrics.csv` | Studies with ≥5 positive and ≥5 negative pairs |
| `results/study_macro_metrics.csv` | Unweighted within-study averages; K metrics omit studies smaller than K |
| `results/uncertainty.json` | 1,000 paired component-bootstrap replicates, AP intervals |
| `results/training_split_audit.json` | Diagnostic grouped fold sizes; no model fitting |
| `results/baseline_comparison.png` / `.svg` | Exportable figures |
| `results/audit.json` | Every retained label/source ID verified; production parity on all pairs |
| `results/precision_sensitivity.csv` | Post-hoc check of rounding the base to composite's four decimals |

Existing published models reviewed: [Smith et al. neoag](https://github.com/vincentlaboratories/neoag)
requires SNV mutant/reference pairs and variant position, unavailable in this
general assay extract. It is not evaluated by fabricating missing inputs; its
weights are not bundled. [DeepNeo](https://academic.oup.com/nar/article/51/W1/W134/7127215)
was identified earlier as a candidate but is not evaluated in this run. Neither
model is assigned a performance number here.
