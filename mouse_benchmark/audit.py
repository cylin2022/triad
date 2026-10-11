"""Final artifact checks, production parity and score precision sensitivity."""
import hashlib
import importlib.metadata
import json
import platform
from pathlib import Path
from unittest.mock import Mock

import numpy as np
import pandas as pd

from curate import POSITIVE, rejection
from evaluate import metrics
from predict import ROOT
from batch_runner import BatchJob


def audit():
    out = ROOT / "results"
    df = pd.read_csv(out / "scored_pairs.csv", dtype={"Reference_IDs": str, "Assay_IDs": str})
    raw = [json.loads(line) for line in (ROOT / "data/iedb_mouse_kb_db.jsonl").open()]
    eligible = {}
    for row in raw:
        if rejection(row) == "eligible":
            key = row["linear_sequence"], row["mhc_allele_name"].replace("H2-", "H-2-")
            eligible.setdefault(key, []).append(row)
    expected = {k: rr for k,rr in eligible.items() if len({r["qualitative_measure"] in POSITIVE for r in rr}) == 1}
    assert len(expected) == len(df)
    assert set(expected) == set(zip(df.Peptide, df.MHC))
    for row in df.itertuples():
        rr = expected[row.Peptide, row.MHC]
        assert row.Label == int(rr[0]["qualitative_measure"] in POSITIVE)
        assert set(row.Assay_IDs.split(";")) == {str(r["tcell_id"]) for r in rr}
        assert set(row.Reference_IDs.split(";")) == {str(r["reference_id"]) for r in rr}
    pred_manifest = json.loads((out / "prediction_manifest.json").read_text())
    for relative, expected_hash in pred_manifest["hashes"].items():
        path = ROOT.parent / relative
        assert hashlib.sha256(path.read_bytes()).hexdigest() == expected_hash, f"Changed prediction input: {relative}"
    # Use the actual current production path for ALL scored pairs, not a copied formula.
    job = BatchJob.__new__(BatchJob)
    job.alleles, job.meta, job.logger, job.job_id = ["H-2-Kb", "H-2-Db"], {}, Mock(), "mouse_benchmark_audit"
    inputs = pd.read_csv(out / "all_predictions.csv")
    inputs = df[["Peptide", "MHC"]].merge(inputs, on=["Peptide", "MHC"], validate="one_to_one")
    job._finalize(inputs)
    actual = job.results.set_index(["Peptide", "MHC"]).sort_index()
    saved = df.set_index(["Peptide", "MHC"]).sort_index()
    np.testing.assert_allclose(actual.Immunogenicity_Score, saved.Calis, rtol=0, atol=1e-12)
    np.testing.assert_allclose(actual.Composite_Score, saved.Mouse_Composite, rtol=0, atol=1e-12)
    # Post-hoc technical sensitivity: precision alone must not drive the interpretation.
    df["Mouse_Base_Rounded"] = df.Mouse_Base.round(4)
    precision = [{"Model": m, **metrics(df, m)} for m in ["Mouse_Base", "Mouse_Base_Rounded", "Mouse_Composite"]]
    pd.DataFrame(precision).to_csv(out / "precision_sensitivity.csv", index=False)
    env = {"python": platform.python_version(), "platform": platform.platform(),
           "packages": {p: importlib.metadata.version(p) for p in ["numpy", "pandas", "scikit-learn", "scipy", "matplotlib"]}}
    (out / "environment.json").write_text(json.dumps(env, indent=2) + "\n")
    result = {"eligible_pairs_before_conflicts": len(eligible), "verified_pairs": len(df),
              "all_assay_labels_and_provenance_match": True, "prediction_input_hashes_match": True,
              "production_Calis_and_composite_parity_all_pairs": True,
              "precision_sensitivity": "post-hoc numerical audit, not a trained model"}
    (out / "audit.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))
    print(pd.DataFrame(precision).to_string(index=False))


if __name__ == "__main__":
    audit()
