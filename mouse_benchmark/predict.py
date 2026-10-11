"""Run unchanged local production scorers, retaining raw netMHCpan output."""
import hashlib
import json
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
sys.path.insert(0, str(REPO / "netmhcpan_platform"))
import netmhcpan_api as api
from immunogenicity_scorer import score_dataframe


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def add_scores(df):
    df = df.copy()
    df["Calis"] = score_dataframe(df)
    pres = np.clip(.5 * df.Score_EL + .5 * np.clip(1 - df.Rank_EL / 2, 0, 1), 0, 1)
    affinity = np.clip(1 - np.log10(df.Affinity_nM.clip(1, 50000)) / np.log10(50000), 0, 1)
    df["Mouse_Base"] = .55 * pres + .45 * affinity
    df["Mouse_Composite"] = np.clip(df.Mouse_Base * np.clip(1 + .35 * df.Calis, .5, 1.5), 0, 1).round(4)
    df["EL_Rank_Score"] = -df.Rank_EL
    df["BA_Rank_Score"] = -df.Rank_BA
    return df


def predict():
    out = ROOT / "results"
    rawdir = out / "netmhcpan_raw"
    rawdir.mkdir(exist_ok=True)
    df = pd.read_csv(out / "curated_pairs.csv")
    peptides = sorted(df.Peptide.unique())
    # Keep scratch files inside the authorized workspace, not a global RAM disk.
    api.TMPDIR = rawdir / "tmp"
    api.TMPDIR.mkdir(exist_ok=True)
    runner = api.NetMHCpanRunner(timeout=1800)
    frames, commands = [], []
    for allele in ["H-2-Kb", "H-2-Db"]:
        pepfile = rawdir / f"{allele}.pep"
        pepfile.write_text("\n".join(peptides) + "\n")
        cmd = runner._build_cmd(str(pepfile), "-p", [allele], True, False, False, False)
        print(f"Predicting {len(peptides)} peptides for {allele}", flush=True)
        raw = runner._run(cmd)
        (rawdir / f"{allele}.txt").write_text(raw)
        frame = api.parse_output(raw, include_ba=True)
        if len(frame) != len(peptides) or set(frame.Peptide) != set(peptides) or set(frame.MHC) != {allele}:
            raise ValueError(f"Missing, duplicate or unexpected predictions for {allele}")
        frames.append(frame)
        commands.append(cmd)
    pred = pd.concat(frames, ignore_index=True)
    pred.to_csv(out / "all_predictions.csv", index=False)
    merged = df.merge(pred, on=["Peptide", "MHC"], how="left", validate="one_to_one")
    required = ["Score_EL", "Rank_EL", "Score_BA", "Rank_BA", "Affinity_nM"]
    if len(merged) != len(df) or not np.isfinite(merged[required].to_numpy()).all():
        raise ValueError("Incomplete or invalid prediction join")
    add_scores(merged).to_csv(out / "scored_pairs.csv", index=False)
    paths = [REPO / "netmhcpan_platform" / n for n in ["immunogenicity_scorer.py", "batch_runner.py", "netmhcpan_api.py"]]
    paths += [api.BINARY, api.SYN_NOCONTEXT, api.HLA_PSEUDO, api.VERSION_FILE, ROOT / "predict.py", out / "curated_pairs.csv"]
    manifest = {
        "utc": datetime.now(timezone.utc).isoformat(), "python": platform.python_version(),
        "numpy": np.__version__, "pandas": pd.__version__, "netmhcpan_version": api.VERSION_FILE.read_text().strip(),
        "commands": commands, "hashes": {str(p.relative_to(REPO)) if p.is_relative_to(REPO) else str(p): sha(p) for p in paths},
        "pair_predictions": len(pred), "evaluated_pairs": len(merged), "missing_predictions": 0,
    }
    (out / "prediction_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"Scored {len(merged)} experimental peptide:allele pairs", flush=True)


if __name__ == "__main__":
    predict()
