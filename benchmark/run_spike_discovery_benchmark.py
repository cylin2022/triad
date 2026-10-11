#!/usr/bin/env python3
"""
run_spike_discovery_benchmark.py
Evaluate TRIAD on a real-world whole-protein epitope discovery task:
SARS-CoV-2 Spike glycoprotein (1,273 aa, 1,265 overlapping 9-mers).
55 experimentally validated HLA-A*02:01 positive CD8+ epitopes from IEDB
vs. 1,210 background non-epitope 9-mers.
"""

import sys
import time
import json
import argparse
import urllib.request
import urllib.parse
import pandas as pd
from pathlib import Path

def get_spike_sequence(fasta_path="multi_fasta.fasta"):
    with open(fasta_path, "r", encoding="utf-8") as f:
        text = f.read()
    for entry in text.split(">"):
        if "P0DTC2" in entry:
            lines = entry.strip().split("\n")
            return "".join(lines[1:])
    raise ValueError("P0DTC2 Spike sequence not found in fasta file.")

def fetch_iedb_spike_epitopes():
    """
    Fetch confirmed HLA-A*02:01 positive 9-mer epitopes for Spike from IEDB.
    """
    params = {
        "parent_source_antigen_name": "ilike.*Spike glycoprotein*",
        "mhc_allele_name": "eq.HLA-A*02:01",
        "qualitative_measure": "eq.Positive",
        "select": "linear_sequence",
        "limit": "1000"
    }
    url = "https://query-api.iedb.org/tcell_search?" + urllib.parse.urlencode(params)
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0 (TRIAD-Benchmark)"})
    with urllib.request.urlopen(req, timeout=30) as resp:
        data = json.loads(resp.read().decode("utf-8"))
    
    seqs = set(d["linear_sequence"].strip().upper() for d in data if d.get("linear_sequence") and len(d.get("linear_sequence").strip()) == 9)
    return seqs

def submit_and_wait(peptides, allele="HLA-A02:01", server_url="http://localhost:5001"):
    url = f"{server_url}/api/predict"
    payload = {
        "mode": "peptide",
        "input": "\n".join(peptides),
        "alleles": [allele],
        "include_ba": True,
        "name": "Spike_Epitope_Discovery"
    }
    req = urllib.request.Request(
        url,
        data=json.dumps(payload).encode("utf-8"),
        headers={"Content-Type": "application/json"}
    )
    with urllib.request.urlopen(req, timeout=30) as resp:
        res = json.loads(resp.read().decode("utf-8"))
        job_id = res["job_id"]
    print(f"[*] Submitted Spike Discovery Job ID: {job_id}")

    # Wait for completion
    status_url = f"{server_url}/api/job/{job_id}"
    start = time.time()
    while time.time() - start < 180:
        req = urllib.request.Request(status_url)
        with urllib.request.urlopen(req, timeout=10) as resp:
            data = json.loads(resp.read().decode("utf-8"))
            status = data.get("status")
            if status in ("completed", "done"):
                print(f"[✓] Job {job_id} finished in {time.time()-start:.1f}s.")
                csv_url = f"{server_url}{data['csv_url']}"
                req_csv = urllib.request.Request(csv_url)
                with urllib.request.urlopen(req_csv, timeout=30) as c_resp:
                    return pd.read_csv(c_resp)
            elif status == "failed":
                raise RuntimeError(f"Job failed: {data.get('error')}")
            time.sleep(1.5)
    raise TimeoutError("Job timed out.")

def main():
    parser = argparse.ArgumentParser(description="Run whole-protein Spike discovery benchmark")
    parser.add_argument("--fasta", default="multi_fasta.fasta", help="Path to multi_fasta.fasta")
    parser.add_argument("--server", default="http://localhost:5001", help="TRIAD server URL")
    parser.add_argument("--outdir", default="benchmark/results", help="Output directory")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    print("[*] Extracting SARS-CoV-2 Spike glycoprotein (P0DTC2)...")
    spike_seq = get_spike_sequence(args.fasta)
    print(f"    Sequence length: {len(spike_seq)} aa")

    print("[*] Generating all 1,265 overlapping 9-mer peptides...")
    peptides = []
    positions = []
    for i in range(len(spike_seq) - 8):
        peptides.append(spike_seq[i:i+9])
        positions.append(i + 1)

    print("[*] Fetching validated HLA-A*02:01 epitopes from IEDB...")
    iedb_epitopes = fetch_iedb_spike_epitopes()
    positives = set(p for p in peptides if p in iedb_epitopes)
    print(f"[✓] Found {len(positives)} experimentally verified epitopes in Spike.")

    df_meta = pd.DataFrame({
        "Peptide": peptides,
        "Protein_Pos": positions,
        "Is_Verified_Epitope": [1 if p in positives else 0 for p in peptides]
    })

    print(f"[*] Submitting {len(peptides)} peptides to TRIAD platform...")
    pred_df = submit_and_wait(peptides, allele="HLA-A02:01", server_url=args.server)

    # Merge predictions with ground-truth epitope annotations
    merged = pd.merge(pred_df, df_meta, on="Peptide", how="inner")
    # De-duplicate if same 9-mer appears in multiple positions (keep first)
    merged = merged.drop_duplicates(subset=["Peptide"]).reset_index(drop=True)

    out_file = outdir / "spike_discovery_predictions.csv"
    merged.to_csv(out_file, index=False)
    print(f"[✓] Successfully saved Spike discovery benchmark to: {out_file}")
    print(f"    Total samples: {len(merged)}")
    print(f"    True Positive Epitopes: {merged['Is_Verified_Epitope'].sum()}")
    print(f"    Background Non-Epitopes: {len(merged) - merged['Is_Verified_Epitope'].sum()}")

if __name__ == "__main__":
    main()
