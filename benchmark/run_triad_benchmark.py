#!/usr/bin/env python3
"""
run_triad_benchmark.py
Submit benchmark peptides to running TRIAD platform, fetch predictions,
and merge with experimental ground-truth labels.
"""

import sys
import time
import json
import argparse
import urllib.request
import urllib.parse
import pandas as pd
from pathlib import Path

def submit_job(peptides, allele="HLA-A02:01", server_url="http://localhost:5001"):
    url = f"{server_url}/api/predict"
    payload = {
        "mode": "peptide",
        "input": "\n".join(peptides),
        "alleles": [allele],
        "include_ba": True,
        "name": f"Benchmark_{allele}"
    }
    req = urllib.request.Request(
        url,
        data=json.dumps(payload).encode("utf-8"),
        headers={"Content-Type": "application/json"}
    )
    with urllib.request.urlopen(req, timeout=30) as resp:
        res = json.loads(resp.read().decode("utf-8"))
        return res["job_id"]

def wait_for_job(job_id, server_url="http://localhost:5001", max_wait_sec=120):
    url = f"{server_url}/api/job/{job_id}"
    start = time.time()
    while time.time() - start < max_wait_sec:
        req = urllib.request.Request(url)
        with urllib.request.urlopen(req, timeout=10) as resp:
            data = json.loads(resp.read().decode("utf-8"))
            status = data.get("status")
            if status in ("completed", "done"):
                print(f"[✓] Job {job_id} completed successfully in {time.time()-start:.1f}s.")
                return data
            elif status == "failed":
                raise RuntimeError(f"Job {job_id} failed: {data.get('error')}")
            else:
                progress = data.get("progress_pct", 0)
                print(f"[*] Job {job_id} in progress ({progress}%)...")
                time.sleep(1.5)
    raise TimeoutError(f"Job {job_id} timed out after {max_wait_sec}s.")

def download_csv(csv_url, server_url="http://localhost:5001"):
    full_url = f"{server_url}{csv_url}"
    req = urllib.request.Request(full_url)
    with urllib.request.urlopen(req, timeout=30) as resp:
        return pd.read_csv(resp)

def run_benchmark(dataset_path, allele="HLA-A02:01", server_url="http://localhost:5001", outdir="benchmark/results"):
    dataset = pd.read_csv(dataset_path)
    print(f"[*] Loaded benchmark dataset from {dataset_path} ({len(dataset)} peptides).")
    
    peptides = dataset["Peptide"].unique().tolist()
    print(f"[*] Submitting {len(peptides)} peptides to TRIAD API ({server_url})...")
    
    job_id = submit_job(peptides, allele=allele, server_url=server_url)
    print(f"[*] Submitted Job ID: {job_id}")
    
    job_info = wait_for_job(job_id, server_url=server_url)
    csv_url = job_info["csv_url"]
    
    print(f"[*] Downloading prediction results from {csv_url}...")
    pred_df = download_csv(csv_url, server_url=server_url)
    print(f"[✓] Predictions received: {len(pred_df)} rows.")

    # Merge predictions with dataset labels
    merged = pd.merge(pred_df, dataset[["Peptide", "Label", "Assay_Outcome", "PubMed_ID"]], on="Peptide", how="inner")
    
    # Save results
    outdir_p = Path(outdir)
    outdir_p.mkdir(parents=True, exist_ok=True)
    out_file = outdir_p / "triad_predictions.csv"
    merged.to_csv(out_file, index=False)
    print(f"[✓] Final merged benchmark predictions saved to: {out_file}")
    print(f"    Total evaluated samples: {len(merged)}")
    print(f"    Columns available: {list(merged.columns)}")
    return merged

def main():
    parser = argparse.ArgumentParser(description="Run TRIAD benchmark on curated dataset")
    parser.add_argument("--data", default="benchmark/data/iedb_HLA-A0201_benchmark.csv", help="Path to input dataset CSV")
    parser.add_argument("--allele", default="HLA-A02:01", help="Target HLA allele")
    parser.add_argument("--server", default="http://localhost:5001", help="TRIAD server base URL")
    parser.add_argument("--outdir", default="benchmark/results", help="Output directory")
    args = parser.parse_args()

    run_benchmark(args.data, allele=args.allele, server_url=args.server, outdir=args.outdir)

if __name__ == "__main__":
    main()
