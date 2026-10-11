#!/usr/bin/env python3
"""
fetch_iedb_benchmark.py
Fetch and curate experimentally validated T-cell assay dataset from IEDB Query API.
Filters for Class I HLA alleles (default HLA-A*02:01), 9-mer canonical peptides,
and balances Positive vs Negative assay outcomes.
"""

import os
import sys
import json
import argparse
import urllib.request
import urllib.parse
import pandas as pd
from pathlib import Path

CANONICAL_AA = set("ACDEFGHIKLMNPQRSTVWY")

def fetch_iedb_assays(allele="HLA-A*02:01", qualitative_measure="Positive", limit=1000):
    """
    Query IEDB Query API (PostgREST) for T-cell assays with specified allele and outcome.
    """
    base_url = "https://query-api.iedb.org/tcell_search"
    params = {
        "mhc_allele_name": f"eq.{allele}",
        "qualitative_measure": f"eq.{qualitative_measure}",
        "select": "tcell_id,linear_sequence,mhc_allele_name,qualitative_measure,pubmed_id",
        "limit": str(limit),
        "order": "tcell_id.desc"
    }
    query_string = urllib.parse.urlencode(params)
    url = f"{base_url}?{query_string}"
    
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0 (TRIAD-Benchmark)"})
    try:
        with urllib.request.urlopen(req, timeout=30) as resp:
            data = json.loads(resp.read().decode("utf-8"))
            return data
    except Exception as e:
        print(f"[ERROR] Failed to fetch {qualitative_measure} assays: {e}", file=sys.stderr)
        return []

def curate_dataset(allele="HLA-A*02:01", target_per_class=250):
    print(f"[*] Fetching positive assays for {allele} from IEDB...")
    pos_data = fetch_iedb_assays(allele, "Positive", limit=target_per_class * 4)
    print(f"    Retrieved {len(pos_data)} raw positive records.")

    print(f"[*] Fetching negative assays for {allele} from IEDB...")
    neg_data = fetch_iedb_assays(allele, "Negative", limit=target_per_class * 4)
    print(f"    Retrieved {len(neg_data)} raw negative records.")

    def process_records(records, label):
        peptides = {}
        for r in records:
            seq = r.get("linear_sequence")
            if not seq or not isinstance(seq, str):
                continue
            seq = seq.strip().upper()
            # 9-mer canonical peptides only
            if len(seq) == 9 and all(aa in CANONICAL_AA for aa in seq):
                if seq not in peptides:
                    peptides[seq] = {
                        "Peptide": seq,
                        "Allele": allele,
                        "Label": label,
                        "Assay_Outcome": r.get("qualitative_measure"),
                        "PubMed_ID": r.get("pubmed_id", "")
                    }
        return peptides

    pos_peptides = process_records(pos_data, 1)
    neg_peptides = process_records(neg_data, 0)

    # Discard any ambiguous peptides that appear in both positive and negative assays
    conflicts = set(pos_peptides.keys()).intersection(set(neg_peptides.keys()))
    if conflicts:
        print(f"[!] Removing {len(conflicts)} conflicting peptides present in both pos and neg assays.")
        for p in conflicts:
            pos_peptides.pop(p, None)
            neg_peptides.pop(p, None)

    # Sample balanced sets
    pos_list = list(pos_peptides.values())[:target_per_class]
    neg_list = list(neg_peptides.values())[:target_per_class]

    print(f"[✓] Curated {len(pos_list)} Positive and {len(neg_list)} Negative 9-mers.")
    combined = pos_list + neg_list
    df = pd.DataFrame(combined)
    return df

def main():
    parser = argparse.ArgumentParser(description="Fetch and curate IEDB benchmark dataset")
    parser.add_argument("--allele", default="HLA-A*02:01", help="Target HLA allele")
    parser.add_argument("--samples-per-class", type=int, default=250, help="Number of samples per class")
    parser.add_argument("--outdir", default="benchmark/data", help="Output directory")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    clean_allele = args.allele.replace("*", "").replace(":", "")
    out_file = outdir / f"iedb_{clean_allele}_benchmark.csv"

    df = curate_dataset(args.allele, args.samples_per_class)
    df.to_csv(out_file, index=False)
    print(f"[✓] Dataset successfully saved to: {out_file}")
    print(f"    Total rows: {len(df)}")
    print(f"    Class balance: {df['Label'].value_counts().to_dict()}")

if __name__ == "__main__":
    main()
