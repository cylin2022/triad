"""
hla_breadth.py
HLA Breadth Analysis Module.
Calculates cross-allele binding capability (HLA-A, HLA-B, HLA-C counts, total binder count, and HLA Breadth Score).
"""

from typing import Dict, List, Any
import pandas as pd
import numpy as np


def calculate_hla_breadth(df: pd.DataFrame) -> pd.DataFrame:
    """
    Calculate HLA Breadth metrics for each unique peptide in the dataset.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing columns ['Peptide', 'MHC', 'Rank_EL', 'BindLevel']

    Returns
    -------
    pd.DataFrame mapping Peptide -> ['n_hla_binders', 'n_hla_sb', 'n_hla_a', 'n_hla_b', 'n_hla_c', 'HLA_Breadth_Score', 'HLA_Breadth_Category']
    """
    if df.empty or "Peptide" not in df.columns or "MHC" not in df.columns or "Rank_EL" not in df.columns:
        return pd.DataFrame(columns=[
            "Peptide", "n_hla_binders", "n_hla_sb", "n_hla_a", "n_hla_b", "n_hla_c",
            "HLA_Breadth_Score", "HLA_Breadth_Category"
        ])

    binders_df = df[df["Rank_EL"] <= 2.0].copy()
    all_peptides = df["Peptide"].unique()
    total_eval_alleles = max(1, df["MHC"].nunique())

    results = []
    for pep in all_peptides:
        pep_binders = binders_df[binders_df["Peptide"] == pep]
        n_binders = len(pep_binders)
        n_sb = len(pep_binders[pep_binders["Rank_EL"] <= 0.5])

        # Gene counts
        n_a = len(pep_binders[pep_binders["MHC"].str.contains("HLA-A", case=False, na=False)])
        n_b = len(pep_binders[pep_binders["MHC"].str.contains("HLA-B", case=False, na=False)])
        n_c = len(pep_binders[pep_binders["MHC"].str.contains("HLA-C", case=False, na=False)])

        # Breadth Score (0.0 to 1.0)
        breadth_score = min(1.0, round((n_sb * 1.0 + (n_binders - n_sb) * 0.5) / max(1, total_eval_alleles), 4))

        if breadth_score >= 0.4 or n_binders >= 3:
            cat = "High"
        elif breadth_score >= 0.15 or n_binders >= 2:
            cat = "Medium"
        else:
            cat = "Low"

        results.append({
            "Peptide": pep,
            "n_hla_binders": n_binders,
            "n_hla_sb": n_sb,
            "n_hla_a": n_a,
            "n_hla_b": n_b,
            "n_hla_c": n_c,
            "HLA_Breadth_Score": breadth_score,
            "HLA_Breadth_Category": cat,
        })

    return pd.DataFrame(results)


if __name__ == "__main__":
    test_data = pd.DataFrame([
        {"Peptide": "AAAWYLWEV", "MHC": "HLA-A*02:01", "Rank_EL": 0.14},
        {"Peptide": "AAAWYLWEV", "MHC": "HLA-A*24:02", "Rank_EL": 0.41},
        {"Peptide": "AAAWYLWEV", "MHC": "HLA-B*07:02", "Rank_EL": 1.20},
        {"Peptide": "FLPSDFFPSV", "MHC": "HLA-A*02:01", "Rank_EL": 0.05},
        {"Peptide": "FLPSDFFPSV", "MHC": "HLA-B*07:02", "Rank_EL": 15.0},
    ])
    res = calculate_hla_breadth(test_data)
    print("=== HLA Breadth Calculation Test ===")
    print(res.to_string(index=False))
