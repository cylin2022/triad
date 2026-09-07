"""
antigen_processing.py
Antigen Processing Module (Proteasomal Cleavage, TAP Transport, and Cellular Presentation).
Computes processing scores for peptide epitope candidates based on Stage 2 of plan.md.
"""

from typing import List, Dict, Any
import pandas as pd
import numpy as np

# C-terminal residues preferred by immunoproteasome & TAP transport
PREFERRED_CTERMINI = {'F': 0.9, 'Y': 0.85, 'W': 0.85, 'L': 0.8, 'I': 0.75, 'V': 0.7, 'M': 0.65, 'K': 0.6, 'R': 0.6}
PREFERRED_NTERMINI = {'A': 0.8, 'L': 0.8, 'R': 0.75, 'K': 0.7, 'I': 0.7, 'V': 0.65}


def calculate_antigen_processing(peptides: List[str], df_el: pd.DataFrame) -> pd.DataFrame:
    """
    Calculate Antigen Processing score for a list of peptides.

    Parameters
    ----------
    peptides : list of str
    df_el : pd.DataFrame with ['Peptide', 'Score_EL']

    Returns
    -------
    pd.DataFrame mapping Peptide -> ['Proteasome_Score', 'TAP_Score', 'Processing_Score']
    """
    if not peptides:
        return pd.DataFrame(columns=["Peptide", "Proteasome_Score", "TAP_Score", "Processing_Score"])

    el_map = {}
    if not df_el.empty and "Peptide" in df_el.columns and "Score_EL" in df_el.columns:
        el_map = df_el.groupby("Peptide")["Score_EL"].max().to_dict()

    rows = []
    for pep in list(dict.fromkeys(peptides)):
        c_term = pep[-1] if pep else ''
        n_term = pep[0] if pep else ''

        proteasome_score = PREFERRED_CTERMINI.get(c_term, 0.3)
        tap_score = (PREFERRED_CTERMINI.get(c_term, 0.3) * 0.6) + (PREFERRED_NTERMINI.get(n_term, 0.3) * 0.4)
        presentation_score = float(el_map.get(pep, 0.0))

        # Composite Stage 2 Processing Score
        processing_score = round(0.50 * presentation_score + 0.30 * proteasome_score + 0.20 * tap_score, 4)

        rows.append({
            "Peptide": pep,
            "Proteasome_Score": round(proteasome_score, 3),
            "TAP_Score": round(tap_score, 3),
            "Processing_Score": processing_score
        })

    return pd.DataFrame(rows)


if __name__ == "__main__":
    test_peps = ["AAAWYLWEV", "FLPSDFFPSV", "KMKGDYFRYF"]
    test_el = pd.DataFrame([
        {"Peptide": "AAAWYLWEV", "Score_EL": 0.364},
        {"Peptide": "FLPSDFFPSV", "Score_EL": 0.955},
        {"Peptide": "KMKGDYFRYF", "Score_EL": 0.006},
    ])
    res = calculate_antigen_processing(test_peps, test_el)
    print("=== Antigen Processing Test ===")
    print(res.to_string(index=False))
