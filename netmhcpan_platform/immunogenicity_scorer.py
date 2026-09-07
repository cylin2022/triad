"""
immunogenicity_scorer.py
IEDB Class I T-cell immunogenicity scoring (Calis et al., PLoS Comput Biol 2013).

Ported directly from the self-contained ImmunogenicityPredictor in IEDB's
ng_tc1-0.1.5-beta (src/immunogenicity_predictor.py), licensed under the
Non-Profit Open Software License 3.0 (full text bundled at
THIRD_PARTY_LICENSES/license-LJI.txt), which explicitly permits
redistribution and network deployment.

Unlike the full ng_tc1 "Next-Gen T Cell Class I" pipeline, this scorer needs
no proteasomal-cleavage (netchop) or MHC-binding (netMHC/netMHCpan/netCTL)
step of its own — those tools are DTU-licensed and cannot be redistributed.
It only needs a peptide sequence (and optionally the presenting HLA allele),
which netMHCpan-4.2 has already produced, so it runs in-process with no
subprocess call and no extra dependency.
"""

from typing import Optional

import pandas as pd

# Position-specific amino acid weights (Calis et al. 2013).
AA_IMMUNOSCALE = {
    "A": 0.127, "C": -0.175, "D": 0.072, "E": 0.325, "F": 0.380,
    "G": 0.110, "H": 0.105, "I": 0.432, "K": -0.700, "L": -0.036,
    "M": -0.570, "N": -0.021, "P": -0.036, "Q": -0.376, "R": 0.168,
    "S": -0.537, "T": 0.126, "V": 0.134, "W": 0.719, "Y": -0.012,
}

# Positional contribution weights for a 9-mer (index 0 = N-terminus).
POSITION_WEIGHTS = [0.00, 0.00, 0.10, 0.31, 0.30, 0.29, 0.26, 0.18, 0.00]

# Allele-specific anchor positions to mask (1-indexed, "cterm" implied last).
ALLELE_POSITION_MASK = {
    "H-2-Db": "2,5,9", "H-2-Dd": "2,3,5", "H-2-Kb": "2,3,9", "H-2-Kd": "2,5,9",
    "H-2-Kk": "2,8,9", "H-2-Ld": "2,5,9", "HLA-A0101": "2,3,9", "HLA-A0201": "1,2,9",
    "HLA-A0202": "1,2,9", "HLA-A0203": "1,2,9", "HLA-A0206": "1,2,9", "HLA-A0211": "1,2,9",
    "HLA-A0301": "1,2,9", "HLA-A1101": "1,2,9", "HLA-A2301": "2,7,9", "HLA-A2402": "2,7,9",
    "HLA-A2601": "1,2,9", "HLA-A2902": "2,7,9", "HLA-A3001": "1,3,9", "HLA-A3002": "2,7,9",
    "HLA-A3101": "1,2,9", "HLA-A3201": "1,2,9", "HLA-A3301": "1,2,9", "HLA-A6801": "1,2,9",
    "HLA-A6802": "1,2,9", "HLA-A6901": "1,2,9", "HLA-B0702": "1,2,9", "HLA-B0801": "2,5,9",
    "HLA-B1501": "1,2,9", "HLA-B1502": "1,2,9", "HLA-B1801": "1,2,9", "HLA-B2705": "2,3,9",
    "HLA-B3501": "1,2,9", "HLA-B3901": "1,2,9", "HLA-B4001": "1,2,9", "HLA-B4002": "1,2,9",
    "HLA-B4402": "2,3,9", "HLA-B4403": "2,3,9", "HLA-B4501": "1,2,9", "HLA-B4601": "1,2,9",
    "HLA-B5101": "1,2,9", "HLA-B5301": "1,2,9", "HLA-B5401": "1,2,9", "HLA-B5701": "1,2,9",
    "HLA-B5801": "1,2,9",
}

VALID_AA = set(AA_IMMUNOSCALE.keys())


def _normalize_allele(allele: Optional[str]) -> Optional[str]:
    if not allele:
        return None
    return allele.replace("*", "").replace(":", "").replace("H2", "H-2")


def score_peptide(peptide: str, allele: Optional[str] = None) -> Optional[float]:
    """
    Score a single peptide's predicted T-cell immunogenicity.

    Uses allele-specific anchor-position masking when the allele is one of
    the alleles Calis et al. characterized; otherwise falls back to the
    generic default mask (positions 1, 2, and the C-terminus).

    Returns None if the peptide contains characters outside the 20 standard
    amino acids.
    """
    peptide = peptide.strip().upper()
    peplen = len(peptide)
    if peplen == 0 or any(aa not in VALID_AA for aa in peptide):
        return None

    cterm = peplen - 1
    mask_str = ALLELE_POSITION_MASK.get(_normalize_allele(allele))
    if mask_str:
        mask_positions = [int(p) - 1 for p in mask_str.split(",")]
    else:
        mask_positions = [0, 1, cterm]

    if peplen > 9:
        weights = POSITION_WEIGHTS[:5] + [0.30] * (peplen - 9) + POSITION_WEIGHTS[5:]
    else:
        weights = POSITION_WEIGHTS

    score = 0.0
    for pos, aa in enumerate(peptide):
        if pos not in mask_positions:
            score += weights[pos] * AA_IMMUNOSCALE[aa]
    return round(score, 5)


def score_dataframe(
    df: pd.DataFrame,
    peptide_col: str = "Peptide",
    allele_col: str = "MHC",
) -> pd.Series:
    """Vectorized T-cell immunogenicity scoring for a netMHCpan-4.2 results DataFrame."""
    if df.empty or peptide_col not in df.columns:
        return pd.Series(dtype="float64", index=df.index)

    alleles = df[allele_col] if allele_col in df.columns else [None] * len(df)
    return pd.Series(
        [score_peptide(pep, allele) for pep, allele in zip(df[peptide_col], alleles)],
        index=df.index,
        dtype="float64",
    )
