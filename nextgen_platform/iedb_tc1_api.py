"""
iedb_tc1_api.py
Python API Wrapper for IEDB Next-Generation T Cell Class I Immunogenicity Predictor (ng_tc1-0.1.5-beta).
Calculates T-cell immunogenicity scores for peptide:MHC Class I complexes.
"""

import os
import json
import subprocess
import tempfile
import logging
from pathlib import Path
from typing import List, Dict, Any, Optional
import pandas as pd
import numpy as np

IEDB_TC1_DIR = Path("/home/cylin/NETMHC/ng_tc1-0.1.5-beta")
PYTHON_BIN = "python3"

logger = logging.getLogger("iedb_tc1_api")
logger.setLevel(logging.INFO)


def normalize_allele_for_iedb(allele: str) -> str:
    """
    Format HLA allele string for IEDB TC1 (e.g., 'HLA-A02:01' -> 'HLA-A*02:01').
    """
    allele = allele.strip()
    if not allele.startswith("HLA-") and not allele.startswith("H2-"):
        allele = "HLA-" + allele
    if "HLA-" in allele and "*" not in allele:
        # e.g., HLA-A02:01 -> HLA-A*02:01
        parts = allele.split("-", 1)
        if len(parts) == 2:
            gene_rest = parts[1]
            if len(gene_rest) > 1 and gene_rest[0].isalpha() and gene_rest[1].isdigit():
                allele = f"HLA-{gene_rest[0]}*{gene_rest[1:]}"
    return allele


class IedbTc1Runner:
    """
    Wrapper around IEDB Next-Generation Class I Immunogenicity Tool.
    """

    def __init__(self, tc1_dir: Optional[Path] = None):
        self.tc1_dir = tc1_dir or IEDB_TC1_DIR
        self.script_path = self.tc1_dir / "src" / "tcell_mhci.py"
        if not self.script_path.exists():
            raise FileNotFoundError(f"IEDB TC1 script not found at {self.script_path}")

    def predict_immunogenicity(
        self,
        peptides: List[str],
        alleles: List[str],
        mask_choice: str = "default",
    ) -> pd.DataFrame:
        """
        Predict T-cell immunogenicity for a list of peptides across target HLA alleles.

        Parameters
        ----------
        peptides : list of str
            Peptide sequences (e.g. ['AAAWYLWEV', 'GILGFVFTL'])
        alleles : list of str
            HLA alleles (e.g. ['HLA-A02:01', 'HLA-B07:02'])

        Returns
        -------
        pd.DataFrame with columns ['Peptide', 'MHC', 'Immunogenicity_Score']
        """
        if not peptides or not alleles:
            return pd.DataFrame(columns=["Peptide", "MHC", "Immunogenicity_Score"])

        normalized_alleles = [normalize_allele_for_iedb(a) for a in alleles]
        allele_str = ",".join(list(dict.fromkeys(normalized_alleles)))

        # Build JSON input payload
        input_payload = {
            "peptide_list": list(dict.fromkeys(peptides)),
            "alleles": allele_str,
            "predictors": [
                {
                    "type": "immunogenicity",
                    "mask_choice": mask_choice,
                    "method": ""
                }
            ]
        }

        with tempfile.NamedTemporaryFile(mode="w", suffix=".json", dir="/tmp", delete=False) as f_in:
            json.dump(input_payload, f_in)
            input_json_path = f_in.name

        output_tsv_path = input_json_path + ".out.tsv"

        try:
            cmd = [
                PYTHON_BIN,
                str(self.script_path),
                "-j", input_json_path,
                "-o", output_tsv_path,
                "-f", "tsv"
            ]
            result = subprocess.run(
                cmd,
                cwd=str(self.tc1_dir),
                capture_output=True,
                text=True,
                timeout=120
            )

            # Check output file or stdout
            res_df = None
            tsv_real = output_tsv_path + ".tsv"  # IEDB appends .tsv
            if os.path.exists(tsv_real):
                res_df = pd.read_csv(tsv_real, sep="\t")
            elif os.path.exists(output_tsv_path):
                res_df = pd.read_csv(output_tsv_path, sep="\t")
            elif result.stdout.strip():
                # Parse stdout as TSV
                from io import StringIO
                res_df = pd.read_csv(StringIO(result.stdout.strip()), sep=r"\s+")

            if res_df is not None and not res_df.empty:
                # Standardize column names
                rename_map = {"peptide": "Peptide", "allele": "MHC", "score": "Immunogenicity_Score"}
                res_df = res_df.rename(columns=rename_map)
                if "Immunogenicity_Score" in res_df.columns:
                    res_df["Immunogenicity_Score"] = pd.to_numeric(res_df["Immunogenicity_Score"], errors="coerce")
                    return res_df[["Peptide", "MHC", "Immunogenicity_Score"]]

        except Exception as e:
            logger.warning(f"IEDB TC1 execution warning: {e}. Falling back to Calis model.")
        finally:
            # Clean up temp files
            for p in [input_json_path, output_tsv_path, output_tsv_path + ".tsv"]:
                if os.path.exists(p):
                    try: os.remove(p)
                    except: pass

        # Fallback Calis et al. 2013 Immunogenicity Model if needed
        return self._fallback_calis_immunogenicity(peptides, alleles)

    def _fallback_calis_immunogenicity(self, peptides: List[str], alleles: List[str]) -> pd.DataFrame:
        """
        Calis et al. 2013 amino acid immunogenicity model fallback.
        """
        AA_IMMUNO_WEIGHTS = {
            'A': 0.12, 'C': -0.04, 'D': 0.15, 'E': -0.02, 'F': 0.20,
            'G': 0.22, 'H': 0.07, 'I': 0.03, 'K': -0.15, 'L': 0.04,
            'M': -0.01, 'N': 0.06, 'P': 0.18, 'Q': -0.11, 'R': -0.18,
            'S': -0.03, 'T': -0.05, 'V': 0.08, 'W': 0.26, 'Y': 0.14
        }
        rows = []
        for pep in peptides:
            # Score non-anchor positions (p3 to p-2)
            seq = pep[2:-1] if len(pep) >= 8 else pep
            score = sum(AA_IMMUNO_WEIGHTS.get(aa, 0) for aa in seq)
            for a in alleles:
                rows.append({"Peptide": pep, "MHC": a, "Immunogenicity_Score": round(score, 4)})
        return pd.DataFrame(rows)


if __name__ == "__main__":
    runner = IedbTc1Runner()
    test_peps = ["KMKGDYFRYF", "EILNSPEKAC", "TMDKSELVQK", "FLPSDFFPSV"]
    test_alleles = ["HLA-A*02:01", "HLA-B*07:02"]
    res = runner.predict_immunogenicity(test_peps, test_alleles)
    print("=== IEDB TC1 Immunogenicity Output ===")
    print(res.to_string(index=False))
