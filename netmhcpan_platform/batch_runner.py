"""
batch_runner.py
Batch processing pipeline for NetMHCpan-4.2.
Handles large peptide lists, multiple FASTA files, and multi-allele jobs.
"""

import os
import math
import json
import logging
from pathlib import Path
from datetime import datetime
from typing import List, Optional, Dict, Any, Callable

import pandas as pd
import numpy as np

from netmhcpan_api import NetMHCpanRunner
from resource_tuner import get_optimal_batch_config, cleanup_ram_disk
from immunogenicity_scorer import score_dataframe
from population_coverage import calculate_population_coverage

# ─── Logging ──────────────────────────────────────────────────────────────────

LOG_DIR = Path(__file__).parent / "logs"
LOG_DIR.mkdir(parents=True, exist_ok=True)


def _get_logger(job_id: str) -> logging.Logger:
    logger = logging.getLogger(job_id)
    logger.setLevel(logging.DEBUG)
    if not logger.handlers:
        fh = logging.FileHandler(LOG_DIR / f"{job_id}.log")
        fh.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s] %(message)s"))
        logger.addHandler(fh)
        ch = logging.StreamHandler()
        ch.setFormatter(logging.Formatter("[%(levelname)s] %(message)s"))
        logger.addHandler(ch)
    return logger


# ─── Batch Job ─────────────────────────────────────────────────────────────────

class BatchJob:
    """
    Represents a batch prediction job.
    Manages a batch prediction job across multiple peptides / FASTA files.
    Auto-tuned for high-performance multi-core execution and RAM disk I/O.
    """

    def __init__(
        self,
        job_id: str,
        alleles: List[str],
        include_ba: bool = False,
        chunk_size: int = 500,
        progress_callback: Optional[Callable[[int, int], None]] = None,
    ):
        self.job_id = job_id
        self.alleles = alleles
        self.include_ba = include_ba
        self.chunk_size = chunk_size
        self.progress_callback = progress_callback
        self.logger = _get_logger(job_id)
        self.runner = NetMHCpanRunner()

        self.results: Optional[pd.DataFrame] = None
        self.status: str = "running"   # pending | running | done | error
        self.error_msg: Optional[str] = None
        self.created_at: str = datetime.utcnow().isoformat()
        self.completed_at: Optional[str] = None
        self.name: Optional[str] = None
        self.raw_input: Optional[str] = None
        self.mode: str = "peptide"
        self.lengths: List[int] = [9]
        self.meta: Dict[str, Any] = {}

    # ── Peptide batch ────────────────────────────────────────────────────────

    def run_peptides(self, peptides: List[str]) -> pd.DataFrame:
        """
        Run prediction on a list of peptide strings in chunks.
        """
        self.status = "running"
        self.meta["input_type"] = "peptides"
        self.meta["n_peptides"] = len(peptides)
        self.meta["n_alleles"] = len(self.alleles)

        cfg = get_optimal_batch_config(len(peptides), len(self.alleles))
        adaptive_chunk = cfg["chunk_size"]
        chunks = _chunk_list(peptides, adaptive_chunk)
        total = len(chunks)
        n_workers = min(total, cfg["allocated_workers"])

        self.meta["progress_current"] = 0
        self.meta["progress_total"] = total
        self.meta["progress_pct"] = 12
        self.meta["progress_stage"] = f"⚡ Allocated {n_workers} CPU Workers + RAM Disk for {len(peptides)} peptides..."
        self.meta["allocated_workers"] = n_workers
        self.meta["hardware_summary"] = cfg["summary"]

        self.logger.info(
            f"Job {self.job_id}: {len(peptides)} peptides × {len(self.alleles)} alleles "
            f"-> {total} chunks across {n_workers} auto-tuned CPU workers (RAM Disk: Active)"
        )

        from concurrent.futures import ThreadPoolExecutor, as_completed

        all_dfs = [None] * total
        completed_count = 0

        def _predict_chunk(idx_chunk):
            idx, chunk = idx_chunk
            runner = NetMHCpanRunner()
            df = runner.predict_peptides(
                chunk, alleles=self.alleles, include_ba=self.include_ba
            )
            return idx, df

        with ThreadPoolExecutor(max_workers=n_workers) as executor:
            future_to_chunk = {
                executor.submit(_predict_chunk, (i, chunk)): i
                for i, chunk in enumerate(chunks)
            }
            for future in as_completed(future_to_chunk):
                i = future_to_chunk[future]
                try:
                    idx, df = future.result()
                    all_dfs[idx] = df
                except Exception as e:
                    self.logger.error(f"  Chunk {i+1} failed: {e}")
                    raise
                completed_count += 1
                pct = int(12 + 76 * (completed_count / total))
                self.meta["progress_current"] = completed_count
                self.meta["progress_total"] = total
                self.meta["progress_pct"] = pct
                self.meta["progress_stage"] = f"⚡ Running NetMHCpan ({completed_count}/{total} chunks, {n_workers} Workers, {pct}%)..."
                if self.progress_callback:
                    self.progress_callback(completed_count, total)

        cleanup_ram_disk()
        self.meta["progress_pct"] = 90
        self.meta["progress_stage"] = "Aggregating peptide:MHC predictions & calculating metrics..."
        result = pd.concat([d for d in all_dfs if d is not None], ignore_index=True) if all_dfs else pd.DataFrame()
        result = result.drop_duplicates(subset=["Peptide", "MHC"]).reset_index(drop=True)
        self._finalize(result)
        return result

    # ── FASTA batch ──────────────────────────────────────────────────────────

    def run_fasta_files(
        self,
        fasta_files: List[str],
        lengths: Optional[List[int]] = None,
    ) -> pd.DataFrame:
        """
        Run prediction on multiple FASTA files.
        """
        self.status = "running"
        self.meta["input_type"] = "fasta"
        self.meta["n_files"] = len(fasta_files)
        lengths = lengths or [9]
        total = len(fasta_files)

        self.meta["progress_current"] = 0
        self.meta["progress_total"] = total
        self.meta["progress_pct"] = 10
        self.meta["progress_stage"] = f"Scanning {total} FASTA files for lengths {lengths}..."

        import multiprocessing
        from concurrent.futures import ThreadPoolExecutor, as_completed

        n_cores = multiprocessing.cpu_count()
        n_workers = min(total, max(1, n_cores // 2))
        all_dfs = [None] * total
        completed_count = 0

        def _predict_fasta_file(idx_path):
            idx, path = idx_path
            runner = NetMHCpanRunner()
            df = runner.predict_fasta(
                path,
                alleles=self.alleles,
                lengths=lengths,
                include_ba=self.include_ba,
            )
            df["source_file"] = Path(path).name
            return idx, df

        with ThreadPoolExecutor(max_workers=n_workers) as executor:
            future_to_file = {
                executor.submit(_predict_fasta_file, (i, p)): i
                for i, p in enumerate(fasta_files)
            }
            for future in as_completed(future_to_file):
                i = future_to_file[future]
                try:
                    idx, df = future.result()
                    all_dfs[idx] = df
                except Exception as e:
                    self.logger.error(f"  File {fasta_files[i]} failed: {e}")
                completed_count += 1
                pct = int(10 + 75 * (completed_count / total))
                self.meta["progress_current"] = completed_count
                self.meta["progress_total"] = total
                self.meta["progress_pct"] = pct
                self.meta["progress_stage"] = f"Processed {completed_count}/{total} FASTA files ({pct}%)..."
                if self.progress_callback:
                    self.progress_callback(completed_count, total)

        self.meta["progress_pct"] = 90
        self.meta["progress_stage"] = "Aggregating FASTA epitope predictions..."
        result = pd.concat([d for d in all_dfs if d is not None], ignore_index=True) if all_dfs else pd.DataFrame()
        self._finalize(result)
        return result

    # ── Finalize ─────────────────────────────────────────────────────────────

    def _finalize(self, df: pd.DataFrame):
        if not df.empty and "Peptide" in df.columns:
            df["Immunogenicity_Score"] = score_dataframe(df)

            # ── Multi-Dimensional Composed Index (Composite_Score) ─────────
            # 1. Presentation Likelihood S_pres (0.0 ~ 1.0)
            if "Score_EL" in df.columns and "Rank_EL" in df.columns:
                s_rank = np.clip(1.0 - (df["Rank_EL"] / 2.0), 0.0, 1.0)
                s_pres = np.clip(0.5 * df["Score_EL"] + 0.5 * s_rank, 0.0, 1.0)
            elif "Rank_EL" in df.columns:
                s_pres = np.clip(1.0 - (df["Rank_EL"] / 2.0), 0.0, 1.0)
            elif "Score_EL" in df.columns:
                s_pres = df["Score_EL"].clip(0.0, 1.0)
            else:
                s_pres = pd.Series(0.0, index=df.index)

            # 2. Binding Affinity S_aff (0.0 ~ 1.0, log10 IC50 scale)
            if "Affinity_nM" in df.columns and df["Affinity_nM"].notna().any():
                aff_clean = pd.to_numeric(df["Affinity_nM"], errors="coerce").fillna(50000.0).clip(lower=1.0, upper=50000.0)
                s_aff = np.clip(1.0 - (np.log10(aff_clean) / np.log10(50000.0)), 0.0, 1.0)
            else:
                s_aff = s_pres.copy()

            # 3. T-Cell Immunogenicity S_imm (0.0 ~ 1.0 based on Calis et al. 2013)
            if "Immunogenicity_Score" in df.columns:
                imm_clean = pd.to_numeric(df["Immunogenicity_Score"], errors="coerce").fillna(0.0)
                s_imm = np.clip((imm_clean + 0.4) / 0.8, 0.0, 1.0)
            else:
                s_imm = pd.Series(0.5, index=df.index)

            # 4. HLA Breadth S_breadth (0.0 ~ 1.0)
            n_target_alleles = max(1, len(self.alleles))
            if "Rank_EL" in df.columns and "MHC" in df.columns:
                binders_mask = df["Rank_EL"] <= 2.0
                allele_counts = df[binders_mask].groupby("Peptide")["MHC"].nunique().to_dict()
                s_breadth = df["Peptide"].map(lambda p: min(1.0, allele_counts.get(p, 0) / n_target_alleles))
            else:
                s_breadth = pd.Series(0.0, index=df.index)

            # Weighted Formula: 40% Pres + 25% Aff + 25% Imm + 10% Breadth
            df["Composite_Score"] = (
                0.40 * s_pres +
                0.25 * s_aff +
                0.25 * s_imm +
                0.10 * s_breadth
            ).round(4)

            # Decision Tier Assignment
            def _assign_tier(row):
                rank_el = row.get("Rank_EL", 999.0)
                imm_score = row.get("Immunogenicity_Score", 0.0)
                comp_score = row.get("Composite_Score", 0.0)
                if (rank_el <= 0.5 and (imm_score > 0.0 or pd.isna(imm_score))) or comp_score >= 0.65:
                    return "Tier 1"
                elif (rank_el <= 2.0 and (imm_score >= -0.15 or pd.isna(imm_score))) or comp_score >= 0.45:
                    return "Tier 2"
                else:
                    return "Tier 3"

            df["Tier"] = df.apply(_assign_tier, axis=1)

            # Ensure Identity column exists and fill nulls/blanks with 'PEPLIST'
            if "Identity" not in df.columns:
                df["Identity"] = "PEPLIST"
            else:
                df["Identity"] = df["Identity"].fillna("PEPLIST").replace("", "PEPLIST")

            # Reorder columns to place Identity (Protein ID / Identifier) as the first column
            desired_order = [
                "Identity", "Pos", "MHC", "Peptide", "Core",
                "Score_EL", "Rank_EL", "Affinity_nM",
                "Immunogenicity_Score", "Composite_Score", "Tier", "BindLevel"
            ]
            existing_cols = [c for c in desired_order if c in df.columns]
            other_cols = [c for c in df.columns if c not in existing_cols]
            df = df[existing_cols + other_cols]

            # Sort by Composite_Score descending by default
            df = df.sort_values(by=["Composite_Score", "Rank_EL"], ascending=[False, True]).reset_index(drop=True)

        self.results = df
        self.status = "done"
        self.completed_at = datetime.utcnow().isoformat()

        n_sb = int((df["Rank_EL"] <= 0.5).sum()) if not df.empty and "Rank_EL" in df.columns else 0
        n_wb = int(((df["Rank_EL"] > 0.5) & (df["Rank_EL"] <= 2.0)).sum()) if not df.empty and "Rank_EL" in df.columns else 0
        n_tier1 = int((df["Tier"] == "Tier 1").sum()) if not df.empty and "Tier" in df.columns else 0
        n_tier2 = int((df["Tier"] == "Tier 2").sum()) if not df.empty and "Tier" in df.columns else 0
        n_tier3 = int((df["Tier"] == "Tier 3").sum()) if not df.empty and "Tier" in df.columns else 0

        cov_world = calculate_population_coverage(self.alleles, "World")
        cov_asia = calculate_population_coverage(self.alleles, "East Asia")

        top_imm_pep, top_imm_mhc, top_imm_score, top_imm_identity, top_imm_pos = "—", "", 0.0, "", None
        if not df.empty and "Immunogenicity_Score" in df.columns and df["Immunogenicity_Score"].notna().any():
            top_row = df.loc[df["Immunogenicity_Score"].idxmax()]
            top_imm_pep = str(top_row["Peptide"])
            top_imm_mhc = str(top_row["MHC"])
            top_imm_score = float(top_row["Immunogenicity_Score"])
            if "Identity" in top_row:
                top_imm_identity = str(top_row["Identity"])
            if "Pos" in top_row:
                try:
                    top_imm_pos = int(top_row["Pos"])
                except Exception:
                    pass

        top_comp_pep, top_comp_mhc, top_comp_score, top_comp_identity, top_comp_pos, top_comp_tier = "—", "", 0.0, "", None, ""
        if not df.empty and "Composite_Score" in df.columns:
            top_comp_row = df.iloc[0]
            top_comp_pep = str(top_comp_row.get("Peptide", "—"))
            top_comp_mhc = str(top_comp_row.get("MHC", ""))
            top_comp_score = float(top_comp_row.get("Composite_Score", 0.0))
            top_comp_identity = str(top_comp_row.get("Identity", ""))
            top_comp_tier = str(top_comp_row.get("Tier", ""))
            if "Pos" in top_comp_row:
                try:
                    top_comp_pos = int(top_comp_row["Pos"])
                except Exception:
                    pass

        self.meta.update({
            "n_strong_binders": n_sb,
            "n_weak_binders": n_wb,
            "n_total": len(df),
            "n_tier1": n_tier1,
            "n_tier2": n_tier2,
            "n_tier3": n_tier3,
            "cov_world_pct": cov_world["cumulative_coverage_pct"],
            "cov_asia_pct": cov_asia["cumulative_coverage_pct"],
            "top_imm_pep": top_imm_pep,
            "top_imm_mhc": top_imm_mhc,
            "top_imm_score": top_imm_score,
            "top_imm_identity": top_imm_identity,
            "top_imm_pos": top_imm_pos,
            "top_comp_pep": top_comp_pep,
            "top_comp_mhc": top_comp_mhc,
            "top_comp_score": top_comp_score,
            "top_comp_identity": top_comp_identity,
            "top_comp_pos": top_comp_pos,
            "top_comp_tier": top_comp_tier,
        })
        self.logger.info(
            f"Job {self.job_id} done: {len(df)} predictions, "
            f"{n_sb} SB, {n_wb} WB"
        )

    def to_summary_dict(self) -> Dict[str, Any]:
        return {
            "job_id": self.job_id,
            "name": self.name or self.meta.get("job_name") or f"Job #{self.job_id}",
            "status": self.status,
            "mode": self.mode,
            "alleles": self.alleles,
            "lengths": self.lengths,
            "include_ba": self.include_ba,
            "raw_input": self.raw_input,
            "created_at": self.created_at,
            "completed_at": self.completed_at,
            "meta": self.meta,
            "error": self.error_msg,
        }

    # ── Export ───────────────────────────────────────────────────────────────

    def save_csv(self, path: str) -> str:
        """Save results to CSV. Returns the path."""
        if self.results is None:
            raise ValueError("No results yet. Run the job first.")
        self.results.to_csv(path, index=False)
        self.logger.info(f"Saved CSV: {path}")
        return path

    def save_excel(self, path: str) -> str:
        """Save results to Excel with summary sheet. Returns the path."""
        if self.results is None:
            raise ValueError("No results yet. Run the job first.")
        with pd.ExcelWriter(path, engine="openpyxl") as writer:
            self.results.to_excel(writer, sheet_name="Predictions", index=False)
            summary_df = pd.DataFrame([{
                "Job ID": self.job_id,
                "Total Predictions": self.meta.get("n_total", 0),
                "Strong Binders (≤0.5%)": self.meta.get("n_strong_binders", 0),
                "Weak Binders (0.5-2%)": self.meta.get("n_weak_binders", 0),
                "Alleles": ", ".join(self.alleles),
                "Completed At": self.completed_at,
            }])
            summary_df.to_excel(writer, sheet_name="Summary", index=False)
        self.logger.info(f"Saved Excel: {path}")
        return path


# ─── Helpers ──────────────────────────────────────────────────────────────────

def _chunk_list(lst: List, size: int) -> List[List]:
    """Split a list into chunks of at most `size` elements."""
    return [lst[i:i + size] for i in range(0, len(lst), size)]


def parse_peptide_file(path: str) -> List[str]:
    """
    Parse a plain-text file of peptides (one per line).
    Lines starting with '#' are ignored. Optional score column is dropped.
    """
    peptides = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            # Take only the first column (ignore score column if present)
            peptides.append(line.split()[0])
    return list(dict.fromkeys(peptides))  # deduplicate preserving order


# ─── Quick self-test ──────────────────────────────────────────────────────────

if __name__ == "__main__":
    test_pep = Path("/home/cylin/NETMHC/netMHCpan-4.2/test/test.pep")
    peptides = parse_peptide_file(str(test_pep))
    print(f"Loaded {len(peptides)} peptides from {test_pep}")

    job = BatchJob(
        job_id="selftest_001",
        alleles=["HLA-A02:01", "HLA-B07:02"],
        include_ba=False,
        chunk_size=5,
    )
    df = job.run_peptides(peptides)
    print("\n=== Results ===")
    print(df[["Peptide", "MHC", "Score_EL", "Rank_EL", "BindLevel"]].to_string(index=False))
    print("\n=== Summary ===")
    print(json.dumps(job.to_summary_dict(), indent=2))
