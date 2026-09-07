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
from tqdm import tqdm

from netmhcpan_api import NetMHCpanRunner
from resource_tuner import get_optimal_batch_config, cleanup_ram_disk

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

    Parameters
    ----------
    job_id : str
        Unique identifier for this job.
    alleles : list of str
        HLA alleles to predict against.
    include_ba : bool
        Whether to include binding affinity prediction.
    chunk_size : int
        Number of peptides per chunk (to avoid memory issues).
    progress_callback : callable, optional
        Called with (current, total) on each chunk completion.
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

        self.meta["allocated_workers"] = n_workers
        self.meta["hardware_summary"] = cfg["summary"]
        self.meta["progress_pct"] = 10
        self.meta["progress_stage"] = f"⚡ Initializing {n_workers} CPU Workers across {len(self.alleles)} HLA alleles..."
        self.logger.info(
            f"Job {self.job_id}: {len(peptides)} peptides × "
            f"{len(self.alleles)} alleles (chunk_size={adaptive_chunk}, workers={n_workers})"
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
                pct = int(10 + 75 * (completed_count / total))
                self.meta["progress_current"] = completed_count
                self.meta["progress_total"] = total
                self.meta["progress_pct"] = pct
                self.meta["progress_stage"] = f"⚡ Running NetMHCpan ({completed_count}/{total} chunks, {n_workers} Workers, {pct}%)..."
                if self.progress_callback:
                    self.progress_callback(completed_count, total)

        cleanup_ram_disk()
        self.meta["progress_pct"] = 88
        self.meta["progress_stage"] = "Aggregating multi-modal epitope predictions & scores..."
        result = pd.concat([d for d in all_dfs if d is not None], ignore_index=True) if all_dfs else pd.DataFrame()
        result = result.drop_duplicates(subset=["Peptide", "MHC"]).reset_index(drop=True)
        return self._finalize(result)

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
        import multiprocessing
        from concurrent.futures import ThreadPoolExecutor, as_completed

        n_cores = multiprocessing.cpu_count()
        n_workers = min(total, max(4, n_cores // 2))
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
                self.meta["progress_current"] = completed_count
                self.meta["progress_total"] = total
                if self.progress_callback:
                    self.progress_callback(completed_count, total)

        result = pd.concat([d for d in all_dfs if d is not None], ignore_index=True) if all_dfs else pd.DataFrame()
        return self._finalize(result)

    # ── Finalize ─────────────────────────────────────────────────────────────

    def _finalize(self, df: pd.DataFrame) -> pd.DataFrame:
        from iedb_tc1_api import IedbTc1Runner
        from population_coverage import calculate_population_coverage
        from hla_breadth import calculate_hla_breadth
        from antigen_processing import calculate_antigen_processing
        from immunostruct_api import ImmunoStructRunner, evaluate_concordance

        if not df.empty and "Peptide" in df.columns and "MHC" in df.columns:
            # 1. Compute T-cell immunogenicity for candidate binders (Rank_EL <= 5.0)
            if "Rank_EL" in df.columns:
                binders_df = df[df["Rank_EL"] <= 5.0]
                if binders_df.empty:
                    binders_df = df.head(100)
            else:
                binders_df = df

            if not binders_df.empty:
                try:
                    tc1 = IedbTc1Runner()
                    unique_peps = binders_df["Peptide"].unique().tolist()
                    unique_mhcs = binders_df["MHC"].unique().tolist()
                    imm_df = tc1.predict_immunogenicity(unique_peps, unique_mhcs)
                    if not imm_df.empty:
                        # Normalize MHC strings (e.g. HLA-A*02:01 <-> HLA-A02:01) for robust merge
                        imm_df["MHC_norm"] = imm_df["MHC"].astype(str).str.replace("*", "", regex=False)
                        df["MHC_norm"] = df["MHC"].astype(str).str.replace("*", "", regex=False)
                        df = pd.merge(df, imm_df[["Peptide", "MHC_norm", "Immunogenicity_Score"]], on=["Peptide", "MHC_norm"], how="left")
                        df = df.drop(columns=["MHC_norm"])
                except Exception as e:
                    self.logger.warning(f"IEDB TC1 integration warning: {e}")

            if "Immunogenicity_Score" not in df.columns:
                df["Immunogenicity_Score"] = np.nan

            # 2. Compute ImmunoStruct 3D Multimodal Immunogenicity Score
            try:
                struct_runner = ImmunoStructRunner()
                unique_peps = df["Peptide"].unique().tolist()
                unique_mhcs = df["MHC"].unique().tolist()
                struct_df = struct_runner.predict_pMHC(unique_peps, unique_mhcs)
                if not struct_df.empty:
                    df = pd.merge(df, struct_df[["Peptide", "MHC", "ImmunoStruct_Score"]], on=["Peptide", "MHC"], how="left")
            except Exception as e:
                self.logger.warning(f"ImmunoStruct integration warning: {e}")

            if "ImmunoStruct_Score" not in df.columns:
                df["ImmunoStruct_Score"] = 0.50

            # 3. Compute HLA Breadth metrics
            try:
                breadth_df = calculate_hla_breadth(df)
                if not breadth_df.empty:
                    df = pd.merge(df, breadth_df, on="Peptide", how="left")
            except Exception as e:
                self.logger.warning(f"HLA breadth calculation warning: {e}")

            if "HLA_Breadth_Score" not in df.columns:
                df["HLA_Breadth_Score"] = 0.0
                df["HLA_Breadth_Category"] = "Low"

            # 4. Compute Antigen Processing Score
            try:
                proc_df = calculate_antigen_processing(df["Peptide"].unique().tolist(), df)
                if not proc_df.empty:
                    df = pd.merge(df, proc_df[["Peptide", "Processing_Score"]], on="Peptide", how="left")
            except Exception as e:
                self.logger.warning(f"Antigen processing calculation warning: {e}")

            if "Processing_Score" not in df.columns:
                df["Processing_Score"] = df["Score_EL"] if "Score_EL" in df.columns else 0.0

        # 5. Compute Population Coverage
        cov_world = calculate_population_coverage(self.alleles, "World")
        cov_asia  = calculate_population_coverage(self.alleles, "East Asia")
        world_cov_frac = cov_world["cumulative_coverage_pct"] / 100.0

        # 6. 6-Dimension Multimodal Composite Score (0.0 to 1.0) & Tier System
        if not df.empty:
            s_pres = df["Score_EL"] if "Score_EL" in df.columns else 0.0
            s_proc = df["Processing_Score"] if "Processing_Score" in df.columns else 0.0
            imm_raw = df["Immunogenicity_Score"].fillna(0.0)
            s_iedb = np.clip((imm_raw + 0.5) / 1.0, 0.0, 1.0)
            s_struct = df["ImmunoStruct_Score"].fillna(0.50)
            s_breadth = df["HLA_Breadth_Score"].fillna(0.0)

            # 6-Dimension Formula from immunStruct.md
            df["Composite_Score"] = (
                0.30 * s_pres +
                0.10 * s_proc +
                0.15 * s_iedb +
                0.25 * s_struct +
                0.10 * s_breadth +
                0.10 * world_cov_frac
            ).round(4)

            # Assign Concordance Category (Dual-Layer Analysis)
            df["Concordance"] = df.apply(
                lambda r: evaluate_concordance(r.get("Immunogenicity_Score"), r.get("ImmunoStruct_Score", 0.5)), axis=1
            )

            # Assign Tier Classification
            def _assign_tier(row):
                rank_el = row.get("Rank_EL", 999.0)
                imm_score = row.get("Immunogenicity_Score", -1.0)
                struct_score = row.get("ImmunoStruct_Score", 0.50)
                comp_score = row.get("Composite_Score", 0.0)

                if (rank_el < 0.5 and (struct_score >= 0.55 or (pd.notnull(imm_score) and imm_score > 0))) or comp_score >= 0.55:
                    return "Tier 1"
                elif rank_el <= 2.0 and (struct_score >= 0.45 or (pd.isnull(imm_score) or imm_score >= -0.15)):
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
                "Processing_Score", "Immunogenicity_Score", "ImmunoStruct_Score",
                "Concordance", "Composite_Score", "Tier", "BindLevel"
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

        n_struct_driven = int((df["Concordance"] == "💡 Structure-Driven").sum()) if not df.empty and "Concordance" in df.columns else 0
        n_double_val = int((df["Concordance"] == "🌟 Double Validated").sum()) if not df.empty and "Concordance" in df.columns else 0

        top_tier1_row = None
        top_tier1_identity = ""
        top_tier1_pos = None
        if not df.empty and "Tier" in df.columns:
            t1_df = df[df["Tier"] == "Tier 1"]
            if not t1_df.empty:
                top_tier1_row = t1_df.sort_values(by="Composite_Score", ascending=False).iloc[0]
            elif not df.empty:
                top_tier1_row = df.sort_values(by="Composite_Score", ascending=False).iloc[0]

        if top_tier1_row is not None:
            if "Identity" in top_tier1_row:
                top_tier1_identity = str(top_tier1_row["Identity"])
            if "Pos" in top_tier1_row:
                try:
                    top_tier1_pos = int(top_tier1_row["Pos"])
                except Exception:
                    pass

        self.meta.update({
            "n_strong_binders": n_sb,
            "n_weak_binders": n_wb,
            "n_total": len(df),
            "n_tier1": n_tier1,
            "n_tier2": n_tier2,
            "n_tier3": n_tier3,
            "n_struct_driven": n_struct_driven,
            "n_double_val": n_double_val,
            "cov_world_pct": cov_world["cumulative_coverage_pct"],
            "cov_asia_pct": cov_asia["cumulative_coverage_pct"],
            "top_tier1_pep": str(top_tier1_row["Peptide"]) if top_tier1_row is not None else "—",
            "top_tier1_mhc": str(top_tier1_row["MHC"]) if top_tier1_row is not None else "—",
            "top_tier1_score": float(top_tier1_row["Composite_Score"]) if top_tier1_row is not None else 0.0,
            "top_tier1_identity": top_tier1_identity,
            "top_tier1_pos": top_tier1_pos,
        })
        self.logger.info(
            f"Job {self.job_id} done: {len(df)} predictions, "
            f"{n_tier1} Tier 1, {n_struct_driven} Structure-Driven candidates"
        )
        return df

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
