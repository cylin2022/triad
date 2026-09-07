"""
netmhcpan_api.py
Python API wrapper for NetMHCpan-4.2 binary.
Directly calls the compiled binary without tcsh dependency.
"""

import os
import re
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import List, Optional, Union
import pandas as pd

from resource_tuner import get_ram_disk_tmpdir

# ─── Configuration ────────────────────────────────────────────────────────────
# NETMHC_HOME lets the Docker image ship without the licensed netMHCpan-4.2
# package: the container points this at an empty mount point (or a directory
# populated by the /setup install wizard in app.py) and defaults to the
# historical on-host path for native, non-container installs.

NETMHC_BASE    = Path(os.environ.get("NETMHC_HOME", "/home/cylin/NETMHC/netMHCpan-4.2"))
BINARY         = NETMHC_BASE / "Linux_x86_64/bin/netMHCpan-4.2"
RDIR           = NETMHC_BASE / "Linux_x86_64"
SYN_NOCONTEXT  = RDIR / "data/synlist_nocontext.bin"
SYN_CONTEXT    = RDIR / "data/synlist_context.bin"
SYN_IEDB       = RDIR / "data/synlist_iedb.bin"
SYN_CEDAR      = RDIR / "data/synlist_cedar.bin"
HLA_PSEUDO     = RDIR / "data/MHC_pseudo.dat"
VERSION_FILE   = RDIR / "data/version"
THR_FMT        = str(RDIR / "data/threshold/%s.thr.%s")
ALLNAMES       = RDIR / "data/allelenames"
TMPDIR         = get_ram_disk_tmpdir()

# ─── Result Parser ─────────────────────────────────────────────────────────────

HEADER_RE = re.compile(
    r"^\s*(\d+)\s+"            # Pos
    r"(HLA-\S+)\s+"           # MHC
    r"(\S+)\s+"               # Peptide
    r"(\S+)\s+"               # Core
    r"(\d+)\s+(\d+)\s+"       # Of, Gp
    r"(\d+)\s+(\d+)\s+"       # Gl, Ip
    r"(\d+)\s+"               # Il
    r"(\S+)\s+"               # Icore
    r"(\S+)\s+"               # Identity
    r"([\d.]+)\s+"            # Score_EL
    r"([\d.]+)"               # %Rank_EL
    r"(?:\s+([\d.]+))?"       # Exp (optional)
    r"(?:\s+<=\s+(\S+))?"     # BindLevel (optional)
)

BA_RE = re.compile(
    r"^\s*(\d+)\s+"
    r"(HLA-\S+)\s+"
    r"(\S+)\s+"
    r"(\S+)\s+"
    r"(\d+)\s+(\d+)\s+"
    r"(\d+)\s+(\d+)\s+"
    r"(\d+)\s+"
    r"(\S+)\s+"
    r"(\S+)\s+"
    r"([\d.]+)\s+"
    r"([\d.]+)\s+"
    r"([\d.]+)\s+"
    r"([\d.]+)"
    r"(?:\s+([\d.]+))?"
    r"(?:\s+<=\s+(\S+))?"
)


def is_installed() -> bool:
    """Check whether the licensed netMHCpan-4.2 package is present under NETMHC_HOME."""
    return all(p.exists() for p in (
        BINARY, SYN_NOCONTEXT, SYN_CONTEXT, HLA_PSEUDO, VERSION_FILE, ALLNAMES,
    ))


def parse_output(raw: str, include_ba: bool = False) -> pd.DataFrame:
    """Parse raw netMHCpan output text into a DataFrame."""
    rows = []
    current_allele = None

    for line in raw.splitlines():
        # Capture current allele
        allele_match = re.match(r"^# Allele:\s+(\S+)", line)
        if allele_match:
            current_allele = allele_match.group(1)
            continue

        # Skip headers / separators
        if line.startswith("#") or line.startswith("-") or line.startswith(" Pos"):
            continue

        stripped = line.strip()
        if not stripped:
            continue

        # Try to parse data line
        if include_ba:
            m = BA_RE.match(line)
            if m:
                g = m.groups()
                rows.append({
                    "Identity": g[10] if g[10] else "PEPLIST",
                    "Pos": int(g[0]),
                    "MHC": g[1],
                    "Peptide": g[2],
                    "Core": g[3],
                    # netMHCpan-4.2 -BA column order is:
                    #   Score_EL  %Rank_EL  Score_BA  %Rank_BA  Aff(nM)  [BindLevel]
                    "Score_EL": float(g[11]),
                    "Rank_EL": float(g[12]),
                    "Score_BA": float(g[13]),
                    "Rank_BA": float(g[14]),
                    "Affinity_nM": float(g[15]) if g[15] else None,
                    "BindLevel": g[16] if g[16] else "",
                })
        else:
            m = HEADER_RE.match(line)
            if m:
                g = m.groups()
                rows.append({
                    "Identity": g[10] if g[10] else "PEPLIST",
                    "Pos": int(g[0]),
                    "MHC": g[1],
                    "Peptide": g[2],
                    "Core": g[3],
                    "Score_EL": float(g[11]),
                    "Rank_EL": float(g[12]),
                    "Exp": float(g[13]) if g[13] else None,
                    "BindLevel": g[14] if g[14] else "",
                })

    df = pd.DataFrame(rows)
    if not df.empty:
        cols = ["Identity", "Pos", "MHC", "Peptide", "Core"] + [c for c in df.columns if c not in ["Identity", "Pos", "MHC", "Peptide", "Core"]]
        df = df[cols]
        if "Rank_EL" in df.columns:
            df = df.sort_values("Rank_EL").reset_index(drop=True)
    return df


# ─── Main Runner ───────────────────────────────────────────────────────────────

class NetMHCpanRunner:
    """
    Python wrapper for the NetMHCpan-4.2 binary.

    Usage
    -----
    runner = NetMHCpanRunner()

    # From a peptide list
    df = runner.predict_peptides(["AAAWYLWEV", "AEFGPWQTV"], alleles=["HLA-A02:01"])

    # From a FASTA string or file
    df = runner.predict_fasta(fasta_str, alleles=["HLA-A02:01", "HLA-B07:02"])
    """

    def __init__(
        self,
        binary: Path = BINARY,
        default_alleles: Optional[List[str]] = None,
        rank_strong: float = 0.5,
        rank_weak: float = 2.0,
        timeout: int = 300,
    ):
        self.binary = str(binary)
        if not os.path.isfile(self.binary):
            raise FileNotFoundError(f"NetMHCpan binary not found: {self.binary}")

        self.default_alleles = default_alleles or ["HLA-A02:01"]
        self.rank_strong = rank_strong
        self.rank_weak = rank_weak
        self.timeout = timeout

    # ── internal ──────────────────────────────────────────────────────────────

    def _build_cmd(
        self,
        input_file: str,
        input_type_flag: str,
        alleles: List[str],
        include_ba: bool,
        include_pathogen: bool,
        include_neo: bool,
        use_context: bool,
        extra_flags: Optional[List[str]] = None,
    ) -> List[str]:
        syn = str(SYN_CONTEXT if use_context else SYN_NOCONTEXT)

        cmd = [
            self.binary,
            input_type_flag, input_file,
            "-rdir", str(RDIR),
            "-syn", syn,
            "-hlapseudo", str(HLA_PSEUDO),
            "-tdir", str(TMPDIR / "netMHCpan_XXXXXX"),
            "-version", str(VERSION_FILE),
            "-thrfmt", THR_FMT,
            "-allname", str(ALLNAMES),
            "-a", ",".join(alleles),
            "-rankS", str(self.rank_strong),
            "-rankW", str(self.rank_weak),
        ]
        if include_ba:
            cmd.append("-BA")
        if include_pathogen:
            cmd.append("-pathogen")
        if include_neo:
            cmd.append("-neo")
        if use_context:
            cmd.append("-context")
        if extra_flags:
            cmd.extend(extra_flags)
        return cmd

    def _run(self, cmd: List[str]) -> str:
        """Execute the binary and return stdout."""
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=self.timeout,
        )
        if result.returncode != 0:
            raise RuntimeError(
                f"NetMHCpan returned exit code {result.returncode}.\n"
                f"STDERR:\n{result.stderr}"
            )
        return result.stdout

    # ── public API ────────────────────────────────────────────────────────────

    def predict_peptides(
        self,
        peptides: List[str],
        alleles: Optional[List[str]] = None,
        include_ba: bool = False,
        include_pathogen: bool = False,
        include_neo: bool = False,
    ) -> pd.DataFrame:
        """
        Predict binding for a list of peptide strings.

        Parameters
        ----------
        peptides : list of str
            Amino acid sequences (e.g. ['AAAWYLWEV', 'AEFGPWQTV']).
        alleles : list of str, optional
            HLA alleles (e.g. ['HLA-A02:01']).  Defaults to self.default_alleles.
        include_ba : bool
            If True, also predict binding affinity (BA) values.

        Returns
        -------
        pd.DataFrame
        """
        alleles = alleles or self.default_alleles
        TMPDIR.mkdir(parents=True, exist_ok=True)
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".pep", dir=str(TMPDIR), delete=False
        ) as f:
            f.write("\n".join(peptides) + "\n")
            pep_file = f.name
        try:
            cmd = self._build_cmd(
                pep_file, "-p", alleles, include_ba,
                include_pathogen, include_neo, False
            )
            raw = self._run(cmd)
            return parse_output(raw, include_ba=include_ba)
        finally:
            os.unlink(pep_file)

    def predict_fasta(
        self,
        fasta: Union[str, Path],
        alleles: Optional[List[str]] = None,
        lengths: Optional[List[int]] = None,
        include_ba: bool = False,
        include_pathogen: bool = False,
        include_neo: bool = False,
        use_context: bool = False,
    ) -> pd.DataFrame:
        """
        Predict binding from a FASTA string or file path.

        Parameters
        ----------
        fasta : str or Path
            Either a FASTA-formatted string or a path to a FASTA file.
        lengths : list of int, optional
            Peptide lengths to enumerate (default: [9]).
        """
        alleles = alleles or self.default_alleles
        lengths = lengths or [9]

        # Write to temp file if it's a string
        if isinstance(fasta, str) and not os.path.isfile(fasta):
            TMPDIR.mkdir(parents=True, exist_ok=True)
            with tempfile.NamedTemporaryFile(
                mode="w", suffix=".fsa", dir=str(TMPDIR), delete=False
            ) as f:
                f.write(fasta)
                fasta_file = f.name
            own_file = True
        else:
            fasta_file = str(fasta)
            own_file = False

        extra = ["-l", ",".join(map(str, lengths))]
        try:
            cmd = self._build_cmd(
                fasta_file, "-f", alleles, include_ba,
                include_pathogen, include_neo, use_context, extra
            )
            raw = self._run(cmd)
            return parse_output(raw, include_ba=include_ba)
        finally:
            if own_file:
                os.unlink(fasta_file)

    def predict_pmhc(
        self,
        pmhc_lines: List[str],
        include_ba: bool = False,
    ) -> pd.DataFrame:
        """
        Predict binding from peptide-MHC pairs.

        Parameters
        ----------
        pmhc_lines : list of str
            Each line: "<peptide> <allele>", e.g. "AAAWYLWEV HLA-A02:01"
        """
        TMPDIR.mkdir(parents=True, exist_ok=True)
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".pmhc", dir=str(TMPDIR), delete=False
        ) as f:
            f.write("\n".join(pmhc_lines) + "\n")
            pmhc_file = f.name

        # For pMHC, alleles are embedded in the file; use a dummy allele
        try:
            cmd = [
                self.binary,
                "-pmhc", pmhc_file,
                "-rdir", str(RDIR),
                "-syn", str(SYN_NOCONTEXT),
                "-hlapseudo", str(HLA_PSEUDO),
                "-tdir", str(TMPDIR / "netMHCpan_XXXXXX"),
                "-version", str(VERSION_FILE),
                "-thrfmt", THR_FMT,
                "-allname", str(ALLNAMES),
            ]
            if include_ba:
                cmd.append("-BA")
            raw = self._run(cmd)
            return parse_output(raw, include_ba=include_ba)
        finally:
            os.unlink(pmhc_file)

    def list_alleles(self) -> List[str]:
        """Return the list of all supported HLA alleles."""
        cmd = [
            self.binary,
            "-listMHC",
            "-rdir", str(RDIR),
            "-allname", str(ALLNAMES),
            "-version", str(VERSION_FILE),
            "-thrfmt", THR_FMT,
            "-hlapseudo", str(HLA_PSEUDO),
            "-syn", str(SYN_NOCONTEXT),
            "-tdir", str(TMPDIR / "netMHCpan_XXXXXX"),
        ]
        raw = self._run(cmd)
        alleles = []
        for line in raw.splitlines():
            line = line.strip()
            if line and not line.startswith("#") and not line.startswith("-"):
                alleles.append(line)
        return alleles

    def get_strong_binders(self, df: pd.DataFrame) -> pd.DataFrame:
        """Filter DataFrame to strong binders (%Rank_EL ≤ rank_strong)."""
        return df[df["Rank_EL"] <= self.rank_strong].copy()

    def get_weak_binders(self, df: pd.DataFrame) -> pd.DataFrame:
        """Filter DataFrame to weak binders (rank_strong < %Rank_EL ≤ rank_weak)."""
        return df[
            (df["Rank_EL"] > self.rank_strong) & (df["Rank_EL"] <= self.rank_weak)
        ].copy()


# ─── Quick self-test ──────────────────────────────────────────────────────────

if __name__ == "__main__":
    runner = NetMHCpanRunner()
    test_peptides = [
        "AAAWYLWEV", "AAGLQDCTM", "AARNIVRRA",
        "AARPDDPTL", "AASCGGAVF", "AASKQQMLM",
        "AASSTHRKV", "AEALLADGL", "AEESLSLEA", "AEFGPWQTV",
    ]
    print("=== NetMHCpan API Self-Test ===")
    df = runner.predict_peptides(test_peptides, alleles=["HLA-A02:01"])
    print(df.to_string(index=False))
    print(f"\nStrong binders: {len(runner.get_strong_binders(df))}")
    print(f"Weak binders  : {len(runner.get_weak_binders(df))}")
