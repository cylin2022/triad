"""
app.py — Flask Web Application for NetMHCpan-4.2
"""

import os
import re
import json
import uuid
import shutil
import tarfile
import threading
from pathlib import Path
from datetime import datetime
from typing import Optional

from flask import (
    Flask, render_template, request, jsonify,
    send_file, abort, url_for, redirect
)

import pandas as pd
import netmhcpan_api
from netmhcpan_api import NetMHCpanRunner, NETMHC_BASE
from batch_runner import BatchJob
from visualizer import build_dashboard_html
from oss_binding_predictors import predict_smm, predict_comblib

# ─── App Setup ────────────────────────────────────────────────────────────────

PLATFORM_VERSION = "v1.1"

BASE_DIR   = Path(__file__).parent
RESULT_DIR = BASE_DIR / "results"
SETUP_TMP_DIR = BASE_DIR / "uploads" / "_setup"   # created on demand by /setup/upload
RESULT_DIR.mkdir(exist_ok=True)

app = Flask(__name__)
app.secret_key = os.urandom(24)
app.config["MAX_CONTENT_LENGTH"] = 300 * 1024 * 1024  # 300 MB, netMHCpan pkg is ~120MB

# In-memory job store (production: use Redis / DB)
JOBS: dict[str, BatchJob] = {}
JOBS_LOCK = threading.Lock()

# ─── netMHCpan install status ──────────────────────────────────────────────────
# The Docker image ships without the licensed netMHCpan-4.2 binary/data (its
# license forbids redistribution). Users provide it themselves — either by
# mounting a volume at NETMHC_HOME, or via the /setup wizard below — so the
# runner must be created lazily instead of crashing the app at import time.

_runner = None

def is_netmhc_installed() -> bool:
    return netmhcpan_api.is_installed()

def get_runner() -> NetMHCpanRunner:
    global _runner
    if _runner is None:
        _runner = NetMHCpanRunner()
    return _runner

# Load allele list once at startup
ALLELE_LIST_PATH = BASE_DIR / "static" / "alleles.json"
with open(ALLELE_LIST_PATH) as f:
    ALL_ALLELES = json.load(f)

# ─── Helpers ──────────────────────────────────────────────────────────────────

VALID_AA = re.compile(r"^[ACDEFGHIKLMNPQRSTVWY]+$")

def _validate_peptides(raw: str):
    """Return (valid_list, error_msg)."""
    lines = [l.strip() for l in raw.strip().splitlines() if l.strip()]
    peptides, errors = [], []
    for i, line in enumerate(lines, 1):
        pep = line.split()[0].upper()
        if not (8 <= len(pep) <= 14):
            errors.append(f"Line {i}: '{pep}' has invalid length ({len(pep)}) — must be 8–14 aa")
            continue
        if not VALID_AA.match(pep):
            errors.append(f"Line {i}: '{pep}' contains invalid amino acids")
            continue
        peptides.append(pep)
    return peptides, errors

def _clean_df_for_json(df_subset):
    """Convert DataFrame to dict records replacing float NaN/Inf with None for valid JSON."""
    import numpy as np
    records = df_subset.to_dict(orient="records")
    for r in records:
        for k, v in r.items():
            if isinstance(v, float) and (np.isnan(v) or np.isinf(v)):
                r[k] = None
    return records

# ─── Routes ───────────────────────────────────────────────────────────────────

@app.route("/")
def index():
    if not is_netmhc_installed():
        return redirect(url_for("setup_page"))
    return render_template("index.html")


@app.route("/api/alleles")
def api_alleles():
    # The allele list stores names without the '*' separator (HLA-A02:01), but
    # WHO nomenclature — and netMHCpan's own output — writes HLA-A*02:01, so
    # search has to match either form.
    q = request.args.get("q", "").lower().replace("*", "")
    if q:
        filtered = [a for a in ALL_ALLELES if q in a.lower().replace("*", "")][:50]
    else:
        filtered = ALL_ALLELES[:50]
    return jsonify(filtered)


@app.route("/api/demo")
def api_demo():
    """Run an instant demonstration prediction with standard benchmark peptides."""
    if not is_netmhc_installed():
        return jsonify({"error": "netMHCpan-4.2 is not installed yet.", "setup_url": url_for("setup_page")}), 503
    demo_peptides = [
        "AAAWYLWEV", "AEFGPWQTV", "YLLPAIVHI",
        "FLPSDFFPSV", "GILGFVFTL", "SLYNTVATLY",
        "RMFPNAPYL", "LLDFVRMGV", "AAGLQDCTM"
    ]
    demo_alleles = ["HLA-A02:01", "HLA-B07:02"]

    job_id = "demo_" + str(uuid.uuid4())[:6]
    job = BatchJob(job_id=job_id, alleles=demo_alleles, include_ba=True)
    with JOBS_LOCK:
        JOBS[job_id] = job

    try:
        df = job.run_peptides(demo_peptides)
        csv_path   = str(RESULT_DIR / f"{job_id}.csv")
        excel_path = str(RESULT_DIR / f"{job_id}.xlsx")
        job.save_csv(csv_path)
        job.save_excel(excel_path)

        resp = job.to_summary_dict()
        charts = build_dashboard_html(df, job.meta)
        resp["charts"]     = charts
        resp["table_data"] = _clean_df_for_json(df)
        resp["columns"]    = list(df.columns)
        resp["csv_url"]    = url_for("download_result", job_id=job_id, fmt="csv")
        resp["excel_url"]  = url_for("download_result", job_id=job_id, fmt="xlsx")
        return jsonify(resp)
    except Exception as e:
        return jsonify({"error": str(e)}), 500


@app.route("/api/predict", methods=["POST"])
def api_predict():
    """
    Body (JSON):
      {
        "mode": "peptide" | "fasta",
        "input": "<peptide lines or FASTA text>",
        "alleles": ["HLA-A02:01", ...],
        "include_ba": false,
        "lengths": [9]        # only for fasta mode
      }
    """
    if not is_netmhc_installed():
        return jsonify({"error": "netMHCpan-4.2 is not installed yet.", "setup_url": url_for("setup_page")}), 503

    data = request.get_json(force=True)
    mode       = data.get("mode", "peptide")
    raw_input  = data.get("input", "").strip()
    alleles    = data.get("alleles", ["HLA-A02:01"])
    include_ba = bool(data.get("include_ba", False))
    lengths    = data.get("lengths", [9])

    if not raw_input:
        return jsonify({"error": "No input provided"}), 400
    if not alleles:
        return jsonify({"error": "No alleles selected"}), 400

    # Normalize allele names (HLA-A*02:01 → HLA-A02:01, A02:01 → HLA-A02:01)
    norm_alleles = []
    for a in alleles:
        a_clean = a.strip().replace("*", "")
        if not a_clean.startswith("HLA-") and not a_clean.startswith("BoLA-") and not a_clean.startswith("SLA-") and not a_clean.startswith("H-2-") and not a_clean.startswith("H2-") and not a_clean.startswith("Gogo-") and not a_clean.startswith("Mamu-") and not a_clean.startswith("Patr-") and not a_clean.startswith("DLA-") and not a_clean.startswith("Eqca-"):
            a_clean = "HLA-" + a_clean
        norm_alleles.append(a_clean)
    # De-duplicate: the same allele can arrive in several notations
    # (e.g. "HLA-A*02:01" and "HLA-A02:01" both normalize to "HLA-A02:01").
    # Without this, FASTA runs emit duplicated rows and double the work.
    alleles = list(dict.fromkeys(norm_alleles))

    # Validate alleles (basic check)
    invalid_alleles = [a for a in alleles if a not in ALL_ALLELES]
    if invalid_alleles:
        return jsonify({"error": f"Unknown alleles: {invalid_alleles[:3]}"}), 400

    # Auto-detect FASTA mode if input contains FASTA headers ('>')
    if raw_input.lstrip().startswith(">") or "\n>" in raw_input:
        mode = "fasta"

    # Create job
    job_id = str(uuid.uuid4())[:8]
    job = BatchJob(job_id=job_id, alleles=alleles, include_ba=include_ba)
    job.raw_input = raw_input
    job.mode = mode
    job.lengths = lengths
    job.name = data.get("name") or f"{mode.upper()} ({len(alleles)} HLA)"
    with JOBS_LOCK:
        job.status = "running"
        job.meta["progress_pct"] = 10
        job.meta["progress_stage"] = "Initializing computation engine & validating inputs..."
        job.meta["job_name"] = job.name
        JOBS[job_id] = job

    def _run():
        try:
            job.status = "running"
            run_mode = job.mode
            if run_mode == "peptide" and (raw_input.lstrip().startswith(">") or "\n>" in raw_input):
                run_mode = "fasta"

            if run_mode == "peptide":
                peptides, errors = _validate_peptides(raw_input)
                if not peptides:
                    job.status = "error"
                    job.error_msg = "; ".join(errors[:5])
                    return
                job.meta["progress_pct"] = 15
                job.meta["progress_stage"] = f"Queuing {len(peptides)} peptides for NetMHCpan-4.2 inference..."
                job.run_peptides(peptides)
            else:  # fasta
                job.meta["input_type"] = "fasta"
                job.meta["progress_pct"] = 25
                job.meta["progress_stage"] = "Scanning FASTA sequence with sliding window..."
                df = get_runner().predict_fasta(
                    raw_input,
                    alleles=alleles,
                    lengths=lengths,
                    include_ba=include_ba,
                )
                job._finalize(df)

            job.meta["progress_pct"] = 92
            job.meta["progress_stage"] = "Generating export reports (CSV & Excel)..."
            # Save CSV & Excel
            csv_path   = str(RESULT_DIR / f"{job_id}.csv")
            excel_path = str(RESULT_DIR / f"{job_id}.xlsx")
            job.save_csv(csv_path)
            job.save_excel(excel_path)
            job.meta["progress_pct"] = 100
            job.meta["progress_stage"] = "Completed successfully!"

        except Exception as e:
            job.status = "error"
            job.error_msg = str(e)

    t = threading.Thread(target=_run, daemon=True)
    t.start()
    return jsonify({"job_id": job_id, "status": "running"})


@app.route("/api/predict_oss", methods=["POST"])
def api_predict_oss():
    """
    Open-source-licensed alternative binding predictors (no DTU dependency):
    SMM, SMM-PMBEC, and comblib_sidney2008 — extracted from IEDB's ng_tc1.

    Body (JSON):
      {
        "peptides": ["GILGFVFTL", ...],
        "allele": "HLA-A*02:01",
        "methods": ["smm", "smmpmbec", "comblib"]   # any subset
      }
    """
    data = request.get_json(force=True) or {}
    peptides = [p.strip().upper() for p in data.get("peptides", []) if p.strip()]
    allele = (data.get("allele") or "").strip()
    methods = data.get("methods") or ["smm", "comblib"]

    if not peptides:
        return jsonify({"error": "No peptides provided"}), 400
    if not allele:
        return jsonify({"error": "No allele provided"}), 400

    rows = []
    for pep in peptides:
        row = {"Peptide": pep, "Allele": allele}
        if "smm" in methods:
            row["SMM_IC50_nM"] = predict_smm([pep], allele, method="smm")[0]
        if "smmpmbec" in methods:
            row["SMMPMBEC_IC50_nM"] = predict_smm([pep], allele, method="smmpmbec")[0]
        if "comblib" in methods:
            row["Comblib_Score"] = predict_comblib([pep], allele)[0]
        rows.append(row)

    return jsonify({"allele": allele, "methods": methods, "results": rows})


@app.route("/api/jobs", methods=["GET"])
def api_all_jobs():
    """Return summary list of all submitted background jobs."""
    with JOBS_LOCK:
        job_list = [j.to_summary_dict() for j in list(JOBS.values())[::-1]]
    return jsonify(job_list)


@app.route("/api/jobs", methods=["DELETE"])
def api_delete_all_jobs():
    """Clear all completed/failed jobs and clean up files."""
    with JOBS_LOCK:
        job_ids = list(JOBS.keys())
        JOBS.clear()

    for j_id in job_ids:
        for ext in ["csv", "xlsx"]:
            p = RESULT_DIR / f"{j_id}.{ext}"
            if p.exists():
                try: p.unlink()
                except Exception: pass
        log_p = Path(__file__).parent / "logs" / f"{j_id}.log"
        if log_p.exists():
            try: log_p.unlink()
            except Exception: pass

    return jsonify({"success": True, "cleared_count": len(job_ids)})


@app.route("/api/job/<job_id>", methods=["GET"])
def api_job_status(job_id):
    with JOBS_LOCK:
        job = JOBS.get(job_id)
    if job is None:
        csv_file = RESULT_DIR / f"{job_id}.csv"
        if csv_file.exists():
            try:
                df = pd.read_csv(csv_file)
                alleles = list(df["MHC"].unique()) if "MHC" in df.columns else []
                job = BatchJob(job_id=job_id, alleles=alleles)
                job.results = df
                job.status = "done"
                with JOBS_LOCK:
                    JOBS[job_id] = job
            except Exception:
                return jsonify({"error": "Job not found"}), 404
        else:
            return jsonify({"error": "Job not found"}), 404

    resp = job.to_summary_dict()
    if job.status == "done" and job.results is not None:
        df = job.results
        # Build charts
        charts = build_dashboard_html(df, job.meta)
        # Convert dataframe to JSON (limit 500 rows in table)
        table_data = _clean_df_for_json(df.head(500))
        resp["charts"]     = charts
        resp["table_data"] = table_data
        resp["columns"]    = list(df.columns)
        resp["csv_url"]    = url_for("download_result", job_id=job_id, fmt="csv")
        resp["excel_url"]  = url_for("download_result", job_id=job_id, fmt="xlsx")
    return jsonify(resp)


@app.route("/api/job/<job_id>", methods=["DELETE"])
def api_delete_job(job_id):
    """Delete a specific job and its generated output files."""
    with JOBS_LOCK:
        job = JOBS.pop(job_id, None)

    # Clean up generated files
    for ext in ["csv", "xlsx"]:
        p = RESULT_DIR / f"{job_id}.{ext}"
        if p.exists():
            try: p.unlink()
            except Exception: pass
    log_p = Path(__file__).parent / "logs" / f"{job_id}.log"
    if log_p.exists():
        try: log_p.unlink()
        except Exception: pass

    if job is None:
        return jsonify({"error": "Job not found"}), 404
    return jsonify({"success": True, "job_id": job_id})


@app.route("/api/job/<job_id>/rename", methods=["POST"])
def api_rename_job(job_id):
    """Rename or add notes to a job."""
    data = request.get_json(silent=True) or {}
    new_name = data.get("name", "").strip()
    if not new_name:
        return jsonify({"error": "Name cannot be empty"}), 400

    with JOBS_LOCK:
        job = JOBS.get(job_id)
        if job is None:
            return jsonify({"error": "Job not found"}), 404
        job.name = new_name
        job.meta["job_name"] = new_name
        summary = job.to_summary_dict()

    return jsonify({"success": True, "job": summary})


@app.route("/api/download/<job_id>/<fmt>")
def download_result(job_id, fmt):
    if fmt not in ("csv", "xlsx"):
        abort(400)
    path = RESULT_DIR / f"{job_id}.{fmt}"
    if not path.exists():
        abort(404)
    mime = "text/csv" if fmt == "csv" else (
        "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
    )
    return send_file(str(path), mimetype=mime,
                     as_attachment=True,
                     download_name=f"netmhcpan_{job_id}.{fmt}")


from resource_tuner import get_system_specs, cleanup_ram_disk

@app.route("/api/system_resources")
def api_system_resources():
    """Return real-time hardware specifications, CPU load, RAM disk, and available memory."""
    specs = get_system_specs()
    return jsonify(specs)


# ─── Setup Wizard (bring-your-own netMHCpan-4.2) ───────────────────────────────
# netMHCpan-4.2's academic license forbids redistribution to third parties, so
# the Docker image ships without it. Users download their own copy from DTU
# under their own license acceptance, then install it here — either by
# mounting a volume at NETMHC_HOME before starting the container, or by
# uploading the tar.gz they downloaded through this wizard.

DTU_DOWNLOAD_URL = "https://services.healthtech.dtu.dk/services/NetMHCpan-4.2/"
ALLOWED_TARBALL_SUFFIXES = (".tar.gz", ".tgz")


@app.route("/api/install_status")
def api_install_status():
    return jsonify({
        "installed": is_netmhc_installed(),
        "netmhc_home": str(NETMHC_BASE),
        "platform_version": PLATFORM_VERSION,
        "download_url": DTU_DOWNLOAD_URL,
    })


@app.route("/setup")
def setup_page():
    return render_template(
        "setup.html",
        installed=is_netmhc_installed(),
        netmhc_home=str(NETMHC_BASE),
        download_url=DTU_DOWNLOAD_URL,
        platform_version=PLATFORM_VERSION,
    )


def _safe_extract(tar: tarfile.TarFile, dest: Path):
    """Extract a tar archive, rejecting absolute paths and path traversal."""
    dest = dest.resolve()
    for member in tar.getmembers():
        member_path = (dest / member.name).resolve()
        if not str(member_path).startswith(str(dest) + os.sep) and member_path != dest:
            raise ValueError(f"Unsafe path in archive: {member.name}")
    if hasattr(tarfile, "data_filter"):
        tar.extractall(dest, filter="data")
    else:
        tar.extractall(dest)


def _find_package_root(staging_dir: Path) -> Optional[Path]:
    """Locate the extracted directory that directly contains Linux_x86_64/bin/netMHCpan-4.2."""
    for candidate in [staging_dir, *staging_dir.iterdir()] if staging_dir.exists() else []:
        if candidate.is_dir() and (candidate / "Linux_x86_64" / "bin" / "netMHCpan-4.2").is_file():
            return candidate
    return None


def _selftest_package(pkg_root: Path) -> tuple[bool, str]:
    """Run a minimal prediction against the staged package to confirm it works before install."""
    import tempfile as _tempfile
    import subprocess as _subprocess

    binary   = pkg_root / "Linux_x86_64" / "bin" / "netMHCpan-4.2"
    rdir     = pkg_root / "Linux_x86_64"
    data     = rdir / "data"
    required = [binary, data / "synlist_nocontext.bin", data / "MHC_pseudo.dat",
                data / "version", data / "allelenames"]
    missing = [str(p) for p in required if not p.exists()]
    if missing:
        return False, f"Package is missing expected files: {missing}"
    if not os.access(binary, os.X_OK):
        try:
            binary.chmod(0o755)
        except Exception:
            return False, f"Binary is not executable and could not be chmod'ed: {binary}"

    with _tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        pep_file = tmp / "selftest.pep"
        pep_file.write_text("GILGFVFTL\n")
        cmd = [
            str(binary), "-p", str(pep_file),
            "-rdir", str(rdir),
            "-syn", str(data / "synlist_nocontext.bin"),
            "-hlapseudo", str(data / "MHC_pseudo.dat"),
            "-tdir", str(tmp / "netMHCpan_XXXXXX"),
            "-version", str(data / "version"),
            "-thrfmt", str(data / "threshold" / "%s.thr.%s"),
            "-allname", str(data / "allelenames"),
            "-a", "HLA-A02:01",
        ]
        try:
            result = _subprocess.run(cmd, capture_output=True, text=True, timeout=60)
        except Exception as e:
            return False, f"Self-test failed to run the binary: {e}"
        if result.returncode != 0:
            return False, f"Self-test exited with code {result.returncode}: {result.stderr[-500:]}"
        if "GILGFVFTL" not in result.stdout:
            return False, "Self-test ran but produced no recognizable prediction output."
    return True, "ok"


@app.route("/setup/upload", methods=["POST"])
def setup_upload():
    if is_netmhc_installed():
        return jsonify({"error": "netMHCpan-4.2 is already installed."}), 400
    if request.form.get("agree_license") != "on":
        return jsonify({"error": "You must confirm you have your own DTU license for netMHCpan."}), 400
    if "file" not in request.files or not request.files["file"].filename:
        return jsonify({"error": "No file uploaded."}), 400

    f = request.files["file"]
    filename = f.filename
    if not filename.lower().endswith(ALLOWED_TARBALL_SUFFIXES):
        return jsonify({"error": "Expected a .tar.gz package as downloaded from DTU."}), 400

    SETUP_TMP_DIR.mkdir(parents=True, exist_ok=True)
    session_dir = SETUP_TMP_DIR / uuid.uuid4().hex
    session_dir.mkdir(parents=True)
    tar_path = session_dir / filename
    staging_dir = session_dir / "staged"

    try:
        f.save(str(tar_path))

        if not tarfile.is_tarfile(tar_path):
            return jsonify({"error": "Uploaded file is not a valid tar.gz archive."}), 400

        staging_dir.mkdir()
        with tarfile.open(tar_path, "r:gz") as tar:
            _safe_extract(tar, staging_dir)

        pkg_root = _find_package_root(staging_dir)
        if pkg_root is None:
            return jsonify({
                "error": "Could not find Linux_x86_64/bin/netMHCpan-4.2 inside the archive. "
                         "Make sure this is the netMHCpan-4.2 Linux package from DTU."
            }), 400

        ok, msg = _selftest_package(pkg_root)
        if not ok:
            return jsonify({"error": f"Self-test failed, install aborted: {msg}"}), 400

        NETMHC_BASE.mkdir(parents=True, exist_ok=True)
        for item in pkg_root.iterdir():
            dest = NETMHC_BASE / item.name
            if dest.exists():
                if dest.is_dir():
                    shutil.rmtree(dest)
                else:
                    dest.unlink()
            shutil.move(str(item), str(dest))

        global _runner
        _runner = None  # force re-init against the newly installed binary

        return jsonify({"success": True, "message": "netMHCpan-4.2 installed successfully.",
                         "redirect": url_for("index")})
    except Exception as e:
        return jsonify({"error": f"Installation failed: {e}"}), 500
    finally:
        shutil.rmtree(session_dir, ignore_errors=True)


# ─── Entry point ──────────────────────────────────────────────────────────────

if __name__ == "__main__":
    specs = get_system_specs()
    port = int(os.environ.get("FLASK_PORT", 5001))
    print("=" * 65)
    print(f"🧬 TRIAD Platform {PLATFORM_VERSION}")
    print(f"👉 Local Access URL: http://localhost:{port}  (or http://127.0.0.1:{port})")
    print(f"   Note: On Windows, use localhost:{port}, do NOT use 0.0.0.0:{port}")
    print(f"🖥️  CPU: {specs['cpu']['model']} ({specs['cpu']['logical_cores']} vCPUs / {specs['cpu']['physical_cores']} Cores)")
    print(f"💾 Memory: {specs['memory']['available_gb']} GB Available / {specs['memory']['total_gb']} GB Total")
    print(f"⚡ RAM Disk: {specs['ram_disk']['path']} ({specs['ram_disk']['available_gb']} GB Free / Active)")
    print(f"🚀 Auto-Tuned Parallel Pool: {specs['cpu']['optimal_workers']} Max Concurrent Workers")
    if is_netmhc_installed():
        print(f"✅ netMHCpan-4.2 detected at {NETMHC_BASE}")
    else:
        print(f"⚠️  netMHCpan-4.2 NOT installed at {NETMHC_BASE} — visit http://localhost:{port}/setup to install it")
    print("=" * 65)
    app.run(debug=False, host="0.0.0.0", port=port, threaded=True, use_reloader=False)
