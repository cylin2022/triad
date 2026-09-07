"""
app.py — Flask Web Application for NetMHCpan-4.2
"""

import os
import re
import json
import uuid
import threading
from pathlib import Path
from datetime import datetime

from flask import (
    Flask, render_template, request, jsonify,
    send_file, abort, url_for
)

from netmhcpan_api import NetMHCpanRunner
from batch_runner import BatchJob, parse_peptide_file
from visualizer import build_dashboard_html

# ─── App Setup ────────────────────────────────────────────────────────────────

BASE_DIR   = Path(__file__).parent
UPLOAD_DIR = BASE_DIR / "uploads"
RESULT_DIR = BASE_DIR / "results"
UPLOAD_DIR.mkdir(exist_ok=True)
RESULT_DIR.mkdir(exist_ok=True)

app = Flask(__name__)
app.secret_key = os.urandom(24)

# In-memory job store (production: use Redis / DB)
JOBS: dict[str, BatchJob] = {}
JOBS_LOCK = threading.Lock()

runner = NetMHCpanRunner()

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

def _fasta_to_peptides_only(text: str):
    """Extract only sequence lines from FASTA (removes headers)."""
    return text  # pass-through: fasta handled by runner

def _allowed_file(filename: str) -> bool:
    return "." in filename and filename.rsplit(".", 1)[1].lower() in {
        "fsa", "fasta", "fa", "pep", "txt"
    }

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
    return render_template("index.html")


@app.route("/api/alleles")
def api_alleles():
    q = request.args.get("q", "").lower()
    if q:
        filtered = [a for a in ALL_ALLELES if q in a.lower()][:50]
    else:
        filtered = ALL_ALLELES[:50]
    return jsonify(filtered)


@app.route("/api/demo")
def api_demo():
    """Run an instant demonstration prediction with standard benchmark peptides."""
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

    # Normalize allele names (HLA-A*02:01 → HLA-A02:01)
    alleles = [a.replace("*", "") for a in alleles]

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
                df = runner.predict_fasta(
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


@app.route("/api/coverage", methods=["POST"])
def api_coverage():
    """Calculate population coverage for selected HLA alleles."""
    from population_coverage import calculate_population_coverage
    data = request.get_json(force=True) or {}
    alleles = data.get("alleles", [])
    region = data.get("region", "World")
    if not alleles:
        return jsonify({"error": "No alleles specified"}), 400
    res = calculate_population_coverage(alleles, region)
    return jsonify(res)


@app.route("/api/upload", methods=["POST"])
def api_upload():
    """Upload a peptide / FASTA file, return its text content."""
    if "file" not in request.files:
        return jsonify({"error": "No file part"}), 400
    f = request.files["file"]
    if not f.filename or not _allowed_file(f.filename):
        return jsonify({"error": "Invalid file type"}), 400

    save_path = UPLOAD_DIR / f"{uuid.uuid4().hex}_{f.filename}"
    f.save(str(save_path))

    text = save_path.read_text(errors="replace")
    return jsonify({"filename": f.filename, "content": text})


from resource_tuner import get_system_specs, cleanup_ram_disk

@app.route("/api/system_resources")
def api_system_resources():
    """Return real-time hardware specifications, CPU load, RAM disk, and available memory."""
    specs = get_system_specs()
    return jsonify(specs)


# ─── Entry point ──────────────────────────────────────────────────────────────

if __name__ == "__main__":
    specs = get_system_specs()
    print("=" * 65)
    print("🧬 Next-Gen Immunoinformatics Portal starting on http://0.0.0.0:5002")
    print(f"🖥️  CPU: {specs['cpu']['model']} ({specs['cpu']['logical_cores']} vCPUs / {specs['cpu']['physical_cores']} Cores)")
    print(f"💾 Memory: {specs['memory']['available_gb']} GB Available / {specs['memory']['total_gb']} GB Total")
    print(f"⚡ RAM Disk: {specs['ram_disk']['path']} ({specs['ram_disk']['available_gb']} GB Free / Active)")
    print(f"🚀 Auto-Tuned Parallel Pool: {specs['cpu']['optimal_workers']} Max Concurrent Workers")
    print("=" * 65)
    app.run(debug=False, host="0.0.0.0", port=5002, threaded=True, use_reloader=False)
