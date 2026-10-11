"""Download a complete, auditable Kb/Db mouse-host T-cell assay snapshot."""
import argparse
import hashlib
import json
import time
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

BASE = "https://query-api.iedb.org/tcell_search"
FIELDS = "tcell_id,structure_id,linear_sequence,structure_type,e_modification,reference_id,pubmed_id,reference_dates,reference_titles,qualitative_measure,mhc_class,mhc_allele_name,mhc_allele_evidence,host_organism_name,host_organism_iri,host_organism_iri_search,assay_names,assay_iris,immunization_description,direct_ex_vivo_bool,source_organism_name,parent_source_antigen_iri,parent_source_antigen_name,curated_source_antigen"


def fetch(outdir):
    outdir.mkdir(parents=True, exist_ok=True)
    target = outdir / "iedb_mouse_kb_db.jsonl"
    if target.exists():
        raise FileExistsError(f"Snapshot already exists: {target}; use a new directory")
    params = {
        "mhc_allele_name": "in.(H2-Kb,H2-Db)",
        "host_organism_iri_search": "cs.{NCBITaxon:10090}",
        "select": FIELDS, "order": "tcell_id.asc", "limit": "1000",
    }
    manifest = {"started_utc": datetime.now(timezone.utc).isoformat(),
                "endpoint": BASE, "filters": params.copy(), "pages": []}
    count, last_id, expected = 0, 0, None
    partial = target.with_suffix(".partial")
    with partial.open("w") as fh:
        while True:
            query = dict(params, tcell_id=f"gt.{last_id}")
            url = BASE + "?" + urllib.parse.urlencode(query)
            for attempt in range(4):
                try:
                    req = urllib.request.Request(url, headers={"Prefer": "count=exact", "User-Agent": "TRIAD-mouse-benchmark/1.0"})
                    with urllib.request.urlopen(req, timeout=90) as response:
                        content_range = response.headers.get("Content-Range", "")
                        rows = json.load(response)
                    break
                except Exception:
                    if attempt == 3:
                        raise
                    time.sleep(2 ** attempt)
            if expected is None:
                expected = int(content_range.split("/")[-1])
            manifest["pages"].append({"url": url, "content_range": content_range, "rows": len(rows)})
            if not rows:
                break
            ids = [r["tcell_id"] for r in rows]
            if ids != sorted(set(ids)) or ids[0] <= last_id:
                raise ValueError("Non-monotonic or duplicate assay IDs")
            for row in rows:
                fh.write(json.dumps(row, ensure_ascii=False) + "\n")
            last_id = ids[-1]
            count += len(rows)
            print(f"Downloaded {count}/{expected} assays", flush=True)
    if count != expected:
        raise ValueError(f"Snapshot changed or download incomplete: {count} != {expected}")
    partial.rename(target)
    manifest.update(rows=count, finished_utc=datetime.now(timezone.utc).isoformat(),
                    sha256=hashlib.sha256(target.read_bytes()).hexdigest())
    (outdir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, default=Path(__file__).parent / "data")
    fetch(parser.parse_args().outdir)
