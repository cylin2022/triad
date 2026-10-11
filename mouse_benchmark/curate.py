"""Deterministic curation and leakage-group inventory; see PROTOCOL.md."""
import collections
import hashlib
import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent
AA = set("ACDEFGHIKLMNPQRSTVWY")
HOSTS = {"Mus musculus C57BL/6", "Mus musculus C57BL/6J", "Mus musculus C57BL/6N"}
POSITIVE = {"Positive", "Positive-Low", "Positive-Intermediate", "Positive-High"}


def rejection(row):
    if row.get("host_organism_name") not in HOSTS:
        return "host_not_explicit_standard_b6"
    seq = row.get("linear_sequence") or ""
    if row.get("structure_type") != "Linear peptide" or not (8 <= len(seq) <= 11) or not set(seq) <= AA:
        return "not_canonical_linear_8_11mer"
    if row.get("e_modification"):
        return "modified_peptide"
    if not str(row.get("assay_names", "")).startswith("IFNg release|"):
        return "not_IFNg_functional_assay"
    if row.get("qualitative_measure") not in POSITIVE | {"Negative"}:
        return "unknown_outcome"
    if not row.get("reference_id"):
        return "missing_reference"
    return "eligible"


def strong_evidence(value):
    return value in {"MHC binding assay", "Single allele present"} or str(value).startswith("T cell assay -")


def years(row):
    return [int(str(v)[:4]) for v in row.get("reference_dates") or [] if str(v)[:4].isdigit()]


def group_components(df):
    """Union studies and exact/one-substitution sequences, even across alleles."""
    parent = list(range(len(df)))

    def root(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    def union(a, b):
        parent[root(a)] = root(b)

    seen = {}
    for i, row in enumerate(df.itertuples()):
        seq = row.Peptide
        keys = [("ref", r) for r in row.Reference_IDs.split(";")]
        keys += [("seq", seq)]
        keys += [("neighbor", len(seq), p, seq[:p] + "*" + seq[p+1:]) for p in range(len(seq))]
        for key in keys:
            if key in seen:
                union(i, seen[key])
            else:
                seen[key] = i
    mapping = {}
    result = []
    for i in range(len(df)):
        r = root(i)
        if r not in mapping:
            mapping[r] = len(mapping)
        result.append(mapping[r])
    return result


def curate():
    data = ROOT / "data"
    out = ROOT / "results"
    out.mkdir(exist_ok=True)
    raw = data / "iedb_mouse_kb_db.jsonl"
    manifest = json.loads((data / "manifest.json").read_text())
    assert hashlib.sha256(raw.read_bytes()).hexdigest() == manifest["sha256"]
    rows = [json.loads(line) for line in raw.open()]
    assert len(rows) == manifest["rows"] == len({r["tcell_id"] for r in rows})
    assert all(r["mhc_allele_name"] in {"H2-Kb", "H2-Db"} and "NCBITaxon:10090" in r["host_organism_iri_search"] for r in rows)
    audit = [{"Assay_ID": r["tcell_id"], "Reason": rejection(r)} for r in rows]
    pd.DataFrame(audit).to_csv(out / "assay_eligibility.csv", index=False)
    eligible = [r for r in rows if rejection(r) == "eligible"]
    groups = collections.defaultdict(list)
    for r in eligible:
        groups[(r["linear_sequence"], r["mhc_allele_name"])].append(r)
    clean, conflicts = [], []
    for (seq, allele), records in sorted(groups.items()):
        labels = {int(r["qualitative_measure"] in POSITIVE) for r in records}
        ref_ids = sorted({r["reference_id"] for r in records})
        yy = [y for r in records for y in years(r)]
        item = {
            "Peptide": seq, "MHC": allele.replace("H2-", "H-2-"), "Length": len(seq),
            "Assay_Count": len(records), "Reference_IDs": ";".join(map(str, ref_ids)),
            "PubMed_IDs": ";".join(sorted({str(r["pubmed_id"]) for r in records if r["pubmed_id"]})),
            "Assay_IDs": ";".join(str(r["tcell_id"]) for r in records),
            "Earliest_Year": min(yy) if yy else None,
            "Latest_Year": max(yy) if yy else None,
            "All_Dates_Known": all(years(r) for r in records),
            "Strong_Restriction": any(strong_evidence(r["mhc_allele_evidence"]) for r in records),
            "Restriction_Evidence": ";".join(sorted({str(r["mhc_allele_evidence"]) for r in records})),
            "Source_Organisms": ";".join(sorted({str(r["source_organism_name"]) for r in records})),
            "Source_Antigens": ";".join(sorted({str(r["parent_source_antigen_iri"]) for r in records})),
        }
        if len(labels) != 1:
            conflicts.append(item)
        else:
            item["Label"] = labels.pop()
            clean.append(item)
    df = pd.DataFrame(clean)
    df["Leakage_Group"] = group_components(df)
    df.to_csv(out / "curated_pairs.csv", index=False)
    pd.DataFrame(conflicts).to_csv(out / "conflicting_pairs.csv", index=False)
    studies = []
    for ref in sorted({r["reference_id"] for r in eligible}):
        rr = [r for r in eligible if r["reference_id"] == ref]
        subset = df[df.Reference_IDs.map(lambda s: str(ref) in s.split(";"))]
        studies.append({"Reference_ID": ref, "Title": rr[0]["reference_titles"],
                        "PubMed_ID": rr[0]["pubmed_id"], "Dates": rr[0]["reference_dates"],
                        "Eligible_Assays": len(rr), "Curated_Pairs": len(subset),
                        "Positive": int(subset.Label.sum()), "Negative": int((subset.Label == 0).sum())})
    pd.DataFrame(studies).sort_values("Curated_Pairs", ascending=False).to_csv(out / "studies.csv", index=False)
    inventory = {
        "raw_assays": len(rows), "raw_unique_pairs": len({(r['linear_sequence'], r['mhc_allele_name']) for r in rows}),
        "raw_references": len({r["reference_id"] for r in rows}),
        "sequential_exclusions": dict(collections.Counter(a["Reason"] for a in audit)),
        "eligible_assays": len(eligible), "eligible_pairs_before_conflicts": len(groups),
        "conflicting_pairs_excluded": len(conflicts), "curated_pairs": len(df),
        "positive": int(df.Label.sum()), "negative": int((df.Label == 0).sum()),
        "eligible_references_before_conflicts": len(studies),
        "curated_references": sum(s["Curated_Pairs"] > 0 for s in studies),
        "unique_peptides": df.Peptide.nunique(),
        "by_allele": df.groupby("MHC").Label.agg(["count", "sum"]).to_dict("index"),
        "raw_hosts": dict(collections.Counter(r["host_organism_name"] for r in rows)),
        "raw_assays_by_type": dict(collections.Counter(r["assay_names"] for r in rows)),
        "eligible_restriction_evidence": dict(collections.Counter(r["mhc_allele_evidence"] for r in eligible)),
        "eligible_source_organisms": dict(collections.Counter(str(r["source_organism_name"]) for r in eligible)),
        "leakage_groups": df.Leakage_Group.nunique(),
        "largest_group_pairs": int(df.Leakage_Group.value_counts().iloc[0]),
        "group_sizes": {str(k): int(v) for k,v in df.Leakage_Group.value_counts().items()},
    }
    (out / "inventory.json").write_text(json.dumps(inventory, indent=2) + "\n")
    print(json.dumps({k:v for k,v in inventory.items() if not isinstance(v,dict)}, indent=2))


if __name__ == "__main__":
    curate()
