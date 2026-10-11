"""
epitope_set.py — v1.2

Epitope-SET population coverage: the share of people for whom a chosen set of
peptides actually produces a presented epitope.

Why this module exists
──────────────────────
Up to v1.1 the platform reported `calculate_population_coverage(job.alleles)` —
the coverage of the ALLELE PANEL the user typed into the search box. That number
is fixed the moment the panel is chosen and does not move when the predictions
come back, so a job in which nothing bound still displayed ~95% coverage. For
vaccine design the question is the other way round:

    Given the N peptides I picked, what fraction of the population has at least
    one HLA allele that presents at least one of them?

That depends on the (peptide x allele) binder matrix, not on the panel.

The model
─────────
Individuals are drawn under Hardy-Weinberg equilibrium with the A / B / C loci
independent — the same assumptions as population_coverage, whose
`absence_probability()` primitive is reused here so the two can never disagree.

Write R for the set of "responsive" alleles (those presenting >= 1 selected
epitope) and P0(T) = absence_probability(T) = P(an individual carries no allele
in T). Then:

    coverage(>=1 epitope) = 1 - P0(R)

For >= 2 distinct epitopes (the redundancy a real vaccine needs, so that one
escape mutation or one down-regulated allele does not leave a person with no
response at all) the counting is exact rather than simulated. For epitope e let

    B_e = alleles presenting at least one epitope OTHER than e

Then "e is the only epitope this individual sees" = (carries >=1 allele of R)
AND (carries no allele of B_e). Because the alleles of R that are not in B_e are
exactly those presenting e alone, that event has probability P0(B_e) - P0(R), so

    P(exactly one epitope) = sum_e [ P0(B_e) - P0(R) ]
    coverage(>=2 epitopes) = 1 - P0(R) - P(exactly one epitope)

which costs one absence_probability() call per epitope — cheap enough to re-run
on every checkbox click, and exact, so no Monte-Carlo noise makes the number
jitter while the user is choosing.
"""

from typing import Any, Dict, Iterable, List, Optional, Sequence

from population_coverage import (
    absence_probability,
    allele_frequency,
    is_human,
    locus_of,
    normalize_allele,
    resolve_region,
)

# A peptide counts as presented by an allele at or below this %Rank_EL.
# 2.0 is netMHCpan's own weak-binder threshold; 0.5 is the strong-binder one.
DEFAULT_RANK_CUT = 2.0
STRONG_RANK_CUT = 0.5


# ─── Binder matrix ────────────────────────────────────────────────────────────

def build_binder_map(
    rows: Iterable[Dict[str, Any]],
    rank_cut: float = DEFAULT_RANK_CUT,
    peptides: Optional[Sequence[str]] = None,
) -> Dict[str, List[str]]:
    """
    Collapse prediction rows into {peptide: [alleles it is presented by]}.

    `rows` are the per-(peptide, allele) records the results table is built from;
    `peptides` optionally restricts the map to a chosen subset, in which case
    peptides with no binder at all are still returned (with an empty list) so the
    caller can tell "not selected" apart from "selected but presented by nobody".
    """
    wanted = set(peptides) if peptides is not None else None
    binder_map: Dict[str, List[str]] = {p: [] for p in (peptides or [])}

    for r in rows:
        pep = r.get("Peptide")
        mhc = r.get("MHC")
        if not pep or not mhc:
            continue
        if wanted is not None and pep not in wanted:
            continue
        try:
            rank = float(r.get("Rank_EL"))
        except (TypeError, ValueError):
            continue
        if rank > rank_cut:
            continue
        lst = binder_map.setdefault(pep, [])
        if mhc not in lst:
            lst.append(mhc)

    return binder_map


def responsive_alleles(binder_map: Dict[str, List[str]]) -> List[str]:
    """Every allele presenting at least one peptide of the set, de-duplicated."""
    seen: Dict[str, None] = {}
    for alleles in binder_map.values():
        for a in alleles:
            seen.setdefault(a, None)
    return list(seen)


# ─── Set coverage ─────────────────────────────────────────────────────────────

def set_coverage(
    binder_map: Dict[str, List[str]],
    region: str = "World",
    panel: Optional[Sequence[str]] = None,
) -> Dict[str, Any]:
    """
    Population coverage of an epitope SET.

    Returns >=1-hit and >=2-hit coverage, plus the alleles of `panel` that the
    set fails to reach (ranked by frequency, i.e. by how much fixing each one
    would be worth).
    """
    resolved = resolve_region(region)
    R = responsive_alleles(binder_map)

    p0 = absence_probability(R, resolved)
    cov1 = 1.0 - p0

    # P(exactly one distinct epitope presented) — see module docstring.
    p_exactly_one = 0.0
    for e in binder_map:
        others = [a for other, alleles in binder_map.items() if other != e for a in alleles]
        p_exactly_one += absence_probability(others, resolved) - p0
    p_exactly_one = max(0.0, p_exactly_one)

    cov2 = max(0.0, min(cov1, 1.0 - p0 - p_exactly_one))

    # Which of the target panel's alleles does the set never reach?
    uncovered = []
    if panel:
        hit = {normalize_allele(a) for a in R}
        for a in dict.fromkeys(panel):
            if normalize_allele(a) in hit:
                continue
            uncovered.append({
                "allele": a,
                "locus": locus_of(a),
                "frequency": round(allele_frequency(a, resolved), 4),
                "human": is_human(a),
            })
        uncovered.sort(key=lambda d: -d["frequency"])

    # Per-allele epitope counts, for the set heatmap and for spotting alleles
    # that are covered by exactly one peptide (a single point of failure).
    per_allele: Dict[str, int] = {}
    for alleles in binder_map.values():
        for a in alleles:
            per_allele[a] = per_allele.get(a, 0) + 1

    return {
        "region": resolved,
        "n_epitopes": len(binder_map),
        "coverage_pct": round(cov1 * 100, 2),
        "coverage2_pct": round(cov2 * 100, 2),
        "n_responsive_alleles": len(R),
        "responsive_alleles": R,
        "uncovered_alleles": uncovered,
        "epitopes_per_allele": per_allele,
        "single_point_alleles": sorted([a for a, n in per_allele.items() if n == 1]),
    }


def marginal_gains(
    selected: Dict[str, List[str]],
    candidates: Dict[str, List[str]],
    region: str = "World",
    objective: str = "coverage_pct",
) -> Dict[str, float]:
    """
    Percentage points each candidate peptide would add to the current set.

    This is what turns a ranking table into a design tool: a Tier 1 peptide
    restricted to an allele the set already covers has a marginal gain near
    zero, while a Tier 2 peptide on an untouched allele can be worth several
    points. Candidates already in `selected` get 0.0.
    """
    resolved = resolve_region(region)
    base = set_coverage(selected, resolved)[objective]

    gains: Dict[str, float] = {}
    for pep, alleles in candidates.items():
        if pep in selected:
            gains[pep] = 0.0
            continue
        trial = dict(selected)
        trial[pep] = alleles
        gains[pep] = round(set_coverage(trial, resolved)[objective] - base, 3)
    return gains


def greedy_select(
    candidates: Dict[str, List[str]],
    region: str = "World",
    goal_pct: float = 90.0,
    max_epitopes: int = 20,
    objective: str = "coverage_pct",
    seed: Optional[Sequence[str]] = None,
    tie_break: Optional[Dict[str, float]] = None,
) -> Dict[str, Any]:
    """
    Smallest epitope set reaching `goal_pct` coverage, added one peptide at a time.

    Greedy is the right algorithm here rather than a convenience: set coverage is
    monotone submodular in the chosen peptides, so the greedy set is within
    (1 - 1/e) of the smallest possible one, and at the scale this runs at (tens
    of candidates) it is normally optimal outright.

    `seed` keeps peptides the user already picked; `tie_break` (e.g. per-peptide
    Composite_Score) decides between candidates whose coverage gain is equal,
    which happens constantly once several peptides share an allele.
    """
    resolved = resolve_region(region)
    tie_break = tie_break or {}

    chosen: Dict[str, List[str]] = {}
    for p in (seed or []):
        if p in candidates:
            chosen[p] = candidates[p]

    trace: List[Dict[str, Any]] = []
    current = set_coverage(chosen, resolved)[objective] if chosen else 0.0

    while len(chosen) < max_epitopes and current < goal_pct:
        best_pep, best_val, best_tb = None, current, -1e9
        for pep, alleles in candidates.items():
            if pep in chosen:
                continue
            trial = dict(chosen)
            trial[pep] = alleles
            val = set_coverage(trial, resolved)[objective]
            tb = tie_break.get(pep, 0.0)
            if val > best_val + 1e-9 or (abs(val - best_val) <= 1e-9 and val > current + 1e-9 and tb > best_tb):
                best_pep, best_val, best_tb = pep, val, tb

        if best_pep is None:
            break  # no remaining candidate adds anything

        chosen[best_pep] = candidates[best_pep]
        trace.append({
            "peptide": best_pep,
            "gain_pct": round(best_val - current, 3),
            "coverage_pct": round(best_val, 2),
            "n_epitopes": len(chosen),
        })
        current = best_val

    result = set_coverage(chosen, resolved)
    result.update({
        "peptides": list(chosen.keys()),
        "goal_pct": goal_pct,
        "goal_met": current >= goal_pct,
        "objective": objective,
        "trace": trace,
    })
    return result


if __name__ == "__main__":
    demo_rows = [
        {"Peptide": "FLPSDFFPSV", "MHC": "HLA-A*02:01", "Rank_EL": 0.05},
        {"Peptide": "KLVALGINAV", "MHC": "HLA-A*02:01", "Rank_EL": 0.12},
        {"Peptide": "AYAQKIFKI",  "MHC": "HLA-A*24:02", "Rank_EL": 0.09},
        {"Peptide": "RTLNAWVKV",  "MHC": "HLA-A*11:01", "Rank_EL": 0.31},
        {"Peptide": "QVPLRPMTYK", "MHC": "HLA-A*11:01", "Rank_EL": 0.44},
        {"Peptide": "QVPLRPMTYK", "MHC": "HLA-B*40:01", "Rank_EL": 1.20},
        {"Peptide": "NPPIPVGEIY", "MHC": "HLA-B*35:01", "Rank_EL": 0.18},
    ]
    panel = ["HLA-A*02:01", "HLA-A*24:02", "HLA-A*11:01", "HLA-B*40:01",
             "HLA-B*35:01", "HLA-C*07:02"]
    cand = build_binder_map(demo_rows)

    print("=== Epitope-set coverage (East Asia) ===")
    for picked in (["FLPSDFFPSV"], ["FLPSDFFPSV", "KLVALGINAV"],
                   ["FLPSDFFPSV", "AYAQKIFKI"], list(cand)):
        sub = {p: cand[p] for p in picked}
        cov = set_coverage(sub, "East Asia", panel)
        print(f"  {len(picked)} epitope(s): >=1 {cov['coverage_pct']:6.2f}%   "
              f">=2 {cov['coverage2_pct']:6.2f}%   "
              f"uncovered: {[u['allele'] for u in cov['uncovered_alleles']] or '-'}")

    print("\n=== Greedy selection to 85% (East Asia) ===")
    g = greedy_select(cand, "East Asia", goal_pct=85.0)
    for step in g["trace"]:
        print(f"  +{step['peptide']:<12} gain {step['gain_pct']:5.2f} pts "
              f"-> {step['coverage_pct']:.2f}%")
    print(f"  goal met: {g['goal_met']}  with {len(g['peptides'])} epitopes")
