"""
population_coverage.py
HLA Population Coverage Calculator.
Estimates the percentage of individuals in target human populations (World, East Asia, Europe)
covered by selected HLA Class I alleles carrying binder epitopes.

Methodology: Bui HH, Sidney J, Dinh K, Southwood S, Newman MJ, Sette A.
"Predicting population coverage of T-cell epitope-based diagnostics and vaccines."
BMC Bioinformatics. 2006;7:153. PMID: 16545123.
(Same Hardy-Weinberg cumulative-coverage approach used by the IEDB Population
Coverage Analysis Resource, https://tools.iedb.org/population/.)

Allele frequency data: Gonzalez-Galarza FF, et al. "Allele frequency net database
(AFND) 2020 update: gold-standard data classification, open access genotype data
and new query tools." Nucleic Acids Research. 2020;48(D1):D783-D788.
(allelefrequencies.net, licensed under CC BY 4.0.)

Frequencies are weighted averages (by reported sample size n) of 2-field
HLA-A/B/C allele records from a snapshot of the AFND dataset assembled by the
community project https://github.com/slowkow/allelefrequencies (MIT license;
full text bundled at THIRD_PARTY_LICENSES/LICENSE-slowkow-allelefrequencies.txt).
They are approximations across the countries/studies present in AFND, not a
demographically-representative census — treat cumulative coverage % as an
illustrative estimate, not a precise epidemiological figure.

The tables themselves live in data/allele_frequencies.json and are generated and
re-verified by tools/build_afnd_tables.py, which also records each region's exact
population set:

    python3 tools/build_afnd_tables.py --verify

Regions
───────
  World        every AFND population.
  East Asia    China, Japan, Taiwan, Mongolia, Hong Kong. Does NOT include Korea.
  Europe       24 European countries.
  Taiwan       island-wide demographic blend, plus every component group below,
               each selectable in its own right.

Taiwan is not one population
────────────────────────────
Its component groups differ far too much to be represented by a single table.
HLA-A*24:02 runs at 0.309 in the Minnan entries and 0.863 in the Paiwan ones;
Saisiyat carry C*07:02 at 0.667 and B*39:01 at 0.549; Tao carry B*15:25 at 0.400,
which no other group does. A vaccine aimed at one community designed against an
island-wide average would misjudge every one of those.

So each group is its own selectable region — Han (all), Minnan, Hakka, Indigenous
(all), the 11 recognised peoples individually, Pingpu (all), Pazeh and Siraya —
and the island-wide "Taiwan" table is an explicit demographic blend rather than a
pooled average. Pooling would not work anyway: AFND holds 11 recognised
Indigenous entries (n = 693) against 6 Han-representative ones (n = 1947), so a
sample-size average would hand Indigenous peoples ~26% of the table against
roughly 2.5% of the population (~580,000 of 23.4 million, Ministry of the
Interior), reporting A*24:02 at 0.29 instead of 0.18. The blend combines the two
group tables at 97.5% / 2.5%; allele frequencies are additive across
subpopulations weighted by population share, so that is the correct combination
rather than an approximation.

Plains Indigenous (Pazeh, Siraya) sit in neither aggregate: they are not
nationally recognised Indigenous peoples and are largely assimilated into the
Han-identifying population, so folding them into either bucket would misstate it.

Several single-people tables rest on 30-55 individuals. REGION_META flags these
(`small_sample`) so the UI can warn before anyone designs against them.

Two data-quality corrections applied in v1.2
────────────────────────────────────────────
  * AFND writes large sample sizes with thousands separators ("3,456,066"), and
    36.6% of HLA-A/B/C records carry one. Parsing n with a bare float() silently
    dropped all of them — 99.0% of all sampled individuals, including the DKMS
    and NMDP registries. v1.1's World and Europe tables were therefore built from
    the ~1% of the data with n < 1000. Fixed in tools/build_afnd_tables.py.
  * With that fixed, one mis-entered record can dominate a pooled region: AFND
    lists C*01:57 at 0.1996 in a single 4,889-person NMDP cohort and nowhere
    else, enough to push a rare allele into the World top-18. Pooled regions now
    require an allele to be reported by >= 3 populations. The curated Taiwan
    groups are exempt, because there a single-population allele is the signal.

────────────────────────────────────────────────────────────────────────────────
v1.2 — locus-aware model
────────────────────────────────────────────────────────────────────────────────
The values stored below are ALLELE frequencies (gene frequencies), not phenotype
frequencies. An individual carries exactly two alleles per locus, so alleles at
the SAME locus compete for the same two slots and must be summed before the
Hardy-Weinberg square is applied; only then are the independent A / B / C loci
combined multiplicatively:

    P(individual carries none of the selected alleles)
        =  prod over loci L of  ( 1 - sum_{a in L} f_a )^2
    coverage = 1 - that

Releases up to v1.1 applied `prod_a (1 - f_a)^2` across every allele regardless
of locus, which treats two alleles at the same locus as independent events and
therefore UNDER-estimates coverage within a locus. `absence_probability()` below
is the shared primitive; `epitope_set.py` reuses it so that panel coverage and
epitope-set coverage can never drift apart.
"""

import json
from pathlib import Path
from typing import List, Dict, Any, Iterable, Tuple

# ─── Frequency tables ─────────────────────────────────────────────────────────
# The tables live in data/allele_frequencies.json, generated and verified by
# tools/build_afnd_tables.py. They are data, not code: 22 regions x 18 alleles is
# too much to keep readable inline, and holding them in one generated file is what
# makes `--verify` a real check rather than a comparison against a hand-edit.

_DATA_PATH = Path(__file__).parent / "data" / "allele_frequencies.json"

try:
    with open(_DATA_PATH, encoding="utf-8") as _fh:
        _DATA: Dict[str, Any] = json.load(_fh)
except FileNotFoundError as exc:  # pragma: no cover - deployment error, not a runtime path
    raise RuntimeError(
        f"Allele frequency tables missing at {_DATA_PATH}. Regenerate them with:\n"
        f"    python3 tools/build_afnd_tables.py --emit --json --out data/allele_frequencies.json"
    ) from exc

DATA_SOURCE: str = _DATA.get("source_citation", "")
DATA_METHOD: str = _DATA.get("method", "")
DATA_GENERATED: str = _DATA.get("generated_at", "")
SMALL_SAMPLE_N: int = int(_DATA.get("small_sample_threshold", 100))

# {region: {allele: frequency}}, each allele stored under both the WHO spelling
# (HLA-A*02:01) and the netMHCpan spelling (HLA-A02:01) so lookups never have to
# care which form the caller holds.
ALLELE_FREQUENCIES: Dict[str, Dict[str, float]] = {}
# {region: {label, group, note, n, small_sample, populations}} — everything the UI
# needs to say where a table came from, next to the selector rather than in a doc.
REGION_META: Dict[str, Dict[str, Any]] = {}

for _region, _meta in _DATA["regions"].items():
    _freqs: Dict[str, float] = {}
    for _allele, _f in _meta["frequencies"].items():
        _freqs[_allele] = _f
        _freqs[_allele.replace("*", "")] = _f
    ALLELE_FREQUENCIES[_region] = _freqs
    REGION_META[_region] = {k: v for k, v in _meta.items() if k != "frequencies"}

REGIONS: List[str] = list(ALLELE_FREQUENCIES.keys())

# Regions grouped for a hierarchical selector, in the order the UI should show them.
REGION_GROUPS: Dict[str, List[str]] = {}
for _region, _meta in REGION_META.items():
    REGION_GROUPS.setdefault(_meta["group"], []).append(_region)

# UI-facing labels that are not themselves tables.
REGION_ALIASES: Dict[str, str] = {
    "Asia": "East Asia",
    "Global": "World",
    "TW": "Taiwan",
    "Taiwan Hoklo": "Taiwan Minnan",
    "Taiwan Yami": "Taiwan Tao",
    "Taiwan Taroko": "Taiwan Truku",
    "Taiwan Ami": "Taiwan Amis",
    "Taiwan Saisiat": "Taiwan Saisiyat",
}

# A locus only has two slots, so the summed allele frequency cannot exceed 1.
# Clamp slightly below 1 so that a hand-edited table can never yield a coverage
# of exactly 100%, which would be an overclaim rather than a computation.
_MAX_LOCUS_SUM = 0.99


def resolve_region(region: str) -> str:
    """Map a UI region label onto a region present in ALLELE_FREQUENCIES."""
    if region in ALLELE_FREQUENCIES:
        return region
    return REGION_ALIASES.get(region, "World")


def normalize_allele(allele: str) -> str:
    """Canonical lookup key: strip whitespace and the '*' separator (HLA-A*02:01 -> HLA-A02:01)."""
    return str(allele).strip().replace("*", "")


def locus_of(allele: str) -> str:
    """
    Return the genetic locus an allele belongs to. Alleles at the same locus
    compete for the same two chromosomal slots and must be summed, not multiplied.

    HLA-A*02:01 -> 'HLA-A'   H-2-Kb -> 'H-2-K'   Mamu-A1*001:01 -> 'Mamu-A1'
    """
    a = normalize_allele(allele)
    if a.startswith("HLA-"):
        # HLA-A02:01 / HLA-DRB1... — locus is everything before the first digit
        rest = a[4:]
        gene = ""
        for ch in rest:
            if ch.isdigit():
                break
            gene += ch
        return "HLA-" + (gene or "?")
    if a.startswith("H-2-") or a.startswith("H2-"):
        prefix, rest = ("H-2-", a[4:]) if a.startswith("H-2-") else ("H2-", a[3:])
        return prefix + (rest[:1].upper() if rest else "?")
    # Other species (BoLA, SLA, Mamu, Patr, ...): split on the first digit or ':'
    gene = ""
    for ch in a:
        if ch.isdigit() or ch == ":":
            break
        gene += ch
    return gene or a


def is_human(allele: str) -> bool:
    return normalize_allele(allele).startswith("HLA-")


def allele_frequency(allele: str, region: str = "World") -> float:
    """Allele (gene) frequency for one allele in one region; 0.0 when unknown."""
    freq_db = ALLELE_FREQUENCIES[resolve_region(region)]
    return freq_db.get(normalize_allele(allele), 0.0)


def split_known(alleles: Iterable[str], region: str = "World") -> Tuple[List[str], List[str], List[str]]:
    """
    Partition alleles into (known_human, unknown_human, non_human).

    'unknown_human' are HLA alleles absent from this region's frequency table.
    They contribute 0 to coverage, which makes every reported figure a LOWER
    BOUND rather than a guess — callers surface the list so the user can see
    exactly which alleles were not counted. (Up to v1.1 these silently received
    a made-up 0.08 frequency, which inflated coverage and also produced a
    meaningless "human population coverage" for mouse H-2 runs.)
    """
    freq_db = ALLELE_FREQUENCIES[resolve_region(region)]
    known, unknown, non_human = [], [], []
    for a in dict.fromkeys(alleles):
        if not is_human(a):
            non_human.append(a)
        elif normalize_allele(a) in freq_db:
            known.append(a)
        else:
            unknown.append(a)
    return known, unknown, non_human


def absence_probability(alleles: Iterable[str], region: str = "World") -> float:
    """
    P(a random individual carries NONE of `alleles`), under Hardy-Weinberg
    equilibrium with loci assumed independent.

        prod over loci L of ( 1 - sum_{a in L} f_a )^2

    This is the single primitive behind every coverage number in the platform:
    panel coverage, epitope-set coverage and marginal gain all call it, so they
    are guaranteed to be mutually consistent.
    """
    freq_db = ALLELE_FREQUENCIES[resolve_region(region)]

    locus_sums: Dict[str, float] = {}
    for a in dict.fromkeys(alleles):
        key = normalize_allele(a)
        f = freq_db.get(key)
        if f is None:
            continue  # unknown or non-human: contributes nothing (lower bound)
        locus_sums[locus_of(a)] = locus_sums.get(locus_of(a), 0.0) + f

    p_absent = 1.0
    for total in locus_sums.values():
        p_absent *= (1.0 - min(total, _MAX_LOCUS_SUM)) ** 2
    return p_absent


def calculate_population_coverage(alleles: List[str], region: str = "World") -> Dict[str, Any]:
    """
    Cumulative population coverage % for a set of HLA Class I alleles — i.e. the
    share of individuals carrying at least one of them.

    NOTE: this measures the ALLELE PANEL, not a vaccine. A panel can cover 95%
    of a population while the epitopes predicted against it cover far fewer
    people. For epitope-set coverage use epitope_set.set_coverage().
    """
    resolved = resolve_region(region)
    known, unknown, non_human = split_known(alleles, resolved)
    unique_alleles = list(dict.fromkeys(alleles))

    cumulative_coverage = round((1.0 - absence_probability(unique_alleles, resolved)) * 100, 2)

    allele_details = []
    for a in unique_alleles:
        f = allele_frequency(a, resolved)
        allele_details.append({
            "allele": a,
            "locus": locus_of(a),
            "frequency": round(f, 4),
            "individual_coverage_pct": round((1.0 - (1.0 - f) ** 2) * 100, 2),
            "known": f > 0.0,
        })
    allele_details.sort(key=lambda d: -d["frequency"])

    return {
        "region": resolved,
        "region_requested": region,
        "cumulative_coverage_pct": cumulative_coverage,
        "n_alleles": len(unique_alleles),
        "n_known": len(known),
        "unknown_alleles": unknown,
        "non_human_alleles": non_human,
        "is_lower_bound": bool(unknown),
        "allele_details": allele_details,
    }


def top_alleles(region: str = "World", limit: int = 30) -> List[Dict[str, Any]]:
    """Alleles of a region ordered by frequency — the pool the greedy panel builder draws from."""
    freq_db = ALLELE_FREQUENCIES[resolve_region(region)]
    # Every allele is stored under both spellings; keep the WHO-style one (with
    # '*') so anything shown in the UI matches netMHCpan's own output format.
    starred = {normalize_allele(n): n for n in freq_db if "*" in n}

    seen, out = set(), []
    for name, f in sorted(freq_db.items(), key=lambda kv: -kv[1]):
        key = normalize_allele(name)
        if key in seen:
            continue
        seen.add(key)
        out.append({"allele": starred.get(key, name), "locus": locus_of(name), "frequency": round(f, 4)})
        if len(out) >= limit:
            break
    return out


def build_panel(region: str = "World", goal_pct: float = 90.0, max_alleles: int = 20) -> Dict[str, Any]:
    """
    Smallest high-frequency allele panel reaching `goal_pct` coverage of `region`.

    Greedy on marginal coverage gain, which is optimal enough here: the objective
    (1 - product of per-locus absence terms) is monotone submodular, so greedy is
    within (1 - 1/e) of the optimum, and with <= 30 candidate alleles the gap is
    in practice zero.
    """
    resolved = resolve_region(region)
    pool = [d["allele"] for d in top_alleles(resolved, limit=40)]

    chosen: List[str] = []
    coverage = 0.0
    while len(chosen) < max_alleles and coverage < goal_pct:
        best, best_cov = None, coverage
        for cand in pool:
            if cand in chosen:
                continue
            cov = (1.0 - absence_probability(chosen + [cand], resolved)) * 100
            if cov > best_cov + 1e-9:
                best, best_cov = cand, cov
        if best is None:
            break
        chosen.append(best)
        coverage = best_cov

    return {
        "region": resolved,
        "goal_pct": goal_pct,
        "alleles": chosen,
        "coverage_pct": round(coverage, 2),
        "goal_met": coverage >= goal_pct,
    }


if __name__ == "__main__":
    test_alleles = ["HLA-A*02:01", "HLA-B*07:02", "HLA-A*24:02"]
    print("=== HLA Population Coverage (v1.2, locus-aware) ===")
    for reg in REGIONS:
        cov = calculate_population_coverage(test_alleles, reg)
        print(f"{reg:10s}: {cov['cumulative_coverage_pct']:6.2f}%  "
              f"({cov['n_known']}/{cov['n_alleles']} alleles with frequency data)")

    print("\n=== Auto-built panels (goal 90%) ===")
    for reg in REGIONS:
        panel = build_panel(reg, 90.0)
        print(f"{reg:10s}: {panel['coverage_pct']:6.2f}% with {len(panel['alleles'])} alleles")
        print(f"            {', '.join(panel['alleles'])}")
