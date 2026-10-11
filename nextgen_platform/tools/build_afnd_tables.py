#!/usr/bin/env python3
"""
build_afnd_tables.py — regenerate the allele-frequency tables in population_coverage.py

Every number in ALLELE_FREQUENCIES is produced by this script from one public
input file, so the tables are auditable and reproducible rather than pasted in by
hand. Run it to verify the shipped tables, or to rebuild them after refreshing
the AFND snapshot:

    python3 tools/build_afnd_tables.py --verify          # check shipped tables
    python3 tools/build_afnd_tables.py --emit            # print the Python dict
    python3 tools/build_afnd_tables.py --emit --region Taiwan

Input
─────
`afnd.tsv` from the community mirror of the Allele Frequency Net Database:
    https://github.com/slowkow/allelefrequencies   (MIT, LICENSE bundled in
    THIRD_PARTY_LICENSES/LICENSE-slowkow-allelefrequencies.txt)
which in turn derives from:
    Gonzalez-Galarza FF, et al. "Allele frequency net database (AFND) 2020
    update." Nucleic Acids Research. 2020;48(D1):D783-D788.  (CC BY 4.0)

The script downloads it on first run and caches it next to this file.

Method
──────
For each region, take every AFND population entry belonging to it, keep HLA-A/B/C
records at 2-field resolution, and average `alleles_over_2n` (the allele/gene
frequency) weighted by the reported sample size `n`. Verified to reproduce the
shipped World and East Asia tables bit-for-bit.

Taiwan is the one region that is NOT a plain weighted average — see TAIWAN_* below.
"""

import argparse
import collections
import csv
import sys
import urllib.request
from pathlib import Path

AFND_URL = "https://raw.githubusercontent.com/slowkow/allelefrequencies/main/afnd.tsv"
AFND_TSV = Path(__file__).parent / "afnd.tsv"

# ─── Region definitions ───────────────────────────────────────────────────────
# A region is a set of prefixes matched against the AFND population label.
# "World" is every population in the database.

EAST_ASIA = ("China", "Japan", "Taiwan", "Mongolia", "Hong Kong")

EUROPE = (
    "Austria", "Belgium", "Bulgaria", "Croatia", "Czech", "England", "Finland",
    "France", "Germany", "Greece", "Ireland", "Italy", "Macedonia", "Netherlands",
    "Poland", "Portugal", "Romania", "Russia", "Scotland", "Serbia", "Spain",
    "Sweden", "Switzerland", "Wales",
)

# ─── Taiwan ───────────────────────────────────────────────────────────────────
# Taiwan is not one population and must not be shipped as one table.
#
# The HLA distributions of its component groups differ enormously — HLA-A*24:02
# runs at 0.309 in the Minnan entries and 0.863 in the Paiwan ones — so a vaccine
# aimed at a specific community cannot be designed against an island-wide average,
# and an island-wide average cannot be built by pooling AFND entries either: AFND
# holds 13 Indigenous entries (n = 799) against 6 Han-representative ones
# (n = 1947), so a plain sample-size average would hand Indigenous peoples ~29% of
# the table against roughly 2.5% of the population.
#
# Every group below is therefore selectable in its own right, and the island-wide
# "Taiwan" table is an explicit demographic blend rather than a pooled average.
# Allele frequencies are additive across subpopulations weighted by population
# share, so the blend is the correct combination, not an approximation.

# Han-majority entries. The four unlabelled/general ones are Han-representative;
# Minnan and Hakka are also selectable separately below.
TAIWAN_HAN = (
    "Taiwan Tzu Chi Cord Blood Bank",  # n=710, general Taiwanese cord blood donors (A/B only)
    "Taiwan Han Chinese",              # n=504
    "Taiwan pop 2",                    # n=364
    "Taiwan pop 3",                    # n=212 (A/B only)
    "Taiwan Minnan pop 1",             # n=102, Minnan / Hoklo
    "Taiwan Hakka",                    # n=55,  Hakka
)

# The 11 nationally recognised Indigenous peoples present in this AFND snapshot.
TAIWAN_INDIGENOUS = (
    "Taiwan Ami", "Taiwan Atayal", "Taiwan Bunun", "Taiwan Paiwan", "Taiwan Puyuma",
    "Taiwan Rukai", "Taiwan Saisiat", "Taiwan Tao", "Taiwan Taroko", "Taiwan Thao",
    "Taiwan Tsou",
)

# Plains Indigenous (Pingpu). Kept out of both the Han and the Indigenous
# aggregates: Pazeh and Siraya are not nationally recognised Indigenous peoples
# and are largely assimilated into the Han-identifying population, so folding them
# into either bucket would misstate that bucket. They remain selectable directly.
TAIWAN_PINGPU = ("Taiwan Pazeh", "Taiwan Siraya")

# Council of Indigenous Peoples / Ministry of the Interior: recognised Indigenous
# peoples are ~580,000 of Taiwan's ~23.4 million residents.
TAIWAN_INDIGENOUS_SHARE = 0.025

# Individually selectable groups: AFND label -> (region key, display label).
TAIWAN_SUBGROUPS = {
    "Taiwan Minnan pop 1": ("Taiwan Minnan", "Minnan / Hoklo 閩南"),
    "Taiwan Hakka":        ("Taiwan Hakka",  "Hakka 客家"),
    "Taiwan Ami":          ("Taiwan Amis",     "Amis 阿美族"),
    "Taiwan Atayal":       ("Taiwan Atayal",   "Atayal 泰雅族"),
    "Taiwan Bunun":        ("Taiwan Bunun",    "Bunun 布農族"),
    "Taiwan Paiwan":       ("Taiwan Paiwan",   "Paiwan 排灣族"),
    "Taiwan Puyuma":       ("Taiwan Puyuma",   "Puyuma 卑南族"),
    "Taiwan Rukai":        ("Taiwan Rukai",    "Rukai 魯凱族"),
    "Taiwan Saisiat":      ("Taiwan Saisiyat", "Saisiyat 賽夏族"),
    "Taiwan Tao":          ("Taiwan Tao",      "Tao / Yami 達悟族"),
    "Taiwan Taroko":       ("Taiwan Truku",    "Truku 太魯閣族"),
    "Taiwan Thao":         ("Taiwan Thao",     "Thao 邵族"),
    "Taiwan Tsou":         ("Taiwan Tsou",     "Tsou 鄒族"),
    "Taiwan Pazeh":        ("Taiwan Pazeh",    "Pazeh 巴宰族 (Pingpu)"),
    "Taiwan Siraya":       ("Taiwan Siraya",   "Siraya 西拉雅族 (Pingpu)"),
}

# Below this sample size an allele-frequency estimate is noisy enough that the UI
# should say so before someone designs against it.
SMALL_SAMPLE_N = 100

REGIONS = (
    "World", "East Asia", "Europe",
    "Taiwan", "Taiwan Han", "Taiwan Indigenous", "Taiwan Pingpu",
    *(key for key, _ in TAIWAN_SUBGROUPS.values()),
)


# ─── Core ─────────────────────────────────────────────────────────────────────

def load_rows(path: Path = AFND_TSV):
    """HLA-A/B/C records from the AFND snapshot, downloading it if absent."""
    if not path.exists():
        print(f"downloading {AFND_URL} -> {path}", file=sys.stderr)
        path.parent.mkdir(parents=True, exist_ok=True)
        urllib.request.urlretrieve(AFND_URL, path)
    with open(path, newline="") as fh:
        return [r for r in csv.DictReader(fh, delimiter="\t")
                if r["group"] == "hla" and r["gene"] in ("A", "B", "C")]


# An allele must be reported by at least this many populations of a region before
# it enters that region's table. Without it a single mis-entered record can
# dominate: AFND lists C*01:57 at 0.1996 in one 4,889-person NMDP cohort and
# nowhere else, which is enough to put a genuinely rare allele into the World
# top-18. Regions built from a single AFND population are exempt (see below).
MIN_POPULATIONS = 3


def _parse_n(raw) -> float:
    """
    Sample size, tolerating the thousands separators AFND uses for large cohorts.

    This matters more than it looks: 36.6% of HLA-A/B/C records carry a comma in
    `n`, and they are the big ones — DKMS German donors (n=3,456,066), NMDP
    European Caucasian (n=1,242,890) and 102 other registries. Parsing `n` with a
    bare float() drops all of them, which is 99.0% of all sampled individuals.
    """
    try:
        return float(str(raw).replace(",", "").strip())
    except (TypeError, ValueError):
        return 0.0


def weighted_average(rows, match=None, exact=None, min_populations=MIN_POPULATIONS):
    """
    Sample-size-weighted mean allele frequency, keyed 'HLA-A*02:01'.

    `match` selects populations by label prefix, `exact` by full label; passing
    neither averages the whole database (the World table).
    """
    num = collections.defaultdict(float)
    den = collections.defaultdict(float)
    seen_pops = collections.defaultdict(set)
    exact = set(exact) if exact else None

    # The support bar guards against one bad record hiding among hundreds of
    # pooled studies, which is a risk for the prefix-matched regions (World,
    # East Asia, Europe). The Taiwan groups are curated lists of named entries
    # where every record has been looked at, and where a single-population allele
    # is usually the point — Tao carries B*15:25 at 0.40 and no other people does.
    # Requiring 3 there would delete exactly the signal the group exists to show.
    if exact is not None:
        min_populations = 1

    for r in rows:
        pop = r["population"]
        if exact is not None and pop not in exact:
            continue
        if match is not None and not pop.startswith(match):
            continue
        allele = r["allele"]
        if allele.count(":") != 1:          # 2-field resolution only
            continue
        try:
            f = float(r["alleles_over_2n"])
        except (TypeError, ValueError):
            continue
        n = _parse_n(r["n"])
        if not n or f != f:                  # missing / NaN
            continue
        key = "HLA-" + allele
        num[key] += f * n
        den[key] += n
        seen_pops[key].add(pop)

    return {k: num[k] / den[k] for k in num
            if den[k] and len(seen_pops[k]) >= min_populations}


def blend(tables_and_weights):
    """Combine subpopulation tables at their demographic shares."""
    out = collections.defaultdict(float)
    for table, w in tables_and_weights:
        for k, v in table.items():
            out[k] += v * w
    return dict(out)


def region_populations(region: str):
    """The AFND population labels a region is built from (None = the whole database)."""
    if region == "World":
        return None
    if region == "East Asia":
        return ("prefix", EAST_ASIA)
    if region == "Europe":
        return ("prefix", EUROPE)
    if region == "Taiwan":
        return ("blend", TAIWAN_HAN + TAIWAN_INDIGENOUS)
    if region == "Taiwan Han":
        return ("exact", TAIWAN_HAN)
    if region == "Taiwan Indigenous":
        return ("exact", TAIWAN_INDIGENOUS)
    if region == "Taiwan Pingpu":
        return ("exact", TAIWAN_PINGPU)
    for afnd_label, (key, _) in TAIWAN_SUBGROUPS.items():
        if key == region:
            return ("exact", (afnd_label,))
    raise ValueError(f"unknown region: {region}")


def build_region(rows, region: str):
    spec = region_populations(region)
    if spec is None:
        return weighted_average(rows)
    kind, pops = spec
    if kind == "prefix":
        return weighted_average(rows, match=pops)
    if kind == "exact":
        return weighted_average(rows, exact=pops)
    # Island-wide Taiwan: average each group on its own, then combine at the
    # groups' demographic shares rather than at their sample sizes.
    han = weighted_average(rows, exact=TAIWAN_HAN)
    ind = weighted_average(rows, exact=TAIWAN_INDIGENOUS)
    return blend([(han, 1.0 - TAIWAN_INDIGENOUS_SHARE),
                  (ind, TAIWAN_INDIGENOUS_SHARE)])


def region_sample_size(rows, region: str) -> int:
    """Largest reported n across the region's populations — the UI shows this."""
    spec = region_populations(region)
    if spec is None:
        pops = {r["population"] for r in rows}
    else:
        kind, sel = spec
        if kind == "prefix":
            pops = {r["population"] for r in rows if r["population"].startswith(sel)}
        else:
            pops = set(sel)
    per_pop = collections.defaultdict(float)
    for r in rows:
        if r["population"] in pops:
            per_pop[r["population"]] = max(per_pop[r["population"]], _parse_n(r["n"]))
    return int(sum(per_pop.values()))


def top_n(table, n=18):
    """The n most frequent alleles — what the shipped tables keep."""
    return dict(sorted(table.items(), key=lambda kv: -kv[1])[:n])


# ─── Entry points ─────────────────────────────────────────────────────────────

# Display metadata. `group` drives the optgroup a region appears under in the UI;
# `note` is shown beneath the selector, so each table's provenance travels with it
# instead of living only in the docs.
REGION_DISPLAY = {
    "World":     ("Global", "🌐 World",
                  "Averaged across every population in AFND — an average of the studies present, "
                  "not a demographically weighted world census."),
    "East Asia": ("Global", "🌏 East Asia",
                  "AFND China, Japan, Taiwan, Mongolia and Hong Kong entries. Does not include Korea."),
    "Europe":    ("Global", "🌍 Europe",
                  "AFND entries from 24 European countries."),
    "Taiwan":    ("Taiwan", "🇹🇼 Taiwan (island-wide)",
                  "Han and recognised Indigenous tables combined at 97.5% / 2.5% population shares. "
                  "To design for one community, select that community instead."),
    "Taiwan Han": ("Taiwan", "Han 漢人 (all)",
                   "Six Han-representative AFND entries: Tzu Chi Cord Blood Bank, Han Chinese, "
                   "pop 2, pop 3, Minnan and Hakka."),
    "Taiwan Indigenous": ("Taiwan", "Indigenous 原住民族 (all)",
                          "The 11 nationally recognised Indigenous peoples in this snapshot, averaged "
                          "by sample size — not weighted by each people's population."),
    "Taiwan Pingpu": ("Taiwan", "Pingpu 平埔族 (all)",
                      "Plains Indigenous (Pazeh, Siraya). Kept out of both the Han and Indigenous "
                      "aggregates: they are not nationally recognised Indigenous peoples and are "
                      "largely assimilated into the Han-identifying population."),
}
for _afnd, (_key, _label) in TAIWAN_SUBGROUPS.items():
    _grp = ("Taiwan Han subgroups" if _key in ("Taiwan Minnan", "Taiwan Hakka")
            else "Taiwan Indigenous peoples")
    REGION_DISPLAY[_key] = (_grp, _label, f"Single AFND population: {_afnd}.")


def build_all(rows, limit):
    """Every region's table, plus the metadata the UI and the docs need."""
    import datetime

    out = {}
    for region in REGIONS:
        spec = region_populations(region)
        if spec is None:
            pops = ["<all AFND populations>"]
        elif spec[0] == "prefix":
            pops = [f"{p}*" for p in spec[1]]
        else:
            pops = list(spec[1])
        group, label, note = REGION_DISPLAY[region]
        n = region_sample_size(rows, region)
        out[region] = {
            "label": label,
            "group": group,
            "note": note,
            "n": n,
            "small_sample": n < SMALL_SAMPLE_N,
            "populations": pops,
            "frequencies": {k: round(v, 4)
                            for k, v in top_n(build_region(rows, region), limit).items()},
        }

    return {
        "generated_by": "tools/build_afnd_tables.py",
        "generated_at": datetime.date.today().isoformat(),
        "source": AFND_URL,
        "source_citation": (
            "Gonzalez-Galarza FF, et al. Allele frequency net database (AFND) 2020 update. "
            "Nucleic Acids Res. 2020;48(D1):D783-D788 (CC BY 4.0), via the MIT-licensed mirror "
            "github.com/slowkow/allelefrequencies"
        ),
        "method": (
            f"2-field HLA-A/B/C records; alleles_over_2n averaged weighted by reported sample "
            f"size n. Each region keeps its {limit} most frequent alleles, so an allele absent "
            f"from a table contributes 0 and coverage is a lower bound."
        ),
        "alleles_per_region": limit,
        "small_sample_threshold": SMALL_SAMPLE_N,
        "regions": out,
    }


def emit(rows, regions, limit, as_json=False, out_path=None):
    import json
    data = build_all(rows, limit)
    if regions:
        data["regions"] = {k: v for k, v in data["regions"].items() if k in regions}

    if as_json:
        text = json.dumps(data, ensure_ascii=False, indent=2) + "\n"
        if out_path:
            out_path = Path(out_path)
            out_path.parent.mkdir(parents=True, exist_ok=True)
            out_path.write_text(text, encoding="utf-8")
            print(f"wrote {out_path} — {len(data['regions'])} regions", file=sys.stderr)
        else:
            sys.stdout.write(text)
        return

    for region, meta in data["regions"].items():
        print(f'    "{region}": {{')
        for allele, f in meta["frequencies"].items():
            print(f'        "{allele.replace("*", "")}": {f}, "{allele}": {f},')
        print("    },")


def verify(rows, limit=18):
    """
    Re-derive every table from the source and compare it to the shipped JSON.

    This is a real test of the data file, not of a hand-edit: if anyone changes
    data/allele_frequencies.json by hand, or refreshes afnd.tsv without
    regenerating, this reports exactly which alleles moved.
    """
    import json

    data_path = Path(__file__).resolve().parent.parent / "data" / "allele_frequencies.json"
    if not data_path.exists():
        print(f"FAIL  no shipped tables at {data_path}", file=sys.stderr)
        return False
    shipped = json.loads(data_path.read_text(encoding="utf-8"))["regions"]
    fresh = build_all(rows, limit)["regions"]

    ok = True
    missing = set(fresh) - set(shipped)
    extra = set(shipped) - set(fresh)
    if missing or extra:
        ok = False
        if missing:
            print(f"FAIL  regions missing from the shipped file: {sorted(missing)}")
        if extra:
            print(f"FAIL  regions in the shipped file but not generated: {sorted(extra)}")

    for region in fresh:
        if region not in shipped:
            continue
        a = shipped[region]["frequencies"]
        b = fresh[region]["frequencies"]
        keys = set(a) | set(b)
        diffs = [(k, abs(a.get(k, 0.0) - b.get(k, 0.0))) for k in keys]
        worst_allele, worst = max(diffs, key=lambda t: t[1])
        n = fresh[region]["n"]
        flag = " small-n" if fresh[region]["small_sample"] else ""
        if worst < 6e-5 and set(a) == set(b):
            print(f"OK    {region:<22} {len(b):>2} alleles  n={n:>9,}{flag}")
        else:
            ok = False
            print(f"FAIL  {region:<22} max |Δ| {worst:.4f} ({worst_allele})")
    return ok


def legacy_delta(rows, limit=18):
    """
    Show how much the v1.2 data-quality fixes moved the pre-existing tables.

    v1.1 parsed `n` with a bare float(), which silently discarded every record
    whose sample size carried a thousands separator. Reproducing that behaviour
    lets us quantify the correction instead of asserting it.
    """
    legacy_rows = [r for r in rows if "," not in (r["n"] or "")]
    print("v1.1 behaviour: records with a comma in `n` were dropped "
          f"({len(rows) - len(legacy_rows):,} of {len(rows):,}).\n")
    for region in ("World", "East Asia", "Europe"):
        old = top_n(weighted_average(legacy_rows,
                                     match=None if region == "World" else
                                     (EAST_ASIA if region == "East Asia" else EUROPE),
                                     min_populations=1), limit)
        new = top_n(build_region(rows, region), limit)
        print(f"=== {region} ===")
        for allele in list(new)[:6]:
            o = old.get(allele, 0.0)
            print(f"   {allele:<14} v1.1 {o:.4f}  ->  v1.2 {new[allele]:.4f}   "
                  f"Δ {new[allele] - o:+.4f}")
        dropped = [a for a in old if a not in new]
        if dropped:
            print(f"   left the top-{limit}: {', '.join(dropped)}")
        print()


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--verify", action="store_true", help="check the shipped tables")
    ap.add_argument("--emit", action="store_true", help="print the tables")
    ap.add_argument("--json", action="store_true", help="emit JSON instead of a Python dict")
    ap.add_argument("--out", type=Path, help="write the JSON here (default: stdout)")
    ap.add_argument("--region", action="append", choices=REGIONS, help="limit --emit to these")
    ap.add_argument("--limit", type=int, default=18, help="alleles kept per region (default 18)")
    ap.add_argument("--legacy-delta", action="store_true",
                    help="show how much the v1.2 data-quality fixes moved the old tables")
    ap.add_argument("--tsv", type=Path, default=AFND_TSV)
    args = ap.parse_args()

    rows = load_rows(args.tsv)
    print(f"# {len(rows):,} HLA-A/B/C records loaded from {args.tsv.name}", file=sys.stderr)

    if args.emit:
        emit(rows, args.region, args.limit, as_json=args.json or bool(args.out), out_path=args.out)
    if args.legacy_delta:
        legacy_delta(rows, args.limit)
    if args.verify or not (args.emit or args.legacy_delta):
        sys.exit(0 if verify(rows, args.limit) else 1)


if __name__ == "__main__":
    main()
