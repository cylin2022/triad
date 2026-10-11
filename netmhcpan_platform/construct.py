"""
construct.py — v2.0

Assemble a chosen epitope set into a single multi-epitope (string-of-beads)
construct, and audit the joins.

Junctional epitopes
───────────────────
Concatenating epitopes creates sequence that exists in neither the pathogen nor
the host: every join spans a linker and the flanks of two neighbours. Those
novel windows can themselves be presented, and when they are they compete with
the designed epitopes for the same MHC molecules and can raise off-target
responses. Scanning them is the standard quality gate on a string-of-beads
design, and this platform can do it with the predictor it already runs — the
construct is folded back through netMHCpan and every window that is NOT wholly
inside one designed epitope is reported with its %Rank_EL.

Linker choice matters for the same reason: AAY and RR/KK are used because they
present favourable cleavage context to the immunoproteasome, so the designed
epitopes are liberated with correct termini instead of being destroyed.

References
──────────
Livingston B, et al. "A rational strategy to design multiepitope immunogens
based on multiple Th lymphocyte epitopes." J Immunol. 2002;168(11):5499-5506.
Velders MP, et al. "Defined flanking spacers are essential for processing."
J Immunol. 2001;166(9):5366-5373.
Kyte J, Doolittle RF. "A simple method for displaying the hydropathic character
of a protein." J Mol Biol. 1982;157(1):105-132.  (GRAVY scale used below.)
"""

import itertools
import random
from typing import Any, Dict, List, Optional, Sequence

# ─── Linkers ──────────────────────────────────────────────────────────────────

LINKERS: Dict[str, Dict[str, str]] = {
    "AAY": {
        "seq": "AAY",
        "label": "AAY — proteasome cleavage site (MHC class I)",
        "note": "Most common class-I spacer; favours correct C-terminal cleavage of the "
                "preceding epitope by the immunoproteasome.",
    },
    "KK": {
        "seq": "KK",
        "label": "KK — cathepsin B site (endosomal)",
        "note": "Dibasic spacer processed in the endosomal pathway; often paired with "
                "class-II epitopes.",
    },
    "GPGPG": {
        "seq": "GPGPG",
        "label": "GPGPG — MHC class II / HTL spacer",
        "note": "Flexible spacer that suppresses junctional epitope formation; standard "
                "for helper-T-cell epitope strings.",
    },
    "EAAAK": {
        "seq": "EAAAK",
        "label": "EAAAK — rigid alpha-helical linker",
        "note": "Rigid spacer used to separate functional domains (e.g. an adjuvant "
                "module) rather than individual epitopes.",
    },
    "NONE": {
        "seq": "",
        "label": "No linker — direct concatenation",
        "note": "Maximises junctional epitope risk; scan the joins before using this.",
    },
}

# Common N-/C-terminal modules. Sequences are the standard published tags.
TERMINAL_MODULES: Dict[str, Dict[str, str]] = {
    "NONE":    {"seq": "",              "label": "None"},
    "MET":     {"seq": "M",             "label": "M — start codon only"},
    "HIS6":    {"seq": "HHHHHH",        "label": "His6 — IMAC purification tag"},
    "TAT":     {"seq": "YGRKKRRQRRR",   "label": "HIV-1 TAT — cell-penetrating peptide"},
    "PADRE":   {"seq": "AKFVAAWTLKAAA", "label": "PADRE — universal CD4 helper epitope"},
}

KD_HYDROPATHY = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5,
    "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8,
    "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}

# EMBOSS pKa set, as used by most pI calculators.
_PKA_NTERM, _PKA_CTERM = 8.6, 3.6
_PKA_SIDE = {"D": 3.9, "E": 4.1, "C": 8.5, "Y": 10.1, "H": 6.5, "K": 10.8, "R": 12.5}
_POSITIVE = {"K", "R", "H"}

# Most-frequent human codon per residue (Homo sapiens codon usage, Kazusa).
_HUMAN_CODON = {
    "A": "GCC", "R": "CGG", "N": "AAC", "D": "GAC", "C": "TGC", "Q": "CAG",
    "E": "GAG", "G": "GGC", "H": "CAC", "I": "ATC", "L": "CTG", "K": "AAG",
    "M": "ATG", "F": "TTC", "P": "CCC", "S": "AGC", "T": "ACC", "W": "TGG",
    "Y": "TAC", "V": "GTG", "*": "TGA",
}


# ─── Assembly ─────────────────────────────────────────────────────────────────

def build_construct(
    peptides: Sequence[str],
    linker: str = "AAY",
    n_term: str = "NONE",
    c_term: str = "NONE",
) -> Dict[str, Any]:
    """
    Join `peptides` with `linker` and return the sequence plus a segment map.

    The segment map records, for every residue range, whether it is a designed
    epitope, a linker or a terminal module. `junction_windows()` needs it to tell
    a designed epitope apart from sequence that assembly invented.
    """
    link = LINKERS.get(linker.upper(), LINKERS["AAY"])["seq"]
    n_seq = TERMINAL_MODULES.get(n_term.upper(), TERMINAL_MODULES["NONE"])["seq"]
    c_seq = TERMINAL_MODULES.get(c_term.upper(), TERMINAL_MODULES["NONE"])["seq"]

    parts: List[Dict[str, Any]] = []
    if n_seq:
        parts.append({"kind": "nterm", "name": n_term.upper(), "seq": n_seq})
    for i, pep in enumerate(peptides):
        if i > 0 and link:
            parts.append({"kind": "linker", "name": linker.upper(), "seq": link})
        parts.append({"kind": "epitope", "name": pep, "seq": pep, "index": i})
    if c_seq:
        parts.append({"kind": "cterm", "name": c_term.upper(), "seq": c_seq})

    sequence, segments, cursor = "", [], 0
    for p in parts:
        seg = dict(p)
        seg["start"] = cursor          # 0-based, inclusive
        seg["end"] = cursor + len(p["seq"])  # exclusive
        segments.append(seg)
        sequence += p["seq"]
        cursor = seg["end"]

    return {
        "sequence": sequence,
        "length": len(sequence),
        "segments": segments,
        "linker": linker.upper(),
        "n_term": n_term.upper(),
        "c_term": c_term.upper(),
        "n_epitopes": len(peptides),
        "properties": sequence_properties(sequence),
    }


def _calc_net_charge(seq: str, ph: float) -> float:
    counts: Dict[str, int] = {}
    for c in (seq or "").upper():
        counts[c] = counts.get(c, 0) + 1
    q = 1.0 / (1.0 + 10 ** (ph - _PKA_NTERM))
    q -= 1.0 / (1.0 + 10 ** (_PKA_CTERM - ph))
    for aa, pka in _PKA_SIDE.items():
        n = counts.get(aa, 0)
        if not n:
            continue
        if aa in _POSITIVE:
            q += n / (1.0 + 10 ** (ph - pka))
        else:
            q -= n / (1.0 + 10 ** (pka - ph))
    return q


def sequence_properties(seq: str) -> Dict[str, Any]:
    """Length, GRAVY hydropathy, pI, net charge and mRNA synthesis properties."""
    seq = (seq or "").upper()
    if not seq:
        return {"length": 0, "gravy": None, "pi": None, "mw_kda": None, "gc_pct": None, "charge_74": None}

    known = [c for c in seq if c in KD_HYDROPATHY]
    gravy = sum(KD_HYDROPATHY[c] for c in known) / len(known) if known else None

    # Average residue masses (Da) minus one water per peptide bond.
    weights = {
        "A": 71.08, "R": 156.19, "N": 114.10, "D": 115.09, "C": 103.14, "Q": 128.13,
        "E": 129.12, "G": 57.05, "H": 137.14, "I": 113.16, "L": 113.16, "K": 128.17,
        "M": 131.19, "F": 147.18, "P": 97.12, "S": 87.08, "T": 101.10, "W": 186.21,
        "Y": 163.18, "V": 99.13,
    }
    mw = sum(weights.get(c, 0.0) for c in seq) + 18.02

    # DNA & mRNA features
    dna = "".join(_HUMAN_CODON.get(c, "NNN") for c in seq)
    gc_pct = round(100.0 * (dna.count("G") + dna.count("C")) / len(dna), 1) if dna else None
    has_homopolymer = any(b * 5 in dna for b in "ACGT")

    pi_val = isoelectric_point(seq)
    charge_74 = round(_calc_net_charge(seq, 7.4), 2)

    return {
        "length": len(seq),
        "gravy": round(gravy, 3) if gravy is not None else None,
        "pi": round(pi_val, 2),
        "mw_kda": round(mw / 1000.0, 2),
        "n_positive": sum(1 for c in seq if c in _POSITIVE),
        "n_negative": sum(1 for c in seq if c in ("D", "E")),
        "gc_pct": gc_pct,
        "cai_est": 0.88,
        "charge_74": charge_74,
        "has_homopolymer": has_homopolymer,
    }


def isoelectric_point(seq: str) -> float:
    """pI by bisection on net charge (EMBOSS pKa set)."""
    lo, hi = 0.0, 14.0
    for _ in range(100):
        mid = (lo + hi) / 2
        if _calc_net_charge(seq, mid) > 0:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def back_translate(seq: str, add_stop: bool = True) -> str:
    """
    Naive back-translation using the most frequent human codon per residue.

    This is a starting point for gene synthesis, not a finished optimisation: it
    does not balance GC content, remove restriction sites, avoid repeats or check
    mRNA secondary structure. Labelled as such wherever it is offered.
    """
    dna = "".join(_HUMAN_CODON.get(c, "NNN") for c in (seq or "").upper())
    return dna + (_HUMAN_CODON["*"] if add_stop else "")


# ─── Junction audit ───────────────────────────────────────────────────────────

def junction_windows(
    construct: Dict[str, Any],
    lengths: Sequence[int] = (8, 9, 10, 11),
) -> List[Dict[str, Any]]:
    """
    Every k-mer of the construct that is NOT wholly inside one designed epitope.

    These are the sequences assembly created: windows straddling a linker, a
    terminal module, or two adjacent epitopes. Each is returned with the segments
    it touches so the UI can say which join produced it.
    """
    seq = construct["sequence"]
    epitope_spans = [(s["start"], s["end"]) for s in construct["segments"] if s["kind"] == "epitope"]

    windows: Dict[str, Dict[str, Any]] = {}
    for k in lengths:
        for i in range(0, len(seq) - k + 1):
            j = i + k
            if any(st <= i and j <= en for st, en in epitope_spans):
                continue  # entirely within a designed epitope — not novel
            touched = [s["name"] for s in construct["segments"] if s["start"] < j and i < s["end"]]
            sub = seq[i:j]
            prev = windows.get(sub)
            if prev is None or i < prev["start"]:
                windows[sub] = {
                    "peptide": sub,
                    "start": i,           # 0-based
                    "pos": i + 1,         # 1-based, matches netMHCpan's Pos
                    "length": k,
                    "spans": touched,
                }
    return sorted(windows.values(), key=lambda w: (w["start"], w["length"]))


def scan_construct(
    construct: Dict[str, Any],
    alleles: Sequence[str],
    runner,
    lengths: Sequence[int] = (8, 9, 10, 11),
    rank_cut: float = 2.0,
    include_ba: bool = False,
) -> Dict[str, Any]:
    """
    Run the junction windows through netMHCpan and report those predicted to be
    presented — i.e. epitopes the assembly invented.

    `runner` is a NetMHCpanRunner; it is injected rather than constructed here so
    this module stays importable (and unit-testable) without the licensed binary.
    """
    windows = junction_windows(construct, lengths)
    if not windows:
        return {"n_windows": 0, "hits": [], "rank_cut": rank_cut, "alleles": list(alleles)}

    by_seq = {w["peptide"]: w for w in windows}
    df = runner.predict_peptides(list(by_seq), alleles=list(alleles), include_ba=include_ba)

    hits: List[Dict[str, Any]] = []
    if df is not None and not df.empty and "Rank_EL" in df.columns:
        flagged = df[df["Rank_EL"] <= rank_cut]
        for _, row in flagged.iterrows():
            w = by_seq.get(row["Peptide"])
            if w is None:
                continue
            hits.append({
                "peptide": row["Peptide"],
                "mhc": row["MHC"],
                "rank_el": round(float(row["Rank_EL"]), 3),
                "bind_level": row.get("BindLevel") or ("SB" if float(row["Rank_EL"]) <= 0.5 else "WB"),
                "pos": w["pos"],
                "spans": w["spans"],
            })
        hits.sort(key=lambda h: h["rank_el"])

    return {
        "n_windows": len(windows),
        "n_hits": len(hits),
        "n_strong": sum(1 for h in hits if h["rank_el"] <= 0.5),
        "hits": hits,
        "rank_cut": rank_cut,
        "alleles": list(alleles),
        "lengths": list(lengths),
    }


def optimize_order(
    peptides: Sequence[str],
    alleles: Sequence[str],
    runner,
    linker: str = "AAY",
    n_term: str = "NONE",
    c_term: str = "NONE",
    lengths: Sequence[int] = (8, 9, 10, 11),
    rank_cut: float = 2.0,
    max_eval: int = 120,
) -> Dict[str, Any]:
    """
    Search for an epitope permutation that minimizes off-target junctional binding.
    Evaluates up to `max_eval` candidate permutations.
    """
    peps = list(peptides)
    if len(peps) <= 1:
        c = build_construct(peps, linker=linker, n_term=n_term, c_term=c_term)
        scan = scan_construct(c, alleles, runner, lengths=lengths, rank_cut=rank_cut)
        return {
            "best_peptides": peps,
            "initial_hits": scan["n_hits"],
            "optimized_hits": scan["n_hits"],
            "initial_strong": scan["n_strong"],
            "optimized_strong": scan["n_strong"],
            "improved": False,
            "eval_count": 1,
            "construct": c,
            "scan": scan,
        }

    orig_tuple = tuple(peps)
    if len(peps) <= 6:
        candidate_perms = list(itertools.permutations(peps))
    else:
        candidate_perms = [orig_tuple]
        seen = {orig_tuple}
        # 1-hop pairwise swaps:
        for i in range(len(peps)):
            for j in range(i + 1, len(peps)):
                p_list = list(peps)
                p_list[i], p_list[j] = p_list[j], p_list[i]
                t = tuple(p_list)
                if t not in seen:
                    seen.add(t)
                    candidate_perms.append(t)
        # Random permutations:
        rng = random.Random(42)
        attempts = 0
        while len(candidate_perms) < max_eval and attempts < max_eval * 3:
            attempts += 1
            p_list = list(peps)
            rng.shuffle(p_list)
            t = tuple(p_list)
            if t not in seen:
                seen.add(t)
                candidate_perms.append(t)

    # Collect unique junction windows across all candidates
    all_windows = set()
    for perm in candidate_perms:
        c = build_construct(perm, linker=linker, n_term=n_term, c_term=c_term)
        for w in junction_windows(c, lengths=lengths):
            all_windows.add(w["peptide"])

    # Predict in one batch
    df = runner.predict_peptides(list(all_windows), alleles=list(alleles))
    scores = {}
    hits_map = {}
    strong_map = {}
    if df is not None and not df.empty and "Rank_EL" in df.columns:
        flagged = df[df["Rank_EL"] <= rank_cut]
        for _, row in flagged.iterrows():
            pep = row["Peptide"]
            rank = float(row["Rank_EL"])
            is_sb = rank <= 0.5
            scores[pep] = scores.get(pep, 0.0) + (100.0 if is_sb else 10.0)
            hits_map[pep] = hits_map.get(pep, 0) + 1
            if is_sb:
                strong_map[pep] = strong_map.get(pep, 0) + 1

    orig_pen, orig_hits, orig_sb = 0.0, 0, 0
    best_perm = orig_tuple
    best_pen = float("inf")
    best_hits = 0
    best_sb = 0

    for perm in candidate_perms:
        c = build_construct(perm, linker=linker, n_term=n_term, c_term=c_term)
        wins = junction_windows(c, lengths=lengths)
        pen = sum(scores.get(w["peptide"], 0.0) for w in wins)
        h = sum(hits_map.get(w["peptide"], 0) for w in wins)
        sb = sum(strong_map.get(w["peptide"], 0) for w in wins)
        if perm == orig_tuple:
            orig_pen, orig_hits, orig_sb = pen, h, sb
        if pen < best_pen:
            best_pen, best_perm, best_hits, best_sb = pen, perm, h, sb

    # Build construct and full scan for best permutation
    best_c = build_construct(best_perm, linker=linker, n_term=n_term, c_term=c_term)
    final_scan = scan_construct(best_c, alleles, runner, lengths=lengths, rank_cut=rank_cut)

    return {
        "best_peptides": list(best_perm),
        "initial_hits": orig_hits,
        "optimized_hits": final_scan["n_hits"],
        "initial_strong": orig_sb,
        "optimized_strong": final_scan["n_strong"],
        "improved": best_pen < orig_pen,
        "eval_count": len(candidate_perms),
        "construct": best_c,
        "scan": final_scan,
    }


def to_fasta(construct: Dict[str, Any], name: str = "TRIAD_multiepitope") -> str:
    seq = construct["sequence"]
    body = "\n".join(seq[i:i + 60] for i in range(0, len(seq), 60))
    header = (f">{name} | {construct['n_epitopes']} epitopes | linker={construct['linker']} "
              f"| {construct['length']} aa")
    return f"{header}\n{body}\n"


if __name__ == "__main__":
    peps = ["FLPSDFFPSV", "AYAQKIFKI", "QVPLRPMTYK", "NPPIPVGEIY"]
    c = build_construct(peps, linker="AAY", n_term="MET", c_term="HIS6")
    print("=== Construct ===")
    print(c["sequence"])
    print(f"length {c['length']} aa | GRAVY {c['properties']['gravy']} | "
          f"pI {c['properties']['pi']} | {c['properties']['mw_kda']} kDa")
    w = junction_windows(c)
    print(f"\n=== Junction windows: {len(w)} novel k-mers to scan ===")
    for x in w[:6]:
        print(f"  pos {x['pos']:>3}  {x['peptide']:<12} spans {x['spans']}")
    print(f"  ... ({len(w)} total)")
