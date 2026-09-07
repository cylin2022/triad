"""
oss_binding_predictors.py
Open-source-licensed MHC-I binding predictors extracted from IEDB's ng_tc1
toolkit: SMM / SMM-PMBEC (Peters lab position-weight matrices) and the
comblib_sidney2008 combinatorial-library matrices.

These are IEDB's own methods, licensed under the Non-Profit Open Software
License 3.0 (see THIRD_PARTY_LICENSES/license-LJI.txt) — unlike netMHCpan-4.1,
netMHC-4.0, netChop, netCTL and netCTLpan (also bundled inside ng_tc1), they
are pure Python + small pre-trained scoring matrices with no DTU-licensed
binary and no external dependency, so they can be redistributed and run
in-process here.

Data (~2.5MB total) is bundled at oss_predictors_data/{smm,comblib}/, copied
directly from ng_tc1-0.1.5-beta/method/{smm-predictor,mhci-comblib-predictor}.

Not included: MHCflurry (also open, Apache 2.0) — its TensorFlow dependency
requires a CPU with AVX support, which this deployment host's CPU lacks
(confirmed by testing: TensorFlow crashes with "Illegal instruction" here).
Not included: MHCnp — ng_tc1's implementation only runs via a private LJI
Docker registry image (harbor.lji.org), not a local/redistributable method.
"""

import re
import pickle
from pathlib import Path
from typing import List, Optional

DATA_DIR = Path(__file__).parent / "oss_predictors_data"
SMM_DIRS = {
    "smm": DATA_DIR / "smm" / "smm",
    "smmpmbec": DATA_DIR / "smm" / "smmpmbec",
}
COMBLIB_PICKLE = DATA_DIR / "comblib" / "comblib_sidney2008" / "dic_pssm_sidney2008.cPickle"

AA_LIST = "ACDEFGHIKLMNPQRSTVWY"

_comblib_cache: Optional[dict] = None
_smm_matrix_cache: dict = {}


def _smm_filename(allele: str, length: int) -> str:
    # "HLA-A*02:01" -> "HLA-A-0201-9.cpickle" (matches ng_tc1's own naming scheme)
    name = allele.replace("*", "-").replace(" ", "-").replace(":", "")
    return f"{name}-{length}.cpickle"


def _load_smm_matrix(allele: str, length: int, method: str):
    cache_key = (method, allele, length)
    if cache_key in _smm_matrix_cache:
        return _smm_matrix_cache[cache_key]

    d = SMM_DIRS.get(method)
    path = d / _smm_filename(allele, length) if d else None
    if path is None or not path.exists():
        _smm_matrix_cache[cache_key] = None
        return None

    with open(path, "rb") as f:
        mat_length = pickle.load(f)
        mat = pickle.load(f)
        offset = pickle.load(f)
    loaded = (mat_length, mat, offset)
    _smm_matrix_cache[cache_key] = loaded
    return loaded


def is_smm_supported(allele: str, length: int, method: str = "smm") -> bool:
    return _load_smm_matrix(allele, length, method) is not None


def predict_smm(peptides: List[str], allele: str, method: str = "smm") -> List[Optional[float]]:
    """IC50 (nM) predictions via SMM or SMM-PMBEC position-weight matrices."""
    results = []
    for pep in peptides:
        pep = pep.strip().upper()
        loaded = _load_smm_matrix(allele, len(pep), method)
        if loaded is None:
            results.append(None)
            continue
        _, mat, offset = loaded
        try:
            score = offset
            for pos, aa in enumerate(pep):
                score += mat[aa][pos]
            results.append(round(10 ** score, 2))
        except (KeyError, IndexError):
            results.append(None)
    return results


def _load_comblib() -> dict:
    global _comblib_cache
    if _comblib_cache is None:
        with open(COMBLIB_PICKLE, "rb") as f:
            raw = pickle.load(f)
        # ng_tc1 applies factor=-1.0 for the comblib_sidney2008 library before scoring
        _comblib_cache = {k: [-1.0 * v for v in vals] for k, vals in raw.items()}
    return _comblib_cache


def _comblib_key(allele: str) -> str:
    # "HLA-B*35:01" -> "HLA_B-3501" (matches ng_tc1's own key convention)
    return allele.replace("-", "_").replace("*", "-").replace(":", "")


def is_comblib_supported(allele: str, length: int) -> bool:
    return (_comblib_key(allele), length) in _load_comblib()


def predict_comblib(peptides: List[str], allele: str) -> List[Optional[float]]:
    """
    comblib_sidney2008 combinatorial-library score predictions.

    Note: unlike SMM, this is NOT an IC50 (nM) value — comblib_sidney2008's
    own get_score_unit() reports plain "Score". Values are small positive
    floats (typically ~1e-6 to ~1e-2); lower = stronger predicted binder,
    same convention as IC50, just a different scale — rank/compare rather
    than reading it as a physical concentration.
    """
    dic_pssm = _load_comblib()
    key_allele = _comblib_key(allele)
    mat_cache: dict = {}
    results = []
    for pep in peptides:
        pep = pep.strip().upper()
        length = len(pep)
        key = (key_allele, length)
        if key not in mat_cache:
            w = dic_pssm.get(key)
            if w is None:
                mat_cache[key] = None
            else:
                offset = w[0]
                mat = {
                    aa: [w[1 + 20 * pos + aa_index] for pos in range(length)]
                    for aa_index, aa in enumerate(AA_LIST)
                }
                mat_cache[key] = (mat, offset)
        loaded = mat_cache[key]
        if loaded is None:
            results.append(None)
            continue
        mat, offset = loaded
        try:
            score = offset
            for pos, aa in enumerate(pep):
                score += mat[aa][pos]
            results.append(float(f"{10 ** score:.6g}"))
        except (KeyError, IndexError):
            results.append(None)
    return results
