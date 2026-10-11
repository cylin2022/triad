"""Frozen baseline comparisons with study-aware uncertainty and sensitivity checks."""
import json
from pathlib import Path

import numpy as np
import pandas as pd
import sklearn
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.model_selection import StratifiedGroupKFold

ROOT = Path(__file__).resolve().parent
MODELS = ["Calis", "EL_Rank_Score", "BA_Rank_Score", "Score_EL", "Mouse_Base", "Mouse_Composite"]
SEED = 20261009


def top_k(y, scores, k):
    """Expected precision under random tie-breaking at the boundary."""
    y, scores = np.asarray(y), np.asarray(scores)
    if len(y) < k:
        return None, None
    boundary = np.sort(scores)[-k]
    above, tied = scores > boundary, scores == boundary
    positives = y[above].sum() + (k - above.sum()) * y[tied].mean()
    return float(positives / k), float(positives / y.sum()) if y.sum() else None


def metrics(df, model):
    y, scores = df.Label.to_numpy(), df[model].to_numpy()
    result = {"n": len(df), "positive": int(y.sum()), "prevalence": float(y.mean()) if len(y) else None}
    both = len(set(y)) == 2
    result["AP"] = float(average_precision_score(y, scores)) if both else None
    result["ROC_AUC"] = float(roc_auc_score(y, scores)) if both else None
    for k in [20, 50]:
        result[f"Precision_at_{k}"], result[f"Recall_at_{k}"] = top_k(y, scores, k)
    return result


def bootstrap(df, repetitions=1000):
    groups, group_index = np.unique(df.Leakage_Group.to_numpy(), return_inverse=True)
    rng = np.random.default_rng(SEED)
    yy = df.Label.to_numpy()
    values = {m: [] for m in MODELS}
    deltas = {m: [] for m in ["Mouse_Base", "EL_Rank_Score", "BA_Rank_Score"]}
    for _ in range(repetitions):
        multiplicities = np.bincount(rng.integers(0, len(groups), len(groups)), minlength=len(groups))
        weights = multiplicities[group_index]
        if len(np.unique(yy[weights > 0])) < 2:
            continue
        for m in MODELS:
            values[m].append(float(average_precision_score(yy, df[m].to_numpy(), sample_weight=weights)))
        for m in deltas:
            deltas[m].append(values["Mouse_Composite"][-1] - values[m][-1])
    def interval(vals):
        return [float(v) for v in np.quantile(vals, [.025, .975])] if vals else None
    return {"method": "paired reference/sequence-component bootstrap", "seed": SEED,
            "requested_replicates": repetitions, "valid_replicates": len(values[MODELS[0]]),
            "groups": len(groups), "largest_group_fraction": float(df.Leakage_Group.value_counts().iloc[0] / len(df)),
            "AP_95_intervals": {m: interval(v) for m,v in values.items()},
            "composite_minus_baseline_AP_95_intervals": {m: interval(v) for m,v in deltas.items()}}


def evaluate():
    out = ROOT / "results"
    df = pd.read_csv(out / "scored_pairs.csv", dtype={"Reference_IDs": str})
    if df.duplicated(["Peptide", "MHC"]).any() or not np.isfinite(df[MODELS].to_numpy()).all():
        raise ValueError("Invalid scored cohort")
    studies = pd.read_csv(out / "studies.csv", dtype={"Reference_ID": str})
    largest_ref = studies.iloc[0].Reference_ID
    has_ref = lambda ref: df.Reference_IDs.map(lambda v: str(ref) in v.split(";"))
    # Established before inspecting scores: direct cross-reactivity/alloreactivity
    # study titles flag a different biological endpoint, not generic immunogenicity.
    cross_refs = set(studies.loc[studies.Title.str.contains(r"cross.react|alloreact", case=False, regex=True), "Reference_ID"])
    cross_mask = df.Reference_IDs.map(lambda v: bool(set(v.split(";")) & cross_refs))
    cohorts = {
        "primary": df,
        "H-2-Kb": df[df.MHC == "H-2-Kb"],
        "H-2-Db": df[df.MHC == "H-2-Db"],
        "9mer": df[df.Length == 9],
        "stronger_restriction": df[df.Strong_Restriction],
        "EL_rank_le_2": df[df.Rank_EL <= 2],
        "post2019_only": df[(df.Earliest_Year >= 2020) & df.All_Dates_Known],
        "without_largest_study": df[~has_ref(largest_ref)],
        "without_crossreactivity_titles": df[~cross_mask],
        "without_largest_leakage_group": df[df.Leakage_Group != df.Leakage_Group.value_counts().index[0]],
    }
    for allele in ["H-2-Kb", "H-2-Db"]:
        cohorts[allele + "_stronger_restriction"] = df[(df.MHC == allele) & df.Strong_Restriction]
    results = [{"Cohort": name, "Model": model, **metrics(part, model)}
               for name,part in cohorts.items() for model in MODELS]
    pd.DataFrame(results).to_csv(out / "metrics.csv", index=False)
    study_metrics = []
    for study in studies.itertuples():
        part = df[has_ref(study.Reference_ID)]
        if (part.Label == 1).sum() < 5 or (part.Label == 0).sum() < 5:
            continue
        for model in MODELS:
            study_metrics.append({"Reference_ID": study.Reference_ID, "Model": model, **metrics(part, model)})
    sm = pd.DataFrame(study_metrics)
    sm.to_csv(out / "within_study_metrics.csv", index=False)
    macro = sm.groupby("Model")[["AP", "ROC_AUC", "Precision_at_20", "Precision_at_50"]].mean()
    macro.to_csv(out / "study_macro_metrics.csv")
    print(pd.DataFrame(results).query("Cohort == 'primary'").to_string(index=False), flush=True)
    print("Computing 1,000 paired component bootstrap replicates...", flush=True)
    uncertainty = bootstrap(df)
    (out / "uncertainty.json").write_text(json.dumps(uncertainty, indent=2) + "\n")
    # Diagnostic only: no models are trained. Exposes giant connected components.
    folds = []
    splitter = StratifiedGroupKFold(n_splits=5, shuffle=True, random_state=SEED)
    for i, (train, test) in enumerate(splitter.split(df, df.Label, df.Leakage_Group)):
        assert not set(df.iloc[train].Leakage_Group) & set(df.iloc[test].Leakage_Group)
        part = df.iloc[test]
        folds.append({"fold": i, "n": len(part), "positive": int(part.Label.sum()),
                      "negative": int((part.Label == 0).sum()), "groups": part.Leakage_Group.nunique()})
    (out / "training_split_audit.json").write_text(json.dumps({"diagnostic_only_no_training": True, "folds": folds}, indent=2) + "\n")
    (out / "evaluation_manifest.json").write_text(json.dumps({
        "sklearn": sklearn.__version__, "seed": SEED, "largest_reference": largest_ref,
        "within_study_eligible_studies": sm.Reference_ID.nunique(), "crossreactivity_title_references": sorted(cross_refs),
        "crossreactivity_title_pairs": int(cross_mask.sum()),
    }, indent=2) + "\n")
    print(json.dumps(uncertainty, indent=2), flush=True)


if __name__ == "__main__":
    evaluate()
