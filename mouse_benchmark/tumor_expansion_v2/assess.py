"""Locked external evaluation across all holdout test cohorts; no model selection or retraining."""
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score

from prepare import ROOT, overlap, sha
from train import verify_lock
REPO = ROOT.parent.parent
sys.path.insert(0, str(ROOT.parent))
from evaluate import top_k

MODELS = [
    'Calis', 'EL_Rank_Score', 'BA_Rank_Score', 'Mouse_Base', 'Mouse_Composite',
    'LR_binding', 'LR_tumor_primary', 'GB_tumor_exploratory', 'LR_tumor_expanded'
]


def metrics(part, model):
    y = part.Label.to_numpy()
    s = part[model].to_numpy()
    both = len(set(y)) == 2
    result = {
        'N': len(part),
        'Positive': int(y.sum()),
        'Prevalence': float(y.mean()) if len(y) else None,
        'AP': float(average_precision_score(y, s)) if both else None,
        'ROC_AUC': float(roc_auc_score(y, s)) if both else None
    }
    for k in [5, 10, 20]:
        prec, rec = top_k(y, s, k)
        result[f'Precision_at_{k}'] = prec
        result[f'Recall_at_{k}'] = rec
    return result


def similarity_groups(part):
    rows = part.to_dict('records')
    parent = list(range(len(rows)))

    def root(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    for i in range(len(rows)):
        for j in range(i):
            if overlap(rows[i], rows[j]):
                parent[root(i)] = root(j)
    return np.unique([root(i) for i in range(len(rows))], return_inverse=True)[1]


def uncertainty(part, primary_col='LR_tumor_primary'):
    # Paired sequence/source-group bootstrap within studies
    group_indices = []
    for _, sub in part.groupby('PMID'):
        groups = similarity_groups(sub)
        group_indices.append([part.index.get_indexer(sub.index[groups == g]) for g in np.unique(groups)])
    rng = np.random.default_rng(20261010)
    comparators = ['Mouse_Base', 'BA_Rank_Score', 'Mouse_Composite', 'LR_binding', 'LR_tumor_expanded']
    deltas = {m: [] for m in comparators if m in part.columns and m != primary_col}
    y = part.Label.to_numpy()
    for _ in range(1000):
        w = np.zeros(len(part))
        for blocks in group_indices:
            for g in rng.integers(0, len(blocks), size=len(blocks)):
                w[blocks[g]] += 1
        if len(set(y[w > 0])) != 2:
            continue
        ap_prim = average_precision_score(y, part[primary_col], sample_weight=w)
        for m in deltas:
            ap_comp = average_precision_score(y, part[m], sample_weight=w)
            deltas[m].append(float(ap_prim - ap_comp))
    return {
        'replicates': 1000,
        'valid_replicates': len(next(iter(deltas.values()))) if deltas else 0,
        'groups_per_study': [len(x) for x in group_indices],
        'method': 'paired sequence/source-group bootstrap, within fixed studies',
        'scope': 'conditional candidate sampling only; not new-study uncertainty',
        'primary_minus_comparator_AP_95': {
            m: np.quantile(v, [0.025, 0.975]).tolist() if len(v) >= 10 else None
            for m, v in deltas.items()
        }
    }


def assess():
    verify_lock()
    out = ROOT / 'run'
    marker = out / 'evaluation_started.json'
    if marker.exists():
        raise FileExistsError('External evaluation already opened; do not tune against it')
    manifest = json.loads((out / 'model_manifest.json').read_text())
    if manifest['lock_sha256'] != sha(ROOT / 'locked/split_lock.json') or manifest['train_code_sha256'] != sha(ROOT / 'train.py'):
        raise ValueError('Training checkpoint provenance mismatch')
    for name, h in manifest['model_files'].items():
        if sha(out / name) != h:
            raise ValueError(f'Model artifact changed: {name}')

    df = pd.read_csv(out / 'scored_cohort.csv', dtype={'PMID': str})
    if not np.isfinite(df[MODELS].to_numpy()).all():
        raise ValueError('Missing prediction values in scored cohort')

    marker.write_text(json.dumps({
        'opened_utc': datetime.now(timezone.utc).isoformat(),
        'model_manifest_sha256': sha(out / 'model_manifest.json'),
        'scored_cohort_sha256': sha(out / 'scored_cohort.csv'),
        'assess_code_sha256': sha(Path(__file__))
    }, indent=2) + '\n')

    results = []

    # Evaluate validation, test_frozen_prior, test_capietto_holdout, and pooled holdouts
    eval_splits = ['validation', 'test_frozen_prior', 'test_capietto_holdout']
    for split in eval_splits:
        part = df[df.Split == split].copy()
        subsets = {'pooled': part}
        subsets.update({a: p for a, p in part.groupby('MHC')})
        subsets.update({f'PMID_{pmid}': p for pmid, p in part.groupby('PMID')})
        for name, p in subsets.items():
            for model in MODELS:
                results.append({'Split': split, 'Cohort': name, 'Model': model, **metrics(p, model)})

    # Also evaluate combined holdout (Prior test + Capietto holdout, N=69+216=285)
    combined_test = df[df.Split.isin(['test_frozen_prior', 'test_capietto_holdout'])].copy()
    for model in MODELS:
        results.append({'Split': 'test_combined', 'Cohort': 'pooled_all_holdouts', 'Model': model, **metrics(combined_test, model)})

    table = pd.DataFrame(results)
    table.to_csv(out / 'external_metrics.csv', index=False)

    # Cross-validation metrics
    cv = pd.read_csv(out / 'development_out_of_fold.csv')
    cv_metrics = []
    for (tumor, model), p in cv.groupby(['Tumor_Model', 'Model']):
        cv_metrics.append({'Tumor_Model': tumor, 'Model': model, **metrics(p.rename(columns={'Score': model}), model)})
    pd.DataFrame(cv_metrics).to_csv(out / 'development_metrics.csv', index=False)

    # Uncertainty calculations
    prior_test = df[df.Split == 'test_frozen_prior'].copy()
    cap_test = df[df.Split == 'test_capietto_holdout'].copy()

    intervals = {
        'test_frozen_prior_pooled': uncertainty(prior_test, 'LR_tumor_primary'),
        'test_frozen_prior_expanded_model': uncertainty(prior_test, 'LR_tumor_expanded'),
        'test_capietto_holdout': uncertainty(cap_test, 'LR_tumor_primary'),
        'test_capietto_holdout_expanded_model': uncertainty(cap_test, 'LR_tumor_expanded'),
    }
    (out / 'uncertainty.json').write_text(json.dumps(intervals, indent=2) + '\n')

    # Print summary tables
    print('\n=== External Test Metrics: Prior Test (LLC1 + CMT93, N=69, Prev=0.710) ===')
    t_prior = table[(table.Split == 'test_frozen_prior') & (table.Cohort == 'pooled')]
    print(t_prior[['Model', 'N', 'Positive', 'Prevalence', 'AP', 'ROC_AUC', 'Precision_at_10', 'Recall_at_10']].to_string(index=False))

    print('\n=== External Test Metrics: Capietto Holdout (MC-38/TRAMP-C1, N=216, Prev=0.102) ===')
    t_cap = table[(table.Split == 'test_capietto_holdout') & (table.Cohort == 'pooled')]
    print(t_cap[['Model', 'N', 'Positive', 'Prevalence', 'AP', 'ROC_AUC', 'Precision_at_10', 'Recall_at_10']].to_string(index=False))

    print('\n=== External Test Metrics: Combined Holdouts (N=285) ===')
    t_comb = table[(table.Split == 'test_combined')]
    print(t_comb[['Model', 'N', 'Positive', 'Prevalence', 'AP', 'ROC_AUC', 'Precision_at_10', 'Recall_at_10']].to_string(index=False))

    (out / 'evaluation_completed.json').write_text(json.dumps({
        'completed_utc': datetime.now(timezone.utc).isoformat(),
        'metrics_sha256': sha(out / 'external_metrics.csv')
    }, indent=2) + '\n')


if __name__ == '__main__':
    assess()
