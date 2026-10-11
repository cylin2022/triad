"""Generate PR and ROC curves for v2 expansion benchmarks."""
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import precision_recall_curve, roc_curve, average_precision_score, roc_auc_score

ROOT = Path(__file__).resolve().parent
RUN = ROOT / 'run'

df = pd.read_csv(RUN / 'scored_cohort.csv')

MODELS_TO_PLOT = [
    ('BA_Rank_Score', 'NetMHCpan BA Rank', '#1f77b4', '-'),
    ('EL_Rank_Score', 'NetMHCpan EL Rank', '#aec7e8', '--'),
    ('Mouse_Composite', 'Current Composite', '#2ca02c', '-.'),
    ('Mouse_Base', 'Mouse Base', '#98df8a', ':'),
    ('LR_binding', 'LR Binding', '#ff7f0e', '-'),
    ('LR_tumor_primary', 'LR Tumor Primary (v2)', '#d62728', '-'),
    ('LR_tumor_expanded', 'LR Tumor Expanded', '#9467bd', '--'),
]

cohorts = [
    ('test_frozen_prior', 'Prior Holdout (LLC1 + CMT93, N=69, Prev=71%)', 'pr_roc_prior_test.png'),
    ('test_capietto_holdout', 'Capietto Holdout (MC-38/TRAMP-C1, N=216, Prev=10%)', 'pr_roc_capietto_holdout.png'),
]

for split_name, title_suffix, out_name in cohorts:
    part = df[df.Split == split_name].copy()
    y = part.Label.to_numpy()
    
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5), dpi=200)
    
    # 1. PR curve
    ax_pr = axes[0]
    for col, label, color, ls in MODELS_TO_PLOT:
        s = part[col].to_numpy()
        p, r, _ = precision_recall_curve(y, s)
        ap = average_precision_score(y, s)
        ax_pr.plot(r, p, label=f'{label} (AP={ap:.3f})', color=color, linestyle=ls, linewidth=1.8)
    
    # Baseline prevalence
    prev = y.mean()
    ax_pr.axhline(prev, color='gray', linestyle=':', label=f'Prevalence ({prev:.3f})')
    ax_pr.set_title(f'Precision-Recall: {title_suffix}', fontsize=11, fontweight='bold')
    ax_pr.set_xlabel('Recall', fontsize=10)
    ax_pr.set_ylabel('Precision', fontsize=10)
    ax_pr.set_xlim([0, 1.02])
    ax_pr.set_ylim([0, 1.02])
    ax_pr.grid(True, alpha=0.3)
    ax_pr.legend(loc='lower left', fontsize=8, framealpha=0.9)
    
    # 2. ROC curve
    ax_roc = axes[1]
    for col, label, color, ls in MODELS_TO_PLOT:
        s = part[col].to_numpy()
        fpr, tpr, _ = roc_curve(y, s)
        auc = roc_auc_score(y, s)
        ax_roc.plot(fpr, tpr, label=f'{label} (AUC={auc:.3f})', color=color, linestyle=ls, linewidth=1.8)
    
    ax_roc.plot([0, 1], [0, 1], 'k--', alpha=0.5, label='Random Chance')
    ax_roc.set_title(f'ROC: {title_suffix}', fontsize=11, fontweight='bold')
    ax_roc.set_xlabel('False Positive Rate', fontsize=10)
    ax_roc.set_ylabel('True Positive Rate', fontsize=10)
    ax_roc.set_xlim([0, 1.02])
    ax_roc.set_ylim([0, 1.02])
    ax_roc.grid(True, alpha=0.3)
    ax_roc.legend(loc='lower right', fontsize=8, framealpha=0.9)
    
    plt.tight_layout()
    plt.savefig(RUN / out_name)
    plt.close()
    print(f'Saved plot to {RUN / out_name}')
