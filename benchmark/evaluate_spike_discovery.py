#!/usr/bin/env python3
"""
evaluate_spike_discovery.py
Statistical evaluation and publication figure generation for the Whole-Protein
SARS-CoV-2 Spike Epitope Discovery Benchmark (1,265 overlapping 9-mers).
"""

import json
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from pathlib import Path
from sklearn.metrics import (
    roc_curve, auc,
    precision_recall_curve, average_precision_score
)

# Styling for Nature/Cell grade aesthetics
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
    "font.size": 11,
    "axes.labelsize": 12,
    "axes.titlesize": 13,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "legend.fontsize": 9.5,
    "figure.titlesize": 14,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "grid.alpha": 0.3,
    "grid.linestyle": "--"
})

COLORS = {
    "triad": "#2563eb",         # Royal Blue
    "presentation": "#d97706",   # Amber Orange
    "affinity": "#9333ea",       # Purple
    "immunogenicity": "#059669"  # Emerald Green
}

def evaluate(df, outdir="benchmark/results"):
    outdir_p = Path(outdir)
    outdir_p.mkdir(parents=True, exist_ok=True)

    y_true = df["Is_Verified_Epitope"].values
    n_pos = int(y_true.sum())
    n_total = len(y_true)
    baseline_prec = n_pos / n_total

    print(f"[*] Total evaluated 9-mers: {n_total}")
    print(f"[*] Verified CD8+ Epitopes: {n_pos} ({baseline_prec*100:.2f}% prevalence)")

    models = [
        {
            "id": "triad_composite",
            "name": "TRIAD 4D Composite Score",
            "scores": df["Composite_Score"].values,
            "color": COLORS["triad"]
        },
        {
            "id": "netmhcpan_el",
            "name": "NetMHCpan-4.2 Presentation (Score_EL)",
            "scores": df["Score_EL"].values,
            "color": COLORS["presentation"]
        },
        {
            "id": "netmhcpan_ba",
            "name": "NetMHCpan-4.2 Binding Affinity (Score_BA)",
            "scores": df["Score_BA"].values,
            "color": COLORS["affinity"]
        },
        {
            "id": "tcr_immunogenicity",
            "name": "TCR Immunogenicity (Calis et al.)",
            "scores": df["Immunogenicity_Score"].values,
            "color": COLORS["immunogenicity"]
        }
    ]

    metrics = {}
    for m in models:
        scores = m["scores"]
        fpr, tpr, _ = roc_curve(y_true, scores)
        roc_auc = auc(fpr, tpr)

        prec, rec, _ = precision_recall_curve(y_true, scores)
        pr_auc = average_precision_score(y_true, scores)

        sorted_idx = np.argsort(scores)[::-1]
        p10 = float(np.mean(y_true[sorted_idx[:10]]))
        r10 = float(np.sum(y_true[sorted_idx[:10]]) / n_pos)

        p20 = float(np.mean(y_true[sorted_idx[:20]]))
        r20 = float(np.sum(y_true[sorted_idx[:20]]) / n_pos)

        p50 = float(np.mean(y_true[sorted_idx[:50]]))
        r50 = float(np.sum(y_true[sorted_idx[:50]]) / n_pos)

        metrics[m["id"]] = {
            "name": m["name"],
            "roc_auc": float(roc_auc),
            "pr_auc": float(pr_auc),
            "top10_precision": p10,
            "top10_recall": r10,
            "top20_precision": p20,
            "top20_recall": r20,
            "top50_precision": p50,
            "top50_recall": r50,
            "fpr": fpr.tolist(),
            "tpr": tpr.tolist(),
            "precision": prec.tolist(),
            "recall": rec.tolist(),
            "color": m["color"]
        }
        print(f"\n[{m['name']}]")
        print(f"   ROC-AUC: {roc_auc:.4f} | PR-AUC: {pr_auc:.4f}")
        print(f"   Top 10: Prec={p10*100:.1f}%, Rec={r10*100:.1f}% ({int(p10*10)}/{min(10, n_pos)} found)")
        print(f"   Top 20: Prec={p20*100:.1f}%, Rec={r20*100:.1f}% ({int(p20*20)}/{n_pos} found)")
        print(f"   Top 50: Prec={p50*100:.1f}%, Rec={r50*100:.1f}% ({int(p50*50)}/{n_pos} found)")

    # Tier Enrichment Analysis
    tier_stats = {}
    for tier in ["Tier 1", "Tier 2", "Tier 3"]:
        sub = df[df["Tier"] == tier]
        total = len(sub)
        pos = int(sub["Is_Verified_Epitope"].sum())
        prec = pos / total if total > 0 else 0.0
        rec = pos / n_pos
        tier_stats[tier] = {
            "total_candidates": total,
            "verified_epitopes": pos,
            "precision_pct": prec * 100,
            "recall_pct": rec * 100
        }
        print(f"\n[{tier}] Total={total} | True Epitopes={pos} | Prec={prec*100:.1f}% | Recall={rec*100:.1f}%")

    # ── Figure Generation ──────────────────────────────────────────────────
    
    # 1. Multi-panel Figure: ROC, PR, Top-K, and Protein Landscape
    fig, axes = plt.subplots(2, 2, figsize=(14, 11), dpi=300)

    # Subplot A: ROC
    ax = axes[0, 0]
    ax.plot([0, 1], [0, 1], linestyle="--", color="#94a3b8", label="Random (AUC = 0.50)", lw=1.2)
    for m_id, res in metrics.items():
        lw = 2.4 if m_id == "triad_composite" else 1.5
        ax.plot(res["fpr"], res["tpr"], color=res["color"], lw=lw,
                label=f"{res['name']} ({res['roc_auc']:.3f})")
    ax.set_title("A. ROC Curves (Spike Epitope Discovery)", loc="left", fontweight="bold")
    ax.set_xlabel("1 - Specificity (False Positive Rate)")
    ax.set_ylabel("Sensitivity (True Positive Rate)")
    ax.legend(loc="lower right")

    # Subplot B: PR Curve
    ax = axes[0, 1]
    ax.axhline(y=baseline_prec, color="#94a3b8", linestyle="--", label=f"Random Prevalence ({baseline_prec*100:.1f}%)", lw=1.2)
    for m_id, res in metrics.items():
        lw = 2.4 if m_id == "triad_composite" else 1.5
        ax.plot(res["recall"], res["precision"], color=res["color"], lw=lw,
                label=f"{res['name']} (PR-AUC = {res['pr_auc']:.3f})")
    ax.set_title("B. Precision-Recall Curves", loc="left", fontweight="bold")
    ax.set_xlabel("Recall")
    ax.set_ylabel("Precision")
    ax.legend(loc="upper right")

    # Subplot C: Precision at Top-K
    ax = axes[1, 0]
    k_labels = ["Top 10", "Top 20", "Top 50"]
    x = np.arange(len(k_labels))
    width = 0.2
    for i, (m_id, res) in enumerate(metrics.items()):
        vals = [res["top10_precision"] * 100, res["top20_precision"] * 100, res["top50_precision"] * 100]
        ax.bar(x + (i - 1.5) * width, vals, width, label=res["name"], color=res["color"])
    ax.set_xticks(x)
    ax.set_xticklabels(k_labels)
    ax.set_ylabel("Precision (%)")
    ax.set_ylim([0, 105])
    ax.set_title("C. Precision at Top-K Candidates", loc="left", fontweight="bold")
    ax.legend(loc="upper right")

    # Subplot D: Tier 1/2/3 Breakdown
    ax = axes[1, 1]
    tiers = list(tier_stats.keys())
    t_precs = [tier_stats[t]["precision_pct"] for t in tiers]
    t_recs = [tier_stats[t]["recall_pct"] for t in tiers]
    t_counts = [tier_stats[t]["total_candidates"] for t in tiers]
    t_pos = [tier_stats[t]["verified_epitopes"] for t in tiers]

    b = ax.bar(tiers, t_precs, color=["#10b981", "#3b82f6", "#94a3b8"], width=0.5, edgecolor="#1e293b")
    for rect, count, pos, rec in zip(b, t_counts, t_pos, t_recs):
        h = rect.get_height()
        ax.annotate(f"{h:.1f}%\n({pos}/{count})\nRec={rec:.0f}%",
                    xy=(rect.get_x() + rect.get_width()/2, h),
                    xytext=(0, 3), textcoords="offset points", ha="center", va="bottom", fontsize=8.5)
    ax.set_ylim([0, 60])
    ax.set_ylabel("Precision / Hit Rate (%)")
    ax.set_title("D. TRIAD Decision Tier Precision & Recall", loc="left", fontweight="bold")

    plt.suptitle("Whole-Protein SARS-CoV-2 Spike Epitope Discovery Benchmark (1,265 K-mers)", fontsize=15, fontweight="bold", y=0.995)
    plt.tight_layout()
    plt.savefig(outdir_p / "fig_spike_discovery_composite.png")
    plt.savefig(outdir_p / "fig_spike_discovery_composite.svg")
    plt.close()
    print("[✓] Saved fig_spike_discovery_composite.png and .svg")

    # 2. Sequence Landscape Track Plot
    plt.figure(figsize=(15, 5), dpi=300)
    plt.plot(df["Protein_Pos"], df["Composite_Score"], color="#cbd5e1", lw=0.9, label="Background Peptides (n=1210)")
    
    # Highlight true positive epitopes
    pos_mask = df["Is_Verified_Epitope"] == 1
    plt.scatter(df.loc[pos_mask, "Protein_Pos"], df.loc[pos_mask, "Composite_Score"],
                color="#dc2626", s=45, zorder=4, label="IEDB Confirmed CD8+ Epitopes (n=55)", edgecolors="#7f1d1d")
    
    # Tier 1 threshold line
    plt.axhline(y=0.65, color="#10b981", linestyle="--", lw=1.5, label="Tier 1 Threshold (Composite >= 0.65)")
    
    # Annotate top famous epitopes
    famous = ["YLQPRTFLL", "FIAGLIAIV", "RLQSLQTYV", "VVFLHVTYV", "KLPDDFTGC"]
    for fam in famous:
        match = df[df["Peptide"] == fam]
        if not match.empty:
            r = match.iloc[0]
            plt.annotate(fam, xy=(r["Protein_Pos"], r["Composite_Score"]),
                         xytext=(0, 10), textcoords="offset points",
                         ha="center", fontsize=8, fontweight="bold", color="#991b1b",
                         arrowprops=dict(arrowstyle="->", color="#991b1b", lw=0.8))

    plt.title("Whole-Protein Epitope Prioritization Landscape: SARS-CoV-2 Spike (HLA-A*02:01)", pad=12, fontweight="bold")
    plt.xlabel("Amino Acid Starting Position (1 ~ 1,265)")
    plt.ylabel("TRIAD Composite Score")
    plt.xlim([0, 1275])
    plt.ylim([0, 1.05])
    plt.legend(loc="upper right", frameon=True, facecolor="white", edgecolor="#e2e8f0")
    plt.tight_layout()
    plt.savefig(outdir_p / "fig_spike_landscape.png")
    plt.savefig(outdir_p / "fig_spike_landscape.svg")
    plt.close()
    print("[✓] Saved fig_spike_landscape.png and .svg")

    # Save summary json
    summary = {
        "dataset": {
            "protein": "SARS-CoV-2 Spike Glycoprotein (P0DTC2)",
            "length_aa": 1273,
            "total_9mers": n_total,
            "verified_epitopes": n_pos,
            "background_peptides": n_total - n_pos
        },
        "models": {k: {
            "name": v["name"],
            "roc_auc": v["roc_auc"],
            "pr_auc": v["pr_auc"],
            "top10_precision": v["top10_precision"],
            "top10_recall": v["top10_recall"],
            "top20_precision": v["top20_precision"],
            "top20_recall": v["top20_recall"],
            "top50_precision": v["top50_precision"],
            "top50_recall": v["top50_recall"]
        } for k, v in metrics.items()},
        "tier_enrichment": tier_stats
    }
    with open(outdir_p / "spike_benchmark_metrics.json", "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)
    print(f"[✓] Spike benchmark summary saved to {outdir_p / 'spike_benchmark_metrics.json'}")

def main():
    parser = argparse.ArgumentParser(description="Evaluate Spike discovery benchmark")
    parser.add_argument("--pred", default="benchmark/results/spike_discovery_predictions.csv", help="Predictions CSV")
    parser.add_argument("--outdir", default="benchmark/results", help="Output directory")
    args = parser.parse_args()

    df = pd.read_csv(args.pred)
    evaluate(df, args.outdir)

if __name__ == "__main__":
    main()
