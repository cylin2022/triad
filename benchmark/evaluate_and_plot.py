#!/usr/bin/env python3
"""
evaluate_and_plot.py
Compute statistical benchmark metrics (ROC-AUC, PR-AUC, Precision@k, Tier enrichment)
and generate publication-grade multi-panel figures for manuscript inclusion.
"""

import os
import json
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from pathlib import Path
from sklearn.metrics import (
    roc_curve, auc,
    precision_recall_curve, average_precision_score,
    confusion_matrix, precision_score, recall_score, f1_score
)

# Publication aesthetic styling
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
    "font.size": 11,
    "axes.labelsize": 12,
    "axes.titlesize": 13,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "legend.fontsize": 10,
    "figure.titlesize": 14,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "grid.alpha": 0.3,
    "grid.linestyle": "--"
})

COLORS = {
    "triad": "#2563eb",       # Royal Blue
    "immunogenicity": "#059669", # Emerald Green
    "presentation": "#d97706", # Amber Orange
    "affinity": "#9333ea"      # Purple
}

def compute_metrics(df):
    y_true = df["Label"].values
    
    # Models to evaluate (name, score_array, color, display_label)
    models = [
        {
            "id": "triad_composite",
            "name": "TRIAD Composite Score",
            "scores": df["Composite_Score"].values,
            "color": COLORS["triad"]
        },
        {
            "id": "tcr_immunogenicity",
            "name": "TCR Immunogenicity (Calis et al.)",
            "scores": df["Immunogenicity_Score"].values,
            "color": COLORS["immunogenicity"]
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
        }
    ]

    results = {}
    for m in models:
        scores = m["scores"]
        fpr, tpr, _ = roc_curve(y_true, scores)
        roc_auc = auc(fpr, tpr)
        
        prec, rec, _ = precision_recall_curve(y_true, scores)
        pr_auc = average_precision_score(y_true, scores)
        
        # Precision @ top K
        sorted_indices = np.argsort(scores)[::-1]
        top10_prec = np.mean(y_true[sorted_indices[:10]])
        top20_prec = np.mean(y_true[sorted_indices[:20]])
        top50_prec = np.mean(y_true[sorted_indices[:50]])

        results[m["id"]] = {
            "name": m["name"],
            "roc_auc": float(roc_auc),
            "pr_auc": float(pr_auc),
            "top10_precision": float(top10_prec),
            "top20_precision": float(top20_prec),
            "top50_precision": float(top50_prec),
            "fpr": fpr.tolist(),
            "tpr": tpr.tolist(),
            "precision": prec.tolist(),
            "recall": rec.tolist(),
            "color": m["color"]
        }
        print(f"[Model: {m['name']}]")
        print(f"   ROC-AUC: {roc_auc:.4f} | PR-AUC: {pr_auc:.4f} | P@10: {top10_prec:.2f} | P@20: {top20_prec:.2f} | P@50: {top50_prec:.2f}")

    # Tier analysis
    tier_stats = {}
    if "Tier" in df.columns:
        for tier in ["Tier 1", "Tier 2", "Tier 3"]:
            sub = df[df["Tier"] == tier]
            total = len(sub)
            pos = int(sub["Label"].sum())
            rate = float(pos / total) if total > 0 else 0.0
            tier_stats[tier] = {
                "total": total,
                "positive": pos,
                "negative": total - pos,
                "positive_rate": rate
            }
        print("\n[Tier Enrichment Analysis]")
        for t, s in tier_stats.items():
            print(f"   {t}: {s['positive']}/{s['total']} ({s['positive_rate']*100:.1f}% Positive)")

    return results, tier_stats

def plot_curves(results, tier_stats, df, outdir):
    outdir_p = Path(outdir)
    outdir_p.mkdir(parents=True, exist_ok=True)

    # 1. ROC Curves
    plt.figure(figsize=(7, 6), dpi=300)
    plt.plot([0, 1], [0, 1], linestyle="--", color="#94a3b8", label="Random Classifier (AUC = 0.5000)", lw=1.2)
    for model_id, res in results.items():
        lw = 2.4 if model_id == "triad_composite" else 1.6
        plt.plot(res["fpr"], res["tpr"], color=res["color"], lw=lw,
                 label=f"{res['name']} (AUC = {res['roc_auc']:.3f})")
    plt.xlim([-0.02, 1.02])
    plt.ylim([-0.02, 1.02])
    plt.xlabel("False Positive Rate (1 - Specificity)")
    plt.ylabel("True Positive Rate (Sensitivity)")
    plt.title("Receiver Operating Characteristic (ROC) — IEDB Benchmark", pad=12)
    plt.legend(loc="lower right", frameon=True, facecolor="white", edgecolor="#e2e8f0")
    plt.tight_layout()
    plt.savefig(outdir_p / "fig_roc_curves.png")
    plt.savefig(outdir_p / "fig_roc_curves.svg")
    plt.close()
    print("[✓] Saved fig_roc_curves.png and .svg")

    # 2. Precision-Recall Curves
    plt.figure(figsize=(7, 6), dpi=300)
    baseline_prec = df["Label"].mean()
    plt.axhline(y=baseline_prec, color="#94a3b8", linestyle="--", label=f"Baseline Prevalence ({baseline_prec:.2f})", lw=1.2)
    for model_id, res in results.items():
        lw = 2.4 if model_id == "triad_composite" else 1.6
        plt.plot(res["recall"], res["precision"], color=res["color"], lw=lw,
                 label=f"{res['name']} (PR-AUC = {res['pr_auc']:.3f})")
    plt.xlim([-0.02, 1.02])
    plt.ylim([-0.02, 1.02])
    plt.xlabel("Recall (Sensitivity)")
    plt.ylabel("Precision (Positive Predictive Value)")
    plt.title("Precision-Recall (PR) Curves — IEDB Benchmark", pad=12)
    plt.legend(loc="lower left", frameon=True, facecolor="white", edgecolor="#e2e8f0")
    plt.tight_layout()
    plt.savefig(outdir_p / "fig_pr_curves.png")
    plt.savefig(outdir_p / "fig_pr_curves.svg")
    plt.close()
    print("[✓] Saved fig_pr_curves.png and .svg")

    # 3. Tier Enrichment Bar Chart
    if tier_stats:
        fig, ax1 = plt.subplots(figsize=(7, 5), dpi=300)
        tiers = list(tier_stats.keys())
        pos_rates = [tier_stats[t]["positive_rate"] * 100 for t in tiers]
        counts = [tier_stats[t]["total"] for t in tiers]

        bar_colors = ["#10b981", "#3b82f6", "#94a3b8"]
        bars = ax1.bar(tiers, pos_rates, color=bar_colors, width=0.55, edgecolor="#1e293b", linewidth=1.1, zorder=3)
        ax1.set_ylabel("True Positive Rate (%)", color="#1e293b", fontweight="bold")
        ax1.set_ylim([0, 100])
        ax1.grid(axis="y", linestyle="--", alpha=0.5, zorder=0)

        # Annotations on bars
        for bar, count, rate in zip(bars, counts, pos_rates):
            height = bar.get_height()
            ax1.annotate(f"{rate:.1f}%\n(n={count})",
                         xy=(bar.get_x() + bar.get_width() / 2, height),
                         xytext=(0, 4), textcoords="offset points",
                         ha="center", va="bottom", fontsize=10, fontweight="semibold")

        ax1.set_title("TRIAD Decision Tier Enrichment: Validated T-Cell Activation", pad=15)
        plt.tight_layout()
        plt.savefig(outdir_p / "fig_tier_enrichment.png")
        plt.savefig(outdir_p / "fig_tier_enrichment.svg")
        plt.close()
        print("[✓] Saved fig_tier_enrichment.png and .svg")

    # 4. Multi-panel Composite Summary Figure
    fig, axes = plt.subplots(2, 2, figsize=(13, 11), dpi=300)
    
    # Subplot A: ROC
    ax = axes[0, 0]
    ax.plot([0, 1], [0, 1], linestyle="--", color="#94a3b8", label="Random", lw=1.2)
    for model_id, res in results.items():
        lw = 2.2 if model_id == "triad_composite" else 1.4
        ax.plot(res["fpr"], res["tpr"], color=res["color"], lw=lw,
                label=f"{res['name']} ({res['roc_auc']:.3f})")
    ax.set_title("A. ROC Curves", loc="left", fontweight="bold")
    ax.set_xlabel("1 - Specificity")
    ax.set_ylabel("Sensitivity")
    ax.legend(fontsize=8.5, loc="lower right")

    # Subplot B: PR Curve
    ax = axes[0, 1]
    ax.axhline(y=baseline_prec, color="#94a3b8", linestyle="--", label="Baseline", lw=1.2)
    for model_id, res in results.items():
        lw = 2.2 if model_id == "triad_composite" else 1.4
        ax.plot(res["recall"], res["precision"], color=res["color"], lw=lw,
                label=f"{res['name']} ({res['pr_auc']:.3f})")
    ax.set_title("B. Precision-Recall Curves", loc="left", fontweight="bold")
    ax.set_xlabel("Recall")
    ax.set_ylabel("Precision")
    ax.legend(fontsize=8.5, loc="lower left")

    # Subplot C: Precision @ Top K
    ax = axes[1, 0]
    top_k_labels = ["Top 10", "Top 20", "Top 50"]
    x = np.arange(len(top_k_labels))
    width = 0.18
    for i, (m_id, res) in enumerate(results.items()):
        vals = [res["top10_precision"] * 100, res["top20_precision"] * 100, res["top50_precision"] * 100]
        ax.bar(x + (i - 1.5) * width, vals, width, label=res["name"], color=res["color"])
    ax.set_xticks(x)
    ax.set_xticklabels(top_k_labels)
    ax.set_ylabel("Precision (%)")
    ax.set_ylim([0, 105])
    ax.set_title("C. Precision at Top-K Candidates", loc="left", fontweight="bold")
    ax.legend(fontsize=8, loc="upper right")

    # Subplot D: Tier Distribution
    ax = axes[1, 1]
    if tier_stats:
        tiers = list(tier_stats.keys())
        rates = [tier_stats[t]["positive_rate"] * 100 for t in tiers]
        t_counts = [tier_stats[t]["total"] for t in tiers]
        b = ax.bar(tiers, rates, color=["#10b981", "#3b82f6", "#94a3b8"], width=0.5)
        for rect, c in zip(b, t_counts):
            h = rect.get_height()
            ax.annotate(f"{h:.1f}%\n(n={c})", xy=(rect.get_x() + rect.get_width()/2, h),
                        xytext=(0, 3), textcoords="offset points", ha="center", va="bottom", fontsize=8.5)
        ax.set_ylim([0, 100])
        ax.set_ylabel("Positive Activation Rate (%)")
        ax.set_title("D. Tier True Positive Enrichment", loc="left", fontweight="bold")

    plt.suptitle("TRIAD Platform Benchmark Validation (IEDB HLA-A*02:01)", fontsize=15, fontweight="bold", y=0.995)
    plt.tight_layout()
    plt.savefig(outdir_p / "fig_benchmark_composite.png")
    plt.savefig(outdir_p / "fig_benchmark_composite.svg")
    plt.close()
    print("[✓] Saved publication composite: fig_benchmark_composite.png and .svg")

def main():
    parser = argparse.ArgumentParser(description="Evaluate benchmark and generate plots")
    parser.add_argument("--pred", default="benchmark/results/triad_predictions.csv", help="Predictions CSV file")
    parser.add_argument("--outdir", default="benchmark/results", help="Output directory for plots")
    args = parser.parse_args()

    df = pd.read_csv(args.pred)
    results, tier_stats = compute_metrics(df)
    plot_curves(results, tier_stats, df, args.outdir)

    # Save summary json
    summary_path = Path(args.outdir) / "benchmark_metrics.json"
    with open(summary_path, "w", encoding="utf-8") as f:
        # Strip large curve arrays for clean json
        clean_res = {}
        for k, v in results.items():
            clean_res[k] = {
                "name": v["name"],
                "roc_auc": v["roc_auc"],
                "pr_auc": v["pr_auc"],
                "top10_precision": v["top10_precision"],
                "top20_precision": v["top20_precision"],
                "top50_precision": v["top50_precision"]
            }
        json.dump({"models": clean_res, "tier_enrichment": tier_stats}, f, indent=2)
    print(f"[✓] Benchmark metrics saved to: {summary_path}")

if __name__ == "__main__":
    main()
