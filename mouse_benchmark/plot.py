"""Standalone retrospective benchmark figures, without inferred probabilities."""
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.metrics import precision_recall_curve

ROOT = Path(__file__).resolve().parent / "results"
NAMES = {"Calis": "Current Calis", "EL_Rank_Score": "netMHCpan EL rank",
         "BA_Rank_Score": "netMHCpan BA rank", "Mouse_Base": "Mouse base (no Calis)",
         "Mouse_Composite": "Current mouse composite"}


def plot():
    df = pd.read_csv(ROOT / "scored_pairs.csv")
    metrics = pd.read_csv(ROOT / "metrics.csv")
    sm = pd.read_csv(ROOT / "within_study_metrics.csv")
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), constrained_layout=True)
    for model, name in NAMES.items():
        p, r, _ = precision_recall_curve(df.Label, df[model])
        ap = metrics.query("Cohort == 'primary' and Model == @model").AP.iloc[0]
        axes[0].step(r, p, where="post", lw=1.5, label=f"{name} (AP {ap:.3f})")
    axes[0].axhline(df.Label.mean(), color="gray", linestyle="--", label=f"Prevalence {df.Label.mean():.3f}")
    axes[0].set(xlabel="Recall", ylabel="Precision", title="Pooled IFN-gamma response benchmark", xlim=(0, 1), ylim=(0, 1))
    axes[0].legend(fontsize=8)
    wide = sm.pivot(index="Reference_ID", columns="Model", values="AP")
    delta = (wide.Mouse_Composite - wide.Mouse_Base).sort_values()
    axes[1].barh(range(len(delta)), delta, color=["#b54b4b" if v < 0 else "#24858b" for v in delta])
    axes[1].set_yticks(range(len(delta)), [str(v) for v in delta.index], fontsize=7)
    axes[1].axvline(0, color="black", lw=.8)
    axes[1].set(xlabel="AP difference: composite minus base", ylabel="IEDB reference ID",
                title="Within-study effect of adding Calis")
    fig.suptitle("Mouse H-2-Kb / H-2-Db: retrospective, heterogeneous public data", fontsize=12)
    for ext in ["png", "svg"]:
        fig.savefig(ROOT / f"baseline_comparison.{ext}", dpi=180)
    plt.close(fig)


if __name__ == "__main__":
    plot()
