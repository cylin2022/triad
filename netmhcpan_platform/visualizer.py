"""
visualizer.py
Interactive visualization for NetMHCpan-4.2 prediction results.
Generates robust, perfectly formatted Plotly charts with auto-scaled viewports.
"""

import io
import base64
from pathlib import Path
from typing import List, Optional, Dict, Any

import numpy as np
import pandas as pd
import plotly.graph_objects as go

# ─── Color Palette ─────────────────────────────────────────────────────────────

COLORS = {
    "strong":    "#00f5c4",   # cyan-green
    "weak":      "#f5a623",   # amber
    "nonbinder": "#6b7280",   # gray
    "bg":        "#0d1117",   # dark background
    "panel":     "#161b22",   # panel
    "panel2":    "#1c2330",   # panel 2
    "border":    "#30363d",   # border
    "text":      "#e6edf3",   # primary text
    "subtext":   "#8b949e",   # secondary text
    "primary":   "#58a6ff",   # blue accent
}


def _add_binder_category(df: pd.DataFrame) -> pd.DataFrame:
    """Add a 'BindCategory' column."""
    df = df.copy()
    if "Rank_EL" not in df.columns:
        df["BindCategory"] = "Non-Binder"
        return df

    conditions = [
        df["Rank_EL"] <= 0.5,
        (df["Rank_EL"] > 0.5) & (df["Rank_EL"] <= 2.0),
    ]
    choices = ["Strong Binder", "Weak Binder"]
    df["BindCategory"] = np.select(conditions, choices, default="Non-Binder")
    return df


def _fig_to_html(fig: go.Figure) -> str:
    """Convert a Plotly figure to an HTML div string."""
    return fig.to_html(full_html=False, include_plotlyjs=False)


# ─── Chart 1: %Rank Distribution (Binned Bar Chart) ───────────────────────────

def rank_distribution_chart(df: pd.DataFrame, title: str = "📊 %Rank EL Binding Affinity Tier Distribution") -> str:
    """Binned Bar chart of %Rank_EL values for 100% reliable rendering."""
    if df.empty or "Rank_EL" not in df.columns:
        return "<p style='color:#8b949e; text-align:center; padding:40px;'>No data for distribution chart.</p>"

    # 7 Standard Bins for %Rank EL
    bins = [0, 0.5, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0]
    labels = ["0–0.5% (SB)", "0.5–2% (WB)", "2–5%", "5–10%", "10–20%", "20–50%", "50–100%"]
    
    binned = pd.cut(df["Rank_EL"], bins=bins, labels=labels, include_lowest=True)
    counts = binned.value_counts().reindex(labels, fill_value=0)

    # Color for each bin matching its biological tier:
    # 0-0.5%: Strong Binder (green)
    # 0.5-2%: Weak Binder (orange)
    # >2%: Non-Binders (slate/gray)
    tier_colors = [
        COLORS["strong"],    # 0–0.5% (SB)
        COLORS["weak"],      # 0.5–2% (WB)
        "#4b5563",           # 2–5% (NB)
        "#6b7280",           # 5–10% (NB)
        "#9ca3af",           # 10–20% (NB)
        "#64748b",           # 20–50% (NB)
        "#475569",           # 50–100% (NB)
    ]

    fig = go.Figure()
    fig.add_trace(go.Bar(
        x=labels,
        y=[int(v) for v in counts.values],
        orientation="v",
        marker_color=tier_colors,
        text=[f"{int(v)}" if v > 0 else "" for v in counts.values],
        textposition="auto",
        textfont=dict(color="#ffffff", size=11, family="JetBrains Mono, monospace"),
        hoverinfo="x+y",
        hovertemplate="<b>%{x}</b><br>Peptides: %{y}<extra></extra>",
    ))

    fig.update_layout(
        paper_bgcolor=COLORS["bg"],
        plot_bgcolor=COLORS["panel"],
        font=dict(color=COLORS["text"], family="Inter, sans-serif"),
        title=dict(text=title, font=dict(size=15, color=COLORS["text"])),
        xaxis=dict(
            title="%Rank EL Interval (Binding Affinity Tiers)",
            type="category",
            categoryorder="array",
            categoryarray=labels,
            tickangle=-20,
            automargin=True,
            gridcolor=COLORS["border"],
        ),
        yaxis=dict(
            title="Peptide Count",
            gridcolor=COLORS["border"],
            zerolinecolor=COLORS["border"],
            rangemode="tozero",
        ),
        showlegend=False,
        height=380,
        margin=dict(t=50, b=75, l=60, r=30),
        autosize=True,
    )
    return _fig_to_html(fig)


# ─── Chart 2: Allele Summary Bar Chart (Horizontal Stacked Bar) ──────────────

def allele_summary_chart(df: pd.DataFrame) -> str:
    """
    Horizontal stacked bar chart: per-allele count of SB / WB / Non-Binder.
    Using complete matrix indexing guarantees all alleles are 100% aligned with exact SB counts.
    """
    if df.empty or "MHC" not in df.columns:
        return "<p style='color:#8b949e; text-align:center; padding:40px;'>No data for allele summary.</p>"

    df = _add_binder_category(df)
    alleles = sorted(df["MHC"].unique())

    # Complete pivot matrix: guarantees every allele has SB, WB, NB rows (with 0 if none)
    matrix = (
        df.groupby(["MHC", "BindCategory"])
        .size()
        .unstack(fill_value=0)
        .reindex(index=alleles, columns=["Strong Binder", "Weak Binder", "Non-Binder"], fill_value=0)
    )
    totals = matrix.sum(axis=1)
    max_total = int(totals.max()) if not totals.empty else 100

    cat_order = ["Strong Binder", "Weak Binder", "Non-Binder"]
    cat_colors_list = [COLORS["strong"], COLORS["weak"], COLORS["nonbinder"]]

    fig = go.Figure()
    for cat, color in zip(cat_order, cat_colors_list):
        vals = [int(v) for v in matrix[cat].values]
        # Only show explicit counts inside Strong Binder and Weak Binder sections
        show_text = [f"{v}" if v > 0 else "" for v in vals] if cat != "Non-Binder" else ["" for _ in vals]
        fig.add_trace(go.Bar(
            y=[str(a) for a in alleles],
            x=vals,
            name=cat,
            orientation="h",
            marker_color=color,
            text=show_text,
            textposition="auto",
            textfont=dict(color="#ffffff", size=11, family="JetBrains Mono, monospace"),
            hovertemplate="<b>%{y}</b><br>" + cat + ": <b>%{x}</b><extra></extra>",
        ))

    # Add unified Total count annotation at the end of each bar
    annotations = []
    for a in alleles:
        tot = int(totals[a])
        annotations.append(dict(
            x=tot,
            y=str(a),
            text=f" <b>Total: {tot}</b>",
            showarrow=False,
            xanchor="left",
            font=dict(color=COLORS["subtext"], size=10, family="Inter, sans-serif")
        ))

    fig.update_layout(
        paper_bgcolor=COLORS["bg"],
        plot_bgcolor=COLORS["panel"],
        font=dict(color=COLORS["text"], family="Inter, sans-serif"),
        title=dict(text="🧬 Candidate Binder Breakdown per HLA Allele", font=dict(size=15, color=COLORS["text"])),
        yaxis=dict(
            title="HLA Allele",
            type="category",
            categoryorder="array",
            categoryarray=[str(a) for a in alleles],
            autorange="reversed",
            automargin=True,
            gridcolor=COLORS["border"],
        ),
        xaxis=dict(
            title="Peptide Count",
            gridcolor=COLORS["border"],
            rangemode="tozero",
            range=[0, max_total * 1.2],
        ),
        barmode="stack",
        annotations=annotations,
        height=380,
        margin=dict(t=50, b=50, l=110, r=40),
        legend=dict(bgcolor=COLORS["panel"], bordercolor=COLORS["border"], borderwidth=1),
        autosize=True,
    )
    return _fig_to_html(fig)


# ─── Chart 3: Volcano Plot (Score EL vs -log10 %Rank) ──────────────────────────

def volcano_chart(df: pd.DataFrame) -> str:
    """
    Volcano scatter plot with auto-calculated Y and X ranges to guarantee
    all data points (SB, WB, Non-Binder) are 100% visible inside the viewport.
    """
    if df.empty or "Rank_EL" not in df.columns or "Score_EL" not in df.columns:
        return "<p style='color:#8b949e; text-align:center; padding:40px;'>No data for volcano plot.</p>"

    df = _add_binder_category(df.copy())
    df = df[df["Rank_EL"] > 0]
    # Transform %Rank (0.001% - 100%) to positive log scale: -log10(%Rank / 100)
    df["log_rank"] = -np.log10(df["Rank_EL"] / 100.0)

    # If dataset is large (>3000 rows), preserve 100% of binders and sample non-binders for ultra-fast WebGL rendering
    if len(df) > 3000:
        binders = df[df["BindCategory"].isin(["Strong Binder", "Weak Binder"])]
        non_binders = df[df["BindCategory"] == "Non-Binder"]
        if len(non_binders) > 1500:
            non_binders = non_binders.sample(n=1500, random_state=42)
        df = pd.concat([binders, non_binders], ignore_index=True)

    cat_colors = {
        "Strong Binder": COLORS["strong"],
        "Weak Binder":   COLORS["weak"],
        "Non-Binder":    COLORS["nonbinder"],
    }

    fig = go.Figure()
    for cat in ["Strong Binder", "Weak Binder", "Non-Binder"]:
        sub = df[df["BindCategory"] == cat]
        if sub.empty:
            continue
        # Use WebGL Accelerated Scattergl for GPU-based sub-10ms rendering
        fig.add_trace(go.Scattergl(
            x=sub["Score_EL"].tolist(),
            y=sub["log_rank"].tolist(),
            mode="markers",
            name=cat,
            marker=dict(
                color=cat_colors[cat],
                size=9,
                opacity=0.85,
            ),
            customdata=sub[["Peptide", "MHC", "Rank_EL"]].values.tolist(),
            hovertemplate=(
                "<b>Peptide:</b> %{customdata[0]}<br>"
                "<b>HLA:</b> %{customdata[1]}<br>"
                "<b>Score EL:</b> %{x:.4f}<br>"
                "<b>%Rank EL:</b> %{customdata[2]:.3f}%<extra></extra>"
            ),
        ))

    # Threshold horizontal lines for SB (0.5%) and WB (2.0%)
    sb_y = -np.log10(0.5 / 100.0)   # ~2.301
    wb_y = -np.log10(2.0 / 100.0)   # ~1.699

    fig.add_hline(
        y=sb_y, line_dash="dash", line_color=COLORS["strong"], line_width=1.5,
        annotation_text="Strong Binder (0.5%)", annotation_position="top right",
        annotation_font_color=COLORS["strong"]
    )
    fig.add_hline(
        y=wb_y, line_dash="dash", line_color=COLORS["weak"], line_width=1.5,
        annotation_text="Weak Binder (2.0%)", annotation_position="bottom right",
        annotation_font_color=COLORS["weak"]
    )

    # Explicitly set Y-axis range to encompass ALL data points + threshold lines
    data_min_y = float(df["log_rank"].min())
    data_max_y = float(df["log_rank"].max())
    min_y = max(0.0, min(data_min_y, wb_y) - 0.3)
    max_y = max(data_max_y, sb_y) + 0.4

    # Explicitly set X-axis range
    min_x = max(0.0, float(df["Score_EL"].min()) - 0.05)
    max_x = min(1.0, float(df["Score_EL"].max()) + 0.08)
    if max_x - min_x < 0.2:
        min_x = 0.0
        max_x = 1.0

    fig.update_layout(
        paper_bgcolor=COLORS["bg"],
        plot_bgcolor=COLORS["panel"],
        font=dict(color=COLORS["text"], family="Inter, sans-serif"),
        title=dict(text="🌋 Presentation Score vs. %Rank EL Affinity (Volcano Plot)", font=dict(size=15, color=COLORS["text"])),
        xaxis=dict(
            title="Score EL (Presentation Score)",
            range=[min_x, max_x],
            gridcolor=COLORS["border"],
            zerolinecolor=COLORS["border"]
        ),
        yaxis=dict(
            title="−log₁₀(%Rank / 100)",
            range=[min_y, max_y],
            gridcolor=COLORS["border"],
            zerolinecolor=COLORS["border"]
        ),
        height=420,
        margin=dict(t=50, b=60, l=70, r=30),
        legend=dict(bgcolor=COLORS["panel"], bordercolor=COLORS["border"], borderwidth=1),
        autosize=True,
    )
    return _fig_to_html(fig)


# ─── Chart 4: Top Peptides × Alleles Heatmap ───────────────────────────────────

def heatmap_chart(df: pd.DataFrame, top_n: int = 25) -> str:
    """
    Heatmap of top-N peptides × alleles (Score_EL).
    Explicitly uses categorical axes and text annotations in cells.
    """
    if df.empty or "Score_EL" not in df.columns or "MHC" not in df.columns:
        return "<p style='color:#8b949e; text-align:center; padding:40px;'>No data for heatmap.</p>"

    # Pivot: rows = peptides, cols = alleles
    pivot = df.pivot_table(
        index="Peptide", columns="MHC", values="Score_EL", aggfunc="max"
    ).fillna(0)

    # Select top_n peptides by max score across alleles
    pivot["_max"] = pivot.max(axis=1)
    pivot = pivot.nlargest(min(top_n, len(pivot)), "_max").drop(columns="_max")

    peptides = [str(p) for p in pivot.index]
    alleles  = [str(a) for a in pivot.columns]
    z_matrix = pivot.values.tolist()
    
    # Text overlay in cells
    text_matrix = [[f"{v:.3f}" if v > 0 else "0.000" for v in row] for row in z_matrix]

    fig = go.Figure(data=go.Heatmap(
        z=z_matrix,
        x=alleles,
        y=peptides,
        text=text_matrix,
        texttemplate="%{text}",
        textfont=dict(size=11, color="#ffffff"),
        colorscale=[
            [0.0,  "#161b22"],
            [0.2,  "#132a4a"],
            [0.5,  "#1d5b96"],
            [0.8,  "#009688"],
            [1.0,  "#00f5c4"],
        ],
        colorbar=dict(
            title=dict(text="Score EL", font=dict(color=COLORS["text"])),
            tickfont=dict(color=COLORS["subtext"]),
        ),
        hoverongaps=False,
    ))

    fig.update_layout(
        paper_bgcolor=COLORS["bg"],
        plot_bgcolor=COLORS["panel"],
        font=dict(color=COLORS["text"], family="Inter, sans-serif"),
        title=dict(
            text=f"🔥 Top-{len(peptides)} Candidate Peptides × HLA Alleles Heatmap (Score EL)",
            font=dict(size=15, color=COLORS["text"])
        ),
        xaxis=dict(title="HLA Allele", type="category", tickangle=-20, automargin=True, gridcolor=COLORS["border"]),
        yaxis=dict(title="Peptide Sequence", type="category", autorange="reversed", gridcolor=COLORS["border"]),
        height=max(360, min(len(peptides) * 28 + 120, 720)),
        margin=dict(t=50, b=70, l=130, r=30),
        autosize=True,
    )
    return _fig_to_html(fig)


# ─── Chart 5: Immunogenicity vs Binding Affinity Quad-Plot ────────────────────

def immunogenicity_quad_chart(df: pd.DataFrame) -> str:
    """
    Quad-Plot: T-Cell Immunogenicity Score vs %Rank EL.
    Identifies Golden Epitope Candidates (High Affinity + High T-Cell Response).
    """
    if df.empty or "Rank_EL" not in df.columns or "Immunogenicity_Score" not in df.columns:
        return "<p style='color:#8b949e; text-align:center; padding:40px;'>IEDB Immunogenicity Scores not available for this run.</p>"

    valid_df = df.dropna(subset=["Immunogenicity_Score"]).copy()
    if valid_df.empty:
        return "<p style='color:#8b949e; text-align:center; padding:40px;'>No T-cell immunogenicity data available.</p>"

    valid_df = _add_binder_category(valid_df)

    cat_colors = {
        "Strong Binder": COLORS["strong"],
        "Weak Binder":   COLORS["weak"],
        "Non-Binder":    COLORS["nonbinder"],
    }

    fig = go.Figure()
    for cat in ["Strong Binder", "Weak Binder", "Non-Binder"]:
        sub = valid_df[valid_df["BindCategory"] == cat]
        if sub.empty:
            continue
        fig.add_trace(go.Scattergl(
            x=sub["Rank_EL"].tolist(),
            y=sub["Immunogenicity_Score"].tolist(),
            mode="markers",
            name=cat,
            marker=dict(color=cat_colors[cat], size=10, opacity=0.85),
            customdata=sub[["Peptide", "MHC", "Score_EL"]].values.tolist(),
            hovertemplate=(
                "<b>Peptide:</b> %{customdata[0]}<br>"
                "<b>HLA:</b> %{customdata[1]}<br>"
                "<b>%Rank EL:</b> %{x:.3f}%<br>"
                "<b>Immunogenicity Score:</b> %{y:.4f}<extra></extra>"
            )
        ))

    # Add threshold lines
    fig.add_vline(x=2.0, line_dash="dash", line_color=COLORS["weak"], annotation_text="Binder Threshold (2.0%)", annotation_position="top right")
    fig.add_hline(y=0.0, line_dash="dash", line_color=COLORS["primary"], annotation_text="Positive Immunogenicity", annotation_position="bottom right")

    fig.update_layout(
        paper_bgcolor=COLORS["bg"],
        plot_bgcolor=COLORS["panel"],
        font=dict(color=COLORS["text"], family="Inter, sans-serif"),
        title=dict(text="🎯 T-Cell Immunogenicity vs. %Rank EL Binding Affinity (Quad-Plot)", font=dict(size=15, color=COLORS["text"])),
        xaxis=dict(title="%Rank EL (Lower = Stronger Binding)", gridcolor=COLORS["border"], zerolinecolor=COLORS["border"]),
        yaxis=dict(title="IEDB T-Cell Immunogenicity Score", gridcolor=COLORS["border"], zerolinecolor=COLORS["border"]),
        height=380,
        margin=dict(t=50, b=50, l=70, r=30),
        legend=dict(bgcolor=COLORS["panel"], bordercolor=COLORS["border"], borderwidth=1),
        autosize=True,
    )
    return _fig_to_html(fig)


# ─── Chart 6: Regional Population Coverage Bar Chart ──────────────────────────

def population_coverage_chart(job_meta: dict) -> str:
    """
    Regional HLA Population Coverage Bar Chart.
    """
    world_pct = job_meta.get("cov_world_pct", 0)
    asia_pct  = job_meta.get("cov_asia_pct", 0)

    regions = ["World Population", "East Asia Population"]
    values = [world_pct, asia_pct]

    fig = go.Figure(data=[go.Bar(
        x=regions,
        y=values,
        marker_color=[COLORS["strong"], COLORS["primary"]],
        text=[f"{v}%" for v in values],
        textposition="auto",
    )])

    fig.update_layout(
        paper_bgcolor=COLORS["bg"],
        plot_bgcolor=COLORS["panel"],
        font=dict(color=COLORS["text"], family="Inter, sans-serif"),
        title=dict(text="HLA Human Population Coverage (%)", font=dict(size=15, color=COLORS["text"])),
        xaxis=dict(title="Target Human Population", gridcolor=COLORS["border"]),
        yaxis=dict(title="Cumulative Coverage %", range=[0, 100], gridcolor=COLORS["border"]),
        height=380,
        margin=dict(t=50, b=50, l=60, r=30),
        autosize=True,
    )
    return _fig_to_html(fig)


# ─── All-in-one Dashboard ─────────────────────────────────────────────────────

def build_dashboard_html(df: pd.DataFrame, job_meta: dict) -> Dict[str, str]:
    """Combine all charts into a dictionary of HTML strings."""
    return {
        "rank_dist":  rank_distribution_chart(df),
        "allele_bar": allele_summary_chart(df),
        "volcano":    volcano_chart(df),
        "heatmap":    heatmap_chart(df),
    }



if __name__ == "__main__":
    test_df = pd.DataFrame([
        {'Pos':1, 'MHC':'HLA-A02:01', 'Peptide':'AAAWYLWEV', 'Core':'AAAWYLWEV', 'Score_EL':0.36, 'Rank_EL':0.68, 'BindLevel':'WB'},
        {'Pos':2, 'MHC':'HLA-B07:02', 'Peptide':'AEFGPWQTV', 'Core':'AEFGPWQTV', 'Score_EL':0.01, 'Rank_EL':5.14, 'BindLevel':''},
    ])
    charts = build_dashboard_html(test_df, {})
    print("Self-test succeeded. Chart keys:", list(charts.keys()))
