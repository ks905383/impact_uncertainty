"""
Scatter plots: country-level pct_reduction vs (a) mean temperature, (b) mean
log GDP per capita.

Reveals where in the climate / income dimensions the disagreement-driving
countries sit. Replaces deck-02's "disagreement vs income / climate" frame
for the new all-pairs averaged metric.

Inputs:
  - dmes/mean_dmes_gdppc.csv         (Kevin's per-country pct_reduction base,
                                      averaged across pairs)
  - dmes/dmes_gdppc.csv              (Kevin's per-PAIR pct_reduction)
  - reg_inputs/gdppc/reginputs_gdppc_025deg.csv  (T + GDP per capita per country-year)

Outputs:
  - tb/figs/loo_scatter.pdf                — averaged metric, 2 panels (T, income)
  - tb/figs/loo_scatter_perpair_T.pdf      — per-pair, 6 panels vs mean T
  - tb/figs/loo_scatter_perpair_income.pdf — per-pair, 6 panels vs mean log income
"""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats as sps


def _ols_overlay(ax, x, y, *, color="0.25") -> None:
    """Fit OLS y ~ x, draw the fitted line over the data x-range, and
    annotate the panel with R^2 + slope + p-value. Robust to <3 points
    (skips). x and y are 1-D arrays (any duplicates / nans dropped)."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    m = np.isfinite(x) & np.isfinite(y)
    x, y = x[m], y[m]
    if x.size < 3:
        return
    res = sps.linregress(x, y)
    xs = np.linspace(x.min(), x.max(), 50)
    ax.plot(xs, res.slope * xs + res.intercept,
            color=color, linewidth=1.0, linestyle="--", zorder=2)
    # Annotation in top-right; small + muted so it doesn't dominate
    ax.text(0.98, 0.97,
            f"$R^2$ = {res.rvalue**2:.2f}\n"
            f"slope = {res.slope:+.2f}\n"
            f"p = {res.pvalue:.2g}",
            transform=ax.transAxes, ha="right", va="top",
            fontsize=7, color=color,
            bbox=dict(boxstyle="round,pad=0.2",
                      facecolor="white", edgecolor="none", alpha=0.75))


DROPBOX = Path(os.environ.get(
    "KS_DROPBOX", "/Users/tombearpark/Library/CloudStorage/Dropbox/ks"))
REPO = Path(os.environ.get(
    "KS_REPO", "/Users/tombearpark/Documents/GitHub/impact_uncertainty"))

DMES_CSV       = DROPBOX / "dmes" / "mean_dmes_gdppc.csv"
DMES_PAIR_CSV  = DROPBOX / "dmes" / "dmes_gdppc.csv"
REG_INPUT      = DROPBOX / "reg_inputs" / "gdppc" / "reginputs_gdppc_025deg.csv"
OUT_PDF        = REPO / "tb" / "figs" / "loo_scatter.pdf"

REANALYSES = ("ERA5-025", "MERRA2", "JRA-3Q", "GMFD")
YEAR_LO, YEAR_HI = 1980, 2010

# Same canonical pair ordering as loo_perpair_maps.py — smaller-baseline
# pairs first so the eye sweeps the disagreement spectrum L→R, T→B.
PAIR_ORDER = [
    ("MERRA2", "GMFD"),
    ("JRA-3Q", "ERA5"),
    ("MERRA2", "ERA5"),
    ("GMFD",   "ERA5"),
    ("MERRA2", "JRA-3Q"),
    ("JRA-3Q", "GMFD"),
]


def country_means() -> pd.DataFrame:
    """Per country: mean temperature (averaged across years and reanalyses),
    mean GDP per capita PPP (averaged across years; constant across reanalyses)."""
    use = ["countryname", "year", "obs_mod", "temp_popweight", "gdpCAPppp_wdi"]
    df = pd.read_csv(REG_INPUT, usecols=use)
    df = df[(df["year"] >= YEAR_LO) & (df["year"] <= YEAR_HI)]
    df = df[df["obs_mod"].isin(REANALYSES)]
    df = df.dropna(subset=["temp_popweight"])
    g = df.groupby("countryname").agg(
        mean_T   = ("temp_popweight", "mean"),
        mean_gdp = ("gdpCAPppp_wdi", "mean"))
    return g.reset_index().rename(columns={"countryname": "country"})


def load_perpair_pct(source: str = "kevin") -> pd.DataFrame:
    """Per (pair, country): pct_reduction = (L1_full − L1_drop) / L1_full * 100.

    source: "kevin"   → Kevin's pop-weighted dmes_gdppc.csv
            "uniform" → recompute via funcs_dme (each obs equal weight)
    """
    if source == "uniform":
        import sys
        from pathlib import Path as _P
        sys.path.insert(0, str(_P(__file__).resolve().parent))
        import funcs_dme  # noqa: E402
        raw = funcs_dme.compute_pair_dme_table(weighting="uniform")
    else:
        raw = pd.read_csv(DMES_PAIR_CSV)
    out_rows = []
    for (a, b), grp in raw.groupby(["proj_base_a", "proj_base_b"]):
        base = float(grp.loc[grp["loo_id"] == "all", "dme"].iloc[0])
        per = grp[grp["loo_id"] != "all"].copy()
        per["pct_reduction"] = (base - per["dme"]) / base * 100.0
        per["pair_a"] = a
        per["pair_b"] = b
        out_rows.append(per.rename(columns={"loo_id": "country"}))
    return pd.concat(out_rows, ignore_index=True)


def load_mean_pct(source: str = "kevin") -> pd.DataFrame:
    """Per-country averaged-across-pairs pct_reduction.

    source: "kevin"   → Kevin's pop-weighted mean_dmes_gdppc.csv
            "uniform" → recompute via funcs_dme (each obs equal weight)
    Returns columns: country, dme_mu, pct_reduction.
    """
    if source == "uniform":
        import sys
        from pathlib import Path as _P
        sys.path.insert(0, str(_P(__file__).resolve().parent))
        import funcs_dme  # noqa: E402
        raw = funcs_dme.compute_mean_dme_table(weighting="uniform")
    else:
        raw = pd.read_csv(DMES_CSV)
        base = float(raw.loc[raw["loo_id"] == "all", "dme_mu"].iloc[0])
        raw["pct_reduction"] = 100.0 * (base - raw["dme_mu"]) / base
    raw = raw.rename(columns={"loo_id": "country"})
    return raw[raw["country"] != "all"].reset_index(drop=True)


def fig_perpair_scatter(perpair: pd.DataFrame, means: pd.DataFrame,
                        x_col: str, x_label: str,
                        title: str) -> plt.Figure:
    """2×3 grid of per-pair scatters of pct_reduction vs `x_col` (a country
    property). Each panel uses its own diverging colour norm centred at 0
    so the within-pair spread is readable. Top-3 / bottom-3 per pair are
    labelled."""
    from matplotlib.colors import TwoSlopeNorm
    cmap = "PRGn_r"

    fig, axes = plt.subplots(2, 3, figsize=(14.5, 8.2),
                             sharex=True, sharey=False)
    for ax, (a, b) in zip(axes.flat, PAIR_ORDER):
        sub = perpair[(perpair["pair_a"] == a) & (perpair["pair_b"] == b)]
        sub = sub.merge(means, on="country", how="inner").dropna(subset=[x_col])
        if sub.empty:
            ax.set_visible(False)
            continue
        vmin = float(sub["pct_reduction"].min())
        vmax = float(sub["pct_reduction"].max())
        norm = TwoSlopeNorm(vmin=min(vmin, -0.1), vcenter=0,
                            vmax=max(vmax, 0.1))
        ax.scatter(sub[x_col], sub["pct_reduction"],
                   c=sub["pct_reduction"], cmap=cmap, norm=norm,
                   s=16, edgecolor="white", linewidth=0.3)
        ax.axhline(0, color="0.4", linewidth=0.5)
        _ols_overlay(ax, sub[x_col].to_numpy(), sub["pct_reduction"].to_numpy())
        ax.set_title(f"{a} vs {b}", fontsize=10, loc="left")
        ax.tick_params(labelsize=8)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        # Label per-pair extremes
        labels = pd.concat([sub.nlargest(3, "pct_reduction"),
                            sub.nsmallest(3, "pct_reduction")])
        for _, row in labels.iterrows():
            ax.annotate(row["country"], (row[x_col], row["pct_reduction"]),
                        fontsize=6.5, ha="left", va="bottom",
                        xytext=(2, 1), textcoords="offset points",
                        color="0.25")

    # Shared x-axis label / y-axis label
    for ax in axes[-1, :]:
        ax.set_xlabel(x_label, fontsize=9)
    for ax in axes[:, 0]:
        ax.set_ylabel("% reduction in pair's L1\nwhen country dropped",
                      fontsize=9)
    fig.suptitle(title, fontsize=12, y=0.995)
    fig.tight_layout()
    return fig


def _build_one(weighting: str, suffix: str, tag: str, means: pd.DataFrame) -> None:
    """Build the 3 scatter outputs for one weighting choice."""
    print(f"\n=========== {tag} ===========")
    print("Loading per-country pct_reduction …")
    source = "kevin" if weighting == "pop" else "uniform"
    pct = load_mean_pct(source=source)
    print(f"  {len(pct)} countries; pct range "
          f"[{pct['pct_reduction'].min():.2f}, {pct['pct_reduction'].max():.2f}]")

    merged = pct.merge(means, on="country", how="inner")
    has_gdp = merged.dropna(subset=["mean_gdp"]).copy()
    print(f"  {len(merged)} with T; {len(has_gdp)} with GDP")

    # Figure: two-panel scatter
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.4))
    cmap = "PRGn_r"
    vmin = float(merged["pct_reduction"].min())
    vmax = float(merged["pct_reduction"].max())
    from matplotlib.colors import TwoSlopeNorm
    norm = TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)

    # Panel (a): pct vs mean T
    ax = axes[0]
    sc = ax.scatter(merged["mean_T"], merged["pct_reduction"],
                    c=merged["pct_reduction"], cmap=cmap, norm=norm,
                    s=22, edgecolor="white", linewidth=0.4)
    ax.axhline(0, color="0.4", linewidth=0.5)
    _ols_overlay(ax, merged["mean_T"].to_numpy(),
                 merged["pct_reduction"].to_numpy())
    ax.set_xlabel("Country-mean annual temperature (°C, 1980–2010)", fontsize=10)
    ax.set_ylabel("% reduction in cross-reanalysis disagreement when country dropped",
                  fontsize=10)
    ax.set_title("(a) vs. country-mean temperature", fontsize=11, loc="left")
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    # Label outliers (top/bottom 5)
    for _, row in pd.concat([merged.nlargest(5, "pct_reduction"),
                              merged.nsmallest(5, "pct_reduction")]).iterrows():
        ax.annotate(row["country"], (row["mean_T"], row["pct_reduction"]),
                    fontsize=7, ha="left", va="bottom",
                    xytext=(3, 2), textcoords="offset points",
                    color="0.25")

    # Panel (b): pct vs mean log GDP
    ax = axes[1]
    has_gdp["log_gdp"] = np.log(has_gdp["mean_gdp"])
    ax.scatter(has_gdp["log_gdp"], has_gdp["pct_reduction"],
               c=has_gdp["pct_reduction"], cmap=cmap, norm=norm,
               s=22, edgecolor="white", linewidth=0.4)
    ax.axhline(0, color="0.4", linewidth=0.5)
    _ols_overlay(ax, has_gdp["log_gdp"].to_numpy(),
                 has_gdp["pct_reduction"].to_numpy())
    ax.set_xlabel("log mean GDP per capita (PPP, 1980–2010)", fontsize=10)
    ax.set_ylabel("% reduction in cross-reanalysis disagreement when country dropped",
                  fontsize=10)
    ax.set_title("(b) vs. country-mean log income", fontsize=11, loc="left")
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    for _, row in pd.concat([has_gdp.nlargest(5, "pct_reduction"),
                              has_gdp.nsmallest(5, "pct_reduction")]).iterrows():
        ax.annotate(row["country"], (row["log_gdp"], row["pct_reduction"]),
                    fontsize=7, ha="left", va="bottom",
                    xytext=(3, 2), textcoords="offset points",
                    color="0.25")

    # Shared colourbar at the bottom
    cax = fig.add_axes([0.30, -0.02, 0.40, 0.022])
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm); sm.set_array([])
    cb = fig.colorbar(sm, cax=cax, orientation="horizontal")
    cb.set_label("% reduction  (purple = drives disagreement; green = drives agreement)",
                 fontsize=8)
    cb.ax.tick_params(labelsize=8)

    fig.suptitle(
        f"Where in climate/income space do the disagreement-driving countries sit?  [{tag}]",
        fontsize=12, y=1.00)
    fig.tight_layout()
    OUT_PDF.parent.mkdir(parents=True, exist_ok=True)
    out_main = OUT_PDF.with_name(OUT_PDF.stem + suffix + ".pdf")
    fig.savefig(out_main, format="pdf", bbox_inches="tight")
    fig.savefig(out_main.with_suffix(".png"), format="png", dpi=180,
                bbox_inches="tight")
    print(f"Saved: {out_main}")

    # ----- Per-pair versions -------------------------------------------------
    print("Loading per-pair pct_reduction …")
    perpair = load_perpair_pct(source=source)
    print(f"  {len(perpair)} (pair, country) rows across "
          f"{perpair.groupby(['pair_a','pair_b']).ngroups} pairs")

    # Merge has_gdp once (adds log_gdp column for income panel)
    means_for_income = means.copy()
    means_for_income["log_gdp"] = np.log(means_for_income["mean_gdp"])

    print("Building per-pair scatter vs mean T …")
    fig_t = fig_perpair_scatter(
        perpair, means, x_col="mean_T",
        x_label="Country-mean annual temperature (°C, 1980–2010)",
        title=f"Per-pair country influence vs. country-mean temperature  [{tag}]")
    out_t = OUT_PDF.with_name(f"loo_scatter_perpair_T{suffix}.pdf")
    fig_t.savefig(out_t, format="pdf", bbox_inches="tight")
    fig_t.savefig(out_t.with_suffix(".png"), format="png", dpi=180,
                  bbox_inches="tight")
    print(f"  saved {out_t}")

    print("Building per-pair scatter vs log mean income …")
    fig_i = fig_perpair_scatter(
        perpair, means_for_income, x_col="log_gdp",
        x_label="log mean GDP per capita (PPP, 1980–2010)",
        title=f"Per-pair country influence vs. country-mean log income  [{tag}]")
    out_i = OUT_PDF.with_name(f"loo_scatter_perpair_income{suffix}.pdf")
    fig_i.savefig(out_i, format="pdf", bbox_inches="tight")
    fig_i.savefig(out_i.with_suffix(".png"), format="png", dpi=180,
                  bbox_inches="tight")
    print(f"  saved {out_i}")


def main():
    print("Loading country-level T + GDP means …")
    means = country_means()
    for weighting, suffix, tag in [
            ("pop",     "",         "pop-weighted (Kevin's CSV)"),
            ("uniform", "_uniform", "uniform (each obs equal weight)")]:
        _build_one(weighting, suffix, tag, means)


if __name__ == "__main__":
    main()
