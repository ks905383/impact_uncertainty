"""
Per-pair residual-temperature disagreement vs country covariates.

For each of the 6 reanalysis pairs and each country, compute the mean
absolute difference of the residualised temperature time series (residuals
from `temp ~ 1 | iso + year + iso[time1] + iso[time2]`, matching the BHM
fixed-effect structure and Tom's earlier R code in
`tb/code/gdp_disagree.R`). Scatter each pair's per-country MAD against
(a) country-mean log GDP per capita and (b) country-mean temperature.

This is the per-PAIR generalisation of the deck-02 figure
`tc_disagr_vs_income.pdf` (which averaged the 10 pairs of the 5-reanalysis
set into one panel). Here we keep the 4-reanalysis set (6 pairs) for
consistency with the new dose-response analysis in deck 06.

Output: tb/figs/tc_disagr_perpair.pdf (2 rows × 6 cols)
  top row    = MAD vs mean T
  bottom row = MAD vs log mean GDP per capita
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats as sps

# Re-use the residualiser already in country_residuals.py — same FE structure
# as Tom's R code (`temp ~ 1 | iso + year + iso[time1] + iso[time2]`).
sys.path.insert(0, str(Path(__file__).resolve().parent))
from country_residuals import residualise_all   # noqa: E402


DROPBOX = Path(os.environ.get(
    "KS_DROPBOX", "/Users/tombearpark/Library/CloudStorage/Dropbox/ks"))
REPO = Path(os.environ.get(
    "KS_REPO", "/Users/tombearpark/Documents/GitHub/impact_uncertainty"))

REG_INPUT = DROPBOX / "reg_inputs" / "gdppc" / "reginputs_gdppc_025deg.csv"
OUT_PDF   = REPO / "tb" / "figs" / "tc_disagr_perpair.pdf"

REANALYSES = ("ERA5-025", "MERRA2", "JRA-3Q", "GMFD")
LABEL = {"ERA5-025": "ERA5", "MERRA2": "MERRA2",
         "JRA-3Q": "JRA-3Q", "GMFD": "GMFD"}
YEAR_LO, YEAR_HI = 1980, 2010

# Same canonical pair ordering as loo_perpair_maps / loo_scatter
PAIR_ORDER = [
    ("MERRA2", "GMFD"),
    ("JRA-3Q", "ERA5-025"),
    ("MERRA2", "ERA5-025"),
    ("GMFD",   "ERA5-025"),
    ("MERRA2", "JRA-3Q"),
    ("JRA-3Q", "GMFD"),
]


def load_panel() -> pd.DataFrame:
    """Read reg_input + force balanced (iso, year) panel across the 4 reanalyses,
    matching the deck-02 R code's balancing step."""
    use = ["countryname", "iso", "year", "obs_mod", "temp_popweight",
           "time", "time2", "gdpCAPppp_wdi"]
    df = pd.read_csv(REG_INPUT, usecols=use)
    df = df[(df["year"] >= YEAR_LO) & (df["year"] <= YEAR_HI)]
    df = df[df["obs_mod"].isin(REANALYSES)]
    df = df.dropna(subset=["temp_popweight"])

    # Balanced (iso, year) — every reanalysis must observe the cell
    n_rea = len(REANALYSES)
    counts = df.groupby(["iso", "year"])["obs_mod"].nunique().reset_index()
    keep = counts[counts["obs_mod"] == n_rea][["iso", "year"]]
    df = df.merge(keep, on=["iso", "year"], how="inner").reset_index(drop=True)
    return df


def country_covariates(df: pd.DataFrame) -> pd.DataFrame:
    """Per-country mean T and mean log GDP/cap. Aggregation matches the
    deck-02 R helper (`make_country_cov`)."""
    # Take one reanalysis to avoid double-counting (T is reanalysis-specific
    # but for ranking we want one representative; deck-02 used mods[1]).
    sub = df[df["obs_mod"] == REANALYSES[0]]
    g = sub.groupby(["iso", "countryname"]).agg(
        mean_T   = ("temp_popweight", "mean"),
        mean_gdp = ("gdpCAPppp_wdi",  "mean"))
    g = g.reset_index()
    g["log_gdp"] = np.log(g["mean_gdp"])
    return g


def per_pair_country_mad(resid: pd.DataFrame) -> pd.DataFrame:
    """For each pair (a, b) and each country, mean over years of
    |resid_a - resid_b|. Returns long DataFrame (pair_a, pair_b, iso,
    countryname, mad)."""
    wide = resid.pivot_table(index=["countryname", "year"],
                             columns="obs_mod", values="resid").reset_index()
    rows = []
    for (a, b) in PAIR_ORDER:
        d = wide.dropna(subset=[a, b]).copy()
        d["d"] = (d[a] - d[b]).abs()
        m = d.groupby("countryname", as_index=False).agg(mad=("d", "mean"))
        m["pair_a"] = a
        m["pair_b"] = b
        rows.append(m)
    return pd.concat(rows, ignore_index=True)


# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------

def _ols_overlay(ax, x, y, color="0.25"):
    x = np.asarray(x, float); y = np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    x, y = x[m], y[m]
    if x.size < 3:
        return
    r = sps.linregress(x, y)
    xs = np.linspace(x.min(), x.max(), 50)
    ax.plot(xs, r.slope * xs + r.intercept,
            color=color, linewidth=1.0, linestyle="--", zorder=2)
    ax.text(0.97, 0.96,
            f"$R^2$={r.rvalue**2:.2f}\np={r.pvalue:.2g}",
            transform=ax.transAxes, ha="right", va="top",
            fontsize=6.5, color=color,
            bbox=dict(boxstyle="round,pad=0.18",
                      facecolor="white", edgecolor="none", alpha=0.75))


def build_figure(pair_mad: pd.DataFrame,
                 cov: pd.DataFrame) -> plt.Figure:
    """2 rows × 6 cols. Row 0 = MAD vs mean T; row 1 = MAD vs log GDP."""
    merged = pair_mad.merge(cov, on="countryname", how="inner")

    fig, axes = plt.subplots(2, 6, figsize=(15.5, 5.6),
                             sharex="row", sharey="row")

    for col, (a, b) in enumerate(PAIR_ORDER):
        sub = merged[(merged["pair_a"] == a) & (merged["pair_b"] == b)]

        # --- top row: vs mean T -----------------------------------------
        ax = axes[0, col]
        ax.scatter(sub["mean_T"], sub["mad"],
                   s=10, color="0.45", alpha=0.7,
                   edgecolor="white", linewidth=0.3)
        _ols_overlay(ax, sub["mean_T"].to_numpy(), sub["mad"].to_numpy())
        ax.set_title(f"{LABEL[a]} vs {LABEL[b]}", fontsize=9, loc="left")
        ax.tick_params(labelsize=7)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        if col == 0:
            ax.set_ylabel("Mean |Δresidual T|  (°C)", fontsize=8)

        # --- bottom row: vs log GDP -------------------------------------
        ax = axes[1, col]
        sub_g = sub.dropna(subset=["log_gdp"])
        ax.scatter(sub_g["log_gdp"], sub_g["mad"],
                   s=10, color="0.45", alpha=0.7,
                   edgecolor="white", linewidth=0.3)
        _ols_overlay(ax, sub_g["log_gdp"].to_numpy(),
                     sub_g["mad"].to_numpy())
        ax.tick_params(labelsize=7)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        if col == 0:
            ax.set_ylabel("Mean |Δresidual T|  (°C)", fontsize=8)

    for ax in axes[0, :]:
        ax.set_xlabel("country-mean T (°C)", fontsize=7.5, labelpad=2)
    for ax in axes[1, :]:
        ax.set_xlabel("log mean GDP/cap (PPP)", fontsize=7.5, labelpad=2)

    fig.suptitle("Per-pair residualised temperature disagreement vs.\\ "
                 "country-mean climate and income",
                 fontsize=12, y=0.99)
    fig.text(0.5, 0.945,
             "Each point = one country.  Residuals from "
             "$T \\sim 1 \\,|\\, \\mathrm{iso} + \\mathrm{year} + "
             "\\mathrm{iso}[t] + \\mathrm{iso}[t^2]$  "
             "(same FE structure as the BHM regression).  "
             "Dashed = OLS fit.",
             fontsize=8, color="0.35", style="italic", ha="center")
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    return fig


def main():
    print("Loading + balancing reg-input panel …")
    df = load_panel()
    n_cy = df.groupby(["iso", "year"]).ngroups
    print(f"  balanced rows: {len(df):,}  ({n_cy:,} country-years × "
          f"{len(REANALYSES)} reanalyses)")

    print("Residualising temperature per reanalysis "
          "(iso FE + year FE + iso×time + iso×time²) …")
    resid = residualise_all(df.rename(columns={"countryname": "countryname"}))
    print(f"  residualised rows: {len(resid):,}")

    print("Computing per-pair country MADs …")
    pair_mad = per_pair_country_mad(resid)
    print(f"  rows: {len(pair_mad):,} "
          f"({pair_mad['countryname'].nunique()} countries × "
          f"{pair_mad.groupby(['pair_a','pair_b']).ngroups} pairs)")
    print(pair_mad.groupby(["pair_a", "pair_b"])["mad"]
          .describe()[["mean", "std", "min", "max"]]
          .round(3).to_string())

    print("Building country covariates …")
    cov = country_covariates(df)
    print(f"  {len(cov)} countries  ({cov['log_gdp'].notna().sum()} with GDP)")

    print("Building figure …")
    fig = build_figure(pair_mad, cov)
    OUT_PDF.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PDF, format="pdf", bbox_inches="tight")
    fig.savefig(OUT_PDF.with_suffix(".png"), format="png", dpi=180,
                bbox_inches="tight")
    print(f"Saved: {OUT_PDF}")


if __name__ == "__main__":
    main()
