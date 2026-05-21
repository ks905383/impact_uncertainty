"""
Bahrain walkthrough: visual explanation of the L1 distance metric.

Produces three figures for the slide deck:

  1. tb/figs/bahrain_me_curves.pdf
       2x2 grid (one panel per reanalysis): ME(T) curves, full sample vs
       LOO-Bahrain, so the reader can see how dropping a single country
       moves each reanalysis's marginal-effect curve.

  2. tb/figs/bahrain_l1_walkthrough_full.pdf
       3-panel step-by-step computation of the L1 distance for ONE chosen
       reanalysis pair (default: GMFD vs ERA5 — Bahrain is the top driver
       for this pair) on the full sample. Shows MEs overlaid, |diff|, and
       the integrand |diff| × pop-weighted T density.

  3. tb/figs/bahrain_l1_walkthrough_noBahrain.pdf
       Same 3-panel layout, recomputed after Bahrain is dropped. The
       reader sees the shaded area shrink — that's where the % reduction
       comes from.

Inputs:
  - Full-sample coeffs: reg_outputs/gdppc/regoutputs_coeffs_gdppc_025deg_{rea}_1980-2010.csv
  - LOO coeffs:         reg_outputs/gdppc/loos/regoutputs_coeffs_gdppc_025deg_{rea}_1980-2010_noBahrain.csv
  - Temperature panel + pop weights: reg_inputs/gdppc/reginputs_gdppc_025deg.csv
"""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec


DROPBOX = Path(os.environ.get(
    "KS_DROPBOX", "/Users/tombearpark/Library/CloudStorage/Dropbox/ks"))
REPO = Path(os.environ.get(
    "KS_REPO", "/Users/tombearpark/Documents/GitHub/impact_uncertainty"))

FULL_DIR  = DROPBOX / "reg_outputs" / "gdppc"
LOO_DIR   = DROPBOX / "reg_outputs" / "gdppc" / "loos"
REG_INPUT = DROPBOX / "reg_inputs"  / "gdppc" / "reginputs_gdppc_025deg.csv"

FIGS = REPO / "tb" / "figs"

# Note: reg_outputs file naming uses "ERA5-025"; we display as "ERA5".
REANALYSES = ("ERA5-025", "MERRA2", "JRA-3Q", "GMFD")
LABEL = {"ERA5-025": "ERA5", "MERRA2": "MERRA2",
         "JRA-3Q": "JRA-3Q", "GMFD": "GMFD"}
COLOUR = {"ERA5-025": "tab:red", "MERRA2": "tab:green",
          "JRA-3Q": "tab:orange", "GMFD": "tab:blue"}

YEAR_LO, YEAR_HI = 1980, 2010

# Which pair to walk through. Default: GMFD vs ERA5, where Bahrain ranks #2
# on the per-pair list of drivers (visible in tb/figs/loo_perpair_maps.png).
WALKTHROUGH_PAIR = ("GMFD", "ERA5-025")
WALKTHROUGH_COUNTRY = "Bahrain"


# ---------------------------------------------------------------------------
# Loaders
# ---------------------------------------------------------------------------

def load_full_coeffs() -> dict[str, tuple[float, float]]:
    out = {}
    for rea in REANALYSES:
        path = FULL_DIR / f"regoutputs_coeffs_gdppc_025deg_{rea}_1980-2010.csv"
        df = pd.read_csv(path)
        out[rea] = (float(df["b1"].iloc[0]), float(df["b2"].iloc[0]))
    return out


def load_loo_coeffs(country: str) -> dict[str, tuple[float, float]]:
    out = {}
    for rea in REANALYSES:
        path = LOO_DIR / f"regoutputs_coeffs_gdppc_025deg_{rea}_1980-2010_no{country}.csv"
        df = pd.read_csv(path)
        out[rea] = (float(df["b1"].iloc[0]), float(df["b2"].iloc[0]))
    return out


def load_reginputs() -> pd.DataFrame:
    use = ["countryname", "year", "obs_mod", "temp_popweight", "growthWDI", "Pop"]
    df = pd.read_csv(REG_INPUT, usecols=use)
    df = df[(df["year"] >= YEAR_LO) & (df["year"] <= YEAR_HI)]
    df = df[df["obs_mod"].isin(REANALYSES)]
    df = df.dropna(subset=["temp_popweight", "growthWDI"])
    return df.reset_index(drop=True)


# ---------------------------------------------------------------------------
# Maths
# ---------------------------------------------------------------------------

def me(b1: float, b2: float, t: np.ndarray) -> np.ndarray:
    """BHM marginal effect: dy/dT = b1 + 2 b2 T."""
    return b1 + 2 * b2 * t


def pair_l1(coeffs: dict[str, tuple[float, float]], pair: tuple[str, str],
            temps: np.ndarray, weights: np.ndarray | None) -> float:
    """L1 distance between ME curves of `pair[0]` and `pair[1]`, integrated
    over the empirical T distribution `temps` (with optional `weights`)."""
    a, b = pair
    ma = me(*coeffs[a], temps)
    mb = me(*coeffs[b], temps)
    diff = np.abs(ma - mb)
    if weights is None:
        return float(np.mean(diff))
    return float(np.sum(weights * diff) / np.sum(weights))


# ---------------------------------------------------------------------------
# Figure 1: ME curves with vs without Bahrain
# ---------------------------------------------------------------------------

def fig_me_curves(full: dict, loo: dict, t_grid: np.ndarray) -> plt.Figure:
    """Two side-by-side panels: full sample (left) vs sample with the country
    dropped (right). Each panel overlays the four reanalyses; the eye
    compares the spread between panels rather than solid-vs-dashed within
    one. Shared y-axis so the vertical drift is directly readable.
    """
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 5.4),
                             sharex=True, sharey=True)
    panel_titles = ("Full sample",
                    f"Sample with {WALKTHROUGH_COUNTRY} dropped")
    for ax, panel, panel_title in zip(axes, (full, loo), panel_titles):
        for rea in REANALYSES:
            ax.plot(t_grid, me(*panel[rea], t_grid),
                    color=COLOUR[rea], linewidth=1.8, label=LABEL[rea])
        ax.axhline(0, color="0.6", linewidth=0.5)
        ax.set_title(panel_title, fontsize=11, loc="left")
        ax.set_xlabel("Annual mean temperature (°C)", fontsize=10)
        ax.tick_params(labelsize=9)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
    axes[0].set_ylabel("Marginal effect of T on GDP-pc growth  (dy/dT)",
                       fontsize=10)
    axes[0].legend(title="reanalysis", fontsize=9, title_fontsize=9,
                   frameon=False, loc="upper right")

    fig.suptitle(
        f"BHM marginal-effect curves: full sample vs.  sample with "
        f"{WALKTHROUGH_COUNTRY} dropped",
        fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    return fig


# ---------------------------------------------------------------------------
# Figure 2 / 3: L1 walkthrough (3 panels)
# ---------------------------------------------------------------------------

def fig_walkthrough(coeffs: dict, pair: tuple[str, str],
                    temps: np.ndarray, weights: np.ndarray,
                    t_grid: np.ndarray, title: str) -> tuple[plt.Figure, float]:
    """3-panel walk-through:  (top) ME_a and ME_b vs T;
    (middle) |ME_a - ME_b|;  (bottom) histogram of T weighted by Pop, and the
    integrand |ME diff| × density (shaded). Returns (fig, L1)."""
    a, b = pair
    ma = me(*coeffs[a], t_grid)
    mb = me(*coeffs[b], t_grid)
    diff = np.abs(ma - mb)
    L1 = pair_l1(coeffs, pair, temps, weights)

    fig = plt.figure(figsize=(10.2, 7.6))
    gs = GridSpec(3, 1, figure=fig, height_ratios=[1.1, 0.9, 1.1],
                  hspace=0.42, left=0.10, right=0.97, top=0.93, bottom=0.08)

    # (a) the two MEs
    ax1 = fig.add_subplot(gs[0])
    ax1.plot(t_grid, ma, color=COLOUR[a], linewidth=1.8, label=f"{LABEL[a]}")
    ax1.plot(t_grid, mb, color=COLOUR[b], linewidth=1.8, label=f"{LABEL[b]}")
    ax1.axhline(0, color="0.6", linewidth=0.5)
    ax1.set_ylabel("ME = $\\beta_1 + 2\\beta_2 T$\n(growth per +1°C)", fontsize=9)
    ax1.set_title("Step 1.  Marginal-effect curves from the two reanalyses",
                  fontsize=10, loc="left")
    ax1.legend(fontsize=9, frameon=False, loc="upper right")
    ax1.tick_params(labelsize=8)
    for sp in ("top", "right"):
        ax1.spines[sp].set_visible(False)

    # (b) absolute difference
    ax2 = fig.add_subplot(gs[1])
    ax2.plot(t_grid, diff, color="black", linewidth=1.6)
    ax2.fill_between(t_grid, 0, diff, alpha=0.10, color="black")
    ax2.set_ylabel("$|\\mathrm{ME}_a - \\mathrm{ME}_b|$", fontsize=9)
    ax2.set_title("Step 2.  Pointwise absolute disagreement between the two MEs",
                  fontsize=10, loc="left")
    ax2.tick_params(labelsize=8)
    for sp in ("top", "right"):
        ax2.spines[sp].set_visible(False)

    # (c) pop-weighted T density + integrand
    ax3 = fig.add_subplot(gs[2])
    # Histogram of the empirical T values (pop-weighted) — normalized to a density
    bins = np.linspace(t_grid.min(), t_grid.max(), 40)
    counts, edges = np.histogram(temps, bins=bins, weights=weights, density=True)
    centres = 0.5 * (edges[1:] + edges[:-1])
    ax3.bar(centres, counts, width=(edges[1] - edges[0]) * 0.95,
            color="0.78", edgecolor="white", linewidth=0.3,
            label="pop-weighted empirical density of T  $f(T)$")

    # Integrand (|diff| × density) on a twin axis, shaded
    # Re-evaluate diff on the bin centres
    diff_b = np.abs(me(*coeffs[a], centres) - me(*coeffs[b], centres))
    integrand = diff_b * counts
    ax3b = ax3.twinx()
    ax3b.fill_between(centres, 0, integrand, color="tab:purple", alpha=0.45,
                      label="$|\\mathrm{ME}_a - \\mathrm{ME}_b| \\cdot f(T)$ (integrand)")
    ax3b.plot(centres, integrand, color="tab:purple", linewidth=1.3)
    ax3b.tick_params(labelsize=8, colors="tab:purple")
    ax3b.set_ylabel("integrand", fontsize=9, color="tab:purple")
    for sp in ("top",):
        ax3b.spines[sp].set_visible(False)

    ax3.set_xlabel("Annual mean temperature (°C)", fontsize=10)
    ax3.set_ylabel("density  $f(T)$", fontsize=9)
    ax3.tick_params(labelsize=8)
    for sp in ("top",):
        ax3.spines[sp].set_visible(False)
    ax3.set_title(
        f"Step 3.  Empirical T density (grey) and the integrand "
        f"$|\\Delta\\mathrm{{ME}}|\\cdot f(T)$ (purple).  "
        f"L1 = $\\int |\\Delta\\mathrm{{ME}}| f(T)\\, dT$ = {L1:.5f}",
        fontsize=10, loc="left")

    # Combined legend
    h1, l1 = ax3.get_legend_handles_labels()
    h2, l2 = ax3b.get_legend_handles_labels()
    ax3.legend(h1 + h2, l1 + l2, fontsize=8, frameon=False, loc="upper left")

    fig.suptitle(title, fontsize=12)
    return fig, L1


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print("Loading coefficients …")
    full = load_full_coeffs()
    loo  = load_loo_coeffs(WALKTHROUGH_COUNTRY)
    print(f"  full sample b1/b2 per reanalysis: ")
    for r in REANALYSES:
        print(f"    {r:10s} b1={full[r][0]:+.5f}  b2={full[r][1]:+.6f}    "
              f"(no {WALKTHROUGH_COUNTRY}: b1={loo[r][0]:+.5f}  b2={loo[r][1]:+.6f})")

    print("Loading regression input panel + pop weights …")
    df = load_reginputs()
    temps_full   = df["temp_popweight"].to_numpy()
    weights_full = df["Pop"].to_numpy()
    print(f"  full sample: N = {len(df):,} (country × year × reanalysis)")

    # The 'no-Bahrain' sample is the same panel with all Bahrain rows removed
    df_noB = df[df["countryname"] != WALKTHROUGH_COUNTRY]
    temps_noB   = df_noB["temp_popweight"].to_numpy()
    weights_noB = df_noB["Pop"].to_numpy()
    print(f"  no-{WALKTHROUGH_COUNTRY}: N = {len(df_noB):,}  "
          f"(removed {len(df) - len(df_noB)} Bahrain rows)")

    # Grid for plotting ME curves over the range of observed T
    t_grid = np.linspace(temps_full.min() - 1, temps_full.max() + 1, 400)

    # --- Figure 1: ME curves ---------------------------------------------
    print("Building figure: ME curves with vs without Bahrain …")
    fig1 = fig_me_curves(full, loo, t_grid)
    out1 = FIGS / "bahrain_me_curves.pdf"
    fig1.savefig(out1, format="pdf", bbox_inches="tight")
    fig1.savefig(out1.with_suffix(".png"), format="png", dpi=180,
                 bbox_inches="tight")
    print(f"  saved {out1}")

    # --- Figure 1 (contrast case): ME curves with vs without Greenland ----
    # Greenland is the largest *negative* country in the averaged metric
    # (drives agreement, i.e. removing it widens the gap). Useful contrast
    # to the Bahrain panel: the dashed lines should fan out, not in.
    print("Building figure: ME curves with vs without Greenland …")
    _orig_country = WALKTHROUGH_COUNTRY
    globals()["WALKTHROUGH_COUNTRY"] = "Greenland"
    try:
        loo_gl = load_loo_coeffs("Greenland")
        fig1g = fig_me_curves(full, loo_gl, t_grid)
        out1g = FIGS / "greenland_me_curves.pdf"
        fig1g.savefig(out1g, format="pdf", bbox_inches="tight")
        fig1g.savefig(out1g.with_suffix(".png"), format="png", dpi=180,
                      bbox_inches="tight")
        print(f"  saved {out1g}")
    finally:
        globals()["WALKTHROUGH_COUNTRY"] = _orig_country

    # --- Figure 2: L1 walkthrough, full sample ----------------------------
    print("Building figure: L1 walkthrough, full sample …")
    fig2, L1_full = fig_walkthrough(
        full, WALKTHROUGH_PAIR, temps_full, weights_full, t_grid,
        title=f"How the L1 distance is computed  ({LABEL[WALKTHROUGH_PAIR[0]]} vs "
              f"{LABEL[WALKTHROUGH_PAIR[1]]}, full sample)")
    out2 = FIGS / "bahrain_l1_walkthrough_full.pdf"
    fig2.savefig(out2, format="pdf", bbox_inches="tight")
    fig2.savefig(out2.with_suffix(".png"), format="png", dpi=180,
                 bbox_inches="tight")
    print(f"  saved {out2}    L1 (full)       = {L1_full:.6f}")

    # --- Figure 3: L1 walkthrough, no-Bahrain -----------------------------
    print("Building figure: L1 walkthrough, no-Bahrain sample …")
    fig3, L1_noB = fig_walkthrough(
        loo, WALKTHROUGH_PAIR, temps_noB, weights_noB, t_grid,
        title=f"Same calculation after dropping {WALKTHROUGH_COUNTRY} "
              f"(refit coeffs + drop Bahrain rows from $f(T)$)")
    out3 = FIGS / "bahrain_l1_walkthrough_noBahrain.pdf"
    fig3.savefig(out3, format="pdf", bbox_inches="tight")
    fig3.savefig(out3.with_suffix(".png"), format="png", dpi=180,
                 bbox_inches="tight")
    pct = (L1_full - L1_noB) / L1_full * 100
    print(f"  saved {out3}    L1 (no-{WALKTHROUGH_COUNTRY}) = {L1_noB:.6f}    "
          f"→  Δ = {pct:+.2f}%")


if __name__ == "__main__":
    main()
