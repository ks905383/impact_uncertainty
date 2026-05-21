"""
Country-influence analysis on the BHM (gdp-pc) dose-response.

For each country, compute % reduction in inter-reanalysis disagreement (mean
pairwise L1 distance between marginal-effect curves of the BHM dose-response)
when that country is dropped from the regression sample. The map / histogram /
country time-series figure shows where observational uncertainty is
concentrated.

Inputs (paths configurable via env vars KS_DROPBOX and KS_REPO):
    - {DROPBOX}/reg_outputs/gdppc/loos/*.csv     LOO coefficient files
    - {DROPBOX}/reg_outputs/gdppc/*.csv          full-sample coefficient files
    - {DROPBOX}/dmes/mean_dmes_gdppc.csv         pre-computed per-LOO mean dme
    - {DROPBOX}/reg_inputs/gdppc/reginputs_gdppc_025deg.csv   T data
    - {REPO}/aux_data/geo_data/country-fixed.shp ADM0 geometry

Output:
    - {REPO}/figures/loo_country_influence.pdf   single A4-portrait figure
    - {REPO}/figures/loo_country_influence.png   same, 200 dpi raster

Run:
    python tb/code/loo_country_influence.py
or:
    KS_DROPBOX=/path/to/ks KS_REPO=/path/to/repo python tb/code/loo_country_influence.py

Dependencies: pandas, numpy, matplotlib, geopandas (>=1.0 for make_valid).
The script does NOT use cartopy.
"""

from __future__ import annotations

import os
import re
import sys
import warnings
from glob import glob
from pathlib import Path

import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.colors import TwoSlopeNorm
from matplotlib.lines import Line2D


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

DROPBOX = Path(os.environ.get(
    "KS_DROPBOX", "/Users/tombearpark/Library/CloudStorage/Dropbox/ks"))
REPO = Path(os.environ.get(
    "KS_REPO", "/Users/tombearpark/Documents/GitHub/impact_uncertainty"))

LOO_DIR = DROPBOX / "reg_outputs" / "gdppc" / "loos"
FULL_DIR = DROPBOX / "reg_outputs" / "gdppc"
MEAN_DME_CSV = DROPBOX / "dmes" / "mean_dmes_gdppc.csv"
REG_INPUTS_CSV = DROPBOX / "reg_inputs" / "gdppc" / "reginputs_gdppc_025deg.csv"
SHP = REPO / "aux_data" / "geo_data" / "country-fixed.shp"

OUT_PDF = REPO / "figures" / "loo_country_influence.pdf"
# Mirror save into tb/figs/ as well — the deck's \graphicspath looks there.
OUT_PDF_MIRROR = REPO / "tb" / "figs" / "loo_country_influence.pdf"

# Four main reanalyses (matches the paper's main analysis and Kevin's dmes).
# Names here are the filename / obs_mod token; we display as REANALYSIS_LABEL.
REANALYSES = ("ERA5-025", "MERRA2", "JRA-3Q", "GMFD")
REANALYSIS_LABEL = {"ERA5-025": "ERA5", "MERRA2": "MERRA2",
                    "JRA-3Q": "JRA-3Q", "GMFD": "GMFD"}

# Reanalysis colours (matches code/funcs_plot.py:11)
REANALYSIS_COLOR = {"ERA5-025": "tab:red", "JRA-3Q": "tab:orange",
                    "MERRA2": "tab:green", "GMFD": "tab:blue"}

YEAR_LO, YEAR_HI = 1980, 2010
N_TOP = 3       # countries per direction in the time-series panel


# ---------------------------------------------------------------------------
# Loaders
# ---------------------------------------------------------------------------

# LOO filename pattern: regoutputs_coeffs_gdppc_025deg_{REA}_1980-2010_no{COUNTRY}.csv
_LOO_RE = re.compile(
    r"^regoutputs_coeffs_gdppc_025deg_(?P<rea>[^_]+)_1980-2010_no(?P<country>.+)\.csv$"
)
# Full-sample pattern (exclude *_noXxx.csv, *_JOINT*, *_ALLREA*, *_ALLYEARS*, *_vcv*)
_FULL_RE = re.compile(
    r"^regoutputs_coeffs_gdppc_025deg_(?P<rea>[^_]+)_1980-2010\.csv$"
)


def load_loo_coeffs() -> pd.DataFrame:
    """One row per (reanalysis, dropped country); columns b1, b2."""
    rows = []
    for fn in sorted(LOO_DIR.glob("regoutputs_coeffs_gdppc_*.csv")):
        m = _LOO_RE.match(fn.name)
        if not m:
            continue
        rea = m.group("rea")
        if rea not in REANALYSES:
            continue
        df = pd.read_csv(fn)
        # File has a one-line header `b1,b2` and a single data row
        b1, b2 = float(df["b1"].iloc[0]), float(df["b2"].iloc[0])
        rows.append({"reanalysis": rea, "country": m.group("country"),
                     "b1": b1, "b2": b2})
    out = pd.DataFrame(rows)
    if out.empty:
        raise RuntimeError(f"No LOO files found under {LOO_DIR}")
    return out


def load_full_coeffs() -> pd.DataFrame:
    """One row per reanalysis; columns b1, b2."""
    rows = []
    for fn in sorted(FULL_DIR.glob("regoutputs_coeffs_gdppc_025deg_*_1980-2010.csv")):
        m = _FULL_RE.match(fn.name)
        if not m:
            continue
        rea = m.group("rea")
        if rea not in REANALYSES:
            continue
        df = pd.read_csv(fn)
        rows.append({"reanalysis": rea,
                     "b1": float(df["b1"].iloc[0]),
                     "b2": float(df["b2"].iloc[0])})
    out = pd.DataFrame(rows).set_index("reanalysis")
    missing = set(REANALYSES) - set(out.index)
    if missing:
        raise RuntimeError(f"Missing full-sample coeffs for reanalyses: {missing}")
    return out


def load_reginputs() -> pd.DataFrame:
    """Country × year × reanalysis temperature panel used in BHM.

    Returns the columns needed downstream: countryname, year, obs_mod,
    temp_popweight, growthWDI.
    """
    use = ["countryname", "year", "obs_mod", "temp_popweight", "growthWDI"]
    df = pd.read_csv(REG_INPUTS_CSV, usecols=use)
    df = df[(df["year"] >= YEAR_LO) & (df["year"] <= YEAR_HI)]
    df = df[df["obs_mod"].isin(REANALYSES)]
    # The BHM regression drops rows with missing growthWDI, so do the same when
    # constructing the empirical T distribution that integrates the dme.
    df = df.dropna(subset=["temp_popweight", "growthWDI"])
    return df.reset_index(drop=True)


def load_kevin_mean_dme() -> pd.DataFrame:
    """Per-LOO mean dme as computed by Kevin (mean across 6 reanalysis pairs)."""
    df = pd.read_csv(MEAN_DME_CSV)
    return df.rename(columns={"loo_id": "country"})


# ---------------------------------------------------------------------------
# Metric: pairwise L1 distance between marginal-effect curves
# ---------------------------------------------------------------------------

def marginal_effect(b1: float, b2: float, t: np.ndarray) -> np.ndarray:
    """dY/dT for the BHM quadratic: y = b1*T + b2*T**2  ⇒  dy/dT = b1 + 2*b2*T."""
    return b1 + 2 * b2 * t


def pairwise_dme(coeffs: dict[str, tuple[float, float]],
                 t_emp: np.ndarray) -> float:
    """Mean over pairs of E_T[|ME_a(T) - ME_b(T)|] for the given coeff dict."""
    names = list(coeffs.keys())
    dists = []
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            ma = marginal_effect(*coeffs[a], t_emp)
            mb = marginal_effect(*coeffs[b], t_emp)
            dists.append(float(np.mean(np.abs(ma - mb))))
    if not dists:
        raise ValueError("Need at least 2 reanalyses to compute pairwise dme")
    return float(np.mean(dists))


# ---------------------------------------------------------------------------
# Country name reconciliation (regression sample ↔ shapefile)
# ---------------------------------------------------------------------------

# The country-fixed.shp file uses World Bank-style country names — the same
# convention as the regression input data. In practice no remapping is needed.
# Add entries here only if the merge below surfaces a real mismatch.
NAME_FIX: dict[str, str] = {}


# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------

def build_figure(pct: pd.DataFrame,
                 reg_inputs: pd.DataFrame,
                 baseline: float,
                 gdf: gpd.GeoDataFrame,
                 top: pd.DataFrame,
                 bottom: pd.DataFrame,
                 unmatched: list[str],
                 weighting_label: str = "pop-weighted") -> plt.Figure:
    """Compose the 3-panel figure.

    Layout (A4 portrait):
        row 0   — map (full width) + horizontal colorbar below it
        row 1   — histogram (left half) + time-series small multiples (right)
    """

    fig = plt.figure(figsize=(8.27, 10.8))
    # Two rows: map on top, histogram + time series below. Big hspace so the
    # colorbar (anchored beneath the map via inset_axes) has room.
    gs = GridSpec(2, 2,
                  height_ratios=[1.0, 1.05],
                  hspace=0.28, wspace=0.22,
                  left=0.06, right=0.97, top=0.95, bottom=0.07)

    # ---- (a) map --------------------------------------------------------
    ax_map = fig.add_subplot(gs[0, :])

    merged = gdf.merge(pct[["country", "pct_reduction"]],
                       left_on="locid", right_on="country", how="left")
    plotted = merged[~merged["pct_reduction"].isna()]
    unplotted = merged[merged["pct_reduction"].isna()]

    # Diverging scale centred at 0. Use the ACTUAL data extrema (not a 95th-pct
    # clip) so the most extreme countries — which are precisely the most
    # interesting ones — are visually distinguishable on the map. The two sides
    # are allowed to be asymmetric: in practice the negative tail (consensus
    # "glue" countries) is much longer than the positive tail.
    actual_min = float(pct["pct_reduction"].min())
    actual_max = float(pct["pct_reduction"].max())
    norm = TwoSlopeNorm(vmin=actual_min, vcenter=0, vmax=actual_max)

    # Diverging purple↔green: PRGn goes purple (low) → green (high); we reverse
    # it so positive = purple (drives disagreement) and negative = green
    # (drives agreement / consensus glue).
    cmap = "PRGn_r"
    unplotted.plot(ax=ax_map, color="lightgray", edgecolor="white",
                   linewidth=0.2)
    plotted.plot(ax=ax_map, column="pct_reduction", cmap=cmap,
                 norm=norm, edgecolor="white", linewidth=0.2, legend=False)
    ax_map.set_xlim(-180, 180)
    ax_map.set_ylim(-60, 85)
    ax_map.set_aspect("equal", adjustable="box")
    ax_map.set_xticks([]); ax_map.set_yticks([])
    for sp in ax_map.spines.values():
        sp.set_visible(False)
    ax_map.set_title(
        "(a) Country influence on cross-reanalysis disagreement in the BHM dose-response"
        + (f"   [{weighting_label}]" if weighting_label else ""),
        fontsize=10.5, loc="left", pad=4)

    # ---- colorbar anchored just below the map -------------------------
    cax = ax_map.inset_axes([0.22, -0.07, 0.56, 0.025])
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cb = fig.colorbar(sm, cax=cax, orientation="horizontal")
    cb.set_label(
        "% reduction in mean pairwise marginal-effect distance when country dropped\n"
        "(purple = country drove disagreement;  green = country drove agreement)",
        fontsize=8)
    cb.ax.tick_params(labelsize=8)

    # ---- (b) histogram --------------------------------------------------
    ax_h = fig.add_subplot(gs[1, 0])
    vals = pct["pct_reduction"].values
    ax_h.hist(vals, bins=30, color="0.55", edgecolor="white")
    ax_h.axvline(0, color="black", linewidth=0.6)
    med = float(np.median(vals))
    mn = float(np.mean(vals))
    ax_h.axvline(med, color="tab:red", linestyle="--", linewidth=0.9,
                 label=f"median = {med:+.2f}%")
    ax_h.axvline(mn, color="tab:blue", linestyle=":", linewidth=0.9,
                 label=f"mean = {mn:+.2f}%")
    ax_h.set_xlabel("% reduction", fontsize=9)
    ax_h.set_ylabel(f"countries (n = {len(vals)})", fontsize=9)
    ax_h.set_title("(b) Distribution across all countries",
                   fontsize=10.5, loc="left")
    ax_h.legend(fontsize=8, frameon=False, loc="upper left")
    ax_h.tick_params(labelsize=8)
    for sp in ("top", "right"):
        ax_h.spines[sp].set_visible(False)

    # ---- (c) time series small multiples -------------------------------
    panel_countries = list(top["country"]) + list(bottom["country"])
    panel_labels = {
        **{c: f"+{v:.1f}%" for c, v in zip(top["country"], top["pct_reduction"])},
        **{c: f"{v:.1f}%"  for c, v in zip(bottom["country"], bottom["pct_reduction"])},
    }

    # 4 rows × 2 cols inside the bottom-right cell: row 0 = title text strip;
    # rows 1–3 = the 6 country panels; bottom = legend.
    sub_gs = gs[1, 1].subgridspec(5, 2,
                                  height_ratios=[0.18, 1, 1, 1, 0.18],
                                  hspace=0.55, wspace=0.18)

    # Title strip for panel (c)
    ax_title = fig.add_subplot(sub_gs[0, :])
    ax_title.set_axis_off()
    ax_title.text(0, 0.4,
                  f"(c) Annual mean temperature [°C], top-{N_TOP} and bottom-{N_TOP} countries",
                  fontsize=10.5, weight="bold", transform=ax_title.transAxes)
    ax_title.text(0, -0.4,
                  "(purple = drives disagreement; green = drives agreement)",
                  fontsize=8, color="0.4", style="italic",
                  transform=ax_title.transAxes)

    # 3×2 panels of country time series. Each panel auto-scales its y-axis
    # independently — the absolute T levels vary by ~30°C across panels
    # (Greenland vs Bahrain) and a shared scale crushes the inter-reanalysis
    # disagreement within each panel into a near-flat line. The point of the
    # panel is to show the spread BETWEEN lines, not the level.
    for k, country in enumerate(panel_countries):
        r, c = divmod(k, 2)
        ax = fig.add_subplot(sub_gs[r + 1, c])   # rows 1, 2, 3
        d_country = reg_inputs[reg_inputs["countryname"] == country]
        for rea in REANALYSES:
            d = d_country[d_country["obs_mod"] == rea].sort_values("year")
            if d.empty:
                continue
            ax.plot(d["year"], d["temp_popweight"],
                    color=REANALYSIS_COLOR[rea], linewidth=0.9)
        ax.set_xlim(YEAR_LO - 0.5, YEAR_HI + 0.5)
        ax.tick_params(labelsize=7)
        if r != 2:
            ax.set_xticklabels([])
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        # Sign-coloured country label: purple for positive (drives
        # disagreement), green for negative (drives agreement). Hex values
        # match the deep ends of the PRGn_r colormap used on the map.
        val = float(pct.loc[pct["country"] == country, "pct_reduction"].iloc[0])
        swatch_color = "#762A83" if val > 0 else "#1B7837"
        ax.set_title(f"{country}  {panel_labels[country]}",
                     fontsize=8.5, loc="left", pad=2, color=swatch_color)

    # Legend strip at the bottom of panel (c)
    ax_leg = fig.add_subplot(sub_gs[4, :])
    ax_leg.set_axis_off()
    legend_handles = [Line2D([0], [0], color=REANALYSIS_COLOR[r], linewidth=1.4,
                              label=REANALYSIS_LABEL[r]) for r in REANALYSES]
    ax_leg.legend(handles=legend_handles, loc="center", ncol=4,
                  fontsize=9, frameon=False)

    # (x-axis label omitted — 1980 / 1990 / 2000 / 2010 tick labels speak for
    # themselves; y-axis unit is included in the panel (c) title)

    # ---- Footer --------------------------------------------------------
    foot = (f"baseline (Kevin's no-LOO 'all' row): {baseline:.4g}    |    "
            f"reanalyses: " + ", ".join(REANALYSIS_LABEL[r] for r in REANALYSES))
    if unmatched:
        foot += f"\n{len(unmatched)} regression-side countries unmatched to shapefile: " \
                + ", ".join(unmatched[:6]) + (" …" if len(unmatched) > 6 else "")
    fig.text(0.5, 0.015, foot, fontsize=7, color="0.4", ha="center")

    return fig


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print("Loading regression coefficients …")
    full = load_full_coeffs()
    loo = load_loo_coeffs()
    print(f"  full sample: {len(full)} reanalyses")
    print(f"  LOO: {loo['country'].nunique()} countries × "
          f"{loo['reanalysis'].nunique()} reanalyses = {len(loo)} rows")

    print("Loading regression input panel for T distribution …")
    reg_inputs = load_reginputs()
    t_emp = reg_inputs["temp_popweight"].to_numpy()
    print(f"  empirical T sample: {t_emp.size} (country × year × reanalysis) "
          f"rows, mean {t_emp.mean():.2f} °C, "
          f"5–95 pct {np.percentile(t_emp, 5):.1f}–{np.percentile(t_emp, 95):.1f} °C")

    print("Loading Kevin's pre-computed per-LOO mean dme …")
    kevin = load_kevin_mean_dme()
    # Kevin's CSV contains a special row loo_id == 'all' which is the no-LOO
    # full-sample baseline (computed with his exact T distribution / weighting).
    # That row is the natural baseline — use it instead of recomputing, so
    # baseline and per-LOO values are on the same footing.
    if "all" not in set(kevin["country"]):
        raise RuntimeError("Expected an 'all' row (no-LOO baseline) in "
                           "mean_dmes_gdppc.csv but didn't find one.")
    baseline = float(kevin.loc[kevin["country"] == "all", "dme_mu"].iloc[0])
    kevin = kevin[kevin["country"] != "all"].reset_index(drop=True)
    print(f"  baseline (Kevin's 'all' row) = {baseline:.6f}")
    print(f"  {len(kevin)} LOO countries; dme_mu range "
          f"[{kevin['dme_mu'].min():.4f}, {kevin['dme_mu'].max():.4f}]")

    # --- Calibration: my own formula vs Kevin's, both on the FULL sample. ----
    # Same offset on baseline and on each LOO ⇒ ratios are consistent.
    full_dict = {r: (full.loc[r, "b1"], full.loc[r, "b2"]) for r in REANALYSES}
    my_baseline = pairwise_dme(full_dict, t_emp)
    ratio = my_baseline / baseline
    print(f"  calibration: my full-sample dme={my_baseline:.6f}, "
          f"Kevin={baseline:.6f}, ratio={ratio:.4f} "
          f"(small offset OK — different T-distribution weighting; we use "
          f"Kevin's value as the baseline so the % reduction is consistent)")
    if not (0.5 < ratio < 2.0):
        raise RuntimeError(
            f"Calibration ratio my_baseline/Kevin = {ratio:.3f} is outside "
            "[0.5, 2.0]. The formula or T distribution disagrees by more "
            "than a small weighting difference — investigate before "
            "trusting the figure.")

    # Cross-check: re-rank countries by Kevin's dme_mu vs my own LOO
    # recomputation. Same ranking ⇒ same map story.
    my_loo_dme = []
    missing_loo = []
    for c in kevin["country"]:
        loo_c = loo[loo["country"] == c]
        if len(loo_c) != len(REANALYSES):
            my_loo_dme.append(np.nan)
            missing_loo.append(c)
            continue
        loo_dict = {r: (loo_c.loc[loo_c["reanalysis"] == r, "b1"].iloc[0],
                        loo_c.loc[loo_c["reanalysis"] == r, "b2"].iloc[0])
                    for r in REANALYSES}
        my_loo_dme.append(pairwise_dme(loo_dict, t_emp))
    rank_corr = float(pd.Series(my_loo_dme).rank().corr(
        kevin["dme_mu"].rank(), method="spearman"))
    print(f"  Spearman rank correlation (my LOO dme vs Kevin's): {rank_corr:.4f}")
    if rank_corr < 0.90:
        raise RuntimeError(
            f"Spearman rank correlation {rank_corr:.3f} between my LOO dme "
            "and Kevin's is below 0.90 — the formulas disagree enough to "
            "change the map story.")
    if missing_loo:
        print(f"  ⚠ {len(missing_loo)} countries in Kevin's CSV have no "
              f"matching LOO files (skipped from calibration; their map "
              f"colour still relies on Kevin's value alone): {missing_loo}")

    # --- Compute % reduction --------------------------------------------------
    pct = kevin.copy()
    pct["pct_reduction"] = (baseline - pct["dme_mu"]) / baseline * 100.0
    pct = pct.sort_values("pct_reduction", ascending=False).reset_index(drop=True)
    print(f"% reduction: median {pct['pct_reduction'].median():+.2f}, "
          f"mean {pct['pct_reduction'].mean():+.2f}, "
          f"min {pct['pct_reduction'].min():+.2f}, "
          f"max {pct['pct_reduction'].max():+.2f}")

    top = pct.head(N_TOP).copy()
    bottom = pct.tail(N_TOP).iloc[::-1].copy()   # most negative first
    print("Top countries (most + % reduction → drove the disagreement):")
    print(top.to_string(index=False))
    print("Bottom countries (most − % reduction → held the consensus):")
    print(bottom.to_string(index=False))

    # --- Load shapefile and merge --------------------------------------------
    print("Loading ADM0 shapefile …")
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        gdf = gpd.read_file(SHP)
    gdf["geometry"] = gdf.make_valid()
    gdf = gdf.rename(columns={"CNTRY_NAME": "locid"})
    gdf = gdf.query("locid != 'Antarctica'")

    # Apply name fixes to the pct table to align with shapefile
    pct["country_orig"] = pct["country"]
    pct["country"] = pct["country"].replace(NAME_FIX)
    assert pct["country"].is_unique, "duplicate country names in pct after NAME_FIX"

    pct_in_gdf = pct["country"].isin(gdf["locid"]).sum()
    unmatched_pct = sorted(set(pct["country"]) - set(gdf["locid"]))
    # The other direction is informational: lots of shapefile polygons (small
    # islands, dependencies) legitimately have no WDI GDP data and so no LOO.
    shp_extras = sorted(set(gdf["locid"]) - set(pct["country"]))
    print(f"  shapefile: {len(gdf)} polygons; "
          f"{pct_in_gdf} / {len(pct)} regression countries matched")
    if unmatched_pct:
        print(f"  ⚠ regression-side unmatched (silently grey on the map): "
              f"{unmatched_pct}")
    print(f"  ({len(shp_extras)} shapefile polygons have no regression entry "
          f"and are drawn grey: e.g. {shp_extras[:6]}…)")

    # --- Build and save pop-weighted figure ----------------------------------
    print("Building pop-weighted figure …")
    fig = build_figure(pct, reg_inputs, baseline, gdf, top, bottom,
                       unmatched_pct,
                       weighting_label="pop-weighted (Kevin's CSV)")
    OUT_PDF.parent.mkdir(parents=True, exist_ok=True)
    OUT_PDF_MIRROR.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PDF, format="pdf", bbox_inches="tight")
    fig.savefig(OUT_PDF_MIRROR, format="pdf", bbox_inches="tight")
    png = OUT_PDF.with_suffix(".png")
    fig.savefig(png, format="png", dpi=200, bbox_inches="tight")
    fig.savefig(OUT_PDF_MIRROR.with_suffix(".png"), format="png", dpi=200,
                bbox_inches="tight")
    print(f"Saved: {OUT_PDF}  (mirror: {OUT_PDF_MIRROR})")

    # --- Build uniform-weighted variant via shared helper --------------------
    # Same figure layout, but pct_reduction is recomputed with each
    # (country × year × reanalysis) row contributing equally instead of being
    # population-weighted.
    print("\n=== Building uniform-weighted variant ===")
    import sys
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import funcs_dme  # noqa: E402

    uni = funcs_dme.compute_mean_dme_table(weighting="uniform")
    uni_baseline = float(uni.loc[uni["loo_id"] == "all", "dme_mu"].iloc[0])
    uni = (uni[uni["loo_id"] != "all"]
              .copy()
              .rename(columns={"loo_id": "country"}))
    uni = uni.sort_values("pct_reduction", ascending=False).reset_index(drop=True)
    print(f"  baseline (uniform) = {uni_baseline:.6f}")
    print(f"  pct range: [{uni['pct_reduction'].min():+.2f}, "
          f"{uni['pct_reduction'].max():+.2f}]")
    uni["country_orig"] = uni["country"]
    uni["country"] = uni["country"].replace(NAME_FIX)
    uni_unmatched = sorted(set(uni["country"]) - set(gdf["locid"]))

    uni_top = uni.head(N_TOP).copy()
    uni_bottom = uni.tail(N_TOP).iloc[::-1].copy()
    print("Uniform top countries (drives disagreement):")
    print(uni_top.to_string(index=False))
    print("Uniform bottom countries (drives agreement):")
    print(uni_bottom.to_string(index=False))

    fig_u = build_figure(uni, reg_inputs, uni_baseline, gdf,
                         uni_top, uni_bottom, uni_unmatched,
                         weighting_label="uniform (each obs equal weight)")
    out_u = OUT_PDF.with_name(OUT_PDF.stem + "_uniform.pdf")
    out_u_m = OUT_PDF_MIRROR.with_name(OUT_PDF_MIRROR.stem + "_uniform.pdf")
    fig_u.savefig(out_u, format="pdf", bbox_inches="tight")
    fig_u.savefig(out_u_m, format="pdf", bbox_inches="tight")
    fig_u.savefig(out_u.with_suffix(".png"), format="png", dpi=200,
                  bbox_inches="tight")
    fig_u.savefig(out_u_m.with_suffix(".png"), format="png", dpi=200,
                  bbox_inches="tight")
    print(f"Saved: {out_u}  (mirror: {out_u_m})")


if __name__ == "__main__":
    main()
