"""
Produce two Extended Data figures for the paper.

ED Fig 1 — Inter-reanalysis disagreement in the temperature SHOCKS used in
            the BHM regression.
   Row 1  (2 × 2 grid)  Residualised T time series for the four countries
                        that agree most (top) and disagree most (bottom)
                        across the 4 reanalyses. Residuals are the BHM
                        regression-equivalent residuals from
                        `T ~ 1 | iso + year + iso[t] + iso[t²]`.
   Row 2  (1 × 2 grid)  Country-level scatter: the per-country mean of the
                        pairwise |Δresidual T| (averaged over 6 pairs)
                        plotted against (a) country-mean T and (b) country-
                        mean log GDP per capita. OLS overlay with R².

ED Fig 2 — Country-by-country drivers of disagreement in the FITTED BHM
            dose-response, in absolute L1 units (uniform weighting:
            each (country × year × reanalysis) observation enters equally).
   Top    (2 × 3 grid)  Per-pair maps of the absolute reduction in the
                        pair's L1 distance when that country is dropped
                        from the regression. Shared diverging colour scale
                        so the six pairs are directly comparable.
   Bottom (1 × 1 wide)  Combined map: mean across the 6 pairs of the same
                        per-country absolute reduction.

Each figure is portrait, ~7" × 9.5", and is written to this folder as both
PDF (vector, for the paper) and PNG (200 dpi, for previewing).

The script is self-contained. It depends only on
  numpy, pandas, scipy, geopandas, matplotlib.
Inputs (paths configurable via env vars KS_DROPBOX / KS_REPO):
  - reg_outputs/gdppc/regoutputs_coeffs_<rea>_1980-2010.csv
  - reg_outputs/gdppc/loos/regoutputs_coeffs_<rea>_1980-2010_no<country>.csv
  - reg_inputs/gdppc/reginputs_gdppc_025deg.csv
  - aux_data/geo_data/country-fixed.shp           (ADM0 polygons)
"""

from __future__ import annotations

import os
import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.colors import TwoSlopeNorm, SymLogNorm, Normalize, LinearSegmentedColormap
from matplotlib.lines import Line2D
from scipy import sparse
from scipy.sparse.linalg import lsqr
from scipy import stats as sps

import geopandas as gpd

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

DROPBOX = Path(os.environ.get(
    "KS_DROPBOX", "/glade/derecho/scratch/schwarzwald/impact_uncertainty/"))
REPO    = Path(os.environ.get(
    "KS_REPO",    "/glade/u/home/schwarzwald/projects/impact_uncertainty"))

REG_INPUT = DROPBOX / "reg_inputs"  / "gdppc" / "reginputs_gdppc_025deg.csv"
FULL_DIR  = DROPBOX / "reg_outputs" / "gdppc"
LOO_DIR   = DROPBOX / "reg_outputs" / "gdppc" / "loos"
SHP       = REPO / "aux_data" / "geo_data" / "country-fixed.shp"

OUT_DIR   = Path(__file__).resolve().parent
OUT_FIG1  = OUT_DIR / "extdata_fig1_temp_disagreement.pdf"
OUT_FIG2  = OUT_DIR / "extdata_fig2_me_disagreement.pdf"

# The 4 main reanalyses; filenames use "ERA5-025" but we display "ERA5".
REANALYSES = ("ERA5-025", "MERRA2", "JRA-3Q", "GMFD")
LABEL  = {"ERA5-025": "ERA5", "MERRA2": "MERRA2",
          "JRA-3Q": "JRA-3Q", "GMFD": "GMFD"}
# Colour-blind-friendly four-tone palette (Wong 2011)
COLOUR = {"ERA5-025": "#E69F00",   # orange
          "MERRA2":   "#0072B2",   # blue
          "JRA-3Q":   "#009E73",   # green
          "GMFD":     "#CC79A7"}   # magenta

YEAR_LO, YEAR_HI = 1980, 2010

# 6 unordered pairs in a canonical order (smallest baseline first)
PAIR_ORDER = [
    ("MERRA2", "GMFD"),
    ("JRA-3Q", "ERA5-025"),
    ("MERRA2", "ERA5-025"),
    ("GMFD",   "ERA5-025"),
    ("MERRA2", "JRA-3Q"),
    ("JRA-3Q", "GMFD"),
]

# Publication-quality matplotlib defaults
plt.rcParams.update({
    "font.family":      "sans-serif",
    "font.sans-serif":  ["Helvetica", "Arial", "DejaVu Sans"],
    "font.size":        8.5,
    "axes.labelsize":   9,
    "axes.titlesize":   9.5,
    "xtick.labelsize":  7.5,
    "ytick.labelsize":  7.5,
    "legend.fontsize":  8,
    "figure.titlesize": 11,
    "pdf.fonttype":     42,   # embed TrueType (editable in Illustrator)
    "ps.fonttype":      42,
    "savefig.bbox":     "tight",
})


# ---------------------------------------------------------------------------
# Loaders
# ---------------------------------------------------------------------------

_FULL_RE = re.compile(
    r"^regoutputs_coeffs_gdppc_025deg_(?P<rea>[^_]+)_1980-2010\.csv$")
_LOO_RE  = re.compile(
    r"^regoutputs_coeffs_gdppc_025deg_(?P<rea>[^_]+)_1980-2010_no(?P<country>.+)\.csv$")


def load_panel(drop_oosamp = True) -> pd.DataFrame:
    """BHM regression-input panel, balanced across the 4 main reanalyses."""
    use = ["countryname", "iso", "year", "obs_mod", "temp_popweight", "Pop",
           "growthWDI", "time", "time2", "gdpCAPppp_wdi","og_insamp"]
    df = pd.read_csv(REG_INPUT, usecols=use)
    df = df[(df["year"] >= YEAR_LO) & (df["year"] <= YEAR_HI)]
    df = df[df["obs_mod"].isin(REANALYSES)]
    df = df.dropna(subset=["temp_popweight"])
     # Drop locs not in the original sample
    df = df.where(df.og_insamp).dropna(how='all')
    # Require every reanalysis to observe each (iso, year)
    counts = df.groupby(["iso", "year"])["obs_mod"].nunique().reset_index()
    keep = counts[counts["obs_mod"] == len(REANALYSES)][["iso", "year"]]
    df = df.merge(keep, on=["iso", "year"], how="inner").reset_index(drop=True)
    return df


def load_full_coeffs() -> dict[str, tuple[float, float]]:
    out: dict[str, tuple[float, float]] = {}
    for fn in sorted(FULL_DIR.glob("regoutputs_coeffs_gdppc_025deg_*_1980-2010.csv")):
        m = _FULL_RE.match(fn.name)
        if not m or m.group("rea") not in REANALYSES:
            continue
        d = pd.read_csv(fn)
        out[m.group("rea")] = (float(d["b1"].iloc[0]), float(d["b2"].iloc[0]))
    missing = set(REANALYSES) - set(out)
    if missing:
        raise RuntimeError(f"missing full-sample coeffs for: {missing}")
    return out


def load_loo_coeffs() -> dict[tuple[str, str], tuple[float, float]]:
    out: dict[tuple[str, str], tuple[float, float]] = {}
    for fn in sorted(LOO_DIR.glob("regoutputs_coeffs_gdppc_*.csv")):
        m = _LOO_RE.match(fn.name)
        if not m or m.group("rea") not in REANALYSES:
            continue
        d = pd.read_csv(fn)
        out[(m.group("rea"), m.group("country"))] = (
            float(d["b1"].iloc[0]), float(d["b2"].iloc[0]))
    return out


# ---------------------------------------------------------------------------
# Residualisation: T ~ 1 | iso + year + iso[time1] + iso[time2]
#   (per reanalysis, using sparse OLS)
# ---------------------------------------------------------------------------

def residualise_one(df_rea: pd.DataFrame) -> np.ndarray:
    iso_c = pd.Categorical(df_rea["iso"])
    yr_c  = pd.Categorical(df_rea["year"])
    t1    = df_rea["time"].to_numpy(dtype=float)
    t2    = df_rea["time2"].to_numpy(dtype=float)
    n     = len(df_rea)
    rows  = np.arange(n)
    iso_codes, yr_codes = iso_c.codes, yr_c.codes
    n_iso, n_yr = len(iso_c.categories), len(yr_c.categories)
    iso_blk = sparse.csc_matrix((np.ones(n), (rows, iso_codes)), shape=(n, n_iso))
    yr_blk  = sparse.csc_matrix((np.ones(n), (rows, yr_codes )), shape=(n, n_yr))
    yr_blk  = yr_blk[:, 1:]   # drop one year for identification
    iso_t1  = sparse.csc_matrix((t1, (rows, iso_codes)), shape=(n, n_iso))
    iso_t2  = sparse.csc_matrix((t2, (rows, iso_codes)), shape=(n, n_iso))
    X = sparse.hstack([iso_blk, yr_blk, iso_t1, iso_t2]).tocsr()
    y = df_rea["temp_popweight"].to_numpy(dtype=float)
    beta = lsqr(X, y, atol=1e-10, btol=1e-10, iter_lim=20000)[0]
    return y - X @ beta


def residualise_all(df: pd.DataFrame) -> pd.DataFrame:
    parts = []
    for rea in REANALYSES:
        sub = df[df["obs_mod"] == rea].copy()
        sub["resid"] = residualise_one(sub)
        parts.append(sub[["countryname", "iso", "year", "obs_mod",
                          "temp_popweight", "resid"]])
    return pd.concat(parts, ignore_index=True)


# ---------------------------------------------------------------------------
# Pairwise quantities
# ---------------------------------------------------------------------------

def per_pair_country_mad(resid: pd.DataFrame) -> pd.DataFrame:
    """Per (pair, country): mean over years of |resid_a - resid_b|."""
    wide = resid.pivot_table(index=["iso", "countryname", "year"],
                             columns="obs_mod", values="resid").reset_index()
    rows = []
    for (a, b) in PAIR_ORDER:
        d = wide.dropna(subset=[a, b]).copy()
        d["d"] = (d[a] - d[b]).abs()
        m = (d.groupby(["iso", "countryname"], as_index=False)
              .agg(mad=("d", "mean")))
        m["pair_a"] = a; m["pair_b"] = b
        rows.append(m)
    return pd.concat(rows, ignore_index=True)


def country_covariates(df: pd.DataFrame) -> pd.DataFrame:
    sub = df[df["obs_mod"] == REANALYSES[0]]
    g = (sub.groupby(["iso", "countryname"])
            .agg(mean_T=("temp_popweight", "mean"),
                 mean_gdp=("gdpCAPppp_wdi", "mean"))
            .reset_index())
    g["log_gdp"] = np.log(g["mean_gdp"])
    return g


# ---------------------------------------------------------------------------
# L1 distance between marginal-effect curves (uniform weighting)
# ---------------------------------------------------------------------------

def _me(b1: float, b2: float, t: np.ndarray) -> np.ndarray:
    """BHM marginal effect dy/dT = b1 + 2 b2 T."""
    return b1 + 2.0 * b2 * t


def pair_l1_uniform(coef_a, coef_b, panel, pair) -> float:
    """L1 between MEs of the two reanalyses in `pair`, integrated UNIFORMLY
    over the empirical T sample from rows of `panel` belonging to either
    reanalysis."""
    a, b = pair
    s = panel[panel["obs_mod"].isin((a, b))]
    if s.empty:
        return float("nan")
    t = s["temp_popweight"].to_numpy()
    return float(np.mean(np.abs(_me(*coef_a, t) - _me(*coef_b, t))))


def compute_perpair_abs_reduction(panel: pd.DataFrame,
                                   full: dict, loo: dict) -> pd.DataFrame:
    """For each (pair, country): abs_reduction = L1_full − L1_drop_country,
    in raw L1 units (same scale across pairs).  Includes loo_id == 'all'
    rows for the per-pair baselines (abs_reduction = 0)."""
    countries = sorted(panel["countryname"].unique())
    rows = []
    for (a, b) in PAIR_ORDER:
        L1_full = pair_l1_uniform(full[a], full[b], panel, (a, b))
        rows.append({"pair_a": a, "pair_b": b, "country": "all",
                     "L1": L1_full, "abs_reduction": 0.0})
        for c in countries:
            if (a, c) not in loo or (b, c) not in loo:
                continue
            sub = panel[panel["countryname"] != c]
            L1_c = pair_l1_uniform(loo[(a, c)], loo[(b, c)], sub, (a, b))
            rows.append({"pair_a": a, "pair_b": b, "country": c,
                         "L1": L1_c,
                         "abs_reduction": L1_full - L1_c})
    return pd.DataFrame(rows)


def compute_mean_abs_reduction(perpair: pd.DataFrame) -> pd.DataFrame:
    """Mean absolute reduction across the 6 pairs, per country."""
    d = perpair[perpair["country"] != "all"]
    return (d.groupby("country", as_index=False)
             .agg(mean_abs_reduction=("abs_reduction", "mean")))


# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------

def _ols_overlay(ax, x, y, color="0.15", show_stats=False, loc="lower right"):
    """OLS fit dashed line. If show_stats=True, also write R² / p in a box
    at the given `loc` corner."""
    x = np.asarray(x, float); y = np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    x, y = x[m], y[m]
    if x.size < 3:
        return
    r = sps.linregress(x, y)
    xs = np.linspace(x.min(), x.max(), 60)
    ax.plot(xs, r.slope * xs + r.intercept,
            color=color, linewidth=1.0, linestyle="--", zorder=2)
    if not show_stats:
        return
    p_txt = f"p < 0.001" if r.pvalue < 1e-3 else f"p = {r.pvalue:.2g}"
    anchors = {
        "upper right": (0.97, 0.96, "right", "top"),
        "upper left":  (0.03, 0.96, "left",  "top"),
        "lower right": (0.97, 0.04, "right", "bottom"),
        "lower left":  (0.03, 0.04, "left",  "bottom"),
    }
    ax_x, ax_y, ha, va = anchors[loc]
    ax.text(ax_x, ax_y,
            f"$R^2$ = {r.rvalue**2:.2f}\n{p_txt}",
            transform=ax.transAxes, ha=ha, va=va,
            fontsize=7, color=color,
            bbox=dict(boxstyle="round,pad=0.22",
                      facecolor="white", edgecolor="0.85", linewidth=0.4))


def _despine(ax, sides=("top", "right")):
    for s in sides:
        ax.spines[s].set_visible(False)


def load_shapefile() -> gpd.GeoDataFrame:
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        gdf = gpd.read_file(SHP)
    gdf["geometry"] = gdf.make_valid()
    gdf = gdf.rename(columns={"CNTRY_NAME": "locid"})
    return gdf.query("locid != 'Antarctica'").reset_index(drop=True)


def _world_map(ax):
    ax.set_xlim(-180, 180)
    ax.set_ylim(-60, 85)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)


# ===========================================================================
# Figure 1 — temperature shock disagreement
# ===========================================================================

def build_fig1(resid: pd.DataFrame,
               pair_mad: pd.DataFrame,
               cov: pd.DataFrame) -> plt.Figure:

    # ---- Pick 2 most-agreeing + 2 most-disagreeing countries by avg MAD ----
    per_country = (pair_mad.groupby(["iso", "countryname"], as_index=False)
                            .agg(mean_mad=("mad", "mean")))
    per_country = per_country.merge(cov[["iso", "log_gdp"]], on="iso",
                                     how="left")
    # Require at least 25 years of residual data for sensible time-series panels
    yr_count = (resid.groupby("countryname")
                     .agg(n_years=("year", "nunique"))
                     .reset_index())
    per_country = per_country.merge(yr_count, on="countryname", how="left")
    per_country = per_country[per_country["n_years"] >= 25].reset_index(drop=True)

    agree = per_country.nsmallest(2, "mean_mad").sort_values("mean_mad")
    disagr = per_country.nlargest(2, "mean_mad").sort_values("mean_mad", ascending=False)
    print("Agreeing 2: " + ", ".join(agree["countryname"].tolist()))
    print("Disagreeing 2: " + ", ".join(disagr["countryname"].tolist()))
    panel_countries = (
        agree.assign(grp="most agreement").to_dict("records") +
        disagr.assign(grp="most disagreement").to_dict("records"))

    # ---- Build figure --------------------------------------------------
    fig = plt.figure(figsize=(7.0, 7.0))
    gs = GridSpec(3, 2, figure=fig,
                  height_ratios=[1.0, 1.0, 1.55],
                  hspace=0.50, wspace=0.28,
                  left=0.11, right=0.97, top=0.93, bottom=0.08)

    # ---- Time series of residuals (top 2x2 grid) -----------------------
    ts_letters = ["a", "b", "c", "d"]
    for k, row in enumerate(panel_countries):
        r, c = divmod(k, 2)
        ax = fig.add_subplot(gs[r, c])
        d_country = resid[resid["countryname"] == row["countryname"]]
        for rea in REANALYSES:
            d = d_country[d_country["obs_mod"] == rea].sort_values("year")
            if d.empty:
                continue
            ax.plot(d["year"], d["resid"],
                    color=COLOUR[rea], linewidth=0.9, alpha=0.95)
        ax.axhline(0, color="0.7", linewidth=0.4, zorder=0)
        ax.tick_params(labelsize=7)
        _despine(ax)
        ax.set_xlim(YEAR_LO - 0.5, YEAR_HI + 0.5)
        tag = "most agreement" if row["grp"] == "most agreement" else "most disagreement"
        ax.set_title(
            f"({ts_letters[k]})  {row['countryname']}  ({tag}; mean |Δresid T| = "
            f"{row['mean_mad']:.2f} °C)",
            fontsize=8.5, loc="left", pad=2)
        if c == 0:
            ax.set_ylabel("Residual T (°C)", fontsize=8.5)
        if r == 0:
            ax.set_xticklabels([])

    # ---- Shared 4-reanalysis legend strip across panel row -------------
    handles = [Line2D([0], [0], color=COLOUR[r], linewidth=1.6,
                       label=LABEL[r]) for r in REANALYSES]
    fig.legend(handles=handles, loc="upper center",
               ncol=4, frameon=False,
               bbox_to_anchor=(0.5, 0.985), fontsize=8.5)

    # ---- Bottom row: 2 scatters ----------------------------------------
    # Per-country mean MAD (across the 6 pairs) is the y-axis of both
    pc_x = per_country.merge(cov[["iso", "mean_T"]], on="iso", how="left")

    panels = [
        ("mean_T",  "Country-mean annual T (°C, 1980–2010)",
         "(e) vs. country-mean temperature"),
        ("log_gdp", "log mean GDP per capita (PPP, 1980–2010)",
         "(f) vs. country-mean log income"),
    ]
    for k, (x_col, x_lab, title) in enumerate(panels):
        ax = fig.add_subplot(gs[2, k])
        sub = pc_x.dropna(subset=[x_col])
        ax.scatter(sub[x_col], sub["mean_mad"],
                   s=12, color="0.35", alpha=0.7,
                   edgecolor="white", linewidth=0.3)
        _ols_overlay(ax, sub[x_col].to_numpy(), sub["mean_mad"].to_numpy())
        # Label top-3 / bottom-3 by mean_mad
        labels = pd.concat([sub.nlargest(3, "mean_mad"),
                            sub.nsmallest(3, "mean_mad")])
        for _, row in labels.iterrows():
            ax.annotate(row["countryname"],
                        (row[x_col], row["mean_mad"]),
                        fontsize=6.5, ha="left", va="bottom",
                        xytext=(2, 2), textcoords="offset points",
                        color="0.2")
        ax.set_xlabel(x_lab, fontsize=8.5)
        if k == 0:
            ax.set_ylabel(
                "Mean |Δresidual T| (°C)\n"
                "(country mean, averaged over 6 pairs)",
                fontsize=8.5)
        ax.set_title(title, fontsize=9, loc="left")
        ax.tick_params(labelsize=7.5)
        _despine(ax)

    # Overall figure title
    fig.suptitle(
        "Extended Data Fig. 1   Inter-reanalysis disagreement in "
        "temperature shocks used to fit the BHM dose-response",
        fontsize=10.5, x=0.5, y=0.995)
    return fig


# ===========================================================================
# Figure 2 — drivers of marginal-effect disagreement
# ===========================================================================

def build_fig2(perpair_abs: pd.DataFrame,
               mean_abs: pd.DataFrame,
               gdf: gpd.GeoDataFrame) -> plt.Figure:

    # Shorter overall figure + a single shared colourbar at the bottom (per-pair
    # values and per-country-mean values are on the same L1 units, so they
    # can share one scale). Tight hspaces so the maps have most of the area.
    fig = plt.figure(figsize=(7.0, 7.2))
    outer = GridSpec(3, 1, figure=fig,
                     height_ratios=[14, 8, 0.6],
                     hspace=0.30,
                     left=0.03, right=0.97, top=0.93, bottom=0.07)

    # ---- Top: 2×3 per-pair maps (sequential purple, positive only) ----
    # Focus on countries that DRIVE DISAGREEMENT (abs_reduction > 0); show
    # negative-reduction (drives-agreement) and no-data countries in light
    # grey. Sequential purple scale shared across all 6 pairs so the
    # absolute scale is comparable.
    perpair = perpair_abs[perpair_abs["country"] != "all"].copy()
    gmin = float(perpair["abs_reduction"].min())
    gmax = float(perpair["abs_reduction"].max())
    # Saturate at p75 of positive values so the top quartile of drivers
    # reach the deepest purple and the rest get a clearly visible shade.
    pos_vals = perpair.loc[perpair["abs_reduction"] > 0, "abs_reduction"].to_numpy()
    pos_sat = float(np.percentile(pos_vals, 75))
    norm = Normalize(vmin=0, vmax=pos_sat)
    # High-contrast purple ramp (truncated `Purples`, skipping the lightest
    # 35%) so even small positive values render visibly.
    _purples_base = plt.get_cmap("Purples")
    cmap = LinearSegmentedColormap.from_list(
        "PurplesHC", _purples_base(np.linspace(0.35, 1.0, 256)))

    top_gs = outer[0].subgridspec(2, 3, hspace=0.18, wspace=0.04)
    for i, (a, b) in enumerate(PAIR_ORDER):
        r, c = divmod(i, 3)
        ax = fig.add_subplot(top_gs[r, c])
        pct = perpair[(perpair["pair_a"] == a) & (perpair["pair_b"] == b)]
        merged = gdf.merge(pct[["country", "abs_reduction"]],
                            left_on="locid", right_on="country", how="left")
        # Positive-only: countries with abs_reduction > 0 painted purple;
        # countries with abs_reduction <= 0 or no data drawn light grey.
        plotted   = merged[merged["abs_reduction"] > 0]
        unplotted = merged[~(merged["abs_reduction"] > 0)]
        unplotted.plot(ax=ax, color="0.92", edgecolor="white", linewidth=0.10)
        plotted.plot(ax=ax, column="abs_reduction", cmap=cmap, norm=norm,
                     edgecolor="white", linewidth=0.10, legend=False)
        _world_map(ax)
        n_pos = int((merged["abs_reduction"] > 0).sum())
        n_tot = int(merged["abs_reduction"].notna().sum())
        ax.set_title(
            f"{LABEL[a]} vs {LABEL[b]}   (n={n_pos}/{n_tot})",
            fontsize=8.5, loc="left", pad=1)

    # Top panel header — anchored to the actual top of the top-row maps
    top_anchor = top_gs[0, 1].get_position(fig)
    fig.text(0.5, top_anchor.y1 + 0.012,
             "(a) Per-pair country influence on cross-reanalysis disagreement  "
             "[uniform, absolute L1]",
             fontsize=10, weight="bold", ha="center")

    # ---- Combined map (mean across pairs) — same colour scale as top ----
    bot_ax = fig.add_subplot(outer[1])
    mean_merged = gdf.merge(
        mean_abs.rename(columns={"country": "locid_join"}),
        left_on="locid", right_on="locid_join", how="left")
    mx = float(mean_merged["mean_abs_reduction"].max())
    plotted   = mean_merged[mean_merged["mean_abs_reduction"] > 0]
    unplotted = mean_merged[~(mean_merged["mean_abs_reduction"] > 0)]
    unplotted.plot(ax=bot_ax, color="0.92", edgecolor="white", linewidth=0.10)
    # Reuse the SAME `norm` (and therefore the same colour scale) as the
    # per-pair maps above. Both quantities are in raw L1 units; the
    # combined-map values are smaller (averaged across 6 pairs) so they
    # naturally render in the mid range of the shared colourbar.
    plotted.plot(ax=bot_ax, column="mean_abs_reduction", cmap=cmap, norm=norm,
                 edgecolor="white", linewidth=0.10, legend=False)
    _world_map(bot_ax)

    bot_anchor = bot_ax.get_position()
    fig.text(0.5, bot_anchor.y1 + 0.012,
             "(b) Combined map: mean across the 6 reanalysis pairs  "
             f"(max observed = {mx:.4f})",
             fontsize=10, weight="bold", ha="center")

    # ---- Shared colourbar at the bottom -------------------------------
    cax = fig.add_subplot(outer[2])
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm); sm.set_array([])
    cb = fig.colorbar(sm, cax=cax, orientation="horizontal")
    cb.set_label(
        "Absolute reduction in L1 distance when country dropped (L1 units, "
        "uniform weighting; positive only; grey = drives agreement or no data)",
        fontsize=8)
    cb.set_ticks([0, pos_sat])
    cb.set_ticklabels([
        "0",
        f"≥ {pos_sat:.4f}\n(strongest drivers; per-pair max = {gmax:.4f})",
    ])
    cb.ax.tick_params(labelsize=7)

    fig.suptitle(
        "Extended Data Fig. 2   Country drivers of cross-reanalysis "
        "disagreement in the fitted BHM dose-response (absolute units)",
        fontsize=10.5, x=0.5, y=0.992)

    return fig


# ===========================================================================
# Main
# ===========================================================================

def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print("Loading regression panel + balancing across 4 reanalyses …")
    panel = load_panel()
    n_cy = panel.groupby(["iso", "year"]).ngroups
    print(f"  balanced rows: {len(panel):,}  ({n_cy:,} country-years × 4)")

    print("Residualising temperature per reanalysis "
          "(iso FE + year FE + iso×t + iso×t²) …")
    resid = residualise_all(panel)

    print("Computing per-pair country MAD of residual T …")
    pair_mad = per_pair_country_mad(resid)

    print("Aggregating country covariates …")
    cov = country_covariates(panel)

    print("Loading regression coefficients (full + LOO) …")
    full = load_full_coeffs()
    loo  = load_loo_coeffs()
    print(f"  full: {len(full)} reanalyses;  LOO: {len(loo):,} (rea × country)")

    print("Computing per-pair absolute reduction in L1 (uniform weighting) …")
    perpair_abs = compute_perpair_abs_reduction(panel, full, loo)
    print("  baselines:")
    for (a, b), grp in (perpair_abs[perpair_abs["country"] == "all"]
                         .groupby(["pair_a", "pair_b"])):
        print(f"    {LABEL[a]:>6s} vs {LABEL[b]:<6s}  L1 = "
              f"{grp['L1'].iloc[0]:.5f}")
    mean_abs = compute_mean_abs_reduction(perpair_abs)

    print("Loading shapefile …")
    gdf = load_shapefile()

    print("Building Extended Data Fig. 1 …")
    fig1 = build_fig1(resid, pair_mad, cov)
    fig1.savefig(OUT_FIG1, format="pdf")
    fig1.savefig(OUT_FIG1.with_suffix(".png"), format="png", dpi=200)
    print(f"  saved {OUT_FIG1}")

    print("Building Extended Data Fig. 2 …")
    fig2 = build_fig2(perpair_abs, mean_abs, gdf)
    fig2.savefig(OUT_FIG2, format="pdf")
    fig2.savefig(OUT_FIG2.with_suffix(".png"), format="png", dpi=200)
    print(f"  saved {OUT_FIG2}")


if __name__ == "__main__":
    main()
