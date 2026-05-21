"""
Residualised temperature time series for selected countries.

Replicates Tom's R residualisation
(tb/code/explore/temp_compare.R lines 97–105):

  feols(temp_1 ~ 1 | iso + year + iso[time1] + iso[time2], data = d)

per reanalysis, where
  - iso         : country fixed effect
  - year        : year fixed effect
  - iso[time1]  : country-specific LINEAR time trend
  - iso[time2]  : country-specific QUADRATIC time trend (time2 = time^2)
  - time        : year - 1960 (pre-computed in the reg-input CSV)

`fixest::feols` partials out all these absorbed effects and returns
residuals. We do the same in Python with statsmodels OLS on the explicit
dummy + interaction design matrix (~200 countries × ~30 years per
reanalysis is small enough for dense least squares).

The residuals isolate the "weather noise" the BHM regression actually sees
once climate + trends are removed. They are the right object for showing
where reanalyses disagree on *weather*, not on *climate*.

Inputs:
  - reg_inputs/gdppc/reginputs_gdppc_025deg.csv
  - dmes/mean_dmes_gdppc.csv  (to pick top-3 / bottom-3 countries)

Outputs:
  - tb/figs/country_temp_raw.pdf       — raw T time series, 6 countries
  - tb/figs/country_temp_residual.pdf  — residualised T, same 6 countries
"""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy import sparse
from scipy.sparse.linalg import lsqr


DROPBOX = Path(os.environ.get(
    "KS_DROPBOX", "/Users/tombearpark/Library/CloudStorage/Dropbox/ks"))
REPO = Path(os.environ.get(
    "KS_REPO", "/Users/tombearpark/Documents/GitHub/impact_uncertainty"))

REG_INPUT = DROPBOX / "reg_inputs" / "gdppc" / "reginputs_gdppc_025deg.csv"
DMES_CSV  = DROPBOX / "dmes" / "mean_dmes_gdppc.csv"
FIGS      = REPO / "tb" / "figs"

REANALYSES = ("ERA5-025", "MERRA2", "JRA-3Q", "GMFD")
LABEL  = {"ERA5-025": "ERA5", "MERRA2": "MERRA2",
          "JRA-3Q": "JRA-3Q", "GMFD": "GMFD"}
COLOUR = {"ERA5-025": "tab:red", "MERRA2": "tab:green",
          "JRA-3Q": "tab:orange", "GMFD": "tab:blue"}

YEAR_LO, YEAR_HI = 1980, 2010
N_PER_SIDE = 3


# ---------------------------------------------------------------------------
# Residualisation
# ---------------------------------------------------------------------------

def residualise_per_reanalysis(df_rea: pd.DataFrame) -> np.ndarray:
    """Project temp_popweight onto [iso FE, year FE, iso×time, iso×time²] and
    return the residuals.

    df_rea has one row per (iso, year). Returns an array aligned with df_rea.
    """
    iso = pd.Categorical(df_rea["iso"])
    yr  = pd.Categorical(df_rea["year"])
    t1  = df_rea["time"].to_numpy(dtype=float)
    t2  = df_rea["time2"].to_numpy(dtype=float)
    iso_codes = iso.codes
    yr_codes  = yr.codes

    n = len(df_rea)
    n_iso, n_yr = len(iso.categories), len(yr.categories)

    # Build the design matrix as 4 sparse blocks:
    #   ISO dummies          (n × n_iso)
    #   YEAR dummies         (n × n_yr)        (drop one for identification)
    #   ISO × time1 dummies  (n × n_iso)
    #   ISO × time2 dummies  (n × n_iso)
    rows = np.arange(n)
    iso_blk = sparse.csc_matrix(
        (np.ones(n), (rows, iso_codes)), shape=(n, n_iso))
    # Drop first year to identify (iso FE absorbs the overall intercept)
    yr_blk = sparse.csc_matrix(
        (np.ones(n), (rows, yr_codes)), shape=(n, n_yr))
    yr_blk = yr_blk[:, 1:]
    iso_t1 = sparse.csc_matrix(
        (t1, (rows, iso_codes)), shape=(n, n_iso))
    iso_t2 = sparse.csc_matrix(
        (t2, (rows, iso_codes)), shape=(n, n_iso))
    X = sparse.hstack([iso_blk, yr_blk, iso_t1, iso_t2]).tocsr()
    y = df_rea["temp_popweight"].to_numpy(dtype=float)

    # Sparse least squares
    sol = lsqr(X, y, atol=1e-10, btol=1e-10, iter_lim=20000)
    beta = sol[0]
    fitted = X @ beta
    return y - fitted


def residualise_all(df: pd.DataFrame) -> pd.DataFrame:
    """For each reanalysis separately, residualise temp_popweight.

    Returns df with new column `resid`."""
    parts = []
    for rea in REANALYSES:
        sub = df[df["obs_mod"] == rea].copy()
        if sub.empty:
            continue
        # Tom's R code subsets to rows non-missing in temp_1 — we do the same.
        sub = sub.dropna(subset=["temp_popweight", "time", "time2", "iso"])
        resid = residualise_per_reanalysis(sub)
        sub["resid"] = resid
        # Diagnostic: |residual| should be small relative to raw T
        rms = float(np.sqrt(np.mean(resid ** 2)))
        print(f"    {rea:10s}  N={len(sub):>5d}   RMS(resid) = {rms:.3f} °C   "
              f"max|r| = {np.max(np.abs(resid)):.2f} °C")
        parts.append(sub[["countryname", "year", "obs_mod", "temp_popweight", "resid"]])
    return pd.concat(parts, ignore_index=True)


# ---------------------------------------------------------------------------
# Selection of panel countries (top-3 + bottom-3 by all-pairs pct_reduction)
# ---------------------------------------------------------------------------

def pick_panel_countries() -> tuple[list[str], dict[str, float]]:
    raw = pd.read_csv(DMES_CSV).rename(columns={"loo_id": "country"})
    base = float(raw.loc[raw["country"] == "all", "dme_mu"].iloc[0])
    pct = raw[raw["country"] != "all"].copy()
    pct["pct_reduction"] = (base - pct["dme_mu"]) / base * 100.0
    pct = pct.sort_values("pct_reduction", ascending=False)
    top = pct.head(N_PER_SIDE)
    bot = pct.tail(N_PER_SIDE).iloc[::-1]
    chosen = list(top["country"]) + list(bot["country"])
    pct_map = {c: float(v) for c, v in zip(pct["country"], pct["pct_reduction"])}
    return chosen, {c: pct_map[c] for c in chosen}


# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------

def plot_panel(df: pd.DataFrame, var: str, pct: dict[str, float],
               title: str, ylabel: str, out: Path) -> None:
    countries = list(pct.keys())
    fig, axes = plt.subplots(3, 2, figsize=(10.5, 7.2), sharex=True)
    for k, c in enumerate(countries):
        r, col = divmod(k, 2)
        ax = axes[r, col]
        d_c = df[df["countryname"] == c]
        for rea in REANALYSES:
            d = d_c[d_c["obs_mod"] == rea].sort_values("year")
            if d.empty:
                continue
            ax.plot(d["year"], d[var], color=COLOUR[rea], linewidth=1.0)
        ax.axhline(0, color="0.6", linewidth=0.4)
        ax.tick_params(labelsize=8)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        # Sign-coloured title
        v = pct[c]
        col_title = "#762A83" if v > 0 else "#1B7837"
        ax.set_title(f"{c}    {v:+.1f}%",
                     fontsize=10, loc="left", pad=2, color=col_title)
    fig.supxlabel("Year", fontsize=10)
    fig.supylabel(ylabel, fontsize=10)
    fig.suptitle(title, fontsize=12)
    legend = [Line2D([0], [0], color=COLOUR[r], linewidth=1.6, label=LABEL[r])
              for r in REANALYSES]
    fig.legend(handles=legend, ncol=4, loc="lower center",
               bbox_to_anchor=(0.5, -0.02), fontsize=9, frameon=False)
    fig.tight_layout(rect=(0, 0.03, 1, 0.96))
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, format="pdf", bbox_inches="tight")
    fig.savefig(out.with_suffix(".png"), format="png", dpi=180,
                bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print("Picking panel countries (top-3 + bottom-3 by all-pairs pct_reduction) …")
    countries, pct = pick_panel_countries()
    for c, v in pct.items():
        print(f"  {c:25s}  {v:+.2f}%")

    print("Loading regression input panel …")
    use = ["countryname", "iso", "year", "obs_mod",
           "temp_popweight", "time", "time2", "growthWDI"]
    df = pd.read_csv(REG_INPUT, usecols=use)
    df = df[(df["year"] >= YEAR_LO) & (df["year"] <= YEAR_HI)]
    df = df[df["obs_mod"].isin(REANALYSES)]
    # The BHM regression sample drops rows with missing growthWDI; the
    # residualisation should be done on the same sample so the FE estimates
    # match the underlying regression.
    df = df.dropna(subset=["temp_popweight", "growthWDI"])
    print(f"  panel: {len(df):,} rows, {df['iso'].nunique()} countries, "
          f"{df['year'].nunique()} years, {df['obs_mod'].nunique()} reanalyses")

    print("Residualising per reanalysis  (temp_popweight ~ iso + year + iso*time + iso*time²) …")
    df_resid = residualise_all(df)

    # --- Raw figure
    plot_panel(
        df_resid, "temp_popweight", pct,
        title="Annual population-weighted mean temperature  (raw)",
        ylabel="T (°C)",
        out=FIGS / "country_temp_raw.pdf")
    print(f"Saved: {FIGS / 'country_temp_raw.pdf'}")

    # --- Residual figure
    plot_panel(
        df_resid, "resid", pct,
        title="Residualised temperature  (after country FE + year FE + "
              "country-specific linear + quadratic time trends)",
        ylabel="T residual (°C)",
        out=FIGS / "country_temp_residual.pdf")
    print(f"Saved: {FIGS / 'country_temp_residual.pdf'}")


if __name__ == "__main__":
    main()
