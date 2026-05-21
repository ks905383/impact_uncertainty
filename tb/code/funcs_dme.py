"""
Shared helpers: load LOO/full coefficients + reg-input panel, recompute L1
distances under configurable weighting, and return tables in the same shape
as Kevin's CSVs so existing figure scripts can swap data source by a single
arg.

Two weightings supported:

  - "pop"     : pair-specific empirical T distribution pooled from the two
                reanalyses in the pair, weighted by population (same recipe
                as Kevin's code/3b_dmgfdiffs_sources.ipynb).
  - "uniform" : same pair-pooled T distribution, but each (country, year,
                reanalysis) observation contributes equally (np.mean).

Helpers:

  - load_full_coeffs(): {rea: (b1, b2)}                              # 4 entries
  - load_loo_coeffs():  {(rea, country): (b1, b2)}                   # ~780
  - load_reg_panel():   DataFrame (countryname, year, obs_mod, temp_popweight, Pop)
  - compute_pair_dme_table(weighting):   columns mirror dmes_gdppc.csv
  - compute_mean_dme_table(weighting):   columns mirror mean_dmes_gdppc.csv

Calibration: run this file directly (`python funcs_dme.py`) to print a
side-by-side check of the "pop" recompute vs Kevin's CSV.
"""

from __future__ import annotations

import os
import re
from pathlib import Path
from typing import Literal

import numpy as np
import pandas as pd


DROPBOX = Path(os.environ.get(
    "KS_DROPBOX", "/Users/tombearpark/Library/CloudStorage/Dropbox/ks"))

FULL_DIR  = DROPBOX / "reg_outputs" / "gdppc"
LOO_DIR   = DROPBOX / "reg_outputs" / "gdppc" / "loos"
REG_INPUT = DROPBOX / "reg_inputs" / "gdppc" / "reginputs_gdppc_025deg.csv"

# Reanalyses in the canonical 4-product set. Kevin's dmes_gdppc.csv uses the
# label "ERA5" for ERA5-025; we keep both. PAIR_KEY normalises any pair to a
# canonical ordered tuple so de-duplicating across files is safe.
REANALYSES = ("ERA5-025", "MERRA2", "JRA-3Q", "GMFD")
# Map between filename-side label (ERA5-025) and Kevin's dme-CSV label (ERA5)
REA_TO_DME = {"ERA5-025": "ERA5", "MERRA2": "MERRA2",
              "JRA-3Q": "JRA-3Q", "GMFD": "GMFD"}
DME_TO_REA = {v: k for k, v in REA_TO_DME.items()}

YEAR_LO, YEAR_HI = 1980, 2010

# All 6 unordered pairs (in filename-naming convention). Used to enumerate.
PAIRS = [(a, b) for i, a in enumerate(REANALYSES) for b in REANALYSES[i+1:]]


# ---------------------------------------------------------------------------
# Loaders
# ---------------------------------------------------------------------------

_FULL_RE = re.compile(
    r"^regoutputs_coeffs_gdppc_025deg_(?P<rea>[^_]+)_1980-2010\.csv$")
_LOO_RE  = re.compile(
    r"^regoutputs_coeffs_gdppc_025deg_(?P<rea>[^_]+)_1980-2010_no(?P<country>.+)\.csv$")


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
    """Returns {(rea, country_name): (b1, b2)} for all 4 main reanalyses."""
    out: dict[tuple[str, str], tuple[float, float]] = {}
    for fn in sorted(LOO_DIR.glob("regoutputs_coeffs_gdppc_*.csv")):
        m = _LOO_RE.match(fn.name)
        if not m or m.group("rea") not in REANALYSES:
            continue
        d = pd.read_csv(fn)
        out[(m.group("rea"), m.group("country"))] = (
            float(d["b1"].iloc[0]), float(d["b2"].iloc[0]))
    return out


def load_reg_panel() -> pd.DataFrame:
    """Country × year × reanalysis panel filtered to BHM sample years and the
    four main reanalyses; rows with missing growthWDI dropped (same screen as
    the regression)."""
    use = ["countryname", "year", "obs_mod", "temp_popweight", "Pop", "growthWDI"]
    df = pd.read_csv(REG_INPUT, usecols=use)
    df = df[(df["year"] >= YEAR_LO) & (df["year"] <= YEAR_HI)]
    df = df[df["obs_mod"].isin(REANALYSES)]
    df = df.dropna(subset=["temp_popweight", "growthWDI", "Pop"])
    return df.reset_index(drop=True)


# ---------------------------------------------------------------------------
# L1 distance under a chosen weighting
# ---------------------------------------------------------------------------

WeightKind = Literal["pop", "uniform"]


def _me(b1: float, b2: float, t: np.ndarray) -> np.ndarray:
    return b1 + 2.0 * b2 * t


def pair_l1(coef_a: tuple[float, float], coef_b: tuple[float, float],
            sub: pd.DataFrame, pair: tuple[str, str],
            weighting: WeightKind) -> float:
    """L1 distance between MEs of `pair[0]` and `pair[1]`, integrated over the
    empirical T distribution from rows of `sub` belonging to either reanalysis.

    `sub` is the (already-filtered) reg_input panel; we further restrict to the
    two reanalyses in `pair` and then evaluate |ME_a(T_i) - ME_b(T_i)| at each
    surviving row, weighted by `Pop` if weighting == "pop", else uniformly.
    """
    a, b = pair
    s = sub[sub["obs_mod"].isin((a, b))]
    if s.empty:
        return float("nan")
    t = s["temp_popweight"].to_numpy()
    diff = np.abs(_me(*coef_a, t) - _me(*coef_b, t))
    if weighting == "pop":
        w = s["Pop"].to_numpy()
        return float(np.sum(w * diff) / np.sum(w))
    return float(np.mean(diff))


# ---------------------------------------------------------------------------
# Tables (in the shape of Kevin's CSVs)
# ---------------------------------------------------------------------------

def compute_pair_dme_table(weighting: WeightKind = "pop") -> pd.DataFrame:
    """Per (pair, loo_id) L1 distance + pct_reduction (vs the pair's own
    baseline). Returns columns: loo_id, proj_base_a, proj_base_b, dme,
    pct_reduction.

    `loo_id == 'all'` row is the full-sample baseline for that pair (pct=0).
    """
    full = load_full_coeffs()
    loo  = load_loo_coeffs()
    panel = load_reg_panel()
    countries = sorted(panel["countryname"].unique())

    rows: list[dict] = []
    for (a, b) in PAIRS:
        # Kevin's CSV convention puts the alphabetically-greater dme-label
        # first in (proj_base_a, proj_base_b); match that here so consumers
        # can swap "kevin" / "uniform" interchangeably without re-mapping.
        a_lbl, b_lbl = REA_TO_DME[a], REA_TO_DME[b]
        if a_lbl < b_lbl:
            pba, pbb = b_lbl, a_lbl
        else:
            pba, pbb = a_lbl, b_lbl
        # Baseline (all countries in)
        L1_all = pair_l1(full[a], full[b], panel, (a, b), weighting)
        rows.append({
            "loo_id": "all",
            "proj_base_a": pba, "proj_base_b": pbb,
            "dme": L1_all, "pct_reduction": 0.0,
        })
        # Per-country LOO
        for c in countries:
            key_a = (a, c); key_b = (b, c)
            if key_a not in loo or key_b not in loo:
                continue  # this country wasn't LOO-fit for one of the rea
            sub_noc = panel[panel["countryname"] != c]
            L1_c = pair_l1(loo[key_a], loo[key_b], sub_noc, (a, b), weighting)
            rows.append({
                "loo_id": c,
                "proj_base_a": pba, "proj_base_b": pbb,
                "dme": L1_c,
                "pct_reduction": 100.0 * (L1_all - L1_c) / L1_all,
            })
    return pd.DataFrame(rows)


def compute_mean_dme_table(weighting: WeightKind = "pop",
                            pair_df: pd.DataFrame | None = None) -> pd.DataFrame:
    """Per loo_id, mean dme across all 6 pairs. Returns columns: loo_id,
    dme_mu, pct_reduction (vs the 'all' row's dme_mu).

    Pass a pre-computed `pair_df` to avoid recomputing.
    """
    if pair_df is None:
        pair_df = compute_pair_dme_table(weighting=weighting)
    mu = pair_df.groupby("loo_id", as_index=False).agg(dme_mu=("dme", "mean"))
    base = float(mu.loc[mu["loo_id"] == "all", "dme_mu"].iloc[0])
    mu["pct_reduction"] = 100.0 * (base - mu["dme_mu"]) / base
    return mu


# ---------------------------------------------------------------------------
# Source switch: read Kevin's CSV when weighting='pop' is wanted AND a
# pre-computed CSV exists (matches the original headline number exactly).
# Useful so the main map can still cite Kevin's value verbatim.
# ---------------------------------------------------------------------------

KEVIN_MEAN_CSV = DROPBOX / "dmes" / "mean_dmes_gdppc.csv"
KEVIN_PAIR_CSV = DROPBOX / "dmes" / "dmes_gdppc.csv"


def kevin_mean_dme_table() -> pd.DataFrame:
    """Read Kevin's CSV and add pct_reduction in the same shape as
    compute_mean_dme_table()."""
    raw = pd.read_csv(KEVIN_MEAN_CSV)
    base = float(raw.loc[raw["loo_id"] == "all", "dme_mu"].iloc[0])
    raw["pct_reduction"] = 100.0 * (base - raw["dme_mu"]) / base
    return raw


def kevin_pair_dme_table() -> pd.DataFrame:
    """Read Kevin's per-pair CSV and add pct_reduction per pair."""
    raw = pd.read_csv(KEVIN_PAIR_CSV)
    out_rows = []
    for (a, b), grp in raw.groupby(["proj_base_a", "proj_base_b"]):
        base = float(grp.loc[grp["loo_id"] == "all", "dme"].iloc[0])
        grp = grp.copy()
        grp["pct_reduction"] = 100.0 * (base - grp["dme"]) / base
        out_rows.append(grp)
    return pd.concat(out_rows, ignore_index=True)


# ---------------------------------------------------------------------------
# Calibration when run as a script
# ---------------------------------------------------------------------------

def _calibration() -> None:
    """Confirm the "pop" recompute lines up with Kevin's CSV. Different
    T-distribution choices (Kevin uses xarray bcd dataset, we use the BHM
    reg-input panel) mean values will differ slightly; we expect the
    *country ranking* to match closely (Spearman > 0.95) and the absolute
    values to be in the same ballpark."""
    print("Computing pop-weighted dme table from BHM reg-input panel …")
    my_pair = compute_pair_dme_table(weighting="pop")
    my_mean = compute_mean_dme_table(weighting="pop", pair_df=my_pair)
    kevin   = kevin_mean_dme_table()

    join = my_mean.merge(kevin, on="loo_id", how="inner",
                         suffixes=("_mine", "_kevin"))
    print(f"\nMatched {len(join)} loo_ids vs Kevin's CSV")
    print(f"  Spearman rank corr (dme_mu, mine vs kevin): "
          f"{join['dme_mu_mine'].rank().corr(join['dme_mu_kevin'].rank(), method='spearman'):.4f}")
    print(f"  Pearson      corr (dme_mu, mine vs kevin): "
          f"{join['dme_mu_mine'].corr(join['dme_mu_kevin']):.4f}")
    print(f"  Mean ratio  (mine / kevin): "
          f"{(join['dme_mu_mine'] / join['dme_mu_kevin']).mean():.4f}")

    # Show top/bottom for sanity
    join_c = join[join["loo_id"] != "all"].copy()
    join_c["pct_mine"]  = join_c["pct_reduction_mine"]
    join_c["pct_kevin"] = join_c["pct_reduction_kevin"]
    print("\nTop-5 by Kevin's pct_reduction:")
    print(join_c.nlargest(5, "pct_kevin")[
        ["loo_id", "pct_kevin", "pct_mine"]].to_string(index=False))
    print("\nBottom-5 by Kevin's pct_reduction:")
    print(join_c.nsmallest(5, "pct_kevin")[
        ["loo_id", "pct_kevin", "pct_mine"]].to_string(index=False))

    print("\n--- Compared to uniform weighting (recompute) ---")
    my_pair_u = compute_pair_dme_table(weighting="uniform")
    my_mean_u = compute_mean_dme_table(weighting="uniform", pair_df=my_pair_u)
    join_u = my_mean_u.merge(kevin, on="loo_id", how="inner",
                              suffixes=("_uniform", "_kevin"))
    join_u_c = join_u[join_u["loo_id"] != "all"]
    print(f"Spearman rank corr (uniform vs Kevin pop): "
          f"{join_u_c['dme_mu_uniform'].rank().corr(join_u_c['dme_mu_kevin'].rank(), method='spearman'):.4f}")
    print("Top-5 by my uniform pct_reduction:")
    print(my_mean_u[my_mean_u['loo_id'] != 'all']
          .nlargest(5, 'pct_reduction')[['loo_id','pct_reduction']]
          .to_string(index=False))
    print("Bottom-5 by my uniform pct_reduction:")
    print(my_mean_u[my_mean_u['loo_id'] != 'all']
          .nsmallest(5, 'pct_reduction')[['loo_id','pct_reduction']]
          .to_string(index=False))


if __name__ == "__main__":
    _calibration()
