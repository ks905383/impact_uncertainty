"""
Per-pair LOO maps + histograms for the BHM (gdp-pc) dose-response.

For each of the six reanalysis pairs, plot a world map and a histogram
underneath showing the % reduction in that pair's marginal-effect distance
when each country is dropped from the regression sample. Same colour scale
across the six maps for direct comparison.

Data source: Kevin's `dmes_gdppc.csv` (per-pair distances, with the special
loo_id == 'all' row giving the full-sample baseline per pair).

Output: tb/figs/loo_perpair_maps.pdf (+ .png).
"""

from __future__ import annotations

import os
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.colors import TwoSlopeNorm, LinearSegmentedColormap, Normalize


# High-contrast purple ramp: truncated `Purples` skipping the lightest 35%,
# so even small positive values render as visible lavender rather than white.
_purples = plt.get_cmap("Purples")
PURPLES_HC = LinearSegmentedColormap.from_list(
    "PurplesHC", _purples(np.linspace(0.35, 1.0, 256)))


DROPBOX = Path(os.environ.get(
    "KS_DROPBOX", "/Users/tombearpark/Library/CloudStorage/Dropbox/ks"))
REPO = Path(os.environ.get(
    "KS_REPO", "/Users/tombearpark/Documents/GitHub/impact_uncertainty"))

DMES_CSV = DROPBOX / "dmes" / "dmes_gdppc.csv"
SHP = REPO / "aux_data" / "geo_data" / "country-fixed.shp"
OUT_PDF = REPO / "tb" / "figs" / "loo_perpair_maps.pdf"

# 4 main reanalyses → 6 unordered pairs. We display in this canonical order
# (smaller-distance pairs first so the eye sweeps through the disagreement
# spectrum left-to-right, top-to-bottom).
# Note: Kevin's dmes_gdppc.csv uses "ERA5" (not "ERA5-025"). Reg-input /
# reg-output files use "ERA5-025". Keep the dmes-side naming here.
PAIR_ORDER = [
    ("MERRA2", "GMFD"),
    ("JRA-3Q", "ERA5"),
    ("MERRA2", "ERA5"),
    ("GMFD", "ERA5"),
    ("MERRA2", "JRA-3Q"),
    ("JRA-3Q", "GMFD"),
]
LABEL = {"ERA5": "ERA5", "MERRA2": "MERRA2",
         "JRA-3Q": "JRA-3Q", "GMFD": "GMFD"}


def load_per_pair(source: str = "kevin") -> tuple[pd.DataFrame, dict[tuple[str, str], float]]:
    """Return (long DataFrame keyed by [pair, country, pct], baseline-per-pair).

    source: "kevin"   → read Kevin's pop-weighted dmes_gdppc.csv as before
            "uniform" → recompute via funcs_dme with uniform weighting
    """
    if source == "uniform":
        import sys
        from pathlib import Path as _P
        sys.path.insert(0, str(_P(__file__).resolve().parent))
        import funcs_dme  # noqa: E402
        raw = funcs_dme.compute_pair_dme_table(weighting="uniform")
    else:
        raw = pd.read_csv(DMES_CSV)
    # Canonical pair key: tuple(sorted) — but we keep the original a/b columns
    # to recover Kevin's ordering for the title strings.
    out_rows = []
    baseline = {}
    for (a, b), grp in raw.groupby(["proj_base_a", "proj_base_b"]):
        baseline_row = grp[grp["loo_id"] == "all"]
        if len(baseline_row) != 1:
            raise RuntimeError(f"Expected 1 'all' row for pair {a}-{b}, got {len(baseline_row)}")
        b_dme = float(baseline_row["dme"].iloc[0])
        per_country = grp[grp["loo_id"] != "all"].copy()
        per_country["pct_reduction"] = (b_dme - per_country["dme"]) / b_dme * 100.0
        per_country["abs_reduction"] = (b_dme - per_country["dme"])  # same units as L1
        per_country["pair"] = f"{LABEL[a]} vs {LABEL[b]}"
        per_country["_a"] = a
        per_country["_b"] = b
        out_rows.append(per_country.rename(columns={"loo_id": "country"}))
        baseline[(a, b)] = b_dme
    out = pd.concat(out_rows, ignore_index=True)
    return out, baseline


def pair_label(a: str, b: str) -> str:
    return f"{LABEL[a]} vs {LABEL[b]}"


def build_figure(perpair: pd.DataFrame,
                 baseline: dict[tuple[str, str], float],
                 gdf: gpd.GeoDataFrame,
                 *,
                 positive_only: bool = False,
                 weighting_tag: str = "",
                 metric: str = "pct") -> plt.Figure:
    """Per-pair maps + histograms.

    When `positive_only=True`, paint only countries with reduction > 0 using
    a sequential purple ramp; countries with reduction <= 0 are drawn in the
    same light grey as countries with no data. The histogram underneath
    still shows the full distribution so the reader can see what has been
    filtered out.

    `metric`:
      - "pct" (default): plot pct_reduction (relative to each pair's
        baseline). Per-pair colour scale (baselines differ ~4.6×).
      - "abs": plot abs_reduction = baseline − dme (same L1 units across all
        pairs). Shared colour scale so the eye can compare which pairs
        produce the biggest *absolute* movement.
    """
    if metric not in ("pct", "abs"):
        raise ValueError(metric)
    val_col = "pct_reduction" if metric == "pct" else "abs_reduction"
    unit_lbl = "% reduction in pair's L1 when country dropped" \
        if metric == "pct" else "absolute reduction in pair's L1  (L1 units)"
    fig = plt.figure(figsize=(15.5, 9.5))
    outer = GridSpec(2, 3, figure=fig,
                     wspace=0.12, hspace=0.38,
                     left=0.03, right=0.97, top=0.91, bottom=0.04)

    cmap = PURPLES_HC if positive_only else plt.get_cmap("PRGn_r")

    # For "abs" metric: shared colour scale + shared histogram x-range
    # across all 6 pairs (units are identical L1 units, so direct comparison
    # is meaningful). For "pct" metric: keep per-pair scaling as before.
    if metric == "abs":
        global_vals = perpair[val_col].to_numpy()
        gmin, gmax = float(global_vals.min()), float(global_vals.max())

    for i, (a, b) in enumerate(PAIR_ORDER):
        r, c = divmod(i, 3)
        # Two-row inner grid: map on top, histogram below. The colour
        # encoding is communicated through coloured histogram bars rather
        # than a separate colourbar (which was crowding the histogram).
        inner = outer[r, c].subgridspec(
            2, 1, height_ratios=[3.2, 1.0], hspace=0.18)
        ax_map = fig.add_subplot(inner[0])
        ax_h   = fig.add_subplot(inner[1])

        pct = perpair[(perpair["_a"] == a) & (perpair["_b"] == b)]
        vals = pct[val_col].to_numpy()
        if metric == "abs":
            vmin, vmax = gmin, gmax
        else:
            vmin, vmax = float(vals.min()), float(vals.max())
        eps = 0.1 if metric == "pct" else max(1e-6, 0.01 * (vmax - vmin))
        if positive_only:
            # Sequential purple ramp, 0 → max(positive). Negative countries
            # are folded into "no data" and drawn light grey.
            norm = Normalize(vmin=0, vmax=max(vmax, eps))
        else:
            # Diverging — make sure the scale spans 0 so colours stay diverging.
            norm = TwoSlopeNorm(vmin=min(vmin, -eps), vcenter=0,
                                vmax=max(vmax, eps))

        merged = gdf.merge(pct[["country", val_col]],
                           left_on="locid", right_on="country", how="left")
        if positive_only:
            plotted = merged[merged[val_col] > 0]
            unplotted = merged[~(merged[val_col] > 0)]
        else:
            plotted = merged[~merged[val_col].isna()]
            unplotted = merged[merged[val_col].isna()]

        unplotted.plot(ax=ax_map, color="lightgray", edgecolor="white",
                       linewidth=0.15)
        plotted.plot(ax=ax_map, column=val_col, cmap=cmap, norm=norm,
                     edgecolor="white", linewidth=0.15, legend=False)
        ax_map.set_xlim(-180, 180)
        ax_map.set_ylim(-60, 85)
        ax_map.set_aspect("equal", adjustable="box")
        ax_map.set_xticks([]); ax_map.set_yticks([])
        for sp in ax_map.spines.values():
            sp.set_visible(False)

        base = baseline[(a, b)]
        ax_map.set_title(
            f"{pair_label(a, b)}    (baseline L1 = {base:.4f})",
            fontsize=10, loc="left", pad=2)

        # Histogram doubles as the colour legend: each bar is filled with the
        # same cmap colour that the map uses for that value range. Always
        # shows the full distribution so the reader sees what's hidden in the
        # positive-only view. For "abs" metric: shared bin edges across all
        # pairs (gmin..gmax) so histograms are directly comparable.
        if metric == "abs":
            hist_lo, hist_hi = min(gmin, -eps), gmax
        elif positive_only:
            hist_lo, hist_hi = min(vals.min(), -0.1), vals.max()
        else:
            hist_lo, hist_hi = vmin, vmax
        full_bins = np.linspace(hist_lo, hist_hi, 26)
        counts, edges = np.histogram(vals, bins=full_bins)
        centers = 0.5 * (edges[:-1] + edges[1:])
        widths = np.diff(edges)
        if positive_only:
            bar_colors = [cmap(norm(c)) if c > 0 else "lightgray" for c in centers]
        else:
            bar_colors = [cmap(norm(c)) for c in centers]
        ax_h.bar(centers, counts, width=widths * 0.95, color=bar_colors,
                 edgecolor="white", linewidth=0.4, align="center")
        ax_h.set_xlim(full_bins[0], full_bins[-1])
        ax_h.axvline(0, color="black", linewidth=0.5)
        ax_h.tick_params(labelsize=7)
        ax_h.set_yticks([])
        for sp in ("top", "right", "left"):
            ax_h.spines[sp].set_visible(False)
        # Annotate the per-pair top-3 drivers (always, ranked by chosen metric)
        top3 = pct.nlargest(3, val_col)
        ax_h.text(0.01, 0.95, "drives disagreement: " +
                  ", ".join(top3["country"].tolist()),
                  transform=ax_h.transAxes, ha="left", va="top",
                  fontsize=6.5, color="#762A83")
        if not positive_only:
            bot3 = pct.nsmallest(3, val_col)
            ax_h.text(0.01, 0.72, "drives agreement: " +
                      ", ".join(bot3["country"].tolist()),
                      transform=ax_h.transAxes, ha="left", va="top",
                      fontsize=6.5, color="#1B7837")
        else:
            n_pos = int((vals > 0).sum()); n_total = int(vals.size)
            ax_h.text(0.01, 0.72,
                      f"{n_pos}/{n_total} countries shown on map",
                      transform=ax_h.transAxes, ha="left", va="top",
                      fontsize=6.5, color="0.35")
        ax_h.set_xlabel(unit_lbl, fontsize=7, labelpad=1)

    scale_note = ("each pair on its own colour scale (baselines differ ~4.6×)"
                  if metric == "pct"
                  else "shared colour scale across pairs (same L1 units)")
    if positive_only:
        title = ("Disagreement-driving countries only "
                 "(per-pair, BHM dose-response)")
        subtitle = ("purple = country drives disagreement   ·   "
                    "grey = country either drives agreement (negative) "
                    f"or has no data   ·   {scale_note}")
    else:
        title = ("Per-pair country influence on cross-reanalysis disagreement "
                 "(BHM dose-response)")
        subtitle = ("purple = country drives disagreement   ·   "
                    f"green = country drives agreement   ·   {scale_note}")
    if metric == "abs":
        title = f"{title} — absolute reduction"
    if weighting_tag:
        title = f"{title}  [{weighting_tag}]"
    fig.suptitle(title, fontsize=13, y=0.965, x=0.5)
    fig.text(0.5, 0.935, subtitle,
             fontsize=8.5, color="0.4", style="italic", ha="center")
    return fig


def main():
    print("Loading shapefile …")
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        gdf = gpd.read_file(SHP)
    gdf["geometry"] = gdf.make_valid()
    gdf = gdf.rename(columns={"CNTRY_NAME": "locid"})
    gdf = gdf.query("locid != 'Antarctica'")

    OUT_PDF.parent.mkdir(parents=True, exist_ok=True)

    # ----- pop-weighted (Kevin's CSV) + uniform-weighted (helper recompute) --
    for source, tag, suffix in [
            ("kevin",   "pop-weighted (Kevin's CSV)", ""),
            ("uniform", "uniform (each obs equal weight)", "_uniform")]:
        print(f"\n=== {tag} ===")
        print("Loading per-pair dme + baselines …")
        perpair, baseline = load_per_pair(source=source)
        print(f"  {len(perpair)} (pair, country) rows across {len(baseline)} pairs")
        for (a, b), v in sorted(baseline.items(), key=lambda kv: kv[1]):
            print(f"    {pair_label(a, b):20s} baseline = {v:.5f}")

        out_main     = OUT_PDF.with_name(OUT_PDF.stem + suffix + ".pdf")
        out_pos      = OUT_PDF.with_name(OUT_PDF.stem + "_positive" + suffix + ".pdf")
        out_main_abs = OUT_PDF.with_name(OUT_PDF.stem + "_abs" + suffix + ".pdf")
        out_pos_abs  = OUT_PDF.with_name(OUT_PDF.stem + "_positive_abs" + suffix + ".pdf")

        for fname, kwargs in [
            (out_main,     dict(positive_only=False, metric="pct")),
            (out_pos,      dict(positive_only=True,  metric="pct")),
            (out_main_abs, dict(positive_only=False, metric="abs")),
            (out_pos_abs,  dict(positive_only=True,  metric="abs")),
        ]:
            print(f"Building {fname.stem} …")
            fig = build_figure(perpair, baseline, gdf,
                               weighting_tag=tag, **kwargs)
            fig.savefig(fname, format="pdf", bbox_inches="tight")
            fig.savefig(fname.with_suffix(".png"), format="png", dpi=180,
                        bbox_inches="tight")
            plt.close(fig)
            print(f"  saved {fname}")


if __name__ == "__main__":
    main()
