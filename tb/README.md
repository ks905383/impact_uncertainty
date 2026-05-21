# Reanalysis comparison for climate-damage estimates

Scripts, figures, and slide decks investigating how the choice of gridded
temperature product affects damage-function regressions, and which
country-years drive the cross-reanalysis disagreement.

The project has two pipelines that share a common dataset family:

1. **R pipeline** (older, decks 01–05) — five reanalyses
   (ERA5-025, GMFD, JRA-3Q, MERRA2, UDel) on a **1° grid**; defines an L1
   distance between marginal-effect curves, attributes the gap, and
   proposes synthesis methods to reconcile the five products.
2. **Python LOO pipeline** (newer, deck 06 + Extended Data) — four
   reanalyses (drops UDel) on a **0.25° grid**; uses leave-one-out
   country attribution to identify which countries drive the cross-
   reanalysis disagreement in the BHM dose-response, partitions the
   disagreement per reanalysis pair, and contrasts pop-weighted vs
   uniform integration.

The Extended Data figures for the paper live in the sibling folder
[`../extended_data/`](../extended_data/) and reuse the Python residualiser
[`code/country_residuals.py`](code/country_residuals.py) from here.

## Folder layout

```
impact_uncertainty/tb/
├── code/
│   ├── *.py                # Python LOO pipeline (deck 06 + Extended Data)
│   ├── disagree_helpers.R  # shared R utilities
│   ├── gdp_disagree.R      # R: disagreement scatters for gdp-pc      (deck 02)
│   ├── mortality_disagree.R# R: same for the mortality panel
│   └── explore/            # older R analyses that produced decks 01–05
├── figs/                   # PDF figures written by all scripts (R + Python)
├── slides/                 # LaTeX source + compiled PDF for 6 decks
└── README.md
```

## R pipeline (decks 01–05)

Run from the Dropbox root `/.../ks` (the R scripts `setwd()` there).
Reads `reg_inputs/gdppc/reginputs_gdppc_1deg.csv` (1°, balanced across
5 reanalyses, 4,592 country-years × 5 reanalyses, 1980–2013). All output
PDFs land in `tb/figs/`.

The original *exploratory* analyses that produced the figures for decks
01, 03, 04, 05 live in `tb/code/explore/`. The two scripts that survived
into the main code/ directory (`gdp_disagree.R`, `mortality_disagree.R`)
were promoted because they are still being consumed and refined.

| Script | What it does | Figures |
|---|---|---|
| `explore/distance.R` | Fits the BHM panel per reanalysis, computes marginal effects with delta-method SEs, defines an L1 distance between ME curves integrated over the empirical T sample, and runs LOO attribution by continent / decade / temp decile. | (interactive) |
| `explore/distance_figs.R` | Sources `distance.R` and writes its plots: ME curves, pairwise distance matrix, LOO bar charts, ERA5↔MERRA2 scatter and bias-by-level, Sierra Leone temperature series. | `me_curves.pdf`, `distance_matrix.pdf`, `loo_{continent,decade,decile}.pdf`, `temp_{scatter,bias}.pdf`, `sle_temp.pdf` |
| `explore/oos_score.R` | Frisch–Waugh–Lovell residualisation per reanalysis, then cross-validated MSE and OOS R² under three CV splits (random 10-fold, leave-country-out, leave-year-out). | `oos_mse.pdf`, `oos_r2.pdf` |
| `explore/temp_compare.R` | Diagnostic comparison of the five temperature products themselves: pairwise correlations and MAD (raw and residualised), country-level disagreement, temporal evolution with a GMFD 2009–2010 drill-down, top outliers. | `tc_*.pdf` |
| `explore/synth_compare.R` | Six synthetic temperature measures (MEAN, MEDIAN, PC1, FACTOR, MIN-VAR, ANCHORED) compared against the five raw datasets on ME precision (β₂, CI width, t-stat), OOS predictive fit, and jackknife stability. Includes EM factor model with per-country loadings. | `sc_*.pdf` |
| `explore/hetero_factor.R` | Heteroskedastic factor model. Five-step: (1) σ_k depends on absolute T (climate axis) but not within-country deviation; (2) σ_k(d) per absolute-T decile; (3) T-dependent MIN-VAR-T synthesis; (4) classical EIV attenuation correction; (5) summary. | `hf_*.pdf` |
| `disagree_helpers.R` | Shared utilities for the two `*_disagree.R` scripts: residualisation, pairwise MAD aggregation, OLS-overlay scatter plotting, heatmap of interaction-model predictions. | — |
| `gdp_disagree.R` | Residualised disagreement scatters for the gdp-pc panel at country / year / country-year granularity; vs log GDP and vs mean T. Produces the deck-02 figures. | `tc_disagr_*.pdf` |
| `mortality_disagree.R` | Analogous to `gdp_disagree.R` but for the Carleton ADM1 mortality panel. | `mt_disagr_*.pdf` |

## Python LOO pipeline (deck 06 + Extended Data)

Run from the repo root: `python3 tb/code/<script>.py`. Reads data from
`$KS_DROPBOX` (env var; default
`/Users/tombearpark/Library/CloudStorage/Dropbox/ks`), specifically:

- `reg_inputs/gdppc/reginputs_gdppc_025deg.csv` — **4-reanalysis** balanced panel at 0.25°
- `reg_outputs/gdppc/regoutputs_coeffs_gdppc_025deg_<rea>_1980-2010.csv` — full-sample BHM coefficients (one per reanalysis)
- `reg_outputs/gdppc/loos/regoutputs_coeffs_gdppc_025deg_<rea>_1980-2010_no<country>.csv` — LOO coefficients per (reanalysis × dropped country); unpack `loos.zip` once if needed
- `dmes/mean_dmes_gdppc.csv` and `dmes/dmes_gdppc.csv` — Kevin's pop-weighted L1 distances (used as canonical "pop" source; the helper recomputes when `weighting="uniform"`)

Python deps: `numpy`, `pandas`, `scipy`, `matplotlib`, `geopandas`.

| Script | What it does | Figures |
|---|---|---|
| `funcs_dme.py` | Shared helper: loaders for coefficients (full + LOO), `pair_l1(..., weighting=)` (pop / uniform), `compute_pair_dme_table(weighting)` and `compute_mean_dme_table(weighting)` returning tables in the same shape as Kevin's CSVs. Run as a script for a calibration check against Kevin's pop-weighted values. | — |
| `country_residuals.py` | Per-reanalysis residualisation of `temp_popweight` against the BHM FE structure (`T ~ 1 | iso + year + iso[t] + iso[t²]`) via sparse OLS. Also produces raw + residual time-series figures for the top-3 / bottom-3 countries by `pct_reduction`. **Re-used by `extended_data/make_extended_data.py`.** | `country_temp_{raw,residual}.pdf` |
| `loo_country_influence.py` | The main "% reduction in cross-reanalysis disagreement when country dropped" multi-panel figure: world map, histogram, and country small-multiples of T time series. Emits both pop-weighted (Kevin's CSV) and uniform (helper-recomputed) variants. Mirror-saves to `tb/figs/` and `figures/` so the deck can find them. | `loo_country_influence{,_uniform}.{pdf,png}` |
| `loo_perpair_maps.py` | Per-pair (6 panels) LOO maps with histograms below each, across (diverging / positive-only) × (pop / uniform) × (% reduction / absolute L1) combinations → 8 variants per run. Histograms double as the colour-bar legend (bars coloured by value). | `loo_perpair_maps[_positive][_abs][_uniform].{pdf,png}` |
| `loo_scatter.py` | Country influence vs mean T and vs log mean GDP. Two layouts (averaged across pairs, and 2×3 per-pair facets), in pop / uniform variants → 6 figures. OLS overlay with R² + p annotations on each panel. | `loo_scatter[_perpair_{T,income}][_uniform].{pdf,png}` |
| `loo_residual_disagree.py` | Per-pair residual-T disagreement vs mean T and log GDP — 6 pairs × 2 covariates = 12 facets. Uses `country_residuals.residualise_all`. | `tc_disagr_perpair.{pdf,png}` |
| `bahrain_walkthrough.py` | Two-facet (full sample / sample with country dropped) overlay of ME curves across 4 reanalyses, for Bahrain (the top driver in the all-pairs metric) and Greenland (the top consensus-glue country). Plus the 3-panel L1-distance walkthrough on the GMFD-vs-ERA5 pair (MEs / |diff| / pop-weighted T density × |diff|), shown for the full sample and after Bahrain is dropped. | `bahrain_me_curves.pdf`, `greenland_me_curves.pdf`, `bahrain_l1_walkthrough_{full,noBahrain}.pdf` |

## Slide decks

Each is a minimal Beamer deck (dove colour, no navigation), self-contained,
sharing the `figs/` directory via `\graphicspath{{../figs/}}`. The filename
prefix `01_`…`06_` indicates intended reading order.

| # | Deck | Role |
|---|---|---|
| 1 | `01_distance_slides` | **Problem statement.** Fit ME curves per reanalysis. Define a scalar L1 distance between them; attribute the distance by continent / decade / temp decile. ERA5 vs MERRA2 is the biggest gap (0.0035), concentrated in the hot tail and driven by tropical countries (Sierra Leone alone ≈ 25%). |
| 2 | `02_temp_compare_slides` | **Diagnose the disagreement.** Raw correlations 0.992+ are deceptively high; residualised 0.74–0.96. Country disagreement is largely *static bias* (CRI, HTI, KGZ…); weather-signal disagreement concentrates in tropical low-income countries (corr with log GDP = −0.45, mean T = +0.33). GMFD 2009–2010 anomaly flagged. |
| 3 | `03_oos_slides` | **Ruling out an easy answer.** Can OOS predictive fit decide which dataset is best? No: all five are statistically indistinguishable (CV-MSE within 0.3%); leave-country-out R² is negative for all. Choice is under-determined by predictive fit. |
| 4 | `04_synth_compare_slides` | **First synthesis attempt.** Six synthetic measures (MEAN, MEDIAN, PC1, FACTOR, MIN-VAR, ANCHORED). Best three on t-stat are all synthetics (FACTOR, PC1, MEAN) — they recover a ~50% less-attenuated β̂₂. OOS-OPT rejected. Factor σ_k estimates: ERA5 ≈ JRA-3Q (~0.06°C) cleanest, GMFD (0.21°C) noisiest. |
| 5 | `05_hetero_factor_slides` | **Refined synthesis + damage-function reconciliation.** Heteroskedasticity is real on the climate axis (MERRA2 1.58×, GMFD 0.71× inverted) and matters for β̂₂. T-dependent MIN-VAR-T synthesis wins (t = 3.86, smallest CI, most stable jackknife). Classical EIV correction collapses the cross-dataset ME spread: raw β̂₂'s are attenuated 10–45%; post-correction, four of five methods agree on ME(25°C) ∈ [−0.018, −0.014]. |
| 6 | `06_loo_country_attribution_slides` | **Country-by-country attribution (Python pipeline, 4 reanalyses, 0.25°).** Re-runs the LOO sweep with a "% reduction in mean pairwise L1 distance" metric on a 0–100% scale; generalises deck 01's pairwise (ERA5-vs-MERRA2) story to all six pairs. Per-pair maps, scatters vs T and income, Bahrain and Greenland worked examples, contrast of pop-weighted vs uniform integration, and per-pair maps in both % and absolute-L1 units. Concludes with a sensitivity section showing the country ranking is robust to the weighting choice (Spearman 0.99 uniform vs Kevin's pop). |

## How to rebuild

```bash
# R pipeline — run from the Dropbox root (the R scripts setwd() there)
cd /Users/tombearpark/Library/CloudStorage/Dropbox/ks
Rscript impact_uncertainty/tb/code/explore/distance.R
Rscript impact_uncertainty/tb/code/explore/distance_figs.R
Rscript impact_uncertainty/tb/code/explore/oos_score.R
Rscript impact_uncertainty/tb/code/explore/temp_compare.R
Rscript impact_uncertainty/tb/code/explore/synth_compare.R
Rscript impact_uncertainty/tb/code/explore/hetero_factor.R
Rscript impact_uncertainty/tb/code/gdp_disagree.R
Rscript impact_uncertainty/tb/code/mortality_disagree.R

# Python pipeline — run from the repo root
cd /Users/tombearpark/Documents/GitHub/impact_uncertainty
python3 tb/code/loo_country_influence.py
python3 tb/code/loo_perpair_maps.py
python3 tb/code/loo_scatter.py
python3 tb/code/loo_residual_disagree.py
python3 tb/code/bahrain_walkthrough.py
python3 tb/code/country_residuals.py
# Optional sanity check: pop-weighted recompute vs Kevin's CSV
python3 tb/code/funcs_dme.py

# Compile slide decks
cd tb/slides
for f in 01_distance_slides 02_temp_compare_slides 03_oos_slides \
         04_synth_compare_slides 05_hetero_factor_slides \
         06_loo_country_attribution_slides; do
  pdflatex -interaction=nonstopmode $f.tex
  pdflatex -interaction=nonstopmode $f.tex   # 2nd pass for TOC
done
rm -f *.aux *.log *.nav *.out *.snm *.toc *.vrb
```

The Extended Data figures (`../extended_data/`) are rebuilt independently
via `python3 extended_data/make_extended_data.py`; see
[`../extended_data/README.md`](../extended_data/README.md) for that
pipeline.
