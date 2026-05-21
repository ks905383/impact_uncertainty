# Extended Data — figures, sources, captions

This folder holds the three Extended Data figures for the paper *Observational
uncertainty dominates climate impact projections* plus the code/source needed to regenerate
them and the proposed figure captions.

## Contents

| File | Type | Description |
|---|---|---|
| `make_extended_data.py` | Python script | Self-contained generator for ED Fig 1 + ED Fig 2. |
| `extdata_fig1_temp_disagreement.{pdf,png}` | Figure | ED Fig 1 — disagreement in temperature shocks. |
| `extdata_fig2_me_disagreement.{pdf,png}` | Figure | ED Fig 2 — country drivers of cross-reanalysis disagreement in the BHM dose-response. |
| `extdata_fig3_mechanism.tex` | LaTeX / TikZ source | Standalone source for the mechanism schematic. |
| `extdata_fig3_mechanism.pdf` | Figure | ED Fig 3 — mechanism schematic (compiled). |

PDFs are the publication-grade outputs; PNGs are 200 dpi previews. ED Fig 3
is a standalone TikZ figure with a tight bounding box (no PNG version
shipped, vector only).

## Requirements

### Python (for ED Fig 1 + ED Fig 2)

- Python ≥ 3.10
- `numpy`, `pandas`, `scipy`, `matplotlib`, `geopandas`

No specific minimum versions enforced; tested with `numpy 1.26`, `pandas
2.1`, `scipy 1.13`, `matplotlib 3.8`, `geopandas 1.0`. `geopandas` is the
only non-trivial dependency — install with
`pip install geopandas` (it will pull `shapely`, `pyproj`, `fiona`).

### LaTeX (for ED Fig 3)

- A TeX distribution with PDF output (TeX Live 2020+ or MacTeX).
- Packages used: `standalone`, `lmodern`, `xcolor`, `amsmath`, `amssymb`,
  `tikz` (with libraries `positioning`, `arrows.meta`, `calc`, `fit`,
  `backgrounds`).

All are in the standard TeX Live full install.

## Input data

### For `make_extended_data.py`

Two root paths, both configurable via environment variables:

| Env var | Default (Tom's setup) | Purpose |
|---|---|---|
| `KS_DROPBOX` | `/Users/tombearpark/Library/CloudStorage/Dropbox/ks` | Root of the shared regression data |
| `KS_REPO` | `/Users/tombearpark/Documents/GitHub/impact_uncertainty` | Root of the GitHub repo (provides the shapefile) |

Specific input files under `KS_DROPBOX/`:

- `reg_inputs/gdppc/reginputs_gdppc_025deg.csv` — balanced BHM
  regression-input panel: country × year × reanalysis rows with
  `temp_popweight`, `Pop`, `gdpCAPppp_wdi`, `growthWDI`, and precomputed
  `time`, `time2` for the trend interactions.
- `reg_outputs/gdppc/regoutputs_coeffs_gdppc_025deg_<rea>_1980-2010.csv` —
  full-sample BHM coefficients (β₁, β₂) for each of the four main
  reanalyses (ERA5-025, MERRA2, JRA-3Q, GMFD). Four files total.
- `reg_outputs/gdppc/loos/regoutputs_coeffs_gdppc_025deg_<rea>_1980-2010_no<country>.csv`
  — leave-one-country-out BHM coefficients for each (reanalysis × dropped
  country). ~780 files total (195 countries × 4 reanalyses).
  - If only the archive `reg_outputs/gdppc/loos.zip` is present, unpack it
    once into `reg_outputs/gdppc/loos/` before running.

And under `KS_REPO/`:

- `aux_data/geo_data/country-fixed.shp` — ADM0 country polygons (already
  in the repo).

### For `extdata_fig3_mechanism.tex`

No data inputs — pure TikZ.

## How to (re)generate

### ED Fig 1 + ED Fig 2

From the repo root:

```bash
python3 extended_data/make_extended_data.py
```

If your data lives elsewhere:

```bash
KS_DROPBOX=/path/to/ks  KS_REPO=/path/to/repo \
  python3 extended_data/make_extended_data.py
```

Outputs are written next to the script as `extdata_fig{1,2}_*.{pdf,png}`.
Runs in ~30 seconds on a modern laptop; the only slow step is the one-time
read of the ~780 LOO coefficient CSVs.

### ED Fig 3

```bash
cd extended_data
pdflatex -interaction=nonstopmode extdata_fig3_mechanism.tex
rm -f extdata_fig3_mechanism.{aux,log,out}
```

The `standalone` document class produces a PDF whose bounding box hugs the
diagram (~205 × 200 mm); no page margin is added.

## Running this on Kevin's setup

The Python script reads from `$KS_DROPBOX` (default
`~/Library/CloudStorage/Dropbox/ks`). To point it at OneDrive instead:

```bash
KS_DROPBOX="/Users/kevin/Library/CloudStorage/OneDrive-Personal/impact_uncertainty" \
KS_REPO="/Users/kevin/path/to/impact_uncertainty" \
  python3 extended_data/make_extended_data.py
```

(Substitute the actual OneDrive path. The script just needs `KS_DROPBOX`
to contain `reg_inputs/gdppc/…` and `reg_outputs/gdppc/…` directly
underneath.)

A few naming conventions the script depends on — verify these match on
Kevin's side before assuming a path-only fix is enough:

- LOO files: `regoutputs_coeffs_gdppc_025deg_<rea>_1980-2010_no<country>.csv`
  exactly (parsed by `_LOO_RE` near the top of the script).
- The `<country>` token must use the full WDI country name including
  spaces and hyphens, e.g. `noPapua New Guinea`, `noGuinea-Bissau`.
- Reanalysis tokens must be one of `ERA5-025 | MERRA2 | JRA-3Q | GMFD`.
  (Kevin's `dmes_gdppc.csv` uses `ERA5` not `ERA5-025`; the helper
  reconciles these via `REA_TO_DME` in `funcs_dme.py` upstream, but the
  Extended Data script reads the coefficient files directly and so uses
  the filename-side spelling.)
- The script skips `UDel` (the current paper uses the 4-reanalysis set).

If any convention differs, edit `REANALYSES`, `_FULL_RE`, and `_LOO_RE`
near the top of `make_extended_data.py`.

ED Fig 3 has no data dependency — recompiling is a one-line `pdflatex`
invocation independent of any Dropbox/OneDrive paths.

---

## Proposed captions

These are first drafts intended for the final paper. They are deliberately
verbose so each figure is self-contained for a reader who has not yet read
the main text; trim as the journal allows.

### Extended Data Fig. 1

**Inter-reanalysis disagreement in the temperature shocks that identify
the BHM dose-response.** Panels **(a)–(d)** show residualised annual
population-weighted mean temperature for the two countries with the
smallest cross-reanalysis disagreement (**a**, China; **b**, Korea) and
the two with the largest (**c**, Maldives; **d**, Seychelles), one line
per reanalysis (ERA5, MERRA2, JRA-3Q, GMFD). Residuals come from
country-by-country regressions of `temp_popweight` on a constant absorbed
by country fixed effects, year fixed effects, and country-specific linear
and quadratic time trends (`temp ~ 1 | iso + year + iso[t] + iso[t²]`,
fit separately on each reanalysis), matching the fixed-effect structure
of the BHM growth regression. The “mean |Δresidual T|” reported in each
panel title is the country-level mean over years (1980–2010) and over the
six unordered reanalysis pairs of |residᵃ − residᵇ|. **(e, f)** Scatter
of the same per-country mean |Δresidual T| (averaged across the six
pairs) against (**e**) the country-mean annual temperature over
1980–2010 and (**f**) the log of the country-mean GDP per capita (PPP,
WDI). Each point is one country (n = 195); dashed line is the OLS fit.
The panel sample is the balanced (country, year, reanalysis) panel
underlying the BHM regression (5,643 country-years × 4 reanalyses =
22,572 rows). Hotter and lower-income countries see systematically larger
inter-reanalysis disagreement in the residualised temperature shocks that
identify β.

### Extended Data Fig. 2

**Country drivers of cross-reanalysis disagreement in the fitted BHM
dose-response, in absolute L1 units.** For each pair of reanalyses (a, b)
from the four-product set {ERA5-025, MERRA2, JRA-3Q, GMFD} (six unordered
pairs), the L1 distance between marginal-effect curves
`MEₐ(T) = β₁ᵃ + 2β₂ᵃT` and `ME_b(T)` is computed as the unweighted mean
of |MEₐ(T) − ME_b(T)| evaluated at every observation in the pair-pooled
empirical temperature sample. For each country *c*, the BHM regression is
refit separately on each reanalysis with *c* removed from the panel, the
ME curves and pair-L1 are recomputed, and the **absolute reduction** is
defined as `L1(full sample) − L1(c dropped)` (positive ⇒ dropping *c*
shrinks the pair’s disagreement; negative ⇒ dropping *c* widens it). All
units are raw L1 (no per-pair normalisation), so the six per-pair panels
share a single colour scale.
**(a)** Per-pair maps: one panel per reanalysis pair. The sequential
purple colour scale is shared across all six panels and is clipped at the
75th percentile of *positive* absolute reductions (≈ 1 × 10⁻⁴); values
above that cut saturate at the deepest purple. The annotation
*n = X / 195* gives the number of countries for which dropping that
country lowers the pair’s L1. Countries with absolute reduction ≤ 0
(consensus-glue countries that hold the pair together) and countries
with no LOO data are drawn in light grey.
**(b)** Combined map: per-country mean across the six pairs of the same
absolute-reduction quantity, on the same colour scale as (a). The
strongest mean drivers are concentrated in tropical Sub-Saharan Africa
and a number of small island states; high-income temperate countries
typically appear pale, consistent with the dose-response disagreement
being driven by the locations where the underlying temperature record
itself disagrees most across reanalyses.

### Extended Data Fig. 3

**Schematic of the pathways by which the choice of reanalysis product
affects projected climate impacts.** The reanalysis-product choice (top,
slate) enters the projection chain through two coupled pathways.
**Pathway A** (left, blue): the chosen reanalysis interacts with seven
analyst-side modelling choices — **(i)** the functional form of the
dose-response D(T) (polynomial, threshold, bin); **(ii)** the fixed-effect
and time-trend specification, which determines what within-unit variation
identifies β and how reanalysis measurement noise is loaded;
**(iii)** the choice of climate exposure metric (mean vs maximum daily
temperature, integrated vs binned); **(iv)** features of the outcome
variable, in particular whether outcome-side measurement errors are
correlated with reanalysis-side errors in temperature; **(v)** outlier
(country, year) cells where extreme T coincides with extreme outcome
realisations and dominate the fitted quadratic; **(vi)** the temporal and
spatial sample (which years and locations enter the panel, jointly
constrained by reanalysis coverage and outcome data availability); and
**(vii)** the aggregation weights used to roll grid-cell climate up to
the regression unit, which act like correlated measurement error and bias
β when those weights co-vary with sub-unit climate heterogeneity. These
interactions produce a fitted dose-response D(T) (terracotta) that
differs across reanalyses, especially in the tails of the temperature
distribution. **Pathway B** (right, blue): the same reanalysis acts as
the bias-correction *target* for the GCM ensemble in the BCD-ME
post-processing step, producing a future-climate trajectory at a given
Global Warming Level (GWL, terracotta) that also varies with the
reanalysis choice. The two pathways meet at the **amplification step**
(centre, blue): the in-sample D(T) is evaluated on out-of-sample future
temperatures, so small reanalysis-driven differences at the extremes of
the temperature distribution become large differences in projected
damages. The terminal box (terracotta) is the projected impact at the
GWL, together with the partitioning of its variance between
reanalysis-choice, GCM model, and internal-variability contributions —
and the implication for the Social Cost of Carbon. Colour encodes role:
**slate** for data inputs, **blue** for analysis steps, **terracotta**
for calculated/derived quantities. Dashed arrows indicate the reanalysis
choice entering each pathway; solid arrows indicate the computational
flow.
