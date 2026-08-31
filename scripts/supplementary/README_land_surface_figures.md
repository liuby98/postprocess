# Supplementary reference-based land evaluation

The supplementary land scripts now evaluate only the manuscript-relevant
variables requested for the gravel experiments:

| Script | CoLM variable(s) | Reference used by default |
|---|---|---|
| `land_evapotranspiration_response.py` | `f_fevpa` (total ET only) | GLEAM v4.2a `E` |
| `land_runoff_response.py` | `f_rnof` (total runoff only) | G-RUN ENSEMBLE multi-model median (MMM) |
| `land_snow_response.py` | `f_scv` (SWE), `f_fsno` (snow-cover fraction) | ERA5-Land SWE benchmark; MODIS Terra MOD10CM C6.1 SCF |
| `land_water_balance_response.py` | `f_qinfl` (infiltration only) | MERRA-2 `QINFIL` reanalysis benchmark |

The former surface-energy group is intentionally removed. ET components,
surface/subsurface runoff, snow depth, recharge, storage, and water-table depth
are no longer plotted by these scripts.

## What every variable produces

Each variable is processed independently, even when two variables share one
entry script. The output set is:

1. `<stem>_off_spatial.pdf` and `<stem>_cpl_spatial.pdf`
   - rows: MAM and JJA;
   - columns: Reference, CTL, EXP, CTL-Reference, EXP-Reference, EXP-CTL;
   - CTL/EXP bias panels report spatial bias, RMSE, and pattern correlation;
   - gray dots on EXP-CTL denote a paired two-sided test with `p < 0.05`.
2. `<stem>_regional_timeseries.pdf` and `.csv`
   - 2001-2017 area-weighted China-land Reference/CTL/EXP time series;
   - 99% grid-sampling confidence bands, following the time-evolution idea in
     Zhang et al. (2024);
   - bias, RMSE, correlation, and KGE, following the multi-metric OpenBench
     evaluation practice.
3. `<stem>_spatial_acc_timeseries.pdf` and `.csv`
   - yearly spatial anomaly correlation (ACC) after removing each source's
     2001-2017 climatology;
   - a companion box plot and paired-test p value, following the evaluation
     structure of Fig. 6 in Zhang et al. (2024).
4. `<stem>_regional_metrics.csv`
   - the exact time-series metrics printed on the figure.

All comparisons use one common time window, the same 0.25-degree target grid,
cosine-latitude area weighting, and a unified valid-cell mask across reference,
CTL, and EXP.

## Color ranges

The fixed wide limits in the original response figures have been removed.
Defaults now use the 2nd-98th percentiles of all Reference/CTL/EXP state fields
and the 95th percentile of absolute bias/change fields. This keeps isolated
extremes from washing out the main spatial signal while retaining out-of-range
values with extended colorbar ends.

The ranges are printed in the log and can be tightened or widened without
editing source code:

```bash
python land_evapotranspiration_response.py \
  --state-quantiles 5 95 \
  --difference-quantile 90
```

## Reference products for 2001-2017

### Total evapotranspiration

**Use GLEAM v4.2a.** It is appropriate for the present ET evaluation and is
consistent with the GLEAM soil-moisture reference already used in this project.
Variable `E` is total terrestrial evaporation/ET in `mm day-1`. The local files
cover every required year:

```text
/share/home/dq135/openbench/Reference/Grid/HigRes/Water/
Evapotranspiration/GLEAM_v4.2a/E_YYYY_GLEAM_v4.2a.nc
```

GLEAM is not an independent direct observation at every grid cell; it is a
satellite/reanalysis-driven land evaporation product. It is nevertheless a
mainstream and suitable primary gridded reference for this analysis.

| Priority | Product | Coverage/resolution | Use and limitation |
|---|---|---|---|
| Primary | GLEAM v4.2a | 1980-2024, daily, 0.1 degrees | Already available locally and consistent with the soil-moisture evaluation. |
| Satellite sensitivity | MOD16A2GF C6.1 | 2000-present, 8-day, 500 m | Easy LP DAAC/AppEEARS access; apply QA and convert 8-day accumulations to rates before using this script. |
| Reanalysis sensitivity | ERA5-Land total evaporation | 1950-present, hourly, 0.1 degrees | Complete coverage and easy CDS access, but land-model-derived. |

Sources: <https://www.gleam.eu/>,
<https://doi.org/10.1038/s41597-025-04610-y>, and
<https://lpdaac.usgs.gov/products/mod16a2gfv061/>.

### Total runoff

Recommended products, all covering 2001-2017:

| Priority | Product | Coverage/resolution | Use and limitation |
|---|---|---|---|
| Primary | G-RUN ENSEMBLE MMM | 1902-2019, monthly, 0.5 degrees | Observation-based machine-learning reconstruction with forcing uncertainty; download the MMM/median field rather than all 525 members when storage is limited. |
| China sensitivity | CNRD v1.0 (`qtot`) | 1961-2018, daily/monthly, 0.25 degrees | Quality-controlled China natural-runoff reconstruction calibrated against 200 natural/near-natural catchments; model-derived and largely natural-flow oriented. |
| Sensitivity | ERA5-Land total runoff | 1950-present, hourly, 0.1 degrees | Very easy CDS access and complete China coverage, but model/reanalysis-derived. |
| Sensitivity | GLDAS-2.1 Noah | 2000-present, 3-hourly, 0.25 degrees | Easy NASA GES DISC access; also land-model-derived. |

G-RUN ENSEMBLE data and method:
<https://doi.org/10.6084/m9.figshare.12794075> and
<https://doi.org/10.1029/2020WR028787>. CNRD data and method:
<https://doi.org/10.6084/m9.figshare.13185410> and
<https://doi.org/10.1175/BAMS-D-20-0094.1>.

Place/preprocess the chosen G-RUN MMM NetCDF under the default directory or
pass its path explicitly. The reference must be a runoff *rate*; if a file is
stored as monthly accumulation, retain a correct `mm month-1` units attribute
so the script can divide by calendar days.

### Snow water equivalent

No single SWE product is ideal over the Tibetan Plateau. Recommended products
covering 2001-2017 are:

| Priority | Product | Coverage/resolution | Use and limitation |
|---|---|---|---|
| Primary full-domain benchmark | ERA5-Land snow depth water equivalent (`sd`) | 1950-present, hourly, 0.1 degrees | Complete spatial coverage and easy CDS access; reanalysis, not direct observation. |
| Observation sensitivity | ESA Snow CCI SWE v3.0 | 1979-2022, daily, 0.1 degrees | Observation-based climate data record; alpine areas are masked, so it must not be the only Tibetan Plateau reference. |
| Multi-product sensitivity | ECCC blended SWE | 1981-2020, daily/monthly, 0.5 degrees | Observation/model blend with full study-period coverage; coarser grid. |

Sources: <https://cds.climate.copernicus.eu/datasets/reanalysis-era5-land>,
<https://catalogue.ceda.ac.uk/uuid/b06c4c5ea7694d30b33e1db04f0ecb6a/>,
and <https://climate-scenarios.canada.ca/?page=blended-snow-nh>.

### Snow-cover fraction

**Use MODIS Terra MOD10CM Collection 6.1 monthly snow cover**, matching the
reference choice in the attached manuscript. It is a 0.05-degree monthly
fractional snow-cover product beginning in 2000 and covers all of 2001-2017.
MOD10C1 daily data can be monthly averaged if MOD10CM is unavailable. ESA Snow
CCI fractional snow cover is a useful second reference for sensitivity.

Sources: <https://modis.gsfc.nasa.gov/data/dataprod/mod10.php> and
<https://climate.esa.int/en/projects/snow/>.

### Infiltration

There is no mainstream global, time-varying, direct observational infiltration
product suitable for gridded 2001-2017 validation. The script therefore uses
**MERRA-2 M2T1NXLND `QINFIL`** as a reanalysis benchmark. `QINFIL` is the soil
water infiltration rate in `kg m-2 s-1`, available globally from 1980 to
present at 0.5 x 0.625 degrees. It must be described as a benchmark/intermodel
comparison, not observational truth.

Source: <https://www.earthdata.nasa.gov/data/catalog/ges-disc-m2t1nxlnd-5.12.4>.

## Reference directory suggestions

The non-GLEAM products are not present in the supplied OpenBench directory
listing. A consistent layout is:

```text
/share/home/dq135/openbench/Reference/Grid/LowRes/Water/
  Total_Runoff/G_RUN_ENSEMBLE/*.nc
  Infiltration/MERRA2/*.nc4

/share/home/dq135/openbench/Reference/Grid/HigRes/Snow/
  Snow_Water_Equivalent/ERA5_Land/*.nc
  Snow_Cover_Fraction/MODIS_MOD10CM/*.nc
```

Other locations work through the per-variable command-line options.
HDF/HDF-EOS products such as native MODIS files must first be QA-filtered and
converted to a regular latitude-longitude CF-NetCDF with an explicit time
coordinate and units attribute.

## Validate inputs

Run from `scripts/supplementary`. ET uses the existing GLEAM directory directly:

```bash
python land_evapotranspiration_response.py --check-only
```

Examples for downloaded references:

```bash
python land_runoff_response.py \
  --runoff-reference-file '/path/to/g-run-mmm*.nc' \
  --runoff-reference-variable runoff \
  --check-only

python land_snow_response.py \
  --swe-reference-file '/path/to/era5_land_swe*.nc' \
  --swe-reference-variable sd \
  --scf-reference-file '/path/to/mod10cm*.nc' \
  --scf-reference-variable Snow_Cover_Monthly_CMG \
  --check-only

python land_water_balance_response.py \
  --infiltration-reference-file '/path/to/MERRA2*tavg1_2d_lnd_Nx*.nc4' \
  --infiltration-reference-variable QINFIL \
  --check-only
```

If a reference file lacks a standard units attribute, supply an explicit
multiplier. For example, a fractional SCF file needs `100`:

```bash
python land_snow_response.py \
  --scf-reference-file '/path/to/scf_fraction*.nc' \
  --scf-reference-factor 100 \
  --swe-reference-file '/path/to/swe*.nc' \
  --check-only
```

The check is intentionally strict: all four model files and every configured
reference must contain March-August data for all years from 2001 through 2017.

## Run all requested figures

```bash
python land_evapotranspiration_response.py
python land_runoff_response.py \
  --runoff-reference-file '/path/to/g-run-mmm*.nc'
python land_snow_response.py \
  --swe-reference-file '/path/to/era5_land_swe*.nc' \
  --scf-reference-file '/path/to/mod10cm*.nc'
python land_water_balance_response.py \
  --infiltration-reference-file '/path/to/MERRA2*tavg1_2d_lnd_Nx*.nc4'
```

Use `--skip-spatial`, `--skip-regional-timeseries`, or
`--skip-acc-timeseries` to redraw only selected products.

## Other defensible supplementary variables

The next most useful additions for the gravel/land-atmosphere mechanism are:

1. precipitation (CN05.1/CMFD as observational reference), to distinguish
   land-parameter effects from water-input differences;
2. land-surface temperature (MODIS LST), directly connected to soil thermal
   and surface-energy responses;
3. surface albedo (MODIS MCD43), directly connected to gravel and snow effects;
4. terrestrial water-storage anomaly (GRACE/GRACE-FO) for 2002-2017 only.

Snowmelt, groundwater recharge, and heat-flux components are not recommended as
the next additions because globally independent references are weak or absent.
The three full-period variables above should be prioritized before expanding
the supplement further.
