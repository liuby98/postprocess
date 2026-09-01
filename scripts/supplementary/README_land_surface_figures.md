# Supplementary land evaluation and process diagnostics

The supplementary land scripts now evaluate only the manuscript-relevant
variables requested for the gravel experiments:

| Script | CoLM variable(s) | Reference used by default |
|---|---|---|
| `land_evapotranspiration_response.py` | `f_fevpa` (total ET only) | GLEAM v4.2a `E` |
| `land_runoff_response.py` | `f_rnof` (total runoff only) | G-RUN ENSEMBLE multi-model median (MMM) |
| `land_snow_response.py` | `f_scv` (SWE), `f_fsno` (snow-cover fraction) | ERA5-Land SWE benchmark; MODIS Terra MOD10CM C6.1 SCF |
| `land_water_balance_response.py` | `f_qinfl` (infiltration only) | No reference; process response diagnostic only |

The former surface-energy group is intentionally removed. ET components,
surface/subsurface runoff, snow depth, recharge, storage, and water-table depth
are no longer plotted by these scripts.

## Reference-evaluated variables

ET, total runoff, SWE, and snow-cover fraction are processed independently.
Their output set is:

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

These comparisons use one common time window, the same 0.25-degree target grid,
cosine-latitude area weighting, and a unified valid-cell mask across reference,
CTL, and EXP.

## Infiltration without reference evaluation

Infiltration is retained only as a model process diagnostic. It does not read
MERRA-2 or any other reference product and does not calculate bias, RMSE, KGE,
or spatial ACC. It produces:

1. `supp_land_infiltration_off_spatial.pdf` and
   `supp_land_infiltration_cpl_spatial.pdf`, with MAM/JJA rows and CTL, EXP,
   and EXP-CTL columns;
2. `supp_land_infiltration_response_timeseries.pdf` and `.csv`, showing the
   yearly China-land EXP-CTL response with a 99% grid-sampling confidence
   interval and a paired test across years.

## Color ranges

The fixed wide limits in the original response figures have been removed.
Defaults now use the 2nd-98th percentiles of all Reference/CTL/EXP state fields
(CTL/EXP for infiltration) and the 95th percentile of absolute bias/change
fields. This keeps isolated
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

On the Figshare file grid, download only the first standalone file,
`G-RUN_ENSEMBLE_MMM.nc` (about 340 MB). Do not use **Download all** and do not
download the multi-GB forcing/member ZIP archives. The standalone MMM file is
the authors' recommended single median estimate and already covers 1902-2019.

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

For seasonal-mean evaluation, use the much smaller **ERA5-Land monthly
averaged** data set rather than the hourly page. Select only `Snow depth water
equivalent`, years 2001-2017, months March-August, monthly averaged reanalysis,
NetCDF, and a China subset. If the web form does not allow or submit all years,
use the CDS API and split the request by year; large ERA5-Land selections are
deliberately restricted to protect the CDS queue.

### Snow-cover fraction

**Use MODIS Terra MOD10CM Collection 6.1 monthly snow cover**, matching the
reference choice in the attached manuscript. It provides monthly mean snow
cover extent on a 0.05-degree CMG beginning in March 2000 and covers all of
2001-2017.
MOD10C1 daily data can be monthly averaged if MOD10CM is unavailable. ESA Snow
CCI fractional snow cover is a useful second reference for sensitivity.

Sources: <https://nsidc.org/data/mod10cm/versions/61>,
<https://modis.gsfc.nasa.gov/data/dataprod/mod10.php>, and
<https://climate.esa.int/en/projects/snow/>.

From the NSIDC product list shown in the download interface, click the Terra
product ID **MOD10CM** (not daily MOD10C1 and not Aqua MYD10CM), open **Data
Access & Tools**, then launch **Earthdata Search** and sign in with a free NASA
Earthdata Login. Apply the temporal filter and download the Version 61
granules. Native files are HDF-EOS2; they must be QA-filtered and converted to
a regular latitude-longitude CF-NetCDF before being passed to the plotting
script.

#### Convert and merge MOD10CM HDF files

`preprocess_mod10cm.py` performs the complete conversion in one pass. It reads
the HDF4 granules one at a time, extracts `Snow_Cover_Monthly_CMG` and
`Snow_Spatial_QA`, keeps provider-valid values in the range 0-100%, masks the
categorical values 211/250/253/254/255, retains only QA=0 by default, subsets
65-145 E and 5-60 N, sorts months parsed from `AYYYYDDD`, and writes one
compressed CF-NetCDF. Missing months are reported but never converted to zero.

Install the HDF/NetCDF readers in the plotting environment:

```bash
conda activate plot
conda install -c conda-forge pyhdf netcdf4 numpy
```

With the directory layout below, first validate all downloaded files:

```bash
cd /share/home/dq135/draw/scripts/supplementary
python preprocess_mod10cm.py --check-only
```

Then convert and merge them:

```bash
python preprocess_mod10cm.py
```

The default output is:

```text
/share/home/dq135/openbench/Reference/Grid/HigRes/Snow/
Snow_Cover_Fraction/MODIS_MOD10CM/MOD10CM_SCF_2001_2017.nc
```

The output stores ascending cell-center latitude, longitude, mid-month `time`,
monthly `time_bnds`, filtered snow-cover percentage, and the original QA field.
If the output exists, use `--overwrite`. To require all 204 calendar months,
add `--strict-coverage`; otherwise documented MODIS outages remain missing. A
different subset can be supplied with, for example, `--bbox 70 10 140 55`.

For the 2001-2017 archive used here, official outages leave June 2001, March
2002, December 2003, and February 2016 unavailable. Only the first two affect
the March-August evaluation. The snow script therefore treats 2001 JJA and
2002 MAM as unavailable SCF season-years and removes the same season-years from
the corresponding CTL and EXP comparisons. It does not substitute zero,
interpolate a missing month, or calculate a two-month mean as a three-month
season. Consequently, SCF MAM and JJA metrics each use 16 common complete
season-years, while ERA5-Land SWE still uses all 17 years.

No intermediate monthly NetCDF files are created. The merged output is
immediately compatible with the default SCF glob in `land_snow_response.py`.

## Reference directory suggestions

The non-GLEAM products are not present in the supplied OpenBench directory
listing. A consistent layout is:

```text
/share/home/dq135/openbench/Reference/Grid/LowRes/Water/
  Total_Runoff/G_RUN_ENSEMBLE/*.nc

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

The check is intentionally strict: all four model files must contain the
requested variables and every configured ET/runoff/snow reference must contain
March-August data for all years from 2001 through 2017. Infiltration requires
only the four model files.

## Run all requested figures

```bash
python land_evapotranspiration_response.py
python land_runoff_response.py \
  --runoff-reference-file '/path/to/g-run-mmm*.nc'
python land_snow_response.py \
  --swe-reference-file '/path/to/era5_land_swe*.nc' \
  --scf-reference-file '/path/to/mod10cm*.nc'
python land_water_balance_response.py
```

Use `--skip-spatial`, `--skip-regional-timeseries`, or (for reference-evaluated
variables) `--skip-acc-timeseries` to redraw only selected products.

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
