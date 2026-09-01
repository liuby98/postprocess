# Supplementary land diagnostics

The supplementary scripts retain only the manuscript variables requested for
the gravel experiments.

| Script | CoLM variable(s) | Default reference | Output mode |
|---|---|---|---|
| `land_evapotranspiration_response.py` | `f_fevpa` | GLEAM v4.2a `E` | spatial only |
| `land_runoff_response.py` | `f_rnof` | G-RUN ENSEMBLE MMM | spatial only |
| `land_snow_response.py` | `f_scv`, `f_fsno` | ERA5-Land SWE, MOD10CM C6.1 SCF | monthly series plus SWE pies |
| `land_water_balance_response.py` | `f_qinfl` | none | spatial process diagnostic only |

The surface-energy group, ET components, runoff components, snow depth,
recharge, storage, and water-table depth are not plotted.

## Snow: monthly gravel-region diagnostics

`land_snow_response.py` no longer creates SWE/SCF spatial maps, seasonal
yearly series, regional metrics, or ACC figures. It creates only:

```text
supp_land_snow_monthly_timeseries_gravel_gt_0p3.pdf
supp_land_snow_monthly_timeseries_gravel_gt_0p3.csv
supp_land_swe_colm_snow_mask_timeseries_gravel_gt_0p3.pdf
supp_land_swe_colm_snow_mask_timeseries_gravel_gt_0p3.csv
supp_land_swe_range_proportions_gravel_gt_0p3.pdf
supp_land_swe_range_proportions_gravel_gt_0p3.csv
```

The monthly figure has four rows: SWE-OFF, SWE-CPL, SCF-OFF, and SCF-CPL.
ERA5-Land and MOD10CM use all available January-December observations. The
model files contain March-August only, so CTL and EXP are placed on the full
monthly axis with January-February and September-December set to missing. This
intentionally breaks the model curves between annual March-August segments.
Every month has a tick; the year is printed only below January.

The additional two-row SWE figure reports OFF and CPL averages over CoLM's
own snow-covered grid cells. For each experiment and month, its common mask is:

```text
gravel > 0.3 AND (CTL f_fsno >= threshold OR EXP f_fsno >= threshold)
```

ERA5-Land `sd`, CTL `f_scv`, and EXP `f_scv` are area-averaged on that
identical dynamic mask. Because CoLM contains March-August only, all three
curves in this mask-conditioned figure contain March-August only. The default
`f_fsno` threshold is 0.01; override it with `--snow-cover-threshold`. No
ERA5-Land snow-cover file is needed. This diagnostic does not replace the
all-gravel-region SWE curve in the four-row main figure.

All snow regional means and SWE pies use only China-land target cells whose
vertically averaged gravel volumetric fraction is strictly greater than 0.3.
The mask comes from:

```text
/share/home/dq013/zhwei/colm/data/CoLMrawdata/soil/vf_gravels_s.nc
```

Variables `vf_gravels_s_l1` through `vf_gravels_s_l8` are combined using the
CoLM ten-layer soil-thickness weights (layers 1-2 and 9-10 are mapped to the
first and eighth source layers). The high-resolution source is sampled to the
common 0.25-degree target grid before applying `gravel > 0.3`. Use
`--gravel-file` or `--gravel-threshold` to override either setting.

The SWE pie figure has MAM/JJA rows and ERA5-Land, OFF CTL, OFF EXP, CPL CTL,
and CPL EXP columns. It reports grid-count proportions in the requested
classes 35-50, 20-35, and 5-20 mm. Percentages are normalized among cells in
5-50 mm so each pie sums to 100%. Counts below 5 mm and above 50 mm are retained
in the CSV and explicitly excluded from the pie denominator.

MOD10CM has four actual missing months in the local 2001-2017 archive
(2001-06, 2002-03, 2003-12, and 2016-02). Those reference points remain
missing; they are never replaced by zero or interpolation.

## ET and total-runoff spatial layout

ET and total runoff create separate OFF and CPL PDFs. Each has MAM/JJA rows and
only four columns:

1. the actual reference product (`GLEAM v4.2a` or `G-RUN ENSEMBLE MMM`);
2. CTL minus that reference;
3. EXP minus that reference;
4. EXP minus CTL.

Longitude and latitude labels are forced horizontal. The state colorbar spans
exactly column 1, the bias colorbar exactly columns 2-3, and the response
colorbar exactly column 4. Tick positions are equally spaced. State limits use
the configured robust percentiles and bias/response limits use symmetric
absolute percentiles.

## Infiltration spatial layout

Infiltration has no reference evaluation. Each OFF/CPL PDF has MAM/JJA rows and
CTL, EXP, and EXP-CTL columns. The state colorbar exactly spans columns 1-2 and
the response colorbar exactly spans column 3. No infiltration time-series,
Bias, RMSE, KGE, or ACC file is generated.

## Default data locations

```text
/share/home/dq135/openbench/Reference/Grid/HigRes/Water/
  Evapotranspiration/GLEAM_v4.2a/E_*_GLEAM_v4.2a.nc

/share/home/dq135/openbench/Reference/Grid/LowRes/Water/
  Total_Runoff/G_RUN_ENSEMBLE/*.nc*

/share/home/dq135/openbench/Reference/Grid/HigRes/Snow/
  Snow_Water_Equivalent/ERA5_Land/*.nc*
  Snow_Cover_Fraction/MODIS_MOD10CM/*.nc*
```

Expected snow files include the 204-month ERA5-Land `sd` NetCDF and the merged
MOD10CM `Snow_Cover_Monthly_CMG` NetCDF produced by `preprocess_mod10cm.py`.
CoLM `f_fsno` supplies the dynamic snow-presence mask; no ERA5-Land `snowc`
file is required. MOD10CM remains the independent reference for evaluating
CoLM `f_fsno`.

## Convert MOD10CM HDF files

Install the HDF4/NetCDF readers once in the plotting environment:

```bash
conda activate plot
conda install -c conda-forge pyhdf netcdf4 numpy
```

Then validate and merge:

```bash
cd /share/home/dq135/draw/scripts/supplementary
python preprocess_mod10cm.py --check-only
python preprocess_mod10cm.py
```

The default merged file is:

```text
/share/home/dq135/openbench/Reference/Grid/HigRes/Snow/
Snow_Cover_Fraction/MODIS_MOD10CM/MOD10CM_SCF_2001_2017.nc
```

The preprocessor extracts `Snow_Cover_Monthly_CMG`, applies the provider value
and QA masks, sorts granules by resolved calendar month, and keeps missing
months missing. Use `--overwrite` to replace an existing output.

## Validate all inputs

Run from `scripts/supplementary`:

```bash
mkdir -p logs figs

python land_snow_response.py --check-only \
  2>&1 | tee logs/land_snow_response_check.log

python land_evapotranspiration_response.py --check-only \
  2>&1 | tee logs/land_et_response_check.log

python land_runoff_response.py --check-only \
  2>&1 | tee logs/land_runoff_response_check.log

python land_water_balance_response.py --check-only \
  2>&1 | tee logs/land_infiltration_response_check.log
```

Paths and variables can be supplied explicitly, for example:

```bash
python land_snow_response.py \
  --swe-reference-file '/path/ERA5_Land_SWE_2001-2017.nc' \
  --swe-reference-variable sd \
  --scf-reference-file '/path/MOD10CM_SCF_2001_2017.nc' \
  --scf-reference-variable Snow_Cover_Monthly_CMG \
  --gravel-file '/path/vf_gravels_s.nc' \
  --check-only
```

## Run

To validate every input and then run snow, ET, total runoff, and infiltration
sequentially, use the batch driver from an activated plotting environment:

```bash
nohup bash run_land_supplementary.sh \
  > logs/run_land_supplementary_driver.log 2>&1 &
echo $! > logs/run_land_supplementary.pid
```

The driver stops on the first failure and writes a separate timestamped log
for every check and figure script. Its default snow paths are listed above;
override `SWE_REFERENCE_FILE`, `MODIS_SCF_REFERENCE_FILE`, or
`SNOW_COVER_THRESHOLD` in the environment when needed.

Individual commands remain available:

```bash
nohup python land_snow_response.py \
  > logs/land_snow_monthly_$(date +%Y%m%d_%H%M%S).log 2>&1 &

python land_evapotranspiration_response.py
python land_runoff_response.py
python land_water_balance_response.py
```

The new entry points do not recreate obsolete snow spatial/ACC files or old
ET/runoff/infiltration time-series files. Existing PDFs from earlier runs are
not automatically deleted; archive or remove them before assembling the final
supplement if their names are no longer listed above.
