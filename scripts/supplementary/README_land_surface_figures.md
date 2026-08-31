# Supplementary land-surface figures

This folder extends the existing coupled soil-temperature and soil-moisture
evaluation scripts with five process-oriented figure groups.  The new figures
use the same China domain, Lambert projection, 0.25-degree CN05.1 land mask,
China/province/South China Sea shapefiles, sans-serif fonts, and PDF output
style as the manuscript scripts already in this repository.

## Figure scripts

| Script | CoLM variables | Scientific role |
|---|---|---|
| `land_snow_response.py` | `f_scv`, `f_snowdp`, `f_fsno` | Snow water equivalent, snow depth, and snow-cover fraction; tests whether snow-related albedo/insulation changes contribute to the gravel response. |
| `land_evapotranspiration_response.py` | `f_fevpa`, `f_fevpl`, `f_etr`, `f_fevpg` | Total ET and its canopy/transpiration/ground components; connects soil-moisture changes to latent-heat partitioning. |
| `land_runoff_response.py` | `f_rsur`, `f_rsub`, `f_rnof` | Surface, subsurface, and total runoff; diagnoses how gravel redistributes hydrological losses. |
| `land_water_balance_response.py` | `f_qinfl`, `f_qcharge`, `f_wat`, `f_zwt` | Infiltration, groundwater recharge, total water storage, and water-table depth; provides a broader water-balance context. |
| `land_energy_partition_response.py` | `f_rnet`, `f_fsena`, `f_lfevpa`, `f_fgrnd`, `f_olrg` | Net radiation, sensible/latent heat, ground heat flux, and outgoing longwave radiation; directly supports the surface-energy mechanism. |

The shared implementation is in `_land_response_common.py`.  Each physical
type still has its own entry script and PDF output.

## Default comparison

- Period: 2001-2017 (17 years).
- Seasons: spring (March-May) and summer (June-August).
- Columns: `(EXP - CTL) OFF` spring, `(EXP - CTL) OFF` summer,
  `(EXP - CTL) CPL` spring, and `(EXP - CTL) CPL` summer.
- Significance: two-sided paired t test across yearly seasonal means.  Gray
  dots mark `p < 0.05`; use `--significance-style mask` to retain only
  significant response values, matching the significant-only treatment in the
  existing coupled soil-temperature/moisture scripts.
- Flux conversion: water fluxes are converted from `mm s-1` to `mm day-1`.
  Snow depth is converted from m to cm, and snow-cover fraction to percentage
  points.

The default input filenames are `colmoff_2001-2017_nogravel.nc`,
`colmoff_2001-2017_gravel.nc`, `colmrun_2001-2017_nogravel.nc`, and
`colmrun_2001-2017_gravel.nc`. Therefore, OFF and CPL use the same 2001-2017
window as the existing supplementary soil-temperature/moisture evaluation.

## Run

Run from `scripts/supplementary`:

```bash
python land_snow_response.py
python land_evapotranspiration_response.py
python land_runoff_response.py
python land_water_balance_response.py
python land_energy_partition_response.py
```

Outputs are written to `scripts/supplementary/figs/`.  Before a long plotting
job, validate all paths, variables, dimensions, and time lengths:

```bash
python land_snow_response.py --check-only
```

All input paths can be overridden without editing source code.  For example:

```bash
python land_runoff_response.py \
  --data-dir /path/to/processed/files \
  --wrfinput /path/to/wrfinput_d01 \
  --mask-file /path/to/CN05_mask.nc \
  --output-dir ./figs
```

Use `python SCRIPT.py --help` for the four individual CTL/EXP overrides and
other plotting options.

## Caption templates

**Snow response.** Spatial distributions of spring and summer differences in
snow water equivalent, snow depth, and snow-cover fraction between the
gravel-inclusive and gravel-free experiments during 2001-2017. Columns show
the offline spring, offline summer, coupled spring, and coupled summer
responses, respectively. Gray dots denote differences significant at the 5%
level based on a two-sided paired Student's t test across yearly seasonal
means.

**Evapotranspiration response.** As in the snow-response figure, but for total
evapotranspiration, leaf evaporation plus transpiration, transpiration, and
ground evaporation. Units are `mm day-1`.

**Runoff response.** As in the snow-response figure, but for surface,
subsurface, and total runoff. Units are `mm day-1`.

**Water-balance response.** As in the snow-response figure, but for
infiltration, groundwater recharge, total water storage, and water-table
depth. Units are `mm day-1`, `mm day-1`, mm, and m, respectively.

**Energy-partition response.** As in the snow-response figure, but for net
radiation, sensible heat flux, latent heat flux, ground heat flux, and outgoing
longwave radiation. Units are `W m-2`.

## Interpretation cautions

- `f_scv` is snow water equivalent (mm), whereas `f_fsno` is fractional snow
  cover.  They should not be described interchangeably.
- `f_fevpa` is the total evapotranspiration flux to the atmosphere.
  `f_fevpl` already contains leaf evaporation plus transpiration, and `f_etr`
  is the transpiration subset; therefore the displayed ET rows must not be
  summed as independent components.
- `f_rnof` is the model's total-runoff diagnostic.  Check its numerical closure
  against `f_rsur + f_rsub` before claiming exact component additivity.
- The supplied NetCDF metadata give `f_qinfl` only the long name `f_qinfl`.
  The figure labels it as infiltration following the established CoLM variable
  convention, but any detailed process attribution should be checked against
  the exact CoLM source version used for the experiments.
- A positive `f_zwt` response means a greater modeled depth to the water table
  (a deeper water table), not more groundwater.
- Snow fields can be very small in summer.  A weak JJA snow response should be
  treated as evidence against, rather than support for, a dominant summer snow
  mechanism.
