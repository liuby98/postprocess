#!/usr/bin/env python3
"""Shared evaluation and response diagnostics for supplementary CoLM figures.

Each public entry script defines one or more variables and calls
``run_category``. For every variable the module produces reference-based
spatial evaluation, China-land seasonal time series, and year-by-year spatial
anomaly-correlation (ACC) time series.

Model files contain March-August monthly means in six consecutive records per
year. Reference products may be daily, monthly, or hourly CF-NetCDF files;
they are reduced to monthly means before MAM/JJA means are calculated.
"""

import argparse
import calendar
import csv
import glob
import os
import re
import warnings
from collections import namedtuple
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
from netCDF4 import Dataset, num2date
from scipy.interpolate import RegularGridInterpolator, griddata
from scipy.stats import ttest_rel

import cartopy.crs as ccrs
import cartopy.io.shapereader as shpreader


warnings.filterwarnings("ignore")
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Helvetica", "Arial", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False

SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_DATA_DIR = Path("/share/home/dq135/draw/data/CRESM_CN15_ALO_GRAVEL")
DEFAULT_WRFINPUT = Path(
    "/share/home/dq135/CRESM_0307S92/CASES/RUN_2001_gravel/wrfinput_d01"
)
DEFAULT_MASK_FILE = Path(
    "/share/home/dq135/reference/CN05.1_Tm_1991_2023_MAM_025x025.nc"
)
DEFAULT_SHAPEFILE_DIR = SCRIPT_DIR.parent.parent / "shapefile_China"
DEFAULT_OUTPUT_DIR = SCRIPT_DIR / "figs"

PLOT_EXTENT = [80, 130, 10, 55]
REFERENCE_SUBSET_EXTENT = [75, 135, 5, 60]
PROJ_CENTRAL_LON = 110.0
PROJ_CENTRAL_LAT = 40.0
PROJ_STD_PARALLELS = (30.0, 60.0)
SEASONS = (("Spring", (3, 4, 5)), ("Summer", (6, 7, 8)))
SimpleDate = namedtuple("SimpleDate", "year month day")


def build_parser(description, variable_specs, reference_mode=True):
    """Create the common CLI plus one reference-data group per variable."""
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("--data-dir", type=Path, default=DEFAULT_DATA_DIR)
    parser.add_argument("--wrfinput", type=Path, default=DEFAULT_WRFINPUT)
    parser.add_argument("--mask-file", type=Path, default=DEFAULT_MASK_FILE)
    parser.add_argument("--mask-variable", default="tm")
    parser.add_argument(
        "--shapefile-dir", type=Path, default=DEFAULT_SHAPEFILE_DIR
    )
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--off-ctl", type=Path, default=None)
    parser.add_argument("--off-exp", type=Path, default=None)
    parser.add_argument("--cpl-ctl", type=Path, default=None)
    parser.add_argument("--cpl-exp", type=Path, default=None)
    parser.add_argument("--nyears", type=int, default=17)
    parser.add_argument("--start-year", type=int, default=2001)
    parser.add_argument("--p-threshold", type=float, default=0.05)
    if reference_mode:
        parser.add_argument("--acc-threshold", type=float, default=0.20)
    parser.add_argument(
        "--significance-style",
        choices=("stipple", "mask", "none"),
        default="stipple",
    )
    parser.add_argument("--sig-stride", type=int, default=5)
    parser.add_argument(
        "--state-quantiles",
        nargs=2,
        type=float,
        default=(2.0, 98.0),
        metavar=("LOW", "HIGH"),
        help=(
            "Robust percentiles for Reference/CTL/EXP color limits."
            if reference_mode
            else "Robust percentiles for CTL/EXP color limits."
        ),
    )
    parser.add_argument(
        "--difference-quantile",
        type=float,
        default=95.0,
        help=(
            "Absolute percentile for symmetric bias and change limits."
            if reference_mode
            else "Absolute percentile for symmetric response limits."
        ),
    )
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument("--skip-spatial", action="store_true")
    parser.add_argument("--skip-regional-timeseries", action="store_true")
    if reference_mode:
        parser.add_argument("--skip-acc-timeseries", action="store_true")
    parser.add_argument(
        "--check-only",
        action="store_true",
        help=(
            "Validate model/reference variables and coverage without plotting."
            if reference_mode
            else "Validate model variables without plotting."
        ),
    )

    if reference_mode:
        for spec in variable_specs:
            key = spec["key"].replace("_", "-")
            dest = spec["key"]
            parser.add_argument(
                f"--{key}-reference-file",
                dest=f"{dest}_reference_files",
                action="append",
                metavar="PATH_OR_GLOB",
                help=(
                    "Reference NetCDF file or quoted glob. Repeat for "
                    "multiple files; overrides the script default."
                ),
            )
            parser.add_argument(
                f"--{key}-reference-variable",
                dest=f"{dest}_reference_variable",
                default=None,
            )
            parser.add_argument(
                f"--{key}-reference-product",
                dest=f"{dest}_reference_product",
                default=None,
                help="Product name printed on figures and CSV files.",
            )
            parser.add_argument(
                f"--{key}-reference-factor",
                dest=f"{dest}_reference_factor",
                type=float,
                default=None,
                help=(
                    "Explicit multiplier from reference values to plotted "
                    "units."
                ),
            )
    return parser


def resolve_model_files(args):
    return {
        "OFF": {
            "ctl": args.off_ctl
            or args.data_dir / "colmoff_2001-2017_nogravel.nc",
            "exp": args.off_exp
            or args.data_dir / "colmoff_2001-2017_gravel.nc",
        },
        "CPL": {
            "ctl": args.cpl_ctl
            or args.data_dir / "colmrun_2001-2017_nogravel.nc",
            "exp": args.cpl_exp
            or args.data_dir / "colmrun_2001-2017_gravel.nc",
        },
    }


def require_file(path, label):
    if not Path(path).is_file():
        raise FileNotFoundError(f"{label} not found: {path}")


def clean_array(data):
    if np.ma.isMaskedArray(data):
        # Convert integer masked arrays (for example MOD10CM uint8) to
        # floating point before replacing their masked cells with NaN.
        data = np.ma.asarray(data, dtype=np.float64).filled(np.nan)
    else:
        data = np.asarray(data, dtype=np.float64)
    data[np.abs(data) > 1.0e30] = np.nan
    data[data == -9999] = np.nan
    return data


def find_coordinate(nc_obj, names):
    lower = {name.lower(): name for name in nc_obj.variables}
    for candidate in names:
        actual = lower.get(candidate.lower())
        if actual is not None:
            return actual, clean_array(nc_obj.variables[actual][:])
    raise KeyError(f"None of the coordinate variables exists: {names}")


def load_plot_grid(wrfinput_file, mask_file, mask_variable):
    require_file(wrfinput_file, "WRF input grid")
    require_file(mask_file, "mask file")
    with Dataset(wrfinput_file) as nc_wrf:
        lat2d = clean_array(nc_wrf.variables["XLAT"][:])
        lon2d = clean_array(nc_wrf.variables["XLONG"][:])
        if lat2d.ndim == 3:
            lat2d = lat2d[0]
        if lon2d.ndim == 3:
            lon2d = lon2d[0]
    with Dataset(mask_file) as nc_mask:
        _, lat1d = find_coordinate(nc_mask, ("lat", "latitude"))
        _, lon1d = find_coordinate(nc_mask, ("lon", "longitude"))
        if mask_variable not in nc_mask.variables:
            raise KeyError(
                f"Mask variable {mask_variable!r} is absent from {mask_file}"
            )
        mask_data = clean_array(nc_mask.variables[mask_variable][:])
        while mask_data.ndim > 2:
            mask_data = mask_data[0]
        land_mask = np.isfinite(mask_data)
    if lat2d.shape != lon2d.shape:
        raise ValueError("XLAT and XLONG shapes do not match")
    if land_mask.shape != (lat1d.size, lon1d.size):
        raise ValueError(
            "Mask shape does not match coordinates: "
            f"{land_mask.shape} vs {(lat1d.size, lon1d.size)}"
        )
    lon_grid, lat_grid = np.meshgrid(lon1d, lat1d)
    return lat2d, lon2d, lat1d, lon1d, lat_grid, lon_grid, land_mask


def read_yearly_season(nc_obj, variable, nyears, months, factor):
    if variable not in nc_obj.variables:
        raise KeyError(f"{variable!r} is absent from {nc_obj.filepath()}")
    nc_var = nc_obj.variables[variable]
    needed_steps = nyears * 6
    if nc_var.ndim != 3:
        raise ValueError(
            f"{variable} in {nc_obj.filepath()} must be (time, y, x); "
            f"found {nc_var.dimensions}"
        )
    if nc_var.shape[0] < needed_steps:
        raise ValueError(
            f"{variable} in {nc_obj.filepath()} has {nc_var.shape[0]} time "
            f"steps; {needed_steps} are required"
        )
    offsets = tuple(month - 3 for month in months)
    indices = [
        6 * year + offset for year in range(nyears) for offset in offsets
    ]
    data = clean_array(nc_var[indices, :, :]) * factor
    return np.nanmean(
        data.reshape(nyears, len(months), data.shape[1], data.shape[2]), axis=1
    )


def regrid_rcm_to_regular(field, lat2d, lon2d, lat1d, lon1d):
    target_lon, target_lat = np.meshgrid(lon1d, lat1d)
    valid = np.isfinite(field) & np.isfinite(lat2d) & np.isfinite(lon2d)
    if np.count_nonzero(valid) < 3:
        return np.full(target_lon.shape, np.nan)
    points = np.column_stack((lon2d[valid], lat2d[valid]))
    output = griddata(
        points, field[valid], (target_lon, target_lat), method="linear"
    )
    if np.isnan(output).any():
        nearest = griddata(
            points, field[valid], (target_lon, target_lat), method="nearest"
        )
        output[np.isnan(output)] = nearest[np.isnan(output)]
    return output


def regrid_regular_to_target(field, lat, lon, lat_grid, lon_grid):
    lat = np.asarray(lat, dtype=float)
    lon = ((np.asarray(lon, dtype=float) + 180.0) % 360.0) - 180.0
    lat_order = np.argsort(lat)
    lon_order = np.argsort(lon)
    lat = lat[lat_order]
    lon = lon[lon_order]
    field = field[np.ix_(lat_order, lon_order)]
    unique_lon, unique_index = np.unique(lon, return_index=True)
    lon = unique_lon
    field = field[:, unique_index]
    interpolator = RegularGridInterpolator(
        (lat, lon),
        field,
        method="linear",
        bounds_error=False,
        fill_value=np.nan,
    )
    points = np.column_stack((lat_grid.ravel(), lon_grid.ravel()))
    return interpolator(points).reshape(lat_grid.shape)


def validate_model_inputs(model_files, variable_specs, nyears):
    needed_steps = nyears * 6
    for experiment, files in model_files.items():
        for role, path in files.items():
            require_file(path, f"{experiment} {role.upper()} file")
            with Dataset(path) as nc_obj:
                for spec in variable_specs:
                    variable = spec["variable"]
                    if variable not in nc_obj.variables:
                        raise KeyError(f"{variable!r} is absent from {path}")
                    nc_var = nc_obj.variables[variable]
                    if nc_var.ndim != 3 or nc_var.shape[0] < needed_steps:
                        raise ValueError(
                            f"{variable} in {path} must be 3-D with at least "
                            f"{needed_steps} time steps; found {nc_var.shape}"
                        )
            print(f"[OK] {experiment} {role.upper()}: {path}")


def reference_cli_value(args, spec, suffix):
    return getattr(args, f"{spec['key']}_reference_{suffix}")


def reference_product(args, spec):
    return (
        reference_cli_value(args, spec, "product")
        or spec["reference_product"]
    )


def resolve_reference_files(args, spec):
    patterns = reference_cli_value(args, spec, "files")
    if not patterns:
        patterns = list(spec.get("reference_globs", ()))
    found = []
    for pattern in patterns:
        expanded = os.path.expandvars(os.path.expanduser(str(pattern)))
        matches = sorted(glob.glob(expanded, recursive=True))
        if not matches and Path(expanded).is_file():
            matches = [expanded]
        found.extend(matches)
    files = [
        Path(item) for item in dict.fromkeys(found) if Path(item).is_file()
    ]
    if not files:
        option = spec["key"].replace("_", "-")
        defaults = "\n  ".join(map(str, patterns)) or "(none)"
        raise FileNotFoundError(
            f"No {reference_product(args, spec)} reference NetCDF files found. "
            f"Checked:\n  {defaults}\nSupply one or more "
            f"--{option}-reference-file 'PATH_OR_GLOB' arguments."
        )
    return files


def resolve_reference_variable(nc_obj, args, spec):
    override = reference_cli_value(args, spec, "variable")
    candidates = [override] if override else spec["reference_variables"]
    lower = {name.lower(): name for name in nc_obj.variables}
    for candidate in candidates:
        if candidate is None:
            continue
        actual = lower.get(candidate.lower())
        if actual is not None:
            return actual
    raise KeyError(
        f"None of {candidates} exists in {nc_obj.filepath()}. Use "
        f"--{spec['key'].replace('_', '-')}-reference-variable."
    )


def find_time_variable(nc_obj, data_variable):
    for name in ("time", "Time", "TIME", "date", "datetime"):
        if name in nc_obj.variables:
            return name
    dims = nc_obj.variables[data_variable].dimensions
    for dim in dims:
        if dim in nc_obj.variables and "time" in dim.lower():
            return dim
    raise KeyError(f"No time coordinate found in {nc_obj.filepath()}")


def infer_dates_without_cf(time_values, path):
    values = np.asarray(time_values)
    if values.size and np.nanmin(values) >= 1.0e7:
        dates = []
        for value in values:
            token = f"{int(value):08d}"
            dates.append(
                SimpleDate(
                    int(token[:4]), int(token[4:6]), int(token[6:8])
                )
            )
        return dates
    match = re.search(r"(?<!\d)((?:19|20)\d{2})(?!\d)", path.name)
    if match is None:
        raise ValueError(
            f"Cannot decode time in {path}: missing CF units and filename year"
        )
    year = int(match.group(1))
    ntime = len(values)
    if ntime == 12:
        return [SimpleDate(year, month, 15) for month in range(1, 13)]
    ndays = 366 if calendar.isleap(year) else 365
    if ntime not in (ndays, ndays * 24):
        raise ValueError(
            f"Cannot infer {ntime} time records in {path}; add CF time metadata"
        )
    repeats = 24 if ntime == ndays * 24 else 1
    dates = []
    for month in range(1, 13):
        for day in range(1, calendar.monthrange(year, month)[1] + 1):
            dates.extend([SimpleDate(year, month, day)] * repeats)
    return dates


def decode_dates(nc_obj, time_name, path):
    time_var = nc_obj.variables[time_name]
    values = clean_array(time_var[:])
    units = getattr(time_var, "units", None)
    if units:
        try:
            return list(
                num2date(
                    values,
                    units=units,
                    calendar=getattr(time_var, "calendar", "standard"),
                    only_use_cftime_datetimes=True,
                )
            )
        except Exception as exc:
            raise ValueError(f"Cannot decode CF time in {path}: {exc}") from exc
    return infer_dates_without_cf(values, path)


def inspect_reference_files(files, args, spec, required_years):
    coverage = {year: set() for year in required_years}
    first_variable = None
    first_units = ""
    for path in files:
        with Dataset(path) as nc_obj:
            variable = resolve_reference_variable(nc_obj, args, spec)
            time_name = find_time_variable(nc_obj, variable)
            dates = decode_dates(nc_obj, time_name, path)
            find_coordinate(nc_obj, ("lat", "latitude", "nav_lat", "y"))
            find_coordinate(nc_obj, ("lon", "longitude", "nav_lon", "x"))
            for date in dates:
                if date.year in coverage and 3 <= date.month <= 8:
                    coverage[date.year].add(date.month)
            if first_variable is None:
                first_variable = variable
                first_units = getattr(
                    nc_obj.variables[variable], "units", ""
                )
    missing = {
        year: sorted(set(range(3, 9)) - months)
        for year, months in coverage.items()
        if months != set(range(3, 9))
    }
    if missing:
        detail = ", ".join(
            f"{year}:{months}" for year, months in missing.items()
        )
        if spec.get("allow_incomplete_seasons", False):
            print(
                f"[WARN] {reference_product(args, spec)} lacks some "
                f"March-August months ({detail}); each affected season-year "
                "will be omitted from Reference, CTL, and EXP together"
            )
        else:
            raise ValueError(
                f"{reference_product(args, spec)} lacks required March-August "
                f"coverage (year:missing months): {detail}"
            )
    print(
        f"[OK] {reference_product(args, spec)}: {len(files)} file(s), "
        f"variable={first_variable}, units={first_units or '(missing)'}, "
        f"coverage={required_years[0]}-{required_years[-1]}"
    )


def coordinate_subset(coordinate, lower, upper, path, name):
    coordinate = np.asarray(coordinate, dtype=float)
    if coordinate.ndim != 1:
        raise ValueError(f"{name} in {path} must be one-dimensional")
    compare = coordinate
    if name == "longitude":
        compare = ((coordinate + 180.0) % 360.0) - 180.0
    indices = np.where((compare >= lower) & (compare <= upper))[0]
    if indices.size < 2:
        raise ValueError(
            f"{path} has insufficient {name} coverage over China"
        )
    return slice(int(indices.min()), int(indices.max()) + 1)


def read_reference_block(
    nc_obj,
    variable_name,
    time_name,
    time_indices,
    lat_name,
    lon_name,
    lat_slice,
    lon_slice,
):
    variable = nc_obj.variables[variable_name]
    time_dim = nc_obj.variables[time_name].dimensions[0]
    lat_dim = nc_obj.variables[lat_name].dimensions[0]
    lon_dim = nc_obj.variables[lon_name].dimensions[0]
    dims = variable.dimensions
    for dim in (time_dim, lat_dim, lon_dim):
        if dim not in dims:
            raise ValueError(
                f"{variable_name} dimensions {dims} do not include {dim!r}"
            )
    time_axis = dims.index(time_dim)
    lat_axis = dims.index(lat_dim)
    lon_axis = dims.index(lon_dim)
    key = [slice(None)] * variable.ndim
    key[time_axis] = list(time_indices)
    key[lat_axis] = lat_slice
    key[lon_axis] = lon_slice
    data = clean_array(variable[tuple(key)])
    order = [time_axis, lat_axis, lon_axis] + [
        axis
        for axis in range(variable.ndim)
        if axis not in (time_axis, lat_axis, lon_axis)
    ]
    data = np.transpose(data, order)
    if data.ndim > 3:
        print(
            f"[WARN] Averaging extra reference dimensions for {variable_name}: "
            f"{data.shape[3:]}"
        )
        data = np.nanmean(data, axis=tuple(range(3, data.ndim)))
    return data


def automatic_reference_factor(units, spec, year, month):
    quantity = spec["reference_quantity"]
    normalized = (
        str(units)
        .lower()
        .replace(" ", "")
        .replace("**", "^")
        .replace("−", "-")
        .replace(".", "")
        .replace("·", "")
    )
    if quantity == "rate_mm_day":
        if any(
            token in normalized
            for token in ("kgm-2s-1", "kg/m2/s", "mms-1", "mm/s")
        ):
            return 86400.0
        if any(
            token in normalized
            for token in ("kgm-2day-1", "mmday-1", "mm/day", "mmd-1")
        ):
            return 1.0
        if any(token in normalized for token in ("ms-1", "m/s")):
            return 86400.0 * 1000.0
        if any(token in normalized for token in ("mday-1", "m/day")):
            return 1000.0
        if "mmmonth-1" in normalized or "mm/month" in normalized:
            return 1.0 / calendar.monthrange(year, month)[1]
    elif quantity == "amount_mm":
        if normalized in ("m", "meter", "metre") or any(
            token in normalized
            for token in ("mofwaterequivalent", "mwe")
        ):
            return 1000.0
        if any(token in normalized for token in ("kgm-2", "kg/m2", "mm")):
            return 1.0
    elif quantity == "percent":
        if "%" in normalized or "percent" in normalized:
            return 1.0
        if normalized in ("1", "fraction", "unitless", ""):
            return 100.0
    default = spec.get("reference_default_factor")
    if default is not None:
        print(
            f"[WARN] Using fallback factor {default:g} for units {units!r} "
            f"({spec['title']})"
        )
        return float(default)
    raise ValueError(
        f"Cannot convert reference units {units!r} to {spec['unit']}. Use "
        f"--{spec['key'].replace('_', '-')}-reference-factor."
    )


def load_reference_yearly(
    files, args, spec, years, lat_grid, lon_grid, land_mask
):
    monthly_sum = {}
    monthly_count = {}
    common_lat = None
    common_lon = None
    override_factor = reference_cli_value(args, spec, "factor")
    required = {
        (year, month) for year in years for month in range(3, 9)
    }
    for path in files:
        with Dataset(path) as nc_obj:
            variable_name = resolve_reference_variable(nc_obj, args, spec)
            time_name = find_time_variable(nc_obj, variable_name)
            dates = decode_dates(nc_obj, time_name, path)
            lat_name, latitude = find_coordinate(
                nc_obj, ("lat", "latitude", "nav_lat", "y")
            )
            lon_name, longitude = find_coordinate(
                nc_obj, ("lon", "longitude", "nav_lon", "x")
            )
            lat_slice = coordinate_subset(
                latitude,
                REFERENCE_SUBSET_EXTENT[2],
                REFERENCE_SUBSET_EXTENT[3],
                path,
                "latitude",
            )
            lon_slice = coordinate_subset(
                longitude,
                REFERENCE_SUBSET_EXTENT[0],
                REFERENCE_SUBSET_EXTENT[1],
                path,
                "longitude",
            )
            lat_sub = latitude[lat_slice]
            lon_sub = longitude[lon_slice]
            if common_lat is None:
                common_lat = np.asarray(lat_sub)
                common_lon = np.asarray(lon_sub)
            elif not (
                common_lat.shape == np.asarray(lat_sub).shape
                and common_lon.shape == np.asarray(lon_sub).shape
                and np.allclose(common_lat, lat_sub)
                and np.allclose(common_lon, lon_sub)
            ):
                raise ValueError(
                    "Reference files use different grids; preprocess them to "
                    "one regular latitude-longitude grid"
                )
            groups = {}
            for index, date in enumerate(dates):
                key = (int(date.year), int(date.month))
                if key in required:
                    groups.setdefault(key, []).append(index)
            units = getattr(nc_obj.variables[variable_name], "units", "")
            for (year, month), indices in groups.items():
                block = read_reference_block(
                    nc_obj,
                    variable_name,
                    time_name,
                    indices,
                    lat_name,
                    lon_name,
                    lat_slice,
                    lon_slice,
                )
                factor = (
                    override_factor
                    if override_factor is not None
                    else automatic_reference_factor(units, spec, year, month)
                )
                block = block * factor
                block_sum = np.nansum(block, axis=0)
                block_count = np.sum(np.isfinite(block), axis=0)
                key = (year, month)
                if key not in monthly_sum:
                    monthly_sum[key] = block_sum
                    monthly_count[key] = block_count
                else:
                    monthly_sum[key] += block_sum
                    monthly_count[key] += block_count
    monthly_mean = {}
    for key in required:
        if key not in monthly_sum:
            if spec.get("allow_incomplete_seasons", False):
                continue
            raise ValueError(f"Reference data missing year/month {key}")
        with np.errstate(invalid="ignore", divide="ignore"):
            monthly_mean[key] = np.where(
                monthly_count[key] > 0,
                monthly_sum[key] / monthly_count[key],
                np.nan,
            )
    output = {}
    for season_name, months in SEASONS:
        yearly = []
        for year in years:
            keys = [(year, month) for month in months]
            missing_months = [
                month for year_month, month in zip(keys, months)
                if year_month not in monthly_mean
            ]
            if missing_months:
                native = np.full(
                    (len(common_lat), len(common_lon)), np.nan, dtype=float
                )
                print(
                    f"[WARN] {spec['title']} {season_name} {year}: "
                    f"missing month(s) {missing_months}; omitting this "
                    "season-year"
                )
            else:
                native = np.nanmean(
                    np.stack([monthly_mean[key] for key in keys]), axis=0
                )
            target = regrid_regular_to_target(
                native, common_lat, common_lon, lat_grid, lon_grid
            )
            yearly.append(np.where(land_mask, target, np.nan))
        output[season_name] = np.stack(yearly)
    return output


def load_model_yearly(
    model_files, spec, args, lat2d, lon2d, lat1d, lon1d, land_mask
):
    output = {}
    for experiment, paths in model_files.items():
        output[experiment] = {"ctl": {}, "exp": {}}
        with Dataset(paths["ctl"]) as nc_ctl, Dataset(paths["exp"]) as nc_exp:
            for season_name, months in SEASONS:
                for role, nc_obj in (("ctl", nc_ctl), ("exp", nc_exp)):
                    native = read_yearly_season(
                        nc_obj,
                        spec["variable"],
                        args.nyears,
                        months,
                        spec["factor"],
                    )
                    target = np.stack(
                        [
                            regrid_rcm_to_regular(
                                field, lat2d, lon2d, lat1d, lon1d
                            )
                            for field in native
                        ]
                    )
                    output[experiment][role][season_name] = np.where(
                        land_mask[None, :, :], target, np.nan
                    )
    return output


def add_boundaries(ax, shapefile_dir):
    for filename, linewidth in (
        ("province.shp", 0.4),
        ("china.shp", 0.6),
        ("south_china_sea.shp", 0.8),
    ):
        path = Path(shapefile_dir) / filename
        if not path.is_file():
            continue
        try:
            reader = shpreader.Reader(str(path))
            ax.add_geometries(
                reader.geometries(),
                crs=ccrs.PlateCarree(),
                facecolor="none",
                edgecolor="black",
                linewidth=linewidth,
                zorder=3,
            )
        except Exception as exc:
            print(f"[WARN] Could not draw {path}: {exc}")


def add_gridlines(ax, show_left, show_bottom):
    gridlines = ax.gridlines(
        crs=ccrs.PlateCarree(),
        draw_labels=True,
        x_inline=False,
        y_inline=False,
        linewidth=0.5,
        color="gray",
        alpha=0.45,
        linestyle="--",
        zorder=2,
    )
    gridlines.xlocator = mticker.FixedLocator(np.arange(70, 135, 10))
    gridlines.ylocator = mticker.FixedLocator(np.arange(10, 60, 10))
    gridlines.top_labels = False
    gridlines.right_labels = False
    gridlines.left_labels = show_left
    gridlines.bottom_labels = show_bottom
    gridlines.xlabel_style = {"size": 7.5}
    gridlines.ylabel_style = {"size": 7.5}


def align_model_years_to_reference(reference, model, years, spec):
    """Apply the reference's season-year availability to all simulations."""
    for season_name, _ in SEASONS:
        available = np.any(
            np.isfinite(reference[season_name]), axis=(1, 2)
        )
        omitted = [
            year for year, keep in zip(years, available) if not keep
        ]
        if omitted:
            print(
                f"[WARN] {spec['title']} {season_name}: omitting years "
                f"{omitted} from Reference/CTL/EXP comparisons"
            )
            for experiment in ("OFF", "CPL"):
                for role in ("ctl", "exp"):
                    model[experiment][role][season_name][~available] = np.nan


def unified_mask(reference, ctl, exp, land_mask):
    available = np.any(np.isfinite(reference), axis=(1, 2))
    if not np.any(available):
        return np.zeros_like(land_mask, dtype=bool)
    return (
        land_mask
        & np.all(np.isfinite(reference[available]), axis=0)
        & np.all(np.isfinite(ctl[available]), axis=0)
        & np.all(np.isfinite(exp[available]), axis=0)
    )


def weighted_mean(field, lat_grid, mask):
    valid = mask & np.isfinite(field)
    if np.count_nonzero(valid) == 0:
        return np.nan
    weights = np.cos(np.deg2rad(lat_grid[valid]))
    return np.average(field[valid], weights=weights)


def weighted_mean_ci(field, lat_grid, mask, z_value=2.575829):
    valid = mask & np.isfinite(field)
    if np.count_nonzero(valid) < 2:
        return np.nan, np.nan, np.nan
    values = field[valid]
    weights = np.cos(np.deg2rad(lat_grid[valid]))
    weights = weights / np.sum(weights)
    mean = np.sum(weights * values)
    denominator = 1.0 - np.sum(weights**2)
    variance = (
        np.sum(weights * (values - mean) ** 2) / denominator
        if denominator > 0
        else np.nan
    )
    effective_n = 1.0 / np.sum(weights**2)
    half_width = z_value * np.sqrt(variance / effective_n)
    return mean, mean - half_width, mean + half_width


def weighted_correlation(first, second, lat_grid, mask):
    valid = mask & np.isfinite(first) & np.isfinite(second)
    if np.count_nonzero(valid) < 3:
        return np.nan
    x = first[valid]
    y = second[valid]
    weights = np.cos(np.deg2rad(lat_grid[valid]))
    weights = weights / np.sum(weights)
    x_anomaly = x - np.sum(weights * x)
    y_anomaly = y - np.sum(weights * y)
    denominator = np.sqrt(
        np.sum(weights * x_anomaly**2)
        * np.sum(weights * y_anomaly**2)
    )
    if denominator == 0:
        return np.nan
    return np.sum(weights * x_anomaly * y_anomaly) / denominator


def field_metrics(simulation, reference, lat_grid, mask):
    valid = mask & np.isfinite(simulation) & np.isfinite(reference)
    if np.count_nonzero(valid) < 3:
        return np.nan, np.nan, np.nan
    difference = simulation - reference
    bias = weighted_mean(difference, lat_grid, valid)
    rmse = np.sqrt(weighted_mean(difference**2, lat_grid, valid))
    correlation = weighted_correlation(
        simulation, reference, lat_grid, valid
    )
    return bias, rmse, correlation


def time_metrics(simulation, reference):
    valid = np.isfinite(simulation) & np.isfinite(reference)
    if np.count_nonzero(valid) < 3:
        return {
            name: np.nan
            for name in ("bias", "rmse", "correlation", "kge")
        }
    sim = simulation[valid]
    ref = reference[valid]
    bias = np.mean(sim - ref)
    rmse = np.sqrt(np.mean((sim - ref) ** 2))
    correlation = np.corrcoef(sim, ref)[0, 1]
    ref_std = np.std(ref, ddof=1)
    ref_mean = np.mean(ref)
    alpha = np.std(sim, ddof=1) / ref_std if ref_std else np.nan
    beta = np.mean(sim) / ref_mean if ref_mean else np.nan
    kge = (
        1.0
        - np.sqrt(
            (correlation - 1.0) ** 2
            + (alpha - 1.0) ** 2
            + (beta - 1.0) ** 2
        )
        if np.all(np.isfinite((correlation, alpha, beta)))
        else np.nan
    )
    return {
        "bias": bias,
        "rmse": rmse,
        "correlation": correlation,
        "kge": kge,
    }


def robust_state_limits(fields, mask, low, high, nonnegative):
    values = np.concatenate(
        [
            field[:, mask].ravel()
            if field.ndim == 3
            else field[mask].ravel()
            for field in fields
        ]
    )
    values = values[np.isfinite(values)]
    if values.size == 0:
        raise ValueError("No finite values available for state color limits")
    minimum, maximum = np.nanpercentile(values, (low, high))
    if nonnegative:
        minimum = 0.0
    if not maximum > minimum:
        maximum = minimum + max(abs(minimum) * 0.1, 1.0e-6)
    return float(minimum), float(maximum)


def robust_symmetric_limit(fields, mask, quantile, fallback):
    values = np.concatenate([np.abs(field[mask]).ravel() for field in fields])
    values = values[np.isfinite(values)]
    limit = np.nanpercentile(values, quantile) if values.size else np.nan
    if not np.isfinite(limit) or limit <= 0:
        limit = fallback
    return float(limit)


def common_plot_data(reference, model, land_mask):
    masks = {}
    for experiment in ("OFF", "CPL"):
        masks[experiment] = {}
        for season_name, _ in SEASONS:
            masks[experiment][season_name] = unified_mask(
                reference[season_name],
                model[experiment]["ctl"][season_name],
                model[experiment]["exp"][season_name],
                land_mask,
            )
    return masks


def draw_map_panel(
    ax,
    field,
    levels,
    cmap,
    lon_grid,
    lat_grid,
    shapefile_dir,
    show_left,
    show_bottom,
):
    ax.set_extent(PLOT_EXTENT, crs=ccrs.PlateCarree())
    add_boundaries(ax, shapefile_dir)
    add_gridlines(ax, show_left, show_bottom)
    return ax.contourf(
        lon_grid,
        lat_grid,
        field,
        levels=levels,
        cmap=cmap,
        extend="both",
        transform=ccrs.PlateCarree(),
        zorder=1,
        antialiased=False,
    )


def plot_spatial_evaluation(
    spec,
    reference,
    model,
    masks,
    args,
    lat_grid,
    lon_grid,
    land_mask,
    product,
):
    state_fields = []
    bias_fields = []
    response_fields = []
    for season_name, _ in SEASONS:
        state_fields.append(reference[season_name])
        for experiment in ("OFF", "CPL"):
            state_fields.extend(
                [
                    model[experiment]["ctl"][season_name],
                    model[experiment]["exp"][season_name],
                ]
            )
            mask = masks[experiment][season_name]
            ctl_mean = np.nanmean(
                model[experiment]["ctl"][season_name], axis=0
            )
            exp_mean = np.nanmean(
                model[experiment]["exp"][season_name], axis=0
            )
            ref_mean = np.nanmean(reference[season_name], axis=0)
            bias_fields.extend(
                [
                    np.where(mask, ctl_mean - ref_mean, np.nan),
                    np.where(mask, exp_mean - ref_mean, np.nan),
                ]
            )
            response_fields.append(
                np.where(mask, exp_mean - ctl_mean, np.nan)
            )
    state_min, state_max = robust_state_limits(
        state_fields,
        land_mask,
        args.state_quantiles[0],
        args.state_quantiles[1],
        spec.get("nonnegative", False),
    )
    bias_limit = robust_symmetric_limit(
        bias_fields,
        land_mask,
        args.difference_quantile,
        spec.get("difference_fallback", 1.0),
    )
    response_limit = robust_symmetric_limit(
        response_fields,
        land_mask,
        args.difference_quantile,
        spec.get("difference_fallback", 1.0),
    )
    state_levels = np.linspace(state_min, state_max, 21)
    bias_levels = np.linspace(-bias_limit, bias_limit, 21)
    response_levels = np.linspace(-response_limit, response_limit, 21)
    print(
        f"Color limits {spec['title']}: state={state_min:.4g}.."
        f"{state_max:.4g}, bias=+/-{bias_limit:.4g}, response=+/-"
        f"{response_limit:.4g} {spec['unit']}"
    )
    projection = ccrs.LambertConformal(
        central_longitude=PROJ_CENTRAL_LON,
        central_latitude=PROJ_CENTRAL_LAT,
        standard_parallels=PROJ_STD_PARALLELS,
    )
    for experiment in ("OFF", "CPL"):
        fig, axes = plt.subplots(
            2,
            6,
            figsize=(19.5, 7.4),
            subplot_kw={"projection": projection},
            squeeze=False,
        )
        plt.subplots_adjust(
            left=0.055,
            right=0.985,
            bottom=0.15,
            top=0.88,
            wspace=0.035,
            hspace=0.10,
        )
        headers = (
            "Reference",
            "CTL",
            "EXP",
            "CTL - Reference",
            "EXP - Reference",
            "EXP - CTL",
        )
        state_mappable = None
        bias_mappable = None
        response_mappable = None
        for row, (season_name, _) in enumerate(SEASONS):
            mask = masks[experiment][season_name]
            ref_yearly = np.where(
                mask[None], reference[season_name], np.nan
            )
            ctl_yearly = np.where(
                mask[None], model[experiment]["ctl"][season_name], np.nan
            )
            exp_yearly = np.where(
                mask[None], model[experiment]["exp"][season_name], np.nan
            )
            ref_mean = np.nanmean(ref_yearly, axis=0)
            ctl_mean = np.nanmean(ctl_yearly, axis=0)
            exp_mean = np.nanmean(exp_yearly, axis=0)
            ctl_bias = ctl_mean - ref_mean
            exp_bias = exp_mean - ref_mean
            response = exp_mean - ctl_mean
            _, p_value = ttest_rel(
                exp_yearly, ctl_yearly, axis=0, nan_policy="omit"
            )
            fields = (
                ref_mean,
                ctl_mean,
                exp_mean,
                ctl_bias,
                exp_bias,
                response,
            )
            for col, field in enumerate(fields):
                plot_field = field
                if col == 5 and args.significance_style == "mask":
                    plot_field = np.where(
                        p_value < args.p_threshold, field, np.nan
                    )
                if col <= 2:
                    levels = state_levels
                    cmap = spec.get("state_cmap", "YlGnBu")
                elif col <= 4:
                    levels = bias_levels
                    cmap = spec.get("difference_cmap", "BrBG")
                else:
                    levels = response_levels
                    cmap = spec.get("difference_cmap", "BrBG")
                mappable = draw_map_panel(
                    axes[row, col],
                    plot_field,
                    levels,
                    cmap,
                    lon_grid,
                    lat_grid,
                    args.shapefile_dir,
                    show_left=(col == 0),
                    show_bottom=(row == 1),
                )
                if col <= 2:
                    state_mappable = mappable
                elif col <= 4:
                    bias_mappable = mappable
                else:
                    response_mappable = mappable
                if col in (3, 4):
                    simulation = ctl_mean if col == 3 else exp_mean
                    bias, rmse, correlation = field_metrics(
                        simulation, ref_mean, lat_grid, mask
                    )
                    axes[row, col].text(
                        0.5,
                        0.02,
                        f"B={bias:.3g}  RMSE={rmse:.3g}  R={correlation:.2f}",
                        transform=axes[row, col].transAxes,
                        ha="center",
                        va="bottom",
                        fontsize=6.8,
                        bbox={
                            "facecolor": "white",
                            "alpha": 0.72,
                            "edgecolor": "none",
                            "pad": 1,
                        },
                        zorder=5,
                    )
                if col == 5 and args.significance_style == "stipple":
                    stride = args.sig_stride
                    significant = (
                        mask
                        & np.isfinite(p_value)
                        & (p_value < args.p_threshold)
                    )
                    sampled = significant[::stride, ::stride]
                    axes[row, col].scatter(
                        lon_grid[::stride, ::stride][sampled],
                        lat_grid[::stride, ::stride][sampled],
                        s=3,
                        color="0.35",
                        edgecolors="none",
                        transform=ccrs.PlateCarree(),
                        rasterized=True,
                        zorder=2,
                    )
            axes[row, 0].text(
                -0.18,
                0.5,
                "MAM" if season_name == "Spring" else "JJA",
                transform=axes[row, 0].transAxes,
                rotation=90,
                fontsize=11,
                va="center",
                ha="center",
            )
        # Cartopy can suppress GeoAxes titles and crop annotations in the
        # left-most column.  Figure-level text keeps all six headers and all
        # twelve panel labels visible in vector PDF output.
        fig.canvas.draw()
        for col, header in enumerate(headers):
            position = axes[0, col].get_position()
            fig.text(
                position.x0 + position.width / 2,
                position.y1 + 0.012,
                header,
                ha="center",
                va="bottom",
                fontsize=10.5,
            )
        for row in range(2):
            for col in range(6):
                position = axes[row, col].get_position()
                fig.text(
                    position.x0 + 0.012 * position.width,
                    position.y1 - 0.035 * position.height,
                    f"({chr(97 + row * 6 + col)})",
                    ha="left",
                    va="top",
                    fontsize=9.5,
                    bbox={
                        "facecolor": "white",
                        "alpha": 0.7,
                        "edgecolor": "none",
                        "pad": 1,
                    },
                    zorder=5,
                )
        period_end = args.start_year + args.nyears - 1
        fig.suptitle(
            f"{spec['title']} - {experiment} evaluation "
            f"({args.start_year}-{period_end})",
            fontsize=14,
            y=0.96,
        )
        state_cax = fig.add_axes([0.08, 0.075, 0.39, 0.022])
        bias_cax = fig.add_axes([0.54, 0.075, 0.25, 0.022])
        response_cax = fig.add_axes([0.84, 0.075, 0.13, 0.022])
        state_cb = fig.colorbar(
            state_mappable, cax=state_cax, orientation="horizontal"
        )
        bias_cb = fig.colorbar(
            bias_mappable, cax=bias_cax, orientation="horizontal"
        )
        response_cb = fig.colorbar(
            response_mappable, cax=response_cax, orientation="horizontal"
        )
        for colorbar, label in (
            (state_cb, f"State ({spec['unit']})"),
            (bias_cb, f"Bias ({spec['unit']})"),
            (response_cb, f"Response ({spec['unit']})"),
        ):
            colorbar.ax.tick_params(labelsize=7.5)
            colorbar.set_label(label, fontsize=8.5, labelpad=2)
        state_cb.locator = mticker.MaxNLocator(nbins=6)
        bias_cb.locator = mticker.MaxNLocator(nbins=5)
        response_cb.locator = mticker.MaxNLocator(nbins=3)
        state_cb.update_ticks()
        bias_cb.update_ticks()
        response_cb.update_ticks()
        fig.text(
            0.5,
            0.018,
            (
                f"Reference: {product}; unified China-land mask. State limits: "
                f"P{args.state_quantiles[0]:g}-P{args.state_quantiles[1]:g}; "
                f"bias/response: P{args.difference_quantile:g}(|x|). "
                f"Dots: paired two-sided t test, p < {args.p_threshold:g}."
            ),
            ha="center",
            fontsize=8,
        )
        output = (
            args.output_dir
            / f"{spec['output_stem']}_{experiment.lower()}_spatial.pdf"
        )
        fig.savefig(output, dpi=args.dpi, bbox_inches="tight")
        plt.close(fig)
        print(f"Saved: {output}")


def regional_series(yearly, lat_grid, mask):
    values, lower, upper = [], [], []
    for field in yearly:
        mean, low, high = weighted_mean_ci(field, lat_grid, mask)
        values.append(mean)
        lower.append(low)
        upper.append(high)
    return np.asarray(values), np.asarray(lower), np.asarray(upper)


def write_csv(path, fieldnames, rows):
    with Path(path).open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    print(f"Saved: {path}")


def plot_regional_timeseries(
    spec, reference, model, masks, args, lat_grid, product, years
):
    fig, axes = plt.subplots(
        2, 2, figsize=(13.2, 7.6), sharex=True, squeeze=False
    )
    colors = {
        "Reference": "black",
        "CTL": "#2166ac",
        "EXP": "#d73027",
    }
    csv_rows = []
    metric_rows = []
    for row, experiment in enumerate(("OFF", "CPL")):
        for col, (season_name, _) in enumerate(SEASONS):
            ax = axes[row, col]
            mask = masks[experiment][season_name]
            data = {
                "Reference": regional_series(
                    reference[season_name], lat_grid, mask
                ),
                "CTL": regional_series(
                    model[experiment]["ctl"][season_name], lat_grid, mask
                ),
                "EXP": regional_series(
                    model[experiment]["exp"][season_name], lat_grid, mask
                ),
            }
            for label, (values, lower, upper) in data.items():
                ax.plot(
                    years,
                    values,
                    color=colors[label],
                    linewidth=1.6,
                    marker="o",
                    markersize=3.2,
                    label=label,
                )
                ax.fill_between(
                    years,
                    lower,
                    upper,
                    color=colors[label],
                    alpha=0.11,
                    linewidth=0,
                )
            ctl_metrics = time_metrics(
                data["CTL"][0], data["Reference"][0]
            )
            exp_metrics = time_metrics(
                data["EXP"][0], data["Reference"][0]
            )
            ax.text(
                0.015,
                0.975,
                (
                    f"CTL  B={ctl_metrics['bias']:.3g}  "
                    f"RMSE={ctl_metrics['rmse']:.3g}  "
                    f"r={ctl_metrics['correlation']:.2f}  "
                    f"KGE={ctl_metrics['kge']:.2f}\n"
                    f"EXP  B={exp_metrics['bias']:.3g}  "
                    f"RMSE={exp_metrics['rmse']:.3g}  "
                    f"r={exp_metrics['correlation']:.2f}  "
                    f"KGE={exp_metrics['kge']:.2f}"
                ),
                transform=ax.transAxes,
                va="top",
                fontsize=8.2,
                bbox={
                    "facecolor": "white",
                    "alpha": 0.78,
                    "edgecolor": "0.8",
                    "pad": 2,
                },
            )
            season_label = "MAM" if season_name == "Spring" else "JJA"
            ax.set_title(f"{experiment} | {season_label}", fontsize=11)
            ax.grid(True, color="0.88", linewidth=0.7)
            ax.set_xlim(years[0], years[-1])
            ax.xaxis.set_major_locator(
                mticker.MaxNLocator(integer=True, nbins=6)
            )
            if col == 0:
                ax.set_ylabel(spec["unit"])
            if row == 1:
                ax.set_xlabel("Year")
            for index, year in enumerate(years):
                csv_rows.append(
                    {
                        "product": product,
                        "experiment": experiment,
                        "season": season_name,
                        "year": year,
                        "reference": data["Reference"][0][index],
                        "ctl": data["CTL"][0][index],
                        "exp": data["EXP"][0][index],
                        "reference_ci99_low": data["Reference"][1][index],
                        "reference_ci99_high": data["Reference"][2][index],
                        "ctl_ci99_low": data["CTL"][1][index],
                        "ctl_ci99_high": data["CTL"][2][index],
                        "exp_ci99_low": data["EXP"][1][index],
                        "exp_ci99_high": data["EXP"][2][index],
                    }
                )
            for role, metrics in (
                ("CTL", ctl_metrics),
                ("EXP", exp_metrics),
            ):
                metric_rows.append(
                    {
                        "product": product,
                        "experiment": experiment,
                        "season": season_name,
                        "simulation": role,
                        **metrics,
                    }
                )
    axes[0, 0].legend(
        frameon=False, ncol=3, loc="lower left", fontsize=9
    )
    fig.suptitle(
        f"{spec['title']}: China-land seasonal time series",
        fontsize=14,
        y=0.98,
    )
    fig.text(
        0.5,
        0.012,
        (
            f"Reference: {product}. Shading is a 99% grid-sampling confidence "
            "interval of the area-weighted mean; spatial autocorrelation is "
            "not corrected."
        ),
        ha="center",
        fontsize=8,
    )
    fig.tight_layout(rect=(0.02, 0.05, 1, 0.95))
    output = (
        args.output_dir / f"{spec['output_stem']}_regional_timeseries.pdf"
    )
    fig.savefig(output, dpi=args.dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {output}")
    write_csv(
        args.output_dir
        / f"{spec['output_stem']}_regional_timeseries.csv",
        list(csv_rows[0]),
        csv_rows,
    )
    write_csv(
        args.output_dir / f"{spec['output_stem']}_regional_metrics.csv",
        list(metric_rows[0]),
        metric_rows,
    )


def spatial_acc_timeseries(simulation, reference, lat_grid, mask):
    simulation_anomaly = simulation - np.nanmean(simulation, axis=0)
    reference_anomaly = reference - np.nanmean(reference, axis=0)
    return np.asarray(
        [
            weighted_correlation(
                simulation_anomaly[index],
                reference_anomaly[index],
                lat_grid,
                mask,
            )
            for index in range(simulation.shape[0])
        ]
    )


def plot_acc_timeseries(
    spec, reference, model, masks, args, lat_grid, product, years
):
    fig = plt.figure(figsize=(11.4, 12.2))
    grid = fig.add_gridspec(
        4, 2, width_ratios=(6.2, 1.35), hspace=0.28, wspace=0.04
    )
    csv_rows = []
    plot_row = 0
    legend_handles = None
    for experiment in ("OFF", "CPL"):
        for season_name, _ in SEASONS:
            mask = masks[experiment][season_name]
            ctl_acc = spatial_acc_timeseries(
                model[experiment]["ctl"][season_name],
                reference[season_name],
                lat_grid,
                mask,
            )
            exp_acc = spatial_acc_timeseries(
                model[experiment]["exp"][season_name],
                reference[season_name],
                lat_grid,
                mask,
            )
            ax = fig.add_subplot(grid[plot_row, 0])
            ctl_line = ax.plot(
                years,
                ctl_acc,
                color="#2166ac",
                marker="o",
                markersize=3.2,
                linewidth=1.5,
                label="CTL",
            )[0]
            exp_line = ax.plot(
                years,
                exp_acc,
                color="#d73027",
                marker="o",
                markersize=3.2,
                linewidth=1.5,
                label="EXP",
            )[0]
            legend_handles = (ctl_line, exp_line)
            ax.axhline(0.0, color="0.25", linewidth=0.8)
            ax.axhline(
                args.acc_threshold,
                color="0.45",
                linestyle="--",
                linewidth=1.0,
            )
            ax.set_ylim(-1.0, 1.0)
            ax.set_xlim(years[0], years[-1])
            ax.xaxis.set_major_locator(
                mticker.MaxNLocator(integer=True, nbins=7)
            )
            ax.grid(True, color="0.9", linewidth=0.6)
            season_label = "MAM" if season_name == "Spring" else "JJA"
            ax.set_ylabel(
                f"Spatial ACC\n{experiment} | {season_label}", fontsize=9.5
            )
            if plot_row < 3:
                ax.tick_params(labelbottom=False)
            else:
                ax.set_xlabel("Year")
            summary = fig.add_subplot(grid[plot_row, 1], sharey=ax)
            box = summary.boxplot(
                [
                    ctl_acc[np.isfinite(ctl_acc)],
                    exp_acc[np.isfinite(exp_acc)],
                ],
                positions=(1, 2),
                widths=0.58,
                whis=(0, 100),
                showmeans=True,
                patch_artist=True,
                meanprops={
                    "marker": "x",
                    "markeredgecolor": "white",
                    "markersize": 5,
                },
                medianprops={"color": "white", "linewidth": 1.2},
            )
            for patch, color in zip(
                box["boxes"], ("#2166ac", "#d73027")
            ):
                patch.set_facecolor(color)
                patch.set_alpha(0.9)
            valid = np.isfinite(ctl_acc) & np.isfinite(exp_acc)
            p_value = (
                ttest_rel(exp_acc[valid], ctl_acc[valid]).pvalue
                if np.count_nonzero(valid) >= 3
                else np.nan
            )
            summary.set_title(f"p={p_value:.3f}", fontsize=8.5, pad=2)
            summary.set_xticks((1, 2), ("CTL", "EXP"), fontsize=7.5)
            summary.tick_params(axis="x", pad=1)
            summary.tick_params(axis="y", labelleft=False)
            summary.grid(True, axis="y", color="0.9", linewidth=0.6)
            for index, year in enumerate(years):
                csv_rows.append(
                    {
                        "product": product,
                        "experiment": experiment,
                        "season": season_name,
                        "year": year,
                        "ctl_spatial_acc": ctl_acc[index],
                        "exp_spatial_acc": exp_acc[index],
                        "paired_p_value": p_value,
                    }
                )
            plot_row += 1
    fig.legend(
        legend_handles,
        ("CTL", "EXP"),
        loc="upper center",
        bbox_to_anchor=(0.46, 0.948),
        ncol=2,
        frameon=False,
    )
    fig.suptitle(
        f"{spec['title']}: interannual spatial-anomaly skill",
        fontsize=14,
        y=0.985,
    )
    period_end = args.start_year + args.nyears - 1
    fig.text(
        0.5,
        0.02,
        (
            f"Reference: {product}. ACC is the area-weighted spatial "
            "correlation of yearly anomalies after removing each source's "
            f"{args.start_year}-{period_end} climatology; dashed line="
            f"{args.acc_threshold:g}. "
            "Boxes show IQR, median, mean (x), and min-max; p values use a "
            "paired two-sided t test."
        ),
        ha="center",
        fontsize=8,
        wrap=True,
    )
    output = (
        args.output_dir
        / f"{spec['output_stem']}_spatial_acc_timeseries.pdf"
    )
    fig.savefig(output, dpi=args.dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {output}")
    write_csv(
        args.output_dir
        / f"{spec['output_stem']}_spatial_acc_timeseries.csv",
        list(csv_rows[0]),
        csv_rows,
    )


def model_response_masks(model, land_mask):
    masks = {}
    for experiment in ("OFF", "CPL"):
        masks[experiment] = {}
        for season_name, _ in SEASONS:
            masks[experiment][season_name] = (
                land_mask
                & np.all(
                    np.isfinite(model[experiment]["ctl"][season_name]),
                    axis=0,
                )
                & np.all(
                    np.isfinite(model[experiment]["exp"][season_name]),
                    axis=0,
                )
            )
    return masks


def plot_response_only_spatial(
    spec, model, masks, args, lat_grid, lon_grid, land_mask
):
    state_fields = []
    response_fields = []
    for experiment in ("OFF", "CPL"):
        for season_name, _ in SEASONS:
            state_fields.extend(
                [
                    model[experiment]["ctl"][season_name],
                    model[experiment]["exp"][season_name],
                ]
            )
            response_fields.append(
                np.nanmean(
                    model[experiment]["exp"][season_name]
                    - model[experiment]["ctl"][season_name],
                    axis=0,
                )
            )
    state_min, state_max = robust_state_limits(
        state_fields,
        land_mask,
        args.state_quantiles[0],
        args.state_quantiles[1],
        spec.get("nonnegative", False),
    )
    response_limit = robust_symmetric_limit(
        response_fields,
        land_mask,
        args.difference_quantile,
        spec.get("difference_fallback", 1.0),
    )
    state_levels = np.linspace(state_min, state_max, 21)
    response_levels = np.linspace(-response_limit, response_limit, 21)
    print(
        f"Color limits {spec['title']}: state={state_min:.4g}.."
        f"{state_max:.4g}, response=+/-{response_limit:.4g} "
        f"{spec['unit']}"
    )
    projection = ccrs.LambertConformal(
        central_longitude=PROJ_CENTRAL_LON,
        central_latitude=PROJ_CENTRAL_LAT,
        standard_parallels=PROJ_STD_PARALLELS,
    )
    for experiment in ("OFF", "CPL"):
        fig, axes = plt.subplots(
            2,
            3,
            figsize=(11.5, 7.4),
            subplot_kw={"projection": projection},
            squeeze=False,
        )
        plt.subplots_adjust(
            left=0.075,
            right=0.985,
            bottom=0.17,
            top=0.88,
            wspace=0.045,
            hspace=0.10,
        )
        headers = ("CTL", "EXP", "EXP - CTL")
        state_mappable = None
        response_mappable = None
        for row, (season_name, _) in enumerate(SEASONS):
            mask = masks[experiment][season_name]
            ctl_yearly = np.where(
                mask[None], model[experiment]["ctl"][season_name], np.nan
            )
            exp_yearly = np.where(
                mask[None], model[experiment]["exp"][season_name], np.nan
            )
            ctl_mean = np.nanmean(ctl_yearly, axis=0)
            exp_mean = np.nanmean(exp_yearly, axis=0)
            response = exp_mean - ctl_mean
            _, p_value = ttest_rel(
                exp_yearly, ctl_yearly, axis=0, nan_policy="omit"
            )
            for col, field in enumerate((ctl_mean, exp_mean, response)):
                plot_field = field
                if col == 2 and args.significance_style == "mask":
                    plot_field = np.where(
                        p_value < args.p_threshold, field, np.nan
                    )
                if col < 2:
                    levels = state_levels
                    cmap = spec.get("state_cmap", "YlGnBu")
                else:
                    levels = response_levels
                    cmap = spec.get("difference_cmap", "BrBG")
                mappable = draw_map_panel(
                    axes[row, col],
                    plot_field,
                    levels,
                    cmap,
                    lon_grid,
                    lat_grid,
                    args.shapefile_dir,
                    show_left=(col == 0),
                    show_bottom=(row == 1),
                )
                if col < 2:
                    state_mappable = mappable
                else:
                    response_mappable = mappable
                    area_change = weighted_mean(
                        response, lat_grid, mask
                    )
                    axes[row, col].text(
                        0.5,
                        0.02,
                        f"Area mean change = {area_change:.3g}",
                        transform=axes[row, col].transAxes,
                        ha="center",
                        va="bottom",
                        fontsize=7.2,
                        bbox={
                            "facecolor": "white",
                            "alpha": 0.72,
                            "edgecolor": "none",
                            "pad": 1,
                        },
                        zorder=5,
                    )
                    if args.significance_style == "stipple":
                        stride = args.sig_stride
                        significant = (
                            mask
                            & np.isfinite(p_value)
                            & (p_value < args.p_threshold)
                        )
                        sampled = significant[::stride, ::stride]
                        axes[row, col].scatter(
                            lon_grid[::stride, ::stride][sampled],
                            lat_grid[::stride, ::stride][sampled],
                            s=3,
                            color="0.35",
                            edgecolors="none",
                            transform=ccrs.PlateCarree(),
                            rasterized=True,
                            zorder=2,
                        )
            axes[row, 0].text(
                -0.18,
                0.5,
                "MAM" if season_name == "Spring" else "JJA",
                transform=axes[row, 0].transAxes,
                rotation=90,
                fontsize=11,
                va="center",
                ha="center",
            )
        fig.canvas.draw()
        for col, header in enumerate(headers):
            position = axes[0, col].get_position()
            fig.text(
                position.x0 + position.width / 2,
                position.y1 + 0.012,
                header,
                ha="center",
                va="bottom",
                fontsize=10.5,
            )
        for row in range(2):
            for col in range(3):
                position = axes[row, col].get_position()
                fig.text(
                    position.x0 + 0.012 * position.width,
                    position.y1 - 0.035 * position.height,
                    f"({chr(97 + row * 3 + col)})",
                    ha="left",
                    va="top",
                    fontsize=9.5,
                    bbox={
                        "facecolor": "white",
                        "alpha": 0.7,
                        "edgecolor": "none",
                        "pad": 1,
                    },
                    zorder=5,
                )
        period_end = args.start_year + args.nyears - 1
        fig.suptitle(
            f"{spec['title']} - {experiment} process response "
            f"({args.start_year}-{period_end})",
            fontsize=14,
            y=0.96,
        )
        state_cax = fig.add_axes([0.12, 0.087, 0.48, 0.022])
        response_cax = fig.add_axes([0.69, 0.087, 0.25, 0.022])
        state_cb = fig.colorbar(
            state_mappable, cax=state_cax, orientation="horizontal"
        )
        response_cb = fig.colorbar(
            response_mappable,
            cax=response_cax,
            orientation="horizontal",
        )
        state_cb.set_label(f"State ({spec['unit']})", fontsize=8.5)
        response_cb.set_label(
            f"Response ({spec['unit']})", fontsize=8.5
        )
        state_cb.set_ticks(np.linspace(state_min, state_max, 5))
        response_cb.set_ticks(
            np.linspace(-response_limit, response_limit, 5)
        )
        for colorbar in (state_cb, response_cb):
            colorbar.ax.tick_params(labelsize=7.5)
        fig.text(
            0.5,
            0.012,
            (
                "Process diagnostic only; no reference evaluation. State "
                f"limits: P{args.state_quantiles[0]:g}-"
                f"P{args.state_quantiles[1]:g}; response: "
                f"P{args.difference_quantile:g}(|x|). Dots: paired "
                f"two-sided t test, p < {args.p_threshold:g}."
            ),
            ha="center",
            fontsize=8,
        )
        output = (
            args.output_dir
            / f"{spec['output_stem']}_{experiment.lower()}_spatial.pdf"
        )
        fig.savefig(output, dpi=args.dpi, bbox_inches="tight")
        plt.close(fig)
        print(f"Saved: {output}")


def plot_response_only_timeseries(
    spec, model, masks, args, lat_grid, years
):
    fig, axes = plt.subplots(
        2, 2, figsize=(13.2, 7.6), sharex=True, squeeze=False
    )
    csv_rows = []
    for row, experiment in enumerate(("OFF", "CPL")):
        for col, (season_name, _) in enumerate(SEASONS):
            ax = axes[row, col]
            mask = masks[experiment][season_name]
            ctl_yearly = np.where(
                mask[None], model[experiment]["ctl"][season_name], np.nan
            )
            exp_yearly = np.where(
                mask[None], model[experiment]["exp"][season_name], np.nan
            )
            response_yearly = exp_yearly - ctl_yearly
            response, lower, upper = regional_series(
                response_yearly, lat_grid, mask
            )
            ctl = regional_series(ctl_yearly, lat_grid, mask)[0]
            exp = regional_series(exp_yearly, lat_grid, mask)[0]
            _, p_value = ttest_rel(exp, ctl, nan_policy="omit")
            ax.fill_between(
                years, lower, upper, color="#d73027", alpha=0.16
            )
            ax.plot(
                years,
                response,
                color="#d73027",
                marker="o",
                markersize=3.2,
                linewidth=1.4,
                label="EXP - CTL",
            )
            ax.axhline(0.0, color="0.25", linewidth=0.8, linestyle="--")
            season_label = "MAM" if season_name == "Spring" else "JJA"
            ax.set_title(f"{experiment} | {season_label}", fontsize=11)
            ax.text(
                0.015,
                0.965,
                f"Mean change={np.nanmean(response):.3g}; p={p_value:.3g}",
                transform=ax.transAxes,
                va="top",
                fontsize=8.8,
                bbox={
                    "facecolor": "white",
                    "alpha": 0.78,
                    "edgecolor": "0.75",
                    "pad": 2,
                },
            )
            ax.grid(True, color="0.88", linewidth=0.6)
            ax.set_xlim(years[0], years[-1])
            ax.xaxis.set_major_locator(
                mticker.MaxNLocator(integer=True, nbins=6)
            )
            if col == 0:
                ax.set_ylabel(spec["unit"])
            if row == 1:
                ax.set_xlabel("Year")
            for index, year in enumerate(years):
                csv_rows.append(
                    {
                        "experiment": experiment,
                        "season": season_name,
                        "year": year,
                        "ctl": ctl[index],
                        "exp": exp[index],
                        "exp_minus_ctl": response[index],
                        "response_ci99_low": lower[index],
                        "response_ci99_high": upper[index],
                    }
                )
    axes[0, 0].legend(frameon=False, loc="best", fontsize=9)
    fig.suptitle(
        f"{spec['title']}: China-land gravel-response time series",
        fontsize=14,
    )
    fig.text(
        0.5,
        0.015,
        (
            "Process diagnostic only; no reference evaluation. Shading is a "
            "99% grid-sampling confidence interval of the area-weighted "
            "EXP-CTL response; p values use a paired test across years."
        ),
        ha="center",
        fontsize=8,
    )
    fig.tight_layout(rect=(0.02, 0.05, 1, 0.95))
    output = (
        args.output_dir
        / f"{spec['output_stem']}_response_timeseries.pdf"
    )
    fig.savefig(output, dpi=args.dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {output}")
    write_csv(
        args.output_dir
        / f"{spec['output_stem']}_response_timeseries.csv",
        list(csv_rows[0]),
        csv_rows,
    )


def validate_arguments(args):
    if args.nyears < 3:
        raise ValueError("--nyears must be at least 3")
    if not 0.0 < args.p_threshold < 1.0:
        raise ValueError("--p-threshold must be between 0 and 1")
    if args.sig_stride < 1:
        raise ValueError("--sig-stride must be at least 1")
    low, high = args.state_quantiles
    if not 0 <= low < high <= 100:
        raise ValueError(
            "--state-quantiles must satisfy 0 <= LOW < HIGH <= 100"
        )
    if not 0 < args.difference_quantile <= 100:
        raise ValueError("--difference-quantile must be in (0, 100]")


def run_category(variable_specs, description):
    parser = build_parser(description, variable_specs)
    args = parser.parse_args()
    validate_arguments(args)
    model_files = resolve_model_files(args)
    validate_model_inputs(model_files, variable_specs, args.nyears)
    years = list(range(args.start_year, args.start_year + args.nyears))
    reference_files = {}
    for spec in variable_specs:
        files = resolve_reference_files(args, spec)
        inspect_reference_files(files, args, spec, years)
        reference_files[spec["key"]] = files
    if args.check_only:
        print(
            "All model and reference inputs passed validation; no figure was "
            "created."
        )
        return
    (
        lat2d,
        lon2d,
        lat1d,
        lon1d,
        lat_grid,
        lon_grid,
        land_mask,
    ) = load_plot_grid(args.wrfinput, args.mask_file, args.mask_variable)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    os.environ["CARTOPY_OFFLINE"] = "true"
    for spec in variable_specs:
        product = reference_product(args, spec)
        print(f"\n=== {spec['title']} | reference: {product} ===")
        reference = load_reference_yearly(
            reference_files[spec["key"]],
            args,
            spec,
            years,
            lat_grid,
            lon_grid,
            land_mask,
        )
        model = load_model_yearly(
            model_files,
            spec,
            args,
            lat2d,
            lon2d,
            lat1d,
            lon1d,
            land_mask,
        )
        align_model_years_to_reference(
            reference, model, years, spec
        )
        masks = common_plot_data(reference, model, land_mask)
        for experiment in ("OFF", "CPL"):
            for season_name, _ in SEASONS:
                count = np.count_nonzero(masks[experiment][season_name])
                if count < 3:
                    raise ValueError(
                        f"Unified mask has only {count} cells for "
                        f"{experiment} {season_name}"
                    )
                print(
                    f"Unified mask {experiment} {season_name}: {count} cells"
                )
        if not args.skip_spatial:
            plot_spatial_evaluation(
                spec,
                reference,
                model,
                masks,
                args,
                lat_grid,
                lon_grid,
                land_mask,
                product,
            )
        if not args.skip_regional_timeseries:
            plot_regional_timeseries(
                spec,
                reference,
                model,
                masks,
                args,
                lat_grid,
                product,
                years,
            )
        if not args.skip_acc_timeseries:
            plot_acc_timeseries(
                spec,
                reference,
                model,
                masks,
                args,
                lat_grid,
                product,
                years,
            )


def run_response_only_category(variable_specs, description):
    """Run model-only process diagnostics without any reference product."""
    parser = build_parser(description, [], reference_mode=False)
    args = parser.parse_args()
    validate_arguments(args)
    model_files = resolve_model_files(args)
    validate_model_inputs(model_files, variable_specs, args.nyears)
    years = list(range(args.start_year, args.start_year + args.nyears))
    if args.check_only:
        print(
            "All model inputs passed validation; no reference data are "
            "required and no figure was created."
        )
        return
    (
        lat2d,
        lon2d,
        lat1d,
        lon1d,
        lat_grid,
        lon_grid,
        land_mask,
    ) = load_plot_grid(args.wrfinput, args.mask_file, args.mask_variable)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    os.environ["CARTOPY_OFFLINE"] = "true"
    for spec in variable_specs:
        print(f"\n=== {spec['title']} | process response only ===")
        model = load_model_yearly(
            model_files,
            spec,
            args,
            lat2d,
            lon2d,
            lat1d,
            lon1d,
            land_mask,
        )
        masks = model_response_masks(model, land_mask)
        for experiment in ("OFF", "CPL"):
            for season_name, _ in SEASONS:
                count = np.count_nonzero(masks[experiment][season_name])
                if count < 3:
                    raise ValueError(
                        f"Unified mask has only {count} cells for "
                        f"{experiment} {season_name}"
                    )
                print(
                    f"Unified mask {experiment} {season_name}: {count} cells"
                )
        if not args.skip_spatial:
            plot_response_only_spatial(
                spec,
                model,
                masks,
                args,
                lat_grid,
                lon_grid,
                land_mask,
            )
        if not args.skip_regional_timeseries:
            plot_response_only_timeseries(
                spec, model, masks, args, lat_grid, years
            )
