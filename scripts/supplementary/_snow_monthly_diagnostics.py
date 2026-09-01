#!/usr/bin/env python3
"""Monthly SWE/SCF diagnostics over grid cells with gravel fraction > 0.3."""

import argparse
import calendar
import csv
import os
from datetime import datetime
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset
from scipy.spatial import cKDTree

from _land_response_common import (
    DEFAULT_DATA_DIR,
    DEFAULT_MASK_FILE,
    DEFAULT_OUTPUT_DIR,
    DEFAULT_WRFINPUT,
    REFERENCE_SUBSET_EXTENT,
    automatic_reference_factor,
    clean_array,
    coordinate_subset,
    decode_dates,
    find_coordinate,
    find_time_variable,
    load_plot_grid,
    read_reference_block,
    reference_cli_value,
    reference_product,
    regrid_regular_to_target,
    resolve_model_files,
    resolve_reference_files,
    resolve_reference_variable,
    validate_model_inputs,
)


DEFAULT_GRAVEL_FILE = Path(
    "/share/home/dq013/zhwei/colm/data/CoLMrawdata/soil/vf_gravels_s.nc"
)
MODEL_MONTHS = (3, 4, 5, 6, 7, 8)
DZ = np.array(
    [0.0175, 0.0276, 0.0455, 0.0750, 0.1236,
     0.2038, 0.3360, 0.5539, 0.9133, 1.5058],
    dtype=float,
)
DZ /= DZ.sum()
GRAVEL_WEIGHTS = np.array(
    [DZ[0] + DZ[1], DZ[2], DZ[3], DZ[4],
     DZ[5], DZ[6], DZ[7], DZ[8] + DZ[9]],
    dtype=float,
)
PIE_RANGES = (
    ("35-50 mm", 35.0, 50.0, True),
    ("20-35 mm", 20.0, 35.0, False),
    ("5-20 mm", 5.0, 20.0, False),
)
PIE_COLORS = ("#2171b5", "#6baed6", "#c6dbef")


def build_snow_parser(description, specs):
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("--data-dir", type=Path, default=DEFAULT_DATA_DIR)
    parser.add_argument("--wrfinput", type=Path, default=DEFAULT_WRFINPUT)
    parser.add_argument("--mask-file", type=Path, default=DEFAULT_MASK_FILE)
    parser.add_argument("--mask-variable", default="tm")
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--off-ctl", type=Path, default=None)
    parser.add_argument("--off-exp", type=Path, default=None)
    parser.add_argument("--cpl-ctl", type=Path, default=None)
    parser.add_argument("--cpl-exp", type=Path, default=None)
    parser.add_argument("--nyears", type=int, default=17)
    parser.add_argument("--start-year", type=int, default=2001)
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument("--gravel-file", type=Path, default=DEFAULT_GRAVEL_FILE)
    parser.add_argument("--gravel-threshold", type=float, default=0.3)
    parser.add_argument(
        "--snow-cover-threshold",
        type=float,
        default=0.01,
        help=(
            "Minimum CoLM f_fsno used to identify snow-covered grid cells "
            "(default: 0.01, i.e. 1 percent)."
        ),
    )
    parser.add_argument("--check-only", action="store_true")
    for spec in specs:
        key = spec["key"].replace("_", "-")
        dest = spec["key"]
        parser.add_argument(
            f"--{key}-reference-file",
            dest=f"{dest}_reference_files",
            action="append",
            metavar="PATH_OR_GLOB",
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
        )
        parser.add_argument(
            f"--{key}-reference-factor",
            dest=f"{dest}_reference_factor",
            type=float,
            default=None,
        )
    return parser


def monthly_dates(start_year, nyears):
    return [
        datetime(year, month, 15)
        for year in range(start_year, start_year + nyears)
        for month in range(1, 13)
    ]


def active_model_dates(start_year, nyears):
    return [
        datetime(year, month, 15)
        for year in range(start_year, start_year + nyears)
        for month in MODEL_MONTHS
    ]


def inspect_monthly_reference(files, args, spec, years):
    required = {(year, month) for year in years for month in range(1, 13)}
    found = set()
    first_variable = None
    first_units = ""
    for path in files:
        with Dataset(path) as nc_obj:
            variable = resolve_reference_variable(nc_obj, args, spec)
            time_name = find_time_variable(nc_obj, variable)
            dates = decode_dates(nc_obj, time_name, path)
            for date in dates:
                key = (int(date.year), int(date.month))
                if key in required:
                    found.add(key)
            if first_variable is None:
                first_variable = variable
                first_units = getattr(nc_obj.variables[variable], "units", "")
    missing = sorted(required - found)
    if missing:
        detail = ", ".join(f"{year}-{month:02d}" for year, month in missing)
        print(
            f"[WARN] {reference_product(args, spec)} lacks {len(missing)} "
            f"monthly record(s): {detail}; the reference curve will contain "
            "gaps at those actual missing months"
        )
    print(
        f"[OK] {reference_product(args, spec)}: {len(files)} file(s), "
        f"variable={first_variable}, units={first_units or '(missing)'}, "
        f"monthly coverage={years[0]}-{years[-1]}"
    )


def _nearest_indices(coordinate, targets):
    coordinate = np.asarray(coordinate, dtype=float)
    targets = np.asarray(targets, dtype=float)
    order = np.argsort(coordinate)
    sorted_coordinate = coordinate[order]
    insertion = np.searchsorted(sorted_coordinate, targets)
    insertion = np.clip(insertion, 1, len(sorted_coordinate) - 1)
    left = insertion - 1
    right = insertion
    choose_right = (
        np.abs(sorted_coordinate[right] - targets)
        < np.abs(sorted_coordinate[left] - targets)
    )
    selected = np.where(choose_right, right, left)
    return order[selected]


def _read_orthogonal(variable, latitude_indices, longitude_indices):
    """Read a small target-grid sample from a very large regular HDF5 grid."""
    lat_order = np.argsort(latitude_indices)
    lon_order = np.argsort(longitude_indices)
    lat_sorted = np.asarray(latitude_indices)[lat_order]
    lon_sorted = np.asarray(longitude_indices)[lon_order]
    data = clean_array(variable[lat_sorted, lon_sorted])
    if data.shape != (len(lat_sorted), len(lon_sorted)):
        raise ValueError(
            f"Unexpected gravel sample shape {data.shape}; expected "
            f"{(len(lat_sorted), len(lon_sorted))}"
        )
    return data[np.argsort(lat_order)][:, np.argsort(lon_order)]


def _gravel_coordinates(nc_obj, path, announce=False):
    """Return explicit coordinates or infer global regular-grid cell centers."""
    lower = {name.lower(): name for name in nc_obj.variables}
    lat_name = next(
        (lower[name] for name in ("lat", "latitude") if name in lower),
        None,
    )
    lon_name = next(
        (lower[name] for name in ("lon", "longitude") if name in lower),
        None,
    )
    if (lat_name is None) != (lon_name is None):
        raise ValueError(
            f"{path} contains only one horizontal coordinate; both latitude "
            "and longitude are required"
        )
    variable = nc_obj.variables["vf_gravels_s_l1"]
    if variable.ndim != 2:
        raise ValueError(
            f"vf_gravels_s_l1 in {path} must be two-dimensional; "
            f"found shape {variable.shape}"
        )
    nlat, nlon = variable.shape
    if lat_name is not None:
        source_lat = clean_array(nc_obj.variables[lat_name][:])
        source_lon = clean_array(nc_obj.variables[lon_name][:])
        if source_lat.ndim != 1 or source_lon.ndim != 1:
            raise ValueError(f"Horizontal coordinates in {path} must be 1-D")
        if source_lat.size != nlat or source_lon.size != nlon:
            raise ValueError(
                f"Coordinate sizes in {path} do not match gravel grid "
                f"{(nlat, nlon)}"
            )
        return source_lat, source_lon

    # The CoLM raw-data gravel file omits coordinate variables.  Its native
    # grid is global, regular latitude-longitude, ordered north-to-south and
    # west-to-east.  Infer cell centers only when the 2:1 global-grid geometry
    # is unambiguous; otherwise stop rather than silently assigning coordinates.
    if nlat < 2 or nlon != 2 * nlat:
        raise ValueError(
            f"{path} has no coordinate variables and shape {(nlat, nlon)} "
            "is not an unambiguous global regular grid (expected nlon=2*nlat)"
        )
    dlat = 180.0 / nlat
    dlon = 360.0 / nlon
    source_lat = 90.0 - (np.arange(nlat, dtype=float) + 0.5) * dlat
    source_lon = -180.0 + (np.arange(nlon, dtype=float) + 0.5) * dlon
    if announce:
        print(
            f"[WARN] {path} has no latitude/longitude variables; inferred "
            f"global cell-center coordinates from shape {(nlat, nlon)} "
            f"({dlon:g} degree resolution)"
        )
    return source_lat, source_lon


def inspect_gravel_file(path):
    if not Path(path).is_file():
        raise FileNotFoundError(f"Gravel file not found: {path}")
    with Dataset(path) as nc_obj:
        missing = [
            f"vf_gravels_s_l{layer}"
            for layer in range(1, 9)
            if f"vf_gravels_s_l{layer}" not in nc_obj.variables
        ]
        if missing:
            raise KeyError(f"Missing gravel variables in {path}: {missing}")
        _gravel_coordinates(nc_obj, path, announce=True)
    print(f"[OK] Gravel mask source: {path}; layers=1-8")


def load_gravel_mask(path, lat1d, lon1d, land_mask, threshold):
    """Thickness-average eight gravel layers at nearest 0.25-degree cells."""
    with Dataset(path) as nc_obj:
        source_lat, source_lon = _gravel_coordinates(nc_obj, path)
        target_lon = ((np.asarray(lon1d) + 180.0) % 360.0) - 180.0
        normalized_source_lon = (
            (np.asarray(source_lon) + 180.0) % 360.0
        ) - 180.0
        lat_indices = _nearest_indices(source_lat, lat1d)
        lon_indices = _nearest_indices(normalized_source_lon, target_lon)
        weighted_sum = np.zeros(land_mask.shape, dtype=float)
        valid_weight = np.zeros(land_mask.shape, dtype=float)
        for layer, weight in enumerate(GRAVEL_WEIGHTS, start=1):
            name = f"vf_gravels_s_l{layer}"
            data = _read_orthogonal(
                nc_obj.variables[name], lat_indices, lon_indices
            )
            data[(data < 0) | (data > 1000)] = np.nan
            finite = np.isfinite(data)
            weighted_sum[finite] += data[finite] * weight
            valid_weight[finite] += weight
    with np.errstate(invalid="ignore", divide="ignore"):
        gravel = np.where(
            valid_weight > 0, weighted_sum / valid_weight, np.nan
        )
    finite = gravel[np.isfinite(gravel)]
    if finite.size and np.nanpercentile(finite, 95) > 1.5:
        print("[WARN] Gravel values appear to be percent; converting to fraction")
        gravel /= 100.0
    mask = land_mask & np.isfinite(gravel) & (gravel > threshold)
    if np.count_nonzero(mask) < 3:
        raise ValueError(
            f"Gravel > {threshold:g} mask contains only "
            f"{np.count_nonzero(mask)} target cells"
        )
    print(
        f"[OK] Gravel > {threshold:g} mask: {np.count_nonzero(mask)} "
        "China-land target cells"
    )
    return gravel, mask


def _source_to_target_index(lat2d, lon2d, lat_grid, lon_grid):
    source_valid = np.isfinite(lat2d) & np.isfinite(lon2d)
    source_flat = np.flatnonzero(source_valid)
    source_points = np.column_stack(
        (lon2d[source_valid], lat2d[source_valid])
    )
    target_points = np.column_stack((lon_grid.ravel(), lat_grid.ravel()))
    tree = cKDTree(source_points)
    _, local_index = tree.query(target_points, k=1)
    source_index = source_flat[local_index]
    inside = (
        (lat_grid >= np.nanmin(lat2d))
        & (lat_grid <= np.nanmax(lat2d))
        & (lon_grid >= np.nanmin(lon2d))
        & (lon_grid <= np.nanmax(lon2d))
    )
    return source_index, inside


def load_model_monthly(
    model_files,
    spec,
    nyears,
    source_index,
    target_shape,
    land_mask,
    inside,
):
    needed = nyears * len(MODEL_MONTHS)
    output = {}
    valid_target = land_mask & inside
    for experiment, paths in model_files.items():
        output[experiment] = {}
        for role, path in paths.items():
            with Dataset(path) as nc_obj:
                variable = nc_obj.variables[spec["variable"]]
                native = clean_array(variable[:needed]) * spec["factor"]
            flat = native.reshape(native.shape[0], -1)
            target = flat[:, source_index].reshape((needed,) + target_shape)
            target[:, ~valid_target] = np.nan
            if spec.get("nonnegative", False):
                target[target < 0] = np.nan
            output[experiment][role] = target
    return output


def load_reference_monthly(
    files, args, spec, years, lat_grid, lon_grid, land_mask
):
    required = {(year, month) for year in years for month in range(1, 13)}
    monthly_sum = {}
    monthly_count = {}
    common_lat = common_lon = None
    override_factor = reference_cli_value(args, spec, "factor")
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
            lat_sub = np.asarray(latitude[lat_slice])
            lon_sub = np.asarray(longitude[lon_slice])
            if common_lat is None:
                common_lat, common_lon = lat_sub, lon_sub
            elif not (
                common_lat.shape == lat_sub.shape
                and common_lon.shape == lon_sub.shape
                and np.allclose(common_lat, lat_sub)
                and np.allclose(common_lon, lon_sub)
            ):
                raise ValueError(
                    "Reference files use different grids; preprocess to one "
                    "regular latitude-longitude grid"
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
                block *= factor
                block_sum = np.nansum(block, axis=0)
                block_count = np.sum(np.isfinite(block), axis=0)
                key = (year, month)
                if key in monthly_sum:
                    monthly_sum[key] += block_sum
                    monthly_count[key] += block_count
                else:
                    monthly_sum[key] = block_sum
                    monthly_count[key] = block_count

    shape = (len(common_lat), len(common_lon))
    output = []
    for year in years:
        for month in range(1, 13):
            key = (year, month)
            if key in monthly_sum:
                with np.errstate(invalid="ignore", divide="ignore"):
                    native = np.where(
                        monthly_count[key] > 0,
                        monthly_sum[key] / monthly_count[key],
                        np.nan,
                    )
            else:
                native = np.full(shape, np.nan, dtype=float)
            target = regrid_regular_to_target(
                native, common_lat, common_lon, lat_grid, lon_grid
            )
            target = np.where(land_mask, target, np.nan)
            if spec.get("nonnegative", False):
                target[target < 0] = np.nan
            output.append(target)
    return np.stack(output)


def area_mean_series(data, lat_grid, mask):
    weights = np.cos(np.deg2rad(lat_grid))
    values = []
    for field in data:
        valid = mask & np.isfinite(field)
        if np.count_nonzero(valid) < 1:
            values.append(np.nan)
        else:
            values.append(np.sum(field[valid] * weights[valid]) / np.sum(weights[valid]))
    return np.asarray(values)


def normalize_snow_fraction(data, label):
    """Normalize a snow-cover field supplied as either 0-1 or 0-100."""
    output = np.asarray(data, dtype=float).copy()
    finite = output[np.isfinite(output)]
    if finite.size == 0:
        raise ValueError(f"{label} contains no finite snow-cover values")
    lower = float(np.nanmin(finite))
    upper = float(np.nanmax(finite))
    if lower < -1.0e-6:
        raise ValueError(f"{label} contains negative snow-cover values: {lower:g}")
    if upper <= 1.0 + 1.0e-6:
        factor = 1.0
    elif upper <= 100.0 + 1.0e-4:
        factor = 0.01
    else:
        raise ValueError(
            f"{label} snow-cover maximum {upper:g} is outside both the "
            "0-1 fraction and 0-100 percent ranges"
        )
    output *= factor
    output[(output < 0.0) | (output > 1.0 + 1.0e-6)] = np.nan
    return np.clip(output, 0.0, 1.0)


def dynamic_mask_area_mean_series(data, dynamic_mask, lat_grid, base_mask):
    """Return area means over a time-varying CoLM snow-presence mask."""
    if data.shape != dynamic_mask.shape:
        raise ValueError(
            "Data and dynamic-mask arrays must have identical shapes; found "
            f"{data.shape} and {dynamic_mask.shape}"
        )
    weights = np.cos(np.deg2rad(lat_grid))
    values = []
    for field, snow_mask in zip(data, dynamic_mask):
        valid = base_mask & snow_mask & np.isfinite(field)
        if np.count_nonzero(valid) < 1:
            values.append(np.nan)
            continue
        values.append(
            np.sum(field[valid] * weights[valid]) / np.sum(weights[valid])
        )
    return np.asarray(values)


def select_active_model_months(data, nyears):
    """Select March-August from a January-December reference array."""
    indices = [
        year_index * 12 + month - 1
        for year_index in range(nyears)
        for month in MODEL_MONTHS
    ]
    return data[indices]


def expand_model_series(active_values, nyears):
    full = np.full(nyears * 12, np.nan, dtype=float)
    cursor = 0
    for year_index in range(nyears):
        for month in MODEL_MONTHS:
            full[year_index * 12 + month - 1] = active_values[cursor]
            cursor += 1
    return full


def seasonal_climatology(data, dates, months):
    selected = [index for index, date in enumerate(dates) if date.month in months]
    return np.nanmean(data[selected], axis=0)


def pie_statistics(field, mask):
    values = field[mask & np.isfinite(field)]
    below = int(np.count_nonzero(values < 5.0))
    above = int(np.count_nonzero(values > 50.0))
    counts = []
    for _, lower, upper, upper_inclusive in PIE_RANGES:
        selected = (values >= lower) & (
            (values <= upper) if upper_inclusive else (values < upper)
        )
        counts.append(int(np.count_nonzero(selected)))
    total = sum(counts)
    percentages = [100.0 * count / total if total else np.nan for count in counts]
    return counts, percentages, total, below, above


def plot_swe_pies(pie_fields, mask, args):
    sources = (
        "ERA5-Land",
        "OFF CTL",
        "OFF EXP",
        "CPL CTL",
        "CPL EXP",
    )
    seasons = (("MAM", (3, 4, 5)), ("JJA", (6, 7, 8)))
    fig, axes = plt.subplots(2, 5, figsize=(15.5, 6.4), squeeze=False)
    csv_rows = []
    for row, (season, _) in enumerate(seasons):
        for col, source in enumerate(sources):
            ax = axes[row, col]
            counts, percentages, total, below, above = pie_statistics(
                pie_fields[season][source], mask
            )
            if total:
                ax.pie(
                    counts,
                    colors=PIE_COLORS,
                    startangle=90,
                    counterclock=False,
                    autopct=lambda value: f"{value:.1f}%" if value >= 2 else "",
                    textprops={"fontsize": 8},
                    wedgeprops={"linewidth": 0.7, "edgecolor": "white"},
                )
            else:
                ax.text(0.5, 0.5, "No 5-50 mm cells", ha="center", va="center")
            if row == 0:
                ax.set_title(source, fontsize=10.5)
            if col == 0:
                ax.text(
                    -0.18,
                    0.5,
                    season,
                    transform=ax.transAxes,
                    rotation=90,
                    va="center",
                    ha="center",
                    fontsize=11,
                    fontweight="bold",
                )
            for (label, lower, upper, _), count, percentage in zip(
                PIE_RANGES, counts, percentages
            ):
                csv_rows.append(
                    {
                        "season": season,
                        "source": source,
                        "range": label,
                        "lower_mm": lower,
                        "upper_mm": upper,
                        "grid_count": count,
                        "percentage_within_5_50_mm": percentage,
                        "total_grids_5_50_mm": total,
                        "grids_below_5_mm": below,
                        "grids_above_50_mm": above,
                        "gravel_threshold": args.gravel_threshold,
                    }
                )
    legend_handles = [
        plt.Line2D(
            [0],
            [0],
            marker="o",
            color="none",
            markerfacecolor=color,
            markeredgecolor="none",
            markersize=10,
        )
        for color in PIE_COLORS
    ]
    fig.legend(
        legend_handles,
        [item[0] for item in PIE_RANGES],
        ncol=3,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.925),
        frameon=False,
    )
    fig.suptitle(
        "SWE grid proportions over gravel-content > 0.3 regions (2001-2017)",
        fontsize=14,
        y=0.985,
    )
    fig.text(
        0.5,
        0.025,
        (
            "Percentages are grid counts normalized within the three requested "
            "SWE classes (5-50 mm); grids below 5 mm or above 50 mm are reported "
            "in the companion CSV but excluded from the pie denominator."
        ),
        ha="center",
        fontsize=8,
    )
    fig.tight_layout(rect=(0.02, 0.06, 0.99, 0.90))
    output = args.output_dir / "supp_land_swe_range_proportions_gravel_gt_0p3.pdf"
    fig.savefig(output, dpi=args.dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {output}")
    csv_path = args.output_dir / "supp_land_swe_range_proportions_gravel_gt_0p3.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(csv_rows[0]))
        writer.writeheader()
        writer.writerows(csv_rows)
    print(f"Saved: {csv_path}")


def _monthly_tick_labels(dates):
    initials = ("J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D")
    labels = []
    for date in dates:
        if date.month == 1:
            labels.append(f"Jan\n{date.year}")
        else:
            labels.append(initials[date.month - 1])
    return labels


def plot_combined_timeseries(series, dates, args):
    panels = (
        ("swe", "OFF", "SWE - OFF", "SWE (mm)"),
        ("swe", "CPL", "SWE - CPL", "SWE (mm)"),
        ("scf", "OFF", "SCF - OFF", "SCF (%)"),
        ("scf", "CPL", "SCF - CPL", "SCF (%)"),
    )
    fig, axes = plt.subplots(4, 1, figsize=(20, 11.5), sharex=True)
    csv_rows = []
    for row, (key, experiment, title, ylabel) in enumerate(panels):
        ax = axes[row]
        block = series[key][experiment]
        ax.plot(
            dates,
            block["reference"],
            color="black",
            linewidth=1.65,
            label=block["reference_label"],
            zorder=3,
        )
        ax.plot(
            dates,
            block["ctl"],
            color="#2166ac",
            linewidth=1.35,
            marker="o",
            markersize=1.8,
            label="CTL (Mar-Aug)",
        )
        ax.plot(
            dates,
            block["exp"],
            color="#d73027",
            linewidth=1.35,
            marker="o",
            markersize=1.8,
            label="EXP (Mar-Aug)",
        )
        for date in dates:
            if date.month == 1:
                ax.axvline(date, color="0.86", linewidth=0.55, zorder=0)
        ax.set_title(title, loc="left", fontsize=10.8, fontweight="bold")
        ax.set_ylabel(ylabel)
        ax.grid(axis="y", color="0.88", linewidth=0.55)
        ax.margins(x=0.002)
        if row == 0:
            ax.legend(ncol=3, frameon=False, loc="upper right", fontsize=8.5)
        for index, date in enumerate(dates):
            csv_rows.append(
                {
                    "date": date.strftime("%Y-%m"),
                    "variable": key.upper(),
                    "experiment": experiment,
                    "reference_product": block["reference_label"],
                    "reference": block["reference"][index],
                    "ctl": block["ctl"][index],
                    "exp": block["exp"][index],
                    "gravel_threshold": args.gravel_threshold,
                }
            )
    axes[-1].set_xticks(dates)
    axes[-1].set_xticklabels(_monthly_tick_labels(dates), fontsize=4.6)
    axes[-1].tick_params(axis="x", length=2.2, pad=2)
    axes[-1].set_xlabel(
        "Month (year is printed only at January; model curves contain Mar-Aug only)"
    )
    fig.suptitle(
        "Monthly SWE and snow-cover fraction over gravel-content > 0.3 regions",
        fontsize=14,
        y=0.992,
    )
    fig.text(
        0.5,
        0.012,
        (
            "Reference curves use all available months. CTL/EXP are intentionally "
            "disconnected between annual March-August segments; missing MOD10CM "
            "granules remain missing rather than being interpolated or set to zero."
        ),
        ha="center",
        fontsize=8.2,
    )
    fig.tight_layout(rect=(0.025, 0.045, 0.995, 0.975))
    output = args.output_dir / "supp_land_snow_monthly_timeseries_gravel_gt_0p3.pdf"
    fig.savefig(output, dpi=args.dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {output}")
    csv_path = args.output_dir / "supp_land_snow_monthly_timeseries_gravel_gt_0p3.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(csv_rows[0]))
        writer.writeheader()
        writer.writerows(csv_rows)
    print(f"Saved: {csv_path}")


def plot_colm_snow_mask_swe_timeseries(series, dates, args):
    """Plot SWE over grid cells diagnosed as snow-covered by CoLM."""
    fig, axes = plt.subplots(2, 1, figsize=(20, 6.6), sharex=True)
    csv_rows = []
    for row, experiment in enumerate(("OFF", "CPL")):
        ax = axes[row]
        block = series[experiment]
        ax.plot(
            dates,
            block["reference"],
            color="black",
            linewidth=1.65,
            label="ERA5-Land on CoLM snow mask (Mar-Aug)",
            zorder=3,
        )
        ax.plot(
            dates,
            block["ctl"],
            color="#2166ac",
            linewidth=1.35,
            marker="o",
            markersize=1.8,
            label="CTL on CoLM snow mask (Mar-Aug)",
        )
        ax.plot(
            dates,
            block["exp"],
            color="#d73027",
            linewidth=1.35,
            marker="o",
            markersize=1.8,
            label="EXP on CoLM snow mask (Mar-Aug)",
        )
        for date in dates:
            if date.month == 1:
                ax.axvline(date, color="0.86", linewidth=0.55, zorder=0)
        ax.set_title(
            f"SWE over CoLM snow-covered grid cells - {experiment}",
            loc="left",
            fontsize=10.8,
            fontweight="bold",
        )
        ax.set_ylabel("SWE (mm)")
        ax.grid(axis="y", color="0.88", linewidth=0.55)
        ax.margins(x=0.002)
        if row == 0:
            ax.legend(ncol=3, frameon=False, loc="upper right", fontsize=8.3)
        for index, date in enumerate(dates):
            csv_rows.append(
                {
                    "date": date.strftime("%Y-%m"),
                    "experiment": experiment,
                    "reference": block["reference"][index],
                    "ctl": block["ctl"][index],
                    "exp": block["exp"][index],
                    "gravel_threshold": args.gravel_threshold,
                    "snow_cover_threshold": args.snow_cover_threshold,
                }
            )
    axes[-1].set_xticks(dates)
    axes[-1].set_xticklabels(_monthly_tick_labels(dates), fontsize=4.6)
    axes[-1].tick_params(axis="x", length=2.2, pad=2)
    axes[-1].set_xlabel(
        "Month (year is printed only at January; model curves contain Mar-Aug only)"
    )
    fig.suptitle(
        "Monthly SWE over CoLM snow-covered, gravel-content > 0.3 grid cells",
        fontsize=14,
        y=0.992,
    )
    fig.text(
        0.5,
        0.012,
        (
            "The monthly mask is the union of CTL and EXP CoLM grid cells with "
            f"f_fsno >= {100.0 * args.snow_cover_threshold:g}%; ERA5-Land, CTL, "
            "and EXP are averaged on that identical mask. This diagnostic does "
            "not replace the all-gravel-region SWE figure."
        ),
        ha="center",
        fontsize=8.2,
    )
    fig.tight_layout(rect=(0.025, 0.055, 0.995, 0.965))
    output = (
        args.output_dir
        / "supp_land_swe_colm_snow_mask_timeseries_gravel_gt_0p3.pdf"
    )
    fig.savefig(output, dpi=args.dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {output}")
    csv_path = (
        args.output_dir
        / "supp_land_swe_colm_snow_mask_timeseries_gravel_gt_0p3.csv"
    )
    with csv_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(csv_rows[0]))
        writer.writeheader()
        writer.writerows(csv_rows)
    print(f"Saved: {csv_path}")


def run_snow_monthly_diagnostics(specs, description):
    parser = build_snow_parser(description, specs)
    args = parser.parse_args()
    if args.nyears < 1:
        raise ValueError("--nyears must be at least 1")
    if not np.isfinite(args.gravel_threshold):
        raise ValueError("--gravel-threshold must be finite")
    if not 0.0 < args.snow_cover_threshold <= 1.0:
        raise ValueError("--snow-cover-threshold must be in the interval (0, 1]")
    model_files = resolve_model_files(args)
    validate_model_inputs(model_files, specs, args.nyears)
    inspect_gravel_file(args.gravel_file)
    years = list(range(args.start_year, args.start_year + args.nyears))
    reference_files = {}
    for spec in specs:
        files = resolve_reference_files(args, spec)
        inspect_monthly_reference(files, args, spec, years)
        reference_files[spec["key"]] = files
    if args.check_only:
        print(
            "All model, SWE/SCF reference, and gravel-mask inputs passed "
            "validation; no figure was created."
        )
        return

    lat2d, lon2d, lat1d, lon1d, lat_grid, lon_grid, land_mask = load_plot_grid(
        args.wrfinput, args.mask_file, args.mask_variable
    )
    _, gravel_mask = load_gravel_mask(
        args.gravel_file,
        lat1d,
        lon1d,
        land_mask,
        args.gravel_threshold,
    )
    source_index, inside = _source_to_target_index(
        lat2d, lon2d, lat_grid, lon_grid
    )
    dates = monthly_dates(args.start_year, args.nyears)
    model_dates = active_model_dates(args.start_year, args.nyears)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    os.environ["CARTOPY_OFFLINE"] = "true"
    series = {}
    colm_snow_mask_series = {}
    pie_fields = {"MAM": {}, "JJA": {}}
    scf_spec = next(spec for spec in specs if spec["key"] == "scf")

    for spec in specs:
        key = spec["key"]
        product = reference_product(args, spec)
        print(f"\n=== Monthly {spec['title']} | reference: {product} ===")
        reference = load_reference_monthly(
            reference_files[key],
            args,
            spec,
            years,
            lat_grid,
            lon_grid,
            land_mask,
        )
        model = load_model_monthly(
            model_files,
            spec,
            args.nyears,
            source_index,
            lat_grid.shape,
            land_mask,
            inside,
        )
        reference_series = area_mean_series(reference, lat_grid, gravel_mask)
        series[key] = {}
        for experiment in ("OFF", "CPL"):
            ctl_active = area_mean_series(
                model[experiment]["ctl"], lat_grid, gravel_mask
            )
            exp_active = area_mean_series(
                model[experiment]["exp"], lat_grid, gravel_mask
            )
            series[key][experiment] = {
                "reference": reference_series,
                "reference_label": spec.get("reference_short_name", product),
                "ctl": expand_model_series(ctl_active, args.nyears),
                "exp": expand_model_series(exp_active, args.nyears),
            }
        if key == "swe":
            model_snow_cover = load_model_monthly(
                model_files,
                scf_spec,
                args.nyears,
                source_index,
                lat_grid.shape,
                land_mask,
                inside,
            )
            reference_active = select_active_model_months(
                reference, args.nyears
            )
            for experiment in ("OFF", "CPL"):
                ctl_fraction = normalize_snow_fraction(
                    model_snow_cover[experiment]["ctl"],
                    f"{experiment} CTL f_fsno",
                )
                exp_fraction = normalize_snow_fraction(
                    model_snow_cover[experiment]["exp"],
                    f"{experiment} EXP f_fsno",
                )
                snow_mask = (
                    (ctl_fraction >= args.snow_cover_threshold)
                    | (exp_fraction >= args.snow_cover_threshold)
                )
                reference_on_mask = dynamic_mask_area_mean_series(
                    reference_active,
                    snow_mask,
                    lat_grid,
                    gravel_mask,
                )
                ctl_on_mask = dynamic_mask_area_mean_series(
                    model[experiment]["ctl"],
                    snow_mask,
                    lat_grid,
                    gravel_mask,
                )
                exp_on_mask = dynamic_mask_area_mean_series(
                    model[experiment]["exp"],
                    snow_mask,
                    lat_grid,
                    gravel_mask,
                )
                colm_snow_mask_series[experiment] = {
                    "reference": expand_model_series(
                        reference_on_mask, args.nyears
                    ),
                    "ctl": expand_model_series(ctl_on_mask, args.nyears),
                    "exp": expand_model_series(exp_on_mask, args.nyears),
                }
            del model_snow_cover
            seasons = {"MAM": (3, 4, 5), "JJA": (6, 7, 8)}
            for season, months in seasons.items():
                pie_fields[season]["ERA5-Land"] = seasonal_climatology(
                    reference, dates, months
                )
                for experiment in ("OFF", "CPL"):
                    for role in ("ctl", "exp"):
                        label = f"{experiment} {role.upper()}"
                        pie_fields[season][label] = seasonal_climatology(
                            model[experiment][role], model_dates, months
                        )
        del reference, model

    plot_combined_timeseries(series, dates, args)
    plot_colm_snow_mask_swe_timeseries(colm_snow_mask_series, dates, args)
    plot_swe_pies(pie_fields, gravel_mask, args)
