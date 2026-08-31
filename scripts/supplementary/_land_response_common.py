#!/usr/bin/env python3
"""Shared plotting utilities for supplementary CoLM land-response figures.

The public entry scripts in this directory each define one physical variable
group and call :func:`run_category`.  All figures compare the gravel-inclusive
experiment (EXP) with the gravel-free control (CTL) over 2001-2017.

Expected processed-file layout
------------------------------
Each model file contains six monthly means per year in this order:
March, April, May, June, July, August.  Variables are expected on a gridded
``(time, south_north, west_east)`` layout, consistent with the existing
``colmoff_*`` and ``colmrun_*`` files used by this repository.
"""

import argparse
import os
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
from netCDF4 import Dataset
from scipy.interpolate import griddata
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
PROJ_CENTRAL_LON = 110.0
PROJ_CENTRAL_LAT = 40.0
PROJ_STD_PARALLELS = (30.0, 60.0)

SEASONS = (
    ("Spring", (0, 1, 2)),
    ("Summer", (3, 4, 5)),
)


def build_parser(description):
    """Create a command-line parser shared by all category scripts."""
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
    parser.add_argument(
        "--nyears",
        type=int,
        default=17,
        help="Number of years from the beginning of each six-month file.",
    )
    parser.add_argument(
        "--start-year",
        type=int,
        default=2001,
        help="Used in labels and checks; data are read from the file start.",
    )
    parser.add_argument(
        "--p-threshold",
        type=float,
        default=0.05,
        help="Two-sided paired t-test threshold.",
    )
    parser.add_argument(
        "--significance-style",
        choices=("stipple", "mask", "none"),
        default="stipple",
        help=(
            "stipple: full response plus gray dots; mask: retain only p-values "
            "below the threshold; none: show the full response without marks."
        ),
    )
    parser.add_argument(
        "--sig-stride",
        type=int,
        default=5,
        help="Subsampling stride for significance dots on the 0.25-degree grid.",
    )
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument(
        "--check-only",
        action="store_true",
        help="Validate files, variables, dimensions, and time length without plotting.",
    )
    return parser


def resolve_model_files(args):
    """Resolve default processed files while allowing command-line overrides."""
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
    if not path.is_file():
        raise FileNotFoundError(f"{label} not found: {path}")


def clean_array(data):
    """Convert a NetCDF slice to floating point and normalize missing values."""
    if np.ma.isMaskedArray(data):
        data = data.filled(np.nan)
    data = np.asarray(data, dtype=np.float64)
    data[np.abs(data) > 1.0e30] = np.nan
    data[data == -9999] = np.nan
    return data


def find_coordinate(nc_obj, names):
    for name in names:
        if name in nc_obj.variables:
            return clean_array(nc_obj.variables[name][:])
    raise KeyError(f"None of the coordinate variables exists: {names}")


def load_plot_grid(wrfinput_file, mask_file, mask_variable):
    """Read the native CRESM grid and the repository's 0.25-degree land mask."""
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
        lat1d = find_coordinate(nc_mask, ("lat", "latitude"))
        lon1d = find_coordinate(nc_mask, ("lon", "longitude"))
        if mask_variable not in nc_mask.variables:
            raise KeyError(
                f"Mask variable {mask_variable!r} is absent from {mask_file}"
            )
        mask_data = clean_array(nc_mask.variables[mask_variable][:])
        while mask_data.ndim > 2:
            mask_data = mask_data[0]
        mask = np.isfinite(mask_data)

    if lat2d.shape != lon2d.shape:
        raise ValueError("XLAT and XLONG shapes do not match")
    if mask.shape != (lat1d.size, lon1d.size):
        raise ValueError(
            "Mask shape does not match its latitude/longitude coordinates: "
            f"{mask.shape} vs {(lat1d.size, lon1d.size)}"
        )

    lon_grid, lat_grid = np.meshgrid(lon1d, lat1d)
    return lat2d, lon2d, lat1d, lon1d, lat_grid, lon_grid, mask


def build_time_indices(nyears, month_offsets):
    return [6 * year + month for year in range(nyears) for month in month_offsets]


def read_yearly_season(nc_obj, variable, nyears, month_offsets, factor):
    """Read and convert one variable to yearly seasonal means."""
    if variable not in nc_obj.variables:
        raise KeyError(f"Variable {variable!r} is absent from {nc_obj.filepath()}")

    nc_var = nc_obj.variables[variable]
    needed_steps = nyears * 6
    if nc_var.ndim != 3:
        raise ValueError(
            f"{variable} in {nc_obj.filepath()} must have 3 dimensions "
            f"(time, y, x); found {nc_var.dimensions}"
        )
    if nc_var.shape[0] < needed_steps:
        raise ValueError(
            f"{variable} in {nc_obj.filepath()} has {nc_var.shape[0]} time steps; "
            f"at least {needed_steps} are required for {nyears} years"
        )

    indices = build_time_indices(nyears, month_offsets)
    data = clean_array(nc_var[indices, :, :]) * factor
    return np.nanmean(
        data.reshape(nyears, len(month_offsets), data.shape[1], data.shape[2]),
        axis=1,
    )


def regrid_rcm_to_regular(field, lat2d, lon2d, lat1d, lon1d, method="linear"):
    """Interpolate a native Lambert-grid field to the 0.25-degree mask grid."""
    target_lon, target_lat = np.meshgrid(lon1d, lat1d)
    valid = np.isfinite(field) & np.isfinite(lat2d) & np.isfinite(lon2d)
    if np.count_nonzero(valid) < 3:
        return np.full(target_lon.shape, np.nan, dtype=np.float64)

    points = np.column_stack((lon2d[valid], lat2d[valid]))
    values = field[valid]
    interpolated = griddata(
        points, values, (target_lon, target_lat), method=method
    )

    if method == "linear" and np.isnan(interpolated).any():
        nearest = griddata(
            points, values, (target_lon, target_lat), method="nearest"
        )
        interpolated[np.isnan(interpolated)] = nearest[np.isnan(interpolated)]
    return interpolated


def calculate_response(
    nc_ctl,
    nc_exp,
    spec,
    nyears,
    month_offsets,
    lat2d,
    lon2d,
    lat1d,
    lon1d,
    land_mask,
):
    """Return regridded EXP-CTL climatology and paired-test p values."""
    ctl_yearly = read_yearly_season(
        nc_ctl, spec["variable"], nyears, month_offsets, spec["factor"]
    )
    exp_yearly = read_yearly_season(
        nc_exp, spec["variable"], nyears, month_offsets, spec["factor"]
    )

    if ctl_yearly.shape != exp_yearly.shape:
        raise ValueError(
            f"CTL and EXP shapes differ for {spec['variable']}: "
            f"{ctl_yearly.shape} vs {exp_yearly.shape}"
        )

    difference = np.nanmean(exp_yearly - ctl_yearly, axis=0)
    _, p_value = ttest_rel(
        exp_yearly, ctl_yearly, axis=0, nan_policy="omit"
    )

    difference_ip = regrid_rcm_to_regular(
        difference, lat2d, lon2d, lat1d, lon1d, method="linear"
    )
    # P values are categorical thresholds for plotting, so nearest-neighbour
    # transfer avoids manufacturing intermediate significance values.
    p_value_ip = regrid_rcm_to_regular(
        p_value, lat2d, lon2d, lat1d, lon1d, method="nearest"
    )

    return (
        np.where(land_mask, difference_ip, np.nan),
        np.where(land_mask, p_value_ip, np.nan),
    )


def validate_inputs(model_files, variable_specs, nyears):
    """Fail early with an actionable message before expensive interpolation."""
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
                    if nc_var.ndim != 3:
                        raise ValueError(
                            f"{variable} in {path} must be (time, y, x); "
                            f"found dimensions {nc_var.dimensions}"
                        )
                    if nc_var.shape[0] < needed_steps:
                        raise ValueError(
                            f"{variable} in {path} has {nc_var.shape[0]} time "
                            f"steps; {needed_steps} are required"
                        )
            print(f"[OK] {experiment} {role.upper()}: {path}")


def add_boundaries(ax, shapefile_dir):
    """Use the same China/province/South China Sea line style as existing plots."""
    for filename, linewidth in (
        ("province.shp", 0.4),
        ("china.shp", 0.6),
        ("south_china_sea.shp", 0.8),
    ):
        shape_path = shapefile_dir / filename
        if not shape_path.is_file():
            continue
        try:
            reader = shpreader.Reader(str(shape_path))
            ax.add_geometries(
                reader.geometries(),
                crs=ccrs.PlateCarree(),
                facecolor="none",
                edgecolor="black",
                linewidth=linewidth,
                zorder=3,
            )
        except Exception as exc:
            print(f"[WARN] Could not draw {shape_path}: {exc}")


def add_gridlines(ax, show_left, show_bottom):
    gridlines = ax.gridlines(
        crs=ccrs.PlateCarree(),
        draw_labels=True,
        x_inline=False,
        y_inline=False,
        linewidth=0.6,
        color="gray",
        alpha=0.5,
        linestyle="--",
        zorder=2,
    )
    gridlines.xlocator = mticker.FixedLocator(np.arange(70, 135, 10))
    gridlines.ylocator = mticker.FixedLocator(np.arange(10, 60, 10))
    gridlines.top_labels = False
    gridlines.right_labels = False
    gridlines.left_labels = show_left
    gridlines.bottom_labels = show_bottom
    gridlines.xlabel_style = {
        "size": 9,
        "rotation": 0,
        "ha": "center",
        "va": "top",
    }
    gridlines.ylabel_style = {
        "size": 9,
        "rotation": 0,
        "ha": "right",
        "va": "center",
    }


def weighted_summary(field, p_value, lat_grid, mask, threshold):
    valid = mask & np.isfinite(field)
    if not np.any(valid):
        return np.nan, np.nan, np.nan, np.nan
    weights = np.cos(np.deg2rad(lat_grid))
    mean_value = np.average(field[valid], weights=weights[valid])
    sig_valid = valid & np.isfinite(p_value)
    sig_fraction = (
        100.0 * np.count_nonzero(sig_valid & (p_value < threshold))
        / np.count_nonzero(sig_valid)
        if np.any(sig_valid)
        else np.nan
    )
    return (
        mean_value,
        np.nanmin(field[valid]),
        np.nanmax(field[valid]),
        sig_fraction,
    )


def run_category(variable_specs, output_filename, description):
    """Validate, calculate, and draw one supplementary variable category."""
    parser = build_parser(description)
    args = parser.parse_args()

    if args.nyears < 2:
        raise ValueError("--nyears must be at least 2 for a paired t test")
    if not 0.0 < args.p_threshold < 1.0:
        raise ValueError("--p-threshold must be between 0 and 1")
    if args.sig_stride < 1:
        raise ValueError("--sig-stride must be at least 1")

    model_files = resolve_model_files(args)
    validate_inputs(model_files, variable_specs, args.nyears)
    if args.check_only:
        print("Input validation completed; no figure was created.")
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

    projection = ccrs.LambertConformal(
        central_longitude=PROJ_CENTRAL_LON,
        central_latitude=PROJ_CENTRAL_LAT,
        standard_parallels=PROJ_STD_PARALLELS,
    )

    nrows = len(variable_specs)
    ncols = 4
    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(17, 3.15 * nrows),
        subplot_kw={"projection": projection},
        squeeze=False,
    )
    plt.subplots_adjust(
        wspace=0.035,
        hspace=0.10,
        left=0.075,
        right=0.89,
        bottom=0.07,
        top=0.93,
    )

    columns = (
        ("OFF", "Spring", SEASONS[0][1]),
        ("OFF", "Summer", SEASONS[1][1]),
        ("CPL", "Spring", SEASONS[0][1]),
        ("CPL", "Summer", SEASONS[1][1]),
    )
    column_titles = (
        "(EXP - CTL) OFF | Spring",
        "(EXP - CTL) OFF | Summer",
        "(EXP - CTL) CPL | Spring",
        "(EXP - CTL) CPL | Summer",
    )

    row_mappables = []
    panel_index = 0

    open_files = {}
    try:
        for experiment, paths in model_files.items():
            open_files[experiment] = {
                "ctl": Dataset(paths["ctl"]),
                "exp": Dataset(paths["exp"]),
            }

        for row_index, spec in enumerate(variable_specs):
            row_mappable = None
            for column_index, (experiment, season_name, offsets) in enumerate(columns):
                print(
                    f"Processing {spec['variable']} | {experiment} | {season_name}..."
                )
                response, p_value = calculate_response(
                    open_files[experiment]["ctl"],
                    open_files[experiment]["exp"],
                    spec,
                    args.nyears,
                    offsets,
                    lat2d,
                    lon2d,
                    lat1d,
                    lon1d,
                    land_mask,
                )

                mean_value, minimum, maximum, sig_fraction = weighted_summary(
                    response,
                    p_value,
                    lat_grid,
                    land_mask,
                    args.p_threshold,
                )
                print(
                    f"  China land: mean={mean_value:.4g}, min={minimum:.4g}, "
                    f"max={maximum:.4g} {spec['unit']}; "
                    f"p<{args.p_threshold:g}: {sig_fraction:.1f}%"
                )

                plot_field = response
                if args.significance_style == "mask":
                    plot_field = np.where(
                        p_value < args.p_threshold, response, np.nan
                    )

                ax = axes[row_index, column_index]
                ax.set_extent(PLOT_EXTENT, crs=ccrs.PlateCarree())
                add_boundaries(ax, args.shapefile_dir)
                add_gridlines(
                    ax,
                    show_left=(column_index == 0),
                    show_bottom=(row_index == nrows - 1),
                )

                row_mappable = ax.contourf(
                    lon_grid,
                    lat_grid,
                    plot_field,
                    levels=spec["levels"],
                    cmap=spec.get("cmap", "RdBu_r"),
                    extend="both",
                    transform=ccrs.PlateCarree(),
                    zorder=1,
                    antialiased=False,
                )

                if args.significance_style == "stipple":
                    stride = args.sig_stride
                    sig = (
                        np.isfinite(p_value)
                        & (p_value < args.p_threshold)
                        & land_mask
                    )
                    sampled_sig = sig[::stride, ::stride]
                    ax.scatter(
                        lon_grid[::stride, ::stride][sampled_sig],
                        lat_grid[::stride, ::stride][sampled_sig],
                        s=3.0,
                        marker="o",
                        color="0.35",
                        edgecolors="none",
                        transform=ccrs.PlateCarree(),
                        zorder=2,
                        rasterized=True,
                    )

                panel_label = chr(97 + panel_index)
                ax.text(
                    0.025,
                    0.965,
                    f"({panel_label})",
                    transform=ax.transAxes,
                    fontsize=10.5,
                    va="top",
                    ha="left",
                    zorder=5,
                    bbox={
                        "facecolor": "white",
                        "alpha": 0.70,
                        "edgecolor": "none",
                        "pad": 1.0,
                    },
                )
                if row_index == 0:
                    ax.set_title(
                        column_titles[column_index],
                        fontsize=11,
                        loc="center",
                        pad=6,
                    )
                panel_index += 1

            axes[row_index, 0].text(
                -0.19,
                0.5,
                spec["title"],
                transform=axes[row_index, 0].transAxes,
                rotation=90,
                fontsize=12,
                va="center",
                ha="center",
            )
            row_mappables.append(row_mappable)
    finally:
        for experiment_files in open_files.values():
            for nc_obj in experiment_files.values():
                nc_obj.close()

    fig.canvas.draw()
    for row_index, (spec, mappable) in enumerate(
        zip(variable_specs, row_mappables)
    ):
        position = axes[row_index, -1].get_position()
        cbar_axis = fig.add_axes(
            [0.905, position.y0, 0.012, position.height]
        )
        colorbar = fig.colorbar(
            mappable,
            cax=cbar_axis,
            orientation="vertical",
            extend="both",
            ticks=spec.get("ticks"),
        )
        colorbar.ax.tick_params(labelsize=8)
        colorbar.set_label(spec["unit"], fontsize=9, labelpad=4)

    period_end = args.start_year + args.nyears - 1
    fig.text(
        0.5,
        0.015,
        (
            f"{args.start_year}-{period_end} climatological response; "
            f"paired two-sided t test, p < {args.p_threshold:g}"
        ),
        ha="center",
        va="bottom",
        fontsize=9,
    )

    output_path = args.output_dir / output_filename
    fig.savefig(output_path, dpi=args.dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {output_path}")
