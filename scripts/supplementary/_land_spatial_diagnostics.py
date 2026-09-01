#!/usr/bin/env python3
"""Spatial-only land-surface diagnostics with compact, aligned layouts."""

import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
from scipy.stats import ttest_rel

import cartopy.crs as ccrs

from _land_response_common import (
    PLOT_EXTENT,
    PROJ_CENTRAL_LAT,
    PROJ_CENTRAL_LON,
    PROJ_STD_PARALLELS,
    SEASONS,
    add_boundaries,
    align_model_years_to_reference,
    build_parser,
    common_plot_data,
    field_metrics,
    inspect_reference_files,
    load_model_yearly,
    load_plot_grid,
    load_reference_yearly,
    model_response_masks,
    reference_product,
    resolve_model_files,
    resolve_reference_files,
    robust_state_limits,
    robust_symmetric_limit,
    validate_arguments,
    validate_model_inputs,
    weighted_mean,
)


def _projection():
    return ccrs.LambertConformal(
        central_longitude=PROJ_CENTRAL_LON,
        central_latitude=PROJ_CENTRAL_LAT,
        standard_parallels=PROJ_STD_PARALLELS,
    )


def _draw_panel(
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
    """Draw one map with explicitly horizontal longitude/latitude labels."""
    ax.set_extent(PLOT_EXTENT, crs=ccrs.PlateCarree())
    add_boundaries(ax, shapefile_dir)
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
    gridlines.xlocator = mticker.FixedLocator(np.arange(80, 131, 10))
    gridlines.ylocator = mticker.FixedLocator(np.arange(10, 60, 10))
    gridlines.top_labels = False
    gridlines.right_labels = False
    gridlines.left_labels = show_left
    gridlines.bottom_labels = show_bottom
    gridlines.rotate_labels = False
    gridlines.xlabel_style = {
        "size": 7.5,
        "rotation": 0,
        "ha": "center",
        "va": "top",
    }
    gridlines.ylabel_style = {
        "size": 7.5,
        "rotation": 0,
        "ha": "right",
        "va": "center",
    }
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


def _panel_text(fig, axes, headers):
    """Use figure-level text so Cartopy does not crop panel labels."""
    fig.canvas.draw()
    for col, header in enumerate(headers):
        position = axes[0, col].get_position()
        fig.text(
            position.x0 + position.width / 2,
            position.y1 + 0.012,
            header,
            ha="center",
            va="bottom",
            fontsize=10.2,
        )
    ncols = axes.shape[1]
    for row in range(axes.shape[0]):
        for col in range(ncols):
            position = axes[row, col].get_position()
            fig.text(
                position.x0 + 0.012 * position.width,
                position.y1 - 0.035 * position.height,
                f"({chr(97 + row * ncols + col)})",
                ha="left",
                va="top",
                fontsize=9.2,
                bbox={
                    "facecolor": "white",
                    "alpha": 0.70,
                    "edgecolor": "none",
                    "pad": 1,
                },
                zorder=5,
            )


def _add_colorbar(fig, mappable, left_axis, right_axis, y, label, ticks):
    """Make a colorbar exactly span the selected map-column boundaries."""
    left = left_axis.get_position().x0
    right = right_axis.get_position().x1
    cax = fig.add_axes([left, y, right - left, 0.022])
    colorbar = fig.colorbar(
        mappable,
        cax=cax,
        orientation="horizontal",
        extendfrac="auto",
    )
    colorbar.set_ticks(ticks)
    colorbar.ax.xaxis.set_major_formatter(mticker.StrMethodFormatter("{x:g}"))
    colorbar.ax.tick_params(labelsize=7.5)
    colorbar.set_label(label, fontsize=8.5, labelpad=2)
    return colorbar


def plot_reference_spatial(
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
    """Plot Reference, two model biases, and EXP-CTL in a 2 x 4 layout."""
    state_fields = [reference[name] for name, _ in SEASONS]
    bias_fields = []
    response_fields = []
    for season_name, _ in SEASONS:
        ref_mean = np.nanmean(reference[season_name], axis=0)
        for experiment in ("OFF", "CPL"):
            mask = masks[experiment][season_name]
            ctl_mean = np.nanmean(model[experiment]["ctl"][season_name], axis=0)
            exp_mean = np.nanmean(model[experiment]["exp"][season_name], axis=0)
            bias_fields.extend(
                [
                    np.where(mask, ctl_mean - ref_mean, np.nan),
                    np.where(mask, exp_mean - ref_mean, np.nan),
                ]
            )
            response_fields.append(np.where(mask, exp_mean - ctl_mean, np.nan))

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
        f"Color limits {spec['title']}: state={state_min:.4g}..{state_max:.4g}, "
        f"bias=+/-{bias_limit:.4g}, response=+/-{response_limit:.4g} "
        f"{spec['unit']}"
    )

    reference_label = spec.get("reference_short_name", product)
    headers = (
        reference_label,
        f"CTL - {reference_label}",
        f"EXP - {reference_label}",
        "EXP - CTL",
    )
    for experiment in ("OFF", "CPL"):
        fig, axes = plt.subplots(
            2,
            4,
            figsize=(14.2, 7.5),
            subplot_kw={"projection": _projection()},
            squeeze=False,
        )
        plt.subplots_adjust(
            left=0.065,
            right=0.988,
            bottom=0.185,
            top=0.88,
            wspace=0.045,
            hspace=0.10,
        )
        state_mappable = bias_mappable = response_mappable = None
        for row, (season_name, _) in enumerate(SEASONS):
            mask = masks[experiment][season_name]
            ref_yearly = np.where(mask[None], reference[season_name], np.nan)
            ctl_yearly = np.where(
                mask[None], model[experiment]["ctl"][season_name], np.nan
            )
            exp_yearly = np.where(
                mask[None], model[experiment]["exp"][season_name], np.nan
            )
            ref_mean = np.nanmean(ref_yearly, axis=0)
            ctl_mean = np.nanmean(ctl_yearly, axis=0)
            exp_mean = np.nanmean(exp_yearly, axis=0)
            fields = (
                ref_mean,
                ctl_mean - ref_mean,
                exp_mean - ref_mean,
                exp_mean - ctl_mean,
            )
            _, p_value = ttest_rel(
                exp_yearly, ctl_yearly, axis=0, nan_policy="omit"
            )
            for col, field in enumerate(fields):
                plot_field = field
                if col == 3 and args.significance_style == "mask":
                    plot_field = np.where(
                        p_value < args.p_threshold, field, np.nan
                    )
                if col == 0:
                    levels = state_levels
                    cmap = spec.get("state_cmap", "YlGnBu")
                elif col < 3:
                    levels = bias_levels
                    cmap = spec.get("difference_cmap", "BrBG")
                else:
                    levels = response_levels
                    cmap = spec.get("difference_cmap", "BrBG")
                mappable = _draw_panel(
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
                if col == 0:
                    state_mappable = mappable
                elif col < 3:
                    bias_mappable = mappable
                    simulation = ctl_mean if col == 1 else exp_mean
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
                else:
                    response_mappable = mappable
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

        _panel_text(fig, axes, headers)
        period_end = args.start_year + args.nyears - 1
        fig.suptitle(
            f"{spec['title']} - {experiment} evaluation "
            f"({args.start_year}-{period_end})",
            fontsize=14,
            y=0.96,
        )
        bar_y = axes[1, 0].get_position().y0 - 0.072
        _add_colorbar(
            fig,
            state_mappable,
            axes[1, 0],
            axes[1, 0],
            bar_y,
            f"{reference_label} ({spec['unit']})",
            np.linspace(state_min, state_max, 5),
        )
        _add_colorbar(
            fig,
            bias_mappable,
            axes[1, 1],
            axes[1, 2],
            bar_y,
            f"Bias ({spec['unit']})",
            np.linspace(-bias_limit, bias_limit, 7),
        )
        _add_colorbar(
            fig,
            response_mappable,
            axes[1, 3],
            axes[1, 3],
            bar_y,
            f"Response ({spec['unit']})",
            np.linspace(-response_limit, response_limit, 5),
        )
        fig.text(
            0.5,
            0.015,
            (
                f"Reference: {product}; unified China-land mask. State: "
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


def plot_response_spatial(
    spec, model, masks, args, lat_grid, lon_grid, land_mask
):
    """Plot CTL, EXP, and EXP-CTL for a variable without a reference."""
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
    for experiment in ("OFF", "CPL"):
        fig, axes = plt.subplots(
            2,
            3,
            figsize=(11.3, 7.5),
            subplot_kw={"projection": _projection()},
            squeeze=False,
        )
        plt.subplots_adjust(
            left=0.075,
            right=0.988,
            bottom=0.185,
            top=0.88,
            wspace=0.045,
            hspace=0.10,
        )
        state_mappable = response_mappable = None
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
                levels = state_levels if col < 2 else response_levels
                cmap = (
                    spec.get("state_cmap", "YlGnBu")
                    if col < 2
                    else spec.get("difference_cmap", "BrBG")
                )
                mappable = _draw_panel(
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
                    area_change = weighted_mean(response, lat_grid, mask)
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

        _panel_text(fig, axes, ("CTL", "EXP", "EXP - CTL"))
        period_end = args.start_year + args.nyears - 1
        fig.suptitle(
            f"{spec['title']} - {experiment} process response "
            f"({args.start_year}-{period_end})",
            fontsize=14,
            y=0.96,
        )
        bar_y = axes[1, 0].get_position().y0 - 0.072
        _add_colorbar(
            fig,
            state_mappable,
            axes[1, 0],
            axes[1, 1],
            bar_y,
            f"State ({spec['unit']})",
            np.linspace(state_min, state_max, 7),
        )
        _add_colorbar(
            fig,
            response_mappable,
            axes[1, 2],
            axes[1, 2],
            bar_y,
            f"Response ({spec['unit']})",
            np.linspace(-response_limit, response_limit, 5),
        )
        fig.text(
            0.5,
            0.015,
            (
                "Process diagnostic only; no reference evaluation. State: "
                f"P{args.state_quantiles[0]:g}-P{args.state_quantiles[1]:g}; "
                f"response: P{args.difference_quantile:g}(|x|). "
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


def _load_common_grid(args):
    return load_plot_grid(args.wrfinput, args.mask_file, args.mask_variable)


def run_reference_spatial_category(variable_specs, description):
    """Validate reference/model inputs and create spatial figures only."""
    parser = build_parser(description, variable_specs, reference_mode=True)
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
            "All model and reference inputs passed validation; spatial-only "
            "mode created no figure."
        )
        return
    lat2d, lon2d, lat1d, lon1d, lat_grid, lon_grid, land_mask = (
        _load_common_grid(args)
    )
    args.output_dir.mkdir(parents=True, exist_ok=True)
    os.environ["CARTOPY_OFFLINE"] = "true"
    for spec in variable_specs:
        product = reference_product(args, spec)
        print(f"\n=== {spec['title']} | spatial reference: {product} ===")
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
        align_model_years_to_reference(reference, model, years, spec)
        masks = common_plot_data(reference, model, land_mask)
        if not args.skip_spatial:
            plot_reference_spatial(
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


def run_response_only_spatial_category(variable_specs, description):
    """Validate model inputs and create reference-free spatial figures only."""
    parser = build_parser(description, [], reference_mode=False)
    args = parser.parse_args()
    validate_arguments(args)
    model_files = resolve_model_files(args)
    validate_model_inputs(model_files, variable_specs, args.nyears)
    if args.check_only:
        print(
            "All model inputs passed validation; no reference data are "
            "required and spatial-only mode created no figure."
        )
        return
    lat2d, lon2d, lat1d, lon1d, lat_grid, lon_grid, land_mask = (
        _load_common_grid(args)
    )
    args.output_dir.mkdir(parents=True, exist_ok=True)
    os.environ["CARTOPY_OFFLINE"] = "true"
    for spec in variable_specs:
        print(f"\n=== {spec['title']} | spatial process response only ===")
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
        if not args.skip_spatial:
            plot_response_spatial(
                spec,
                model,
                masks,
                args,
                lat_grid,
                lon_grid,
                land_mask,
            )
