#!/usr/bin/env python3
"""Plot spatial total-evapotranspiration evaluation."""

from _land_spatial_diagnostics import run_reference_spatial_category


VARIABLES = [
    {
        "key": "et",
        "variable": "f_fevpa",
        "title": "Total evapotranspiration",
        "unit": r"mm day$^{-1}$",
        "factor": 86400.0,
        "output_stem": "supp_land_et",
        "nonnegative": True,
        "state_cmap": "YlGnBu",
        "difference_cmap": "BrBG",
        "difference_fallback": 0.05,
        "reference_product": "GLEAM v4.2a",
        "reference_short_name": "GLEAM v4.2a",
        "reference_globs": (
            "/share/home/dq135/openbench/Reference/Grid/HigRes/Water/"
            "Evapotranspiration/GLEAM_v4.2a/E_*_GLEAM_v4.2a.nc",
        ),
        "reference_variables": ("E", "ET", "evapotranspiration"),
        "reference_quantity": "rate_mm_day",
    }
]


if __name__ == "__main__":
    run_reference_spatial_category(
        VARIABLES,
        "Spatial-only total-ET evaluation and gravel response.",
    )
