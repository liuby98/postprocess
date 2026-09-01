#!/usr/bin/env python3
"""Plot spatial total-runoff evaluation."""

from _land_spatial_diagnostics import run_reference_spatial_category


VARIABLES = [
    {
        "key": "runoff",
        "variable": "f_rnof",
        "title": "Total runoff",
        "unit": r"mm day$^{-1}$",
        "factor": 86400.0,
        "output_stem": "supp_land_total_runoff",
        "nonnegative": True,
        "state_cmap": "YlGnBu",
        "difference_cmap": "BrBG",
        "difference_fallback": 0.05,
        "reference_product": "G-RUN ENSEMBLE MMM",
        "reference_short_name": "G-RUN ENSEMBLE MMM",
        "reference_globs": (
            "/share/home/dq135/openbench/Reference/Grid/LowRes/Water/"
            "Total_Runoff/G_RUN_ENSEMBLE/*.nc*",
        ),
        "reference_variables": (
            "runoff",
            "Runoff",
            "RUNOFF",
            "mrro",
            "runoff_mean",
            "runoff_median",
            "qtot",
            "GRUN",
            "ensemble_mean",
            "G_RUN_ENSEMBLE_MMM",
        ),
        "reference_quantity": "rate_mm_day",
    }
]


if __name__ == "__main__":
    run_reference_spatial_category(
        VARIABLES,
        "Spatial-only total-runoff evaluation and gravel response.",
    )
