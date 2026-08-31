#!/usr/bin/env python3
"""Evaluate infiltration in OFF and CPL experiments."""

from _land_response_common import run_category


VARIABLES = [
    {
        "key": "infiltration",
        "variable": "f_qinfl",
        "title": "Infiltration",
        "unit": r"mm day$^{-1}$",
        "factor": 86400.0,
        "output_stem": "supp_land_infiltration",
        "nonnegative": True,
        "state_cmap": "YlGnBu",
        "difference_cmap": "BrBG",
        "difference_fallback": 0.05,
        "reference_product": "MERRA-2 QINFIL reanalysis benchmark",
        "reference_globs": (
            "/share/home/dq135/openbench/Reference/Grid/LowRes/Water/"
            "Infiltration/MERRA2/*.nc*",
        ),
        "reference_variables": ("QINFIL", "qinfil", "infiltration"),
        "reference_quantity": "rate_mm_day",
    }
]


if __name__ == "__main__":
    run_category(
        VARIABLES,
        "Reference-based infiltration evaluation and gravel-response figures.",
    )
