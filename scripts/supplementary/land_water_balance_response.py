#!/usr/bin/env python3
"""Evaluate infiltration in OFF and CPL experiments."""

from _land_response_common import run_response_only_category


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
    }
]


if __name__ == "__main__":
    run_response_only_category(
        VARIABLES,
        "Reference-free infiltration response diagnostics for gravel experiments.",
    )
