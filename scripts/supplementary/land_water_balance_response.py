#!/usr/bin/env python3
"""Plot spatial infiltration diagnostics without a reference product."""

from _land_spatial_diagnostics import run_response_only_spatial_category


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
    run_response_only_spatial_category(
        VARIABLES,
        "Spatial-only infiltration response diagnostics; no reference product.",
    )
