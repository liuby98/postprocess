#!/usr/bin/env python3
"""Evaluate SWE and snow-cover fraction in OFF and CPL experiments."""

from _land_response_common import run_category


VARIABLES = [
    {
        "key": "swe",
        "variable": "f_scv",
        "title": "Snow water equivalent",
        "unit": "mm",
        "factor": 1.0,
        "output_stem": "supp_land_swe",
        "nonnegative": True,
        "state_cmap": "Blues",
        "difference_cmap": "RdBu_r",
        "difference_fallback": 1.0,
        "reference_product": "ERA5-Land SWE benchmark",
        "reference_globs": (
            "/share/home/dq135/openbench/Reference/Grid/HigRes/Snow/"
            "Snow_Water_Equivalent/ERA5_Land/*.nc*",
        ),
        "reference_variables": (
            "sd",
            "SWE",
            "swe",
            "SNOMAS",
            "snow_water_equivalent",
        ),
        "reference_quantity": "amount_mm",
    },
    {
        "key": "scf",
        "variable": "f_fsno",
        "title": "Snow-cover fraction",
        "unit": "percentage points",
        "factor": 100.0,
        "output_stem": "supp_land_scf",
        "nonnegative": True,
        "state_cmap": "Blues",
        "difference_cmap": "RdBu_r",
        "difference_fallback": 1.0,
        "reference_product": "MODIS Terra MOD10CM Collection 6.1",
        "reference_globs": (
            "/share/home/dq135/openbench/Reference/Grid/HigRes/Snow/"
            "Snow_Cover_Fraction/MODIS_MOD10CM/*.nc*",
        ),
        "reference_variables": (
            "Snow_Cover_Monthly_CMG",
            "NDSI_Snow_Cover",
            "snow_cover_fraction",
            "SCF",
            "FRSNO",
        ),
        "reference_quantity": "percent",
    },
]


if __name__ == "__main__":
    run_category(
        VARIABLES,
        "Reference-based SWE/SCF evaluation and gravel-response figures.",
    )
