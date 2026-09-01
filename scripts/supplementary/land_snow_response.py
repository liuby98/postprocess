#!/usr/bin/env python3
"""Plot grid-box and CoLM-snow-mask SWE/SCF diagnostics."""

from _snow_monthly_diagnostics import run_snow_monthly_diagnostics


VARIABLES = [
    {
        "key": "swe",
        "variable": "f_scv",
        "title": "Snow water equivalent",
        "unit": "mm",
        "factor": 1.0,
        "nonnegative": True,
        "reference_product": "ERA5-Land SWE benchmark",
        "reference_short_name": "ERA5-Land",
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
        "nonnegative": True,
        "reference_product": "MODIS Terra MOD10CM Collection 6.1",
        "reference_short_name": "MOD10CM C6.1",
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
    run_snow_monthly_diagnostics(
        VARIABLES,
        (
            "Monthly grid-box/CoLM-snow-mask SWE, SCF, and SWE-range "
            "diagnostics for gravel > 0.3."
        ),
    )
