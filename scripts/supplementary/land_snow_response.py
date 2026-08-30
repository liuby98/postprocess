#!/usr/bin/env python3
"""Plot offline and coupled gravel responses of snow state variables."""

import numpy as np

from _land_response_common import run_category


VARIABLES = [
    {
        "variable": "f_scv",
        "title": "Snow water equivalent",
        "unit": "mm",
        "factor": 1.0,
        "levels": np.linspace(-20, 20, 33),
        "ticks": np.arange(-20, 21, 10),
        "cmap": "RdBu_r",
    },
    {
        "variable": "f_snowdp",
        "title": "Snow depth",
        "unit": "cm",
        "factor": 100.0,
        "levels": np.linspace(-6, 6, 25),
        "ticks": np.arange(-6, 7, 3),
        "cmap": "RdBu_r",
    },
    {
        "variable": "f_fsno",
        "title": "Snow-cover fraction",
        "unit": "percentage points",
        "factor": 100.0,
        "levels": np.linspace(-12, 12, 25),
        "ticks": np.arange(-12, 13, 6),
        "cmap": "RdBu_r",
    },
]


if __name__ == "__main__":
    run_category(
        VARIABLES,
        "supp_land_snow_response.pdf",
        "Snow-state response to gravel in offline and coupled simulations.",
    )
