#!/usr/bin/env python3
"""Plot offline and coupled gravel responses of water-balance diagnostics."""

import numpy as np

from _land_response_common import run_category


MM_S_TO_MM_DAY = 86400.0

VARIABLES = [
    {
        "variable": "f_qinfl",
        "title": "Infiltration",
        "unit": r"mm day$^{-1}$",
        "factor": MM_S_TO_MM_DAY,
        "levels": np.linspace(-1.0, 1.0, 25),
        "ticks": np.arange(-1.0, 1.01, 0.5),
        "cmap": "BrBG",
    },
    {
        "variable": "f_qcharge",
        "title": "Groundwater recharge",
        "unit": r"mm day$^{-1}$",
        "factor": MM_S_TO_MM_DAY,
        "levels": np.linspace(-0.4, 0.4, 25),
        "ticks": np.arange(-0.4, 0.41, 0.2),
        "cmap": "BrBG",
    },
    {
        "variable": "f_wat",
        "title": "Total water storage",
        "unit": "mm",
        "factor": 1.0,
        "levels": np.linspace(-40, 40, 25),
        "ticks": np.arange(-40, 41, 20),
        "cmap": "BrBG",
    },
    {
        "variable": "f_zwt",
        "title": "Water-table depth",
        "unit": "m",
        "factor": 1.0,
        "levels": np.linspace(-0.8, 0.8, 25),
        "ticks": np.arange(-0.8, 0.81, 0.4),
        "cmap": "BrBG",
    },
]


if __name__ == "__main__":
    run_category(
        VARIABLES,
        "supp_land_water_balance_response.pdf",
        "Water-balance response to gravel in offline and coupled simulations.",
    )
