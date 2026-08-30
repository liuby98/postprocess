#!/usr/bin/env python3
"""Plot offline and coupled gravel responses of evapotranspiration terms."""

import numpy as np

from _land_response_common import run_category


MM_S_TO_MM_DAY = 86400.0

VARIABLES = [
    {
        "variable": "f_fevpa",
        "title": "Total evapotranspiration",
        "unit": r"mm day$^{-1}$",
        "factor": MM_S_TO_MM_DAY,
        "levels": np.linspace(-0.6, 0.6, 25),
        "ticks": np.arange(-0.6, 0.61, 0.3),
        "cmap": "BrBG",
    },
    {
        "variable": "f_fevpl",
        "title": "Leaf evaporation + transpiration",
        "unit": r"mm day$^{-1}$",
        "factor": MM_S_TO_MM_DAY,
        "levels": np.linspace(-0.5, 0.5, 25),
        "ticks": np.arange(-0.5, 0.51, 0.25),
        "cmap": "BrBG",
    },
    {
        "variable": "f_etr",
        "title": "Transpiration",
        "unit": r"mm day$^{-1}$",
        "factor": MM_S_TO_MM_DAY,
        "levels": np.linspace(-0.4, 0.4, 25),
        "ticks": np.arange(-0.4, 0.41, 0.2),
        "cmap": "BrBG",
    },
    {
        "variable": "f_fevpg",
        "title": "Ground evaporation",
        "unit": r"mm day$^{-1}$",
        "factor": MM_S_TO_MM_DAY,
        "levels": np.linspace(-0.4, 0.4, 25),
        "ticks": np.arange(-0.4, 0.41, 0.2),
        "cmap": "BrBG",
    },
]


if __name__ == "__main__":
    run_category(
        VARIABLES,
        "supp_land_evapotranspiration_response.pdf",
        "Evapotranspiration response to gravel in offline and coupled simulations.",
    )
