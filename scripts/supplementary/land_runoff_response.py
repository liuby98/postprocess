#!/usr/bin/env python3
"""Plot offline and coupled gravel responses of runoff components."""

import numpy as np

from _land_response_common import run_category


MM_S_TO_MM_DAY = 86400.0

VARIABLES = [
    {
        "variable": "f_rsur",
        "title": "Surface runoff",
        "unit": r"mm day$^{-1}$",
        "factor": MM_S_TO_MM_DAY,
        "levels": np.linspace(-0.8, 0.8, 25),
        "ticks": np.arange(-0.8, 0.81, 0.4),
        "cmap": "BrBG",
    },
    {
        "variable": "f_rsub",
        "title": "Subsurface runoff",
        "unit": r"mm day$^{-1}$",
        "factor": MM_S_TO_MM_DAY,
        "levels": np.linspace(-0.4, 0.4, 25),
        "ticks": np.arange(-0.4, 0.41, 0.2),
        "cmap": "BrBG",
    },
    {
        "variable": "f_rnof",
        "title": "Total runoff",
        "unit": r"mm day$^{-1}$",
        "factor": MM_S_TO_MM_DAY,
        "levels": np.linspace(-0.8, 0.8, 25),
        "ticks": np.arange(-0.8, 0.81, 0.4),
        "cmap": "BrBG",
    },
]


if __name__ == "__main__":
    run_category(
        VARIABLES,
        "supp_land_runoff_response.pdf",
        "Runoff response to gravel in offline and coupled simulations.",
    )
