#!/usr/bin/env python3
"""Plot offline and coupled gravel responses of surface energy partitioning."""

import numpy as np

from _land_response_common import run_category


VARIABLES = [
    {
        "variable": "f_rnet",
        "title": "Net radiation",
        "unit": r"W m$^{-2}$",
        "factor": 1.0,
        "levels": np.linspace(-6, 6, 25),
        "ticks": np.arange(-6, 7, 3),
        "cmap": "PRGn",
    },
    {
        "variable": "f_fsena",
        "title": "Sensible heat flux",
        "unit": r"W m$^{-2}$",
        "factor": 1.0,
        "levels": np.linspace(-6, 6, 25),
        "ticks": np.arange(-6, 7, 3),
        "cmap": "PRGn",
    },
    {
        "variable": "f_lfevpa",
        "title": "Latent heat flux",
        "unit": r"W m$^{-2}$",
        "factor": 1.0,
        "levels": np.linspace(-6, 6, 25),
        "ticks": np.arange(-6, 7, 3),
        "cmap": "PRGn",
    },
    {
        "variable": "f_fgrnd",
        "title": "Ground heat flux",
        "unit": r"W m$^{-2}$",
        "factor": 1.0,
        "levels": np.linspace(-6, 6, 25),
        "ticks": np.arange(-6, 7, 3),
        "cmap": "PRGn",
    },
    {
        "variable": "f_olrg",
        "title": "Outgoing longwave radiation",
        "unit": r"W m$^{-2}$",
        "factor": 1.0,
        "levels": np.linspace(-6, 6, 25),
        "ticks": np.arange(-6, 7, 3),
        "cmap": "PRGn",
    },
]


if __name__ == "__main__":
    run_category(
        VARIABLES,
        "supp_land_energy_partition_response.pdf",
        "Surface-energy response to gravel in offline and coupled simulations.",
    )
