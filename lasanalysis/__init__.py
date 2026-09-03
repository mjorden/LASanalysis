"""Quick-look petrophysics for LAS well logs.

Modules
-------
load       read a LAS file, mask undeclared sentinels, normalise mnemonics
petro      pure-numpy petrophysics: Vshale, density porosity, crossover, Archie Sw,
           Pickett water-line fit
plot       track plots, neutron-density crossplot, Pickett plot
kgs        search the Kansas Geological Survey LAS index, fetch LAS files by KID
multiwell  run the workflow over many wells (search -> fetch -> analyse -> summary.csv)
"""

from .load import ALIASES, clean_frame, curves, find_curve, read_las, read_log_csv, standardize
from .petro import (
    MATRIX_DENSITY,
    NEUTRON_MATRIX_OFFSET,
    SW_MODELS,
    archie_sw,
    density_porosity,
    fit_water_line,
    matrix_density,
    neutron_density_crossover,
    neutron_lithology_correction,
    pick_rsh,
    pick_rw_from_rwa,
    rwa,
    sw_indonesia,
    sw_simandoux,
    vshale_larionov,
    vshale_linear,
    water_saturation,
)
from .plot import crossplot_neutron_density, pickett_plot, plot_tracks

__all__ = [
    "ALIASES",
    "MATRIX_DENSITY",
    "NEUTRON_MATRIX_OFFSET",
    "SW_MODELS",
    "archie_sw",
    "matrix_density",
    "neutron_lithology_correction",
    "pick_rsh",
    "pick_rw_from_rwa",
    "rwa",
    "sw_indonesia",
    "sw_simandoux",
    "water_saturation",
    "clean_frame",
    "crossplot_neutron_density",
    "curves",
    "density_porosity",
    "find_curve",
    "fit_water_line",
    "neutron_density_crossover",
    "pickett_plot",
    "plot_tracks",
    "read_las",
    "read_log_csv",
    "standardize",
    "vshale_larionov",
    "vshale_linear",
]

__version__ = "0.1.0"
