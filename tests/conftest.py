"""
This provides standard configurations for tomplot tests. This takes the form
of some analytic profiles that are useful to consider when plotting.
"""

import numpy as np
from collections import namedtuple
import pytest
import matplotlib.pyplot as plt
from os.path import abspath, dirname


# ---------------------------------------------------------------------------- #
# Set up config object
# ---------------------------------------------------------------------------- #
# Config option just holds a bunch of properties
opts = ('coords_X', 'coords_Y', 'field_data', 'colour_scheme',
        'tracer_background', 'tracer_max', 'show_plots', 'save_plots',
        'overwrite_plots')
PlotSetup = namedtuple('PlotSetup', opts)
PlotSetup.__new__.__defaults__ = (None,)*len(opts)


# ---------------------------------------------------------------------------- #
# Add routine for making plots to config object
# ---------------------------------------------------------------------------- #

def make_plots(self, plot_name):
    path_to_here = abspath(dirname(__file__))

    if self.show_plots:
        plt.show()
    if self.save_plots:
        plot_path = f'{path_to_here}/tmp_figures/{plot_name}'
        plt.savefig(plot_path, bbox_inches='tight')
    elif self.overwrite_plots:
        plot_path = f'{path_to_here}/figures/{plot_name}'
        plt.savefig(plot_path, bbox_inches='tight')
    return None


setattr(PlotSetup, "make_plots", make_plots)


# ---------------------------------------------------------------------------- #
# Initial conditions
# ---------------------------------------------------------------------------- #

# An initial condition to be plotted with a divergent colour map
def dipole_initial_condition(x, y, Lx, Ly, tracer_background, tracer_max):
    r = np.sqrt(((x - Lx/2) / (Lx/5))**2 + ((y - Ly/2) / (Ly/5))**2)
    return 2*((x - Lx/2) / (Lx/5))*(tracer_max - tracer_background) * np.exp(-r**2)


# An initial condition to be plotted with a linear colour map, but featuring
# a small perturbation that we may wish not to plot
def two_gaussian_initial_condition(x, y, Lx, Ly, tracer_background, tracer_max):
    r1 = np.sqrt(((x - Lx/2) / (Lx/5))**2 + ((y - Ly/3) / (Ly/5))**2)
    r2 = np.sqrt(((x - 3*Lx/4) / (Lx/10))**2 + ((y - 3*Ly/4) / (Ly/10))**2)
    return (tracer_background + (tracer_max - tracer_background) * np.exp(-r1**2)
            + (tracer_max - tracer_background) / 20.0 * np.exp(-r2**2))


# ---------------------------------------------------------------------------- #
# Specific PlotSetup objects
# ---------------------------------------------------------------------------- #

def dipole(tracer_background, tracer_max, show_plots,
           save_plots, overwrite_plots, high_res=False):

    npoints_1d = 40 if high_res else 15
    Lx = Ly = 10.0
    x_1d = np.linspace(0, Lx, npoints_1d)
    y_1d = np.linspace(0, Ly, npoints_1d)
    coords_X, coords_Y = np.meshgrid(x_1d, y_1d, indexing='ij')

    field_data = dipole_initial_condition(coords_X, coords_Y, Lx, Ly,
                                          tracer_background, tracer_max)

    colour_scheme = 'RdBu_r'

    return PlotSetup(coords_X, coords_Y, field_data, colour_scheme,
                     tracer_background, tracer_max, show_plots, save_plots,
                     overwrite_plots)


def two_gaussian(tracer_background, tracer_max,
                 show_plots, save_plots, overwrite_plots, high_res=False):

    npoints_1d = 30 if high_res else 15
    Lx = Ly = 10.0
    x_1d = np.linspace(0, Lx, npoints_1d)
    y_1d = np.linspace(0, Ly, npoints_1d)
    coords_X, coords_Y = np.meshgrid(x_1d, y_1d, indexing='ij')

    field_data = two_gaussian_initial_condition(coords_X, coords_Y, Lx, Ly,
                                                tracer_background, tracer_max)

    colour_scheme = 'Blues'

    return PlotSetup(coords_X, coords_Y, field_data, colour_scheme,
                     tracer_background, tracer_max, show_plots,
                     save_plots, overwrite_plots)


# ---------------------------------------------------------------------------- #
# Routines for allowing plots to be shown from tests
# ---------------------------------------------------------------------------- #

def pytest_addoption(parser):
    parser.addoption("--show_plots", action="store_true", default=False)
    parser.addoption("--save_plots", action="store_true", default=False)
    parser.addoption("--overwrite_plots", action="store_true", default=False)


# ---------------------------------------------------------------------------- #
# The pytest fixture
# ---------------------------------------------------------------------------- #

@pytest.fixture(scope="session")
def plot_setup(pytestconfig):

    show_plots = pytestconfig.getoption("show_plots")
    save_plots = pytestconfig.getoption("save_plots")
    overwrite_plots = pytestconfig.getoption("overwrite_plots")
    if save_plots and overwrite_plots:
        raise ValueError('Cannot save plots and overwrite plots. '
                         + 'Supply only one of these options')

    def _plot_setup(initial_condition, tracer_background, tracer_max, high_res=False):
        if initial_condition == "dipole":
            return dipole(tracer_background, tracer_max, show_plots,
                          save_plots, overwrite_plots, high_res)
        elif initial_condition == "two_gaussian":
            return two_gaussian(tracer_background, tracer_max, show_plots,
                                save_plots, overwrite_plots, high_res)
        else:
            raise ValueError(
                f'Config initial condition {initial_condition} not recognised')

    return _plot_setup
