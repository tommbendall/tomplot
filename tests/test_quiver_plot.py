"""
This tests the plotting of 2D vector-valued fields using quivers.
"""

# TODO: have unstructured test using actual unstructured data!!
# TODO: test slicing / offset options!

import matplotlib.pyplot as plt
from tomplot import automake_field_title, plot_field_quivers
import numpy as np
import pytest

@pytest.mark.parametrize("data_layout", ["structured", "unstructured"])
def test_quiver_plot(data_layout, plot_setup):

    # ------------------------------------------------------------------------ #
    # Settings for test
    # ------------------------------------------------------------------------ #

    data_background = 1.0
    data_max = 2.0
    dipole_setup = plot_setup('dipole', data_background, data_max, npoints_1d=30)
    gaussian_setup = plot_setup('two_gaussian', data_background, data_max, npoints_1d=30)

    # ------------------------------------------------------------------------ #
    # Get coordinates and data
    # ------------------------------------------------------------------------ #

    if data_layout == "structured":
        # Coords should be the same from either setup
        coords_X, coords_Y = dipole_setup.coords_X, dipole_setup.coords_Y
        field_data_X = dipole_setup.field_data
        field_data_Y = gaussian_setup.field_data

    elif data_layout == "unstructured":
        # Coords should be the same from either setup
        coords_X = dipole_setup.coords_X.flatten()
        coords_Y = dipole_setup.coords_Y.flatten()
        field_data_X = dipole_setup.field_data.flatten()
        field_data_Y = gaussian_setup.field_data.flatten()

    title = f'Quiver plot with {data_layout} data'

    # ------------------------------------------------------------------------ #
    # Make plot
    # ------------------------------------------------------------------------ #

    plt.close()

    _, ax = plt.subplots(1, 1, figsize=(5, 5))

    _ = plot_field_quivers(ax, coords_X, coords_Y, field_data_X, field_data_Y)

    automake_field_title(ax, title)

    plot_name = f'quiver_contour_plot_{data_layout}.png'
    dipole_setup.make_plots(plot_name)
