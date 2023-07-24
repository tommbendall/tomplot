"""
This tests the plotting of 2D fields through coloured contours.
"""

import matplotlib.pyplot as plt
from tomplot import (add_colorbar_ax, tomplot_field_title,
                     tomplot_cmap, plot_contoured_field)
import numpy as np
import pytest


@pytest.mark.parametrize("method", ["contour", "tricontour", "scatter"])
def test_field_contour_plot(method, plot_setup):

    # ------------------------------------------------------------------------ #
    # Settings for test
    # ------------------------------------------------------------------------ #

    data_background = 20.0
    data_max = 40.0
    data_diff = data_max - data_background
    setup = plot_setup('dipole', data_background, data_max, npoints_1d=30)
    contours = np.linspace(-data_diff, data_diff, 11)
    title = f'2D field contour plot with {method} method'

    coords_X, coords_Y = setup.coords_X, setup.coords_Y
    field_data = setup.field_data

    cmap, _ = tomplot_cmap(contours, setup.colour_scheme)

    # ------------------------------------------------------------------------ #
    # Make plot
    # ------------------------------------------------------------------------ #

    plt.close()

    _, ax = plt.subplots(1, 1, figsize=(5, 5))

    cf, _ = plot_contoured_field(ax, coords_X, coords_Y, field_data,
                                 method, contours, cmap=cmap,
                                 plot_contour_lines=False)

    add_colorbar_ax(ax, cf, '')
    tomplot_field_title(ax, title)

    plot_name = f'field_contour_plot_{method}.png'
    setup.make_plots(plot_name)
