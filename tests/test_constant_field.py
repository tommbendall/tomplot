"""
This tests plotting filled contours for a constant field.
"""

import matplotlib.pyplot as plt
from tomplot import (add_colorbar_ax, tomplot_field_title,
                     tomplot_cmap, plot_contoured_field)
import numpy as np
import pytest


@pytest.mark.parametrize("method", ["contour", "tricontour"])
def test_field_contour_plot(method, plot_setup):

    # ------------------------------------------------------------------------ #
    # Settings for test
    # ------------------------------------------------------------------------ #

    data_background = 21.9
    data_max = 21.9
    setup = plot_setup('two_gaussian', data_background, data_max, npoints_1d=30)
    contours = np.linspace(21.89, 21.91, 4)
    title = f'Constant field with {method} method'

    coords_X, coords_Y = setup.coords_X, setup.coords_Y
    field_data = setup.field_data

    cmap, _ = tomplot_cmap(contours, 'inferno')

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

    plot_name = f'constant_field_plot_{method}.png'
    setup.make_plots(plot_name)
