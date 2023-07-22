"""
This tests the extension of a cmap to values beyond the contour range.
"""

import matplotlib.pyplot as plt
from tomplot import (add_colorbar_fig, tomplot_field_title,
                     tomplot_cmap, plot_contoured_field)
import numpy as np


def test_extend_cmap(plot_setup):

    # ------------------------------------------------------------------------ #
    # Settings for test
    # ------------------------------------------------------------------------ #

    data_background = 20.0
    data_max = 40.0
    data_diff = data_max - data_background
    setup = plot_setup('dipole', data_background, data_max, npoints_1d=30)
    contours = np.linspace(-data_diff + 4.0, data_diff - 4.0, 11)
    title = 'cmap extension'

    coords_X, coords_Y = setup.coords_X, setup.coords_Y
    field_data = setup.field_data

    cmap, _ = tomplot_cmap(contours, setup.colour_scheme, extend_cmap=True)

    # ------------------------------------------------------------------------ #
    # Make plot
    # ------------------------------------------------------------------------ #

    plt.close()

    fig, ax = plt.subplots(1, 1, figsize=(5, 5))

    cf, _ = plot_contoured_field(ax, coords_X, coords_Y, field_data,
                                 'contour', contours, cmap=cmap,
                                 plot_contour_lines=False)

    add_colorbar_fig(fig, cf)
    tomplot_field_title(ax, title)

    plot_name = 'extend_cmap.png'
    setup.make_plots(plot_name)
