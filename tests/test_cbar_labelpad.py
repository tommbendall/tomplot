"""
This tests the automatic padding of cbar labels.
"""

import matplotlib
import matplotlib.pyplot as plt
from tomplot import (add_colorbar_ax, plot_contoured_field, set_tomplot_style,
                     tomplot_field_title, tomplot_cmap)
import numpy as np
import pytest


@pytest.mark.parametrize("cbar_format", ['.1f', '.3f', '.4e'])
@pytest.mark.parametrize("fontsize", [12, 24, 32])
def test_cbar_labelpad(cbar_format, fontsize, plot_setup):

    plt.close()

    set_tomplot_style(fontsize)

    title = f'cbar labelpad with fontsize {fontsize} and format {cbar_format}'

    # ------------------------------------------------------------------------ #
    # Settings for data
    # ------------------------------------------------------------------------ #

    data_background = 20.0
    data_max = 40.0
    data_diff = data_max - data_background
    setup = plot_setup('dipole', data_background, data_max, npoints_1d=30)
    contours = np.linspace(-data_diff, data_diff, 11)

    coords_X, coords_Y = setup.coords_X, setup.coords_Y
    field_data = setup.field_data

    cmap, _ = tomplot_cmap(contours, setup.colour_scheme)

    # ------------------------------------------------------------------------ #
    # Make figure
    # ------------------------------------------------------------------------ #

    _, ax = plt.subplots(1, 1, figsize=(5, 5))

    cf, _ = plot_contoured_field(ax, coords_X, coords_Y, field_data,
                                 'tricontour', contours, cmap=cmap,
                                 plot_contour_lines=False)

    add_colorbar_ax(ax, cf, 'dipole field', cbar_format=cbar_format)
    tomplot_field_title(ax, title)

    plot_name = f'cbar_labelpad_{cbar_format}_{fontsize}.png'
    setup.make_plots(plot_name)

    # Restore the default style
    matplotlib.style.use('default')
