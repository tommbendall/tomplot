"""
This tests the plotting of 2D fields through coloured contours.
"""

import matplotlib.pyplot as plt
from tomplot import add_colorbar, automake_field_title, plot_contoured_field
import numpy as np
import pytest

@pytest.mark.parametrize("npoints_1d", [20, 40, 60, 100])
@pytest.mark.parametrize("figsize", [(5,5), (10,5), (5,10), (10,10)])
def test_field_markersize(npoints_1d, figsize, plot_setup):


    # ------------------------------------------------------------------------ #
    # Settings for test
    # ------------------------------------------------------------------------ #

    data_background = -10.0
    data_max = -5
    setup = plot_setup('two_gaussian', data_background, data_max, npoints_1d=npoints_1d)
    contours = np.linspace(data_background, data_max, 21)
    title = f'scatter field plot with {npoints_1d} points and figsize {figsize}'

    coords_X, coords_Y = setup.coords_X, setup.coords_Y
    field_data = setup.field_data

    # ------------------------------------------------------------------------ #
    # Make plot
    # ------------------------------------------------------------------------ #

    plt.close()

    _, ax = plt.subplots(1, 1, figsize=figsize)

    cf, _ = plot_contoured_field(ax, coords_X, coords_Y, field_data,
                                 "scatter", contours, plot_contour_lines=False)

    add_colorbar(ax, cf, '')
    automake_field_title(ax, title)

    plot_name = f'field_contour_plot_{npoints_1d}_{figsize[0]}_{figsize[1]}.png'
    setup.make_plots(plot_name)

    assert False, 'Need to find a way of testing if point size is reasonable'
