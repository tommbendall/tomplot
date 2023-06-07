"""
This tests the plotting of 2D fields through coloured contours.
"""

import matplotlib.pyplot as plt
from tomplot import add_colorbar, automake_field_title, plot_contoured_field
import numpy as np
import pytest


# Common settings for two tests
def field_markersize_test_settings(npoints_1d, setup):

    data_background = -10.0
    data_max = -5
    test_setup = setup('two_gaussian', data_background, data_max, npoints_1d=npoints_1d)
    contours = np.linspace(data_background, data_max, 21)

    return test_setup, contours

@pytest.mark.parametrize("figsize", [(5, 10), (10, 5)])
def test_field_markersize_subplots(figsize, plot_setup):

    answers = {510: [16, 8], 1005: [16, 8]}

    plt.close()

    subplots_x = int(figsize[1] / 5)
    subplots_y = int(figsize[0] / 5)
    fig, axarray = plt.subplots(subplots_x, subplots_y, figsize=figsize)

    npoints_1d_array = [50*(i+1) for i, _ in enumerate(axarray)]

    title = f'scatter field subplots with figsize {figsize}'

    markersizes = []

    for ax, npoints_1d in zip(axarray, npoints_1d_array):

        # -------------------------------------------------------------------- #
        # Settings for test
        # -------------------------------------------------------------------- #

        setup, contours = field_markersize_test_settings(npoints_1d, plot_setup)

        coords_X, coords_Y = setup.coords_X, setup.coords_Y
        field_data = setup.field_data

        # ------------------------------------------------------------------------ #
        # Make plot
        # ------------------------------------------------------------------------ #

        cf, _ = plot_contoured_field(ax, coords_X, coords_Y, field_data,
                                    "scatter", contours, plot_contour_lines=False)
        markersizes.append(cf.get_sizes()[0])
        add_colorbar(ax, cf, '')

    fig.suptitle(title)

    plot_name = f'field_contour_subplots_{figsize[0]}_{figsize[1]}.png'
    setup.make_plots(plot_name)

    tol = 1
    for i, markersize in enumerate(markersizes):
        answer = answers[100*figsize[0]+figsize[1]][i]
        assert abs(markersize - answer) / answer < tol, f'Markersize {i} for ' \
            + f'{npoints_1d} points in 1D for figsize {figsize} is not correct'


@pytest.mark.parametrize("npoints_1d", [20, 50, 100])
@pytest.mark.parametrize("figsize", [(5, 5), (10, 10)])
def test_field_markersize_one_plot(npoints_1d, figsize, plot_setup):

    answers = {505: {20: 100, 50: 24, 100: 8},
               1010: {20: 310, 50: 72, 100: 24}}

    # ------------------------------------------------------------------------ #
    # Settings for test
    # ------------------------------------------------------------------------ #

    setup, contours = field_markersize_test_settings(npoints_1d, plot_setup)

    coords_X, coords_Y = setup.coords_X, setup.coords_Y
    field_data = setup.field_data

    title = f'scatter field plot with {npoints_1d} points and figsize {figsize}'

    # ------------------------------------------------------------------------ #
    # Make plot
    # ------------------------------------------------------------------------ #

    plt.close()

    _, ax = plt.subplots(1, 1, figsize=figsize)

    cf, _ = plot_contoured_field(ax, coords_X, coords_Y, field_data,
                                 "scatter", contours, plot_contour_lines=False)

    markersize = cf.get_sizes()[0]
    answer = answers[100*figsize[0]+figsize[1]][npoints_1d]

    add_colorbar(ax, cf, '')
    automake_field_title(ax, title)

    plot_name = f'field_contour_one_plot_{npoints_1d}_{figsize[0]}_{figsize[1]}.png'
    setup.make_plots(plot_name)

    tol = 0.1
    assert abs(markersize - answer) / answer < tol, f'Markersize for ' \
        + f'{npoints_1d} points in 1D for figsize {figsize} is not correct'
