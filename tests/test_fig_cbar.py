"""
This tests the plotting of 2D fields through coloured contours.
"""

import matplotlib.pyplot as plt
from tomplot import add_colorbar_fig, plot_contoured_field
import numpy as np
import pytest

# Generate all combos of different options
def generate_combos(arrangements, locations):
    import itertools
    return list(itertools.product(arrangements, locations))

@pytest.mark.parametrize("subplot_arrangement, cbar_location",
                         # All combinations of square arrangements with location
                         generate_combos([(1, 1), (2, 2)],
                                         ["top", "bottom", "left", "right"])
                         # For non-square arrangements be more specific
                         + generate_combos([(1, 2)], ["left", "right"])
                         + generate_combos([(2, 1)], ["top", "bottom"]))
def test_fig_cbar(subplot_arrangement, cbar_location, plot_setup):

    plt.close()

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

    # ------------------------------------------------------------------------ #
    # Make figure
    # ------------------------------------------------------------------------ #

    subplots_x, subplots_y = subplot_arrangement
    figsize = (5*subplots_y, 5*subplots_x)

    fig, axarray = plt.subplots(subplots_x, subplots_y, figsize=figsize)

    if (subplots_x == 1 and subplots_y == 1):
        # Make axarray into an array, as it will just be a single Axes
        axarray = [axarray]
    else:
        # Flatten in case this is a 2D array
        axarray = axarray.flatten()

    title = f'fig_cbar {cbar_location} with axarray {subplot_arrangement}'

    for ax in axarray:

        # ------------------------------------------------------------------------ #
        # Make plot
        # ------------------------------------------------------------------------ #

        cf, _ = plot_contoured_field(ax, coords_X, coords_Y, field_data,
                                     "contour", contours)

    fig.suptitle(title)
    add_colorbar_fig(fig, cf, location=cbar_location)

    plot_name = f'fig_cbar_{cbar_location}_{subplots_x}_{subplots_y}.png'
    setup.make_plots(plot_name)
