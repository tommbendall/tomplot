"""
In this file we test the function for making a 2D field plot with a removed contour
"""

from tomplot import *
import numpy as np
import pytest


@pytest.mark.parametrize("num_contours", ["odd", "even"])
def test_remove_contour(num_contours):

    if num_contours == "odd":
        remove_contour = 0.0
        restricted_cmaps = [None, 'both']
        colour_levels_scalings = [None, (1.3,1.3)]
        colour_scheme = 'RdBu_r'
    else:
        remove_contour = 0.0
        restricted_cmaps = [None, 'top']
        colour_levels_scalings = [None, 1.2]
        colour_scheme = 'OrRd'


    x = np.linspace(0, 10, 11)
    y = np.linspace(-5,5,15) if num_contours == 'odd' else np.linspace(-1,10,12)

    coords_X, coords_Y = np.meshgrid(x, y, indexing='ij')
    field = coords_Y

    for restricted_cmap, colour_levels_scaling in zip(restricted_cmaps, colour_levels_scalings):

        plt.close()

        fig, ax = plt.subplots(1,1,figsize=(8,8))

        cf = individual_field_contour_plot(coords_X, coords_Y, field,
                                          ax=ax, colour_scheme=colour_scheme,
                                          contours=y, no_cbar=True, title=None,
                                          title_method=None,
                                          remove_contour=remove_contour,
                                          colour_levels_scaling=colour_levels_scaling,
                                          restricted_cmap=restricted_cmap)

        # Add colorbar in its own axis
        cbar_ax = fig.add_axes([0.92, 0.11, 0.02, 0.77])
        cb = fig.colorbar(cf, cax=cbar_ax)

        plt.show()
