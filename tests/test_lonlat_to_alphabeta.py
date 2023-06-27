"""
This tests the conversion to cubed sphere coordinates.
"""

import matplotlib.pyplot as plt
import matplotlib.colors as mat_col
from tomplot import (add_colorbar, tomplot_field_title, plot_contoured_field,
                     lonlat_to_alphabeta)
import numpy as np
import pytest

panel_colours = [mat_col.to_rgba_array('red'),
                 mat_col.to_rgba_array('orange'),
                 mat_col.to_rgba_array('gold'),
                 mat_col.to_rgba_array('forestgreen'),
                 mat_col.to_rgba_array('skyblue'),
                 mat_col.to_rgba_array('mediumpurple'),]


@pytest.mark.parametrize("aspect", ['panel_layout', 'panel_orientation'])
def test_field_markersize_subplots(aspect, plot_setup):

    plt.close()

    # ------------------------------------------------------------------------ #
    # Generate coordinate data
    # ------------------------------------------------------------------------ #

    coords_lon_1d = np.linspace(-180, 180, 51)
    coords_lat_1d = np.linspace(-90, 90, 51)
    coords_lon, coords_lat = np.meshgrid(coords_lon_1d, coords_lat_1d, indexing='ij')

    alpha, beta, panel = lonlat_to_alphabeta(coords_lon, coords_lat)

    if aspect == 'panel_orientation':
        field_data = alpha + 5*beta
        contours = [-6*np.pi/4, 6*np.pi/4]
        cmap = 'RdBu_r'
    elif aspect == 'panel_layout':
        field_data = panel
        contours = [0.5, 6.5]
        cmap = mat_col.ListedColormap(panel_colours)

    # ------------------------------------------------------------------------ #
    # Make plot
    # ------------------------------------------------------------------------ #

    _, ax = plt.subplots(1, 1, figsize=(5, 5))

    cf, _ = plot_contoured_field(ax, coords_lon, coords_lat, field_data,
                                 "scatter", contours, cmap=cmap,
                                 plot_contour_lines=False)

    add_colorbar(ax, cf, '')
    title = f'lonlat_to_alphabeta: {aspect}'
    tomplot_field_title(ax, title)

    plot_name = f'lonlat_to_alphabeta_{aspect}.png'
    setup = plot_setup('none', 0.0, 0.0)
    setup.make_plots(plot_name)
