"""
This tests the plotting of scatter points using an orthographic projection.
"""

import matplotlib.pyplot as plt
from tomplot import (add_colorbar_ax, tomplot_field_title,
                     tomplot_cmap, plot_contoured_field,
                     extract_gusto_coords, extract_gusto_field,
                     tomplot_contours)
from os.path import abspath, dirname
import cartopy.crs as ccrs
from netCDF4 import Dataset

def test_orthographic_scatter_plot(plot_setup):

    # ------------------------------------------------------------------------ #
    # Settings for test
    # ------------------------------------------------------------------------ #

    # Use example data
    results_dir = f'{abspath(dirname(__file__))}/data'
    results_file_name = f'{results_dir}/gusto_sphere_2d_field_output.nc'

    field_name = 'D'
    title = f'Orthographic scatter'
    colour_scheme = 'RdYlBu_r'
    field_label = r'$D \ / $ m'
    method = 'scatter'
    spherical_centre = (90, 60)
    time_idx = 0

    # ------------------------------------------------------------------------ #
    # Extract data
    # ------------------------------------------------------------------------ #

    data_file = Dataset(results_file_name, 'r')

    field_data = extract_gusto_field(data_file, field_name, time_idx=time_idx)
    coords_lon, coords_lat = extract_gusto_coords(data_file, field_name)

    contours = tomplot_contours(field_data)
    cmap, _ = tomplot_cmap(contours, colour_scheme)

    # ------------------------------------------------------------------------ #
    # Make plot
    # ------------------------------------------------------------------------ #

    plt.close()

    # When using cartopy we need to specify projections for each axes separately
    fig = plt.figure(figsize=(5, 5))
    projection = ccrs.Orthographic(spherical_centre[0], spherical_centre[1])
    ax = fig.add_subplot(1, 1, 1, projection=projection)

    cf, _ = plot_contoured_field(ax, coords_lon, coords_lat, field_data,
                                 method, contours, cmap=cmap,
                                 plot_contour_lines=False,
                                 projection=projection)

    add_colorbar_ax(ax, cf, field_label)
    tomplot_field_title(ax, title)

    plot_name = f'orthographic_scatter.png'
    setup = plot_setup('none', 0.0, 0.0)
    setup.make_plots(plot_name)
