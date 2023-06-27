"""
This tests the plotting of 2D vector-valued fields using quivers.
"""

# TODO: have unstructured test using actual unstructured data!!
# TODO: test slicing / offset options!

import matplotlib.pyplot as plt
from tomplot import tomplot_field_title, plot_field_quivers
import pytest


@pytest.mark.parametrize("data_layout", ["structured", "unstructured"])
@pytest.mark.parametrize("spatial_filter", ["with", "without"])
@pytest.mark.parametrize("magnitude_filter", [None, 1.5])
def test_quiver_plot(data_layout, spatial_filter, magnitude_filter, plot_setup):

    # ------------------------------------------------------------------------ #
    # Settings for test
    # ------------------------------------------------------------------------ #

    data_background = 1.0
    data_max = 2.0
    dipole_setup = plot_setup('dipole', data_background, data_max, npoints_1d=30)
    gaussian_setup = plot_setup('two_gaussian', data_background, data_max, npoints_1d=30)

    # ------------------------------------------------------------------------ #
    # Get coordinates and data
    # ------------------------------------------------------------------------ #

    if data_layout == "structured":
        # Coords should be the same from either setup
        coords_X, coords_Y = dipole_setup.coords_X, dipole_setup.coords_Y
        field_data_X = dipole_setup.field_data
        field_data_Y = gaussian_setup.field_data

        if spatial_filter == "with":
            spatial_filter_step = (5, 5)
            spatial_filter_offset = (1, 1)
        else:
            spatial_filter_step = None
            spatial_filter_offset = None

    elif data_layout == "unstructured":
        # Coords should be the same from either setup
        coords_X = dipole_setup.coords_X.flatten()
        coords_Y = dipole_setup.coords_Y.flatten()
        field_data_X = dipole_setup.field_data.flatten()
        field_data_Y = gaussian_setup.field_data.flatten()

        if spatial_filter == "with":
            spatial_filter_step = 5
            spatial_filter_offset = 1
        else:
            spatial_filter_step = None
            spatial_filter_offset = None

    title = f'Quivers {data_layout} data, {spatial_filter} spatial filter'
    if magnitude_filter is not None:
        title += ' and mag filter'

    # ------------------------------------------------------------------------ #
    # Make plot
    # ------------------------------------------------------------------------ #

    plt.close()

    _, ax = plt.subplots(1, 1, figsize=(5, 5))

    _ = plot_field_quivers(ax, coords_X, coords_Y, field_data_X, field_data_Y,
                           magnitude_filter=magnitude_filter,
                           spatial_filter_step=spatial_filter_step,
                           spatial_filter_offset=spatial_filter_offset)

    ax.set_xlim([-0.5, 10.5])
    ax.set_ylim([-0.5, 10.5])

    tomplot_field_title(ax, title)

    spat_fil_txt = 'with_spat_fil' if spatial_filter == 'with' else 'no_spat_fil'
    mag_fil_txt = 'with_mag_fil' if magnitude_filter is not None else 'no_mag_fil'
    plot_name = f'quiver_contour_plot_{data_layout}_{spat_fil_txt}_{mag_fil_txt}.png'
    dipole_setup.make_plots(plot_name)
