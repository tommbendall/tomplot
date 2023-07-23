"""
A tomplot example, plotting 6 different slices with LFRic data.
"""

import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from os.path import abspath, dirname
from tomplot import (set_tomplot_style, tomplot_contours, tomplot_cmap,
                     plot_contoured_field, add_colorbar_ax,
                     regrid_vertical_slice, tomplot_field_title,
                     extract_gusto_vertical_slice, apply_gusto_domain)

# ---------------------------------------------------------------------------- #
# Directory for results and plots
# ---------------------------------------------------------------------------- #
# When copying this example these should not be relative to this file
results_dir = f'{abspath(dirname(__file__))}/../tests/data'
plot_dir = f'{abspath(dirname(__file__))}/../tests/tmp_figures'
results_file_name = f'{results_dir}/gusto_sphere_3d_field_output.nc'
plot_name = f'{plot_dir}/example_gusto_sphere_3x2.png'

# ---------------------------------------------------------------------------- #
# Things that should be altered based on the plot
# ---------------------------------------------------------------------------- #
# Specify lists for variables that are different between subplots
field_names = ['u_zonal', 'u_meridional', 'u_radial',
               'u_zonal', 'u_meridional', 'u_radial']
slice_along_values = ['lat', 'lat', 'lat',
                      'lon', 'lon', 'lon']
field_labels = [r'$u \ / $ m s$^{-1}$', r'$v \ / $ m s$^{-1}$', r'$w \ / $ m s$^{-1}$',
                r'$u \ / $ m s$^{-1}$', r'$v \ / $ m s$^{-1}$', r'$w \ / $ m s$^{-1}$']
colour_schemes = ['OrRd', 'RdBu_r', 'RdBu_r',
                  'OrRd', 'RdBu_r', 'RdBu_r']
# Things that are the same for all subplots
time_idx = -1
contour_method = 'tricontour'
slice_at = 0.0
# 1D grids for vertical regridding
coords_lon_1d = np.linspace(-180, 180, 50)
coords_lat_1d = np.linspace(-90, 90, 50)
# Dictionary to hold plotting grids -- keys are "slice_along" values
plotting_grids = {'lat': coords_lon_1d, 'lon': coords_lat_1d}
# ---------------------------------------------------------------------------- #
# Things that are likely the same for all plots
# ---------------------------------------------------------------------------- #
set_tomplot_style(fontsize=12)
data_file = Dataset(results_file_name, 'r')
fig, axarray = plt.subplots(2, 3, figsize=(12, 8), sharey='row')

# Loop through subplots
for i, (ax, field_name, field_label, colour_scheme, slice_along) in \
    enumerate(zip(axarray.flatten(), field_names, field_labels,
                  colour_schemes, slice_along_values)):
    # ------------------------------------------------------------------------ #
    # Data extraction
    # ------------------------------------------------------------------------ #
    orig_field_data, orig_coords_X, orig_coords_Y, orig_coords_Z = \
        extract_gusto_vertical_slice(data_file, field_name, time_idx,
                                     slice_along=slice_along, slice_at=slice_at)
    # Slices need regridding as points don't cleanly live along lon or lat = 0.0
    field_data, coords_X, coords_Z = regrid_vertical_slice(plotting_grids[slice_along],
                                                           slice_along, slice_at,
                                                           orig_coords_X, orig_coords_Y,
                                                           orig_coords_Z, orig_field_data)
    # ------------------------------------------------------------------------ #
    # Plot data
    # ------------------------------------------------------------------------ #
    contours = tomplot_contours(field_data)
    cmap, lines = tomplot_cmap(contours, colour_scheme)
    cf, _ = plot_contoured_field(ax, coords_X, coords_Z, field_data,
                                 contour_method, contours, cmap=cmap,
                                 line_contours=lines)
    add_colorbar_ax(ax, cf, field_label, location='bottom', cbar_labelpad=-10)
    # Don't add ylabels unless left-most subplots
    ylabel = True if i % 3 == 0 else None
    apply_gusto_domain(ax, data_file, slice_along=slice_along, ylabel=ylabel, ylabelpad=-30)
    tomplot_field_title(ax, None, minmax=True, field_data=field_data)

# These subplots tend to be quite clustered together, so move them apart a bit
fig.subplots_adjust(wspace=0.3, hspace=0.15)

# ---------------------------------------------------------------------------- #
# Save figure
# ---------------------------------------------------------------------------- #
print(f'Saving figure to {plot_name}')
fig.savefig(plot_name, bbox_inches='tight')
plt.close()


