"""
A tomplot example, for making a quick single plot from Gusto data.
"""

import matplotlib.pyplot as plt
from netCDF4 import Dataset
from os.path import abspath, dirname
from tomplot import (set_tomplot_style, tomplot_contours, tomplot_cmap,
                     plot_contoured_field, add_colorbar_ax,
                     tomplot_field_title, extract_gusto_coords,
                     extract_gusto_field)

# ---------------------------------------------------------------------------- #
# Directory for results and plots
# ---------------------------------------------------------------------------- #
# When copying this example these should not be relative to this file
results_dir = f'{abspath(dirname(__file__))}/../tests/data'
plot_dir = f'{abspath(dirname(__file__))}/../tests/tmp_figures'
results_file_name = f'{results_dir}/gusto_sphere_2d_field_output.nc'
plot_name = f'{plot_dir}/example_gusto_quick.png'
# ---------------------------------------------------------------------------- #
# Things that should be altered based on the plot
# ---------------------------------------------------------------------------- #
field_name = 'RelativeVorticity'
colour_scheme = 'RdBu_r'
time_idx = -1
field_label = r'$\zeta \ / $ s$^{-1}$'
contour_method = 'tricontour'
# ---------------------------------------------------------------------------- #
# Things that are likely the same for all plots
# ---------------------------------------------------------------------------- #
set_tomplot_style()
data_file = Dataset(results_file_name, 'r')
# ---------------------------------------------------------------------------- #
# Data extraction
# ---------------------------------------------------------------------------- #
field_data = extract_gusto_field(data_file, field_name, time_idx=time_idx)
coords_X, coords_Y = extract_gusto_coords(data_file, field_name)
time = data_file['time'][time_idx]
# ---------------------------------------------------------------------------- #
# Plot data
# ---------------------------------------------------------------------------- #
fig, ax = plt.subplots(1, 1, figsize=(5, 5))
contours = tomplot_contours(field_data)
cmap, lines = tomplot_cmap(contours, colour_scheme, remove_contour=0.0)
cf, _ = plot_contoured_field(ax, coords_X, coords_Y, field_data, contour_method,
                             contours, cmap=cmap, line_contours=lines)
add_colorbar_ax(ax, cf, field_label)
tomplot_field_title(ax, f't = {time:.1f}', minmax=True, field_data=field_data)
# ---------------------------------------------------------------------------- #
# Save figure
# ---------------------------------------------------------------------------- #
print(f'Saving figure to {plot_name}')
fig.savefig(plot_name, bbox_inches='tight')
plt.close()


