"""
A tomplot example, for making side-by-side plots from Gusto vertical slice data.
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
results_file_name = f'{results_dir}/gusto_slice_2d_field_output.nc'
plot_name = f'{plot_dir}/example_gusto_slice.png'
# ---------------------------------------------------------------------------- #
# Things that should be altered based on the plot
# ---------------------------------------------------------------------------- #
# Specify lists for variables that are different between subplots
field_names = ['rho', 'theta_perturbation']
colour_schemes = ['Blues', 'OrRd']
field_labels = [r'$\rho \ / $ kg m$^{-3}$', r"$\theta' \ / $ K"]
# Things that are the same for both subplots
time_idx = -1
contour_method = 'tricontour'
# ---------------------------------------------------------------------------- #
# Things that are likely the same for all plots
# ---------------------------------------------------------------------------- #
set_tomplot_style()
data_file = Dataset(results_file_name, 'r')
fig, axarray = plt.subplots(1, 2, figsize=(10, 5), sharey='row')

# Loop through subplots
for ax, field_name, field_label, colour_scheme in \
    zip(axarray, field_names, field_labels, colour_schemes):
    # ------------------------------------------------------------------------ #
    # Data extraction
    # ------------------------------------------------------------------------ #
    field_data = extract_gusto_field(data_file, field_name, time_idx=time_idx)
    coords_X, coords_Z = extract_gusto_coords(data_file, field_name)
    time = data_file['time'][time_idx]
    # ------------------------------------------------------------------------ #
    # Plot data
    # ------------------------------------------------------------------------ #
    contours = tomplot_contours(field_data)
    cmap, lines = tomplot_cmap(contours, colour_scheme)
    cf, _ = plot_contoured_field(ax, coords_X, coords_Z, field_data, contour_method,
                                contours, cmap=cmap, line_contours=lines)
    add_colorbar_ax(ax, cf, field_label)
    tomplot_field_title(ax, f't = {time:.1f}', minmax=True, field_data=field_data)
# ---------------------------------------------------------------------------- #
# Save figure
# ---------------------------------------------------------------------------- #
print(f'Saving figure to {plot_name}')
fig.savefig(plot_name, bbox_inches='tight')
plt.close()


