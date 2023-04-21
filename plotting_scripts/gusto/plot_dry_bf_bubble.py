"""
Plots final states of two Bryan and Fritsch bubble tests.
"""

import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np
from tomplot import plot_contoured_field, extract_2D_data, add_colorbar, \
                    automake_field_axis_labels, automake_cmap, \
                    automake_field_title, label_contour_lines

# ---------------------------------------------------------------------------- #
# Variables to alter based on the desired plot and test case
# ---------------------------------------------------------------------------- #

# Parameters for data extraction and figure creation that must be specified
results_filename = '/home/thomas/firedrake/src/gusto/results/dry_bryan_fritsch/field_output.nc'
plotname = '/home/thomas/firedrake/src/gusto/tomplot/test_plots/test_dry_bf.png'
time_idx = -1
figsize = (12,6)

# Things that will be the same for each subplot

# Things that differ with each subplot
field_names = ['rho', 'theta_perturbation']
all_contours = [np.arange(0.6, 1.3, step=0.1), np.linspace(-2.0, 2.0, 11)]
cbar_labels = [r'$\rho \ / $ kg m$^3$', r"$\theta' \ / $ K"]


# ---------------------------------------------------------------------------- #
# Things that are likely the same for all scripts
# ---------------------------------------------------------------------------- #

# This is declared BEFORE figure and ax are initialised
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
font = {'size':24}
plt.rc('font',**font)

# ---------------------------------------------------------------------------- #
# Field plots
# ---------------------------------------------------------------------------- #

fig, axarray = plt.subplots(1,2,figsize=figsize,sharey='row')

for ax, field_name, contours, cbar_label in \
    zip(axarray, field_names, all_contours, cbar_labels):

    data_file = Dataset(results_filename, 'r')

    coords_X, coords_Y, field_data, data_metadata = \
        extract_2D_data(data_file, field_name, time_idx, 'xz',
                        num_points=50)

    data_file.close()

    # Scale coordinate fields to km
    coords_X *= 0.001
    coords_Y *= 0.001

    remove_contour = 0.0 if field_name == 'theta_perturbation' else None
    cmap, contours = automake_cmap(contours, remove_contour=remove_contour)

    cf, cl = plot_contoured_field(ax, coords_X, coords_Y, field_data,
                                 contours, method='contour', cmap=cmap)

    add_colorbar(ax, cf, cbar_label)
    if field_name == 'rho':
        label_contour_lines(ax, cl)
    automake_field_axis_labels(ax, data_metadata)
    automake_field_title(ax, field_name.replace('_', ' '), minmax=True,
                         field_data=field_data)


print(f'Plotting to {plotname}')
fig.savefig(plotname, bbox_inches='tight')
plt.close()

