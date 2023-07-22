"""
Plots final states of two Skamarock/Klemp gravity-wave tests.
"""

import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np
import pandas as pd
from tomplot import plot_contoured_field, add_colorbar, \
                    tomplot_field_axis_labels, tomplot_cmap, \
                    tomplot_field_title, set_tomplot_style

# ---------------------------------------------------------------------------- #
# Variables to alter based on the desired plot and test case
# ---------------------------------------------------------------------------- #

# Parameters for data extraction and figure creation that must be specified
results_stem = '/data/users/tbendall/cylc-run/dispersion_fix-gungho_model-meto-spice-gw/share/output/intel_64-bit_fast-debug/gungho_model/'
mesh = 'BiP300x4-1000x2000_dt-10p0'
plotname = '/home/h01/tbendall/results/adv_vs_cons_gw.png'
time_idx = -1
figsize = (12,6)

# Things that will be the same for each subplot
colour_scheme = 'RdBu_r'
field_name = 'theta'
contours = np.linspace(-5e-3, 5e-3, 21)
cbar_ticks = [-5e-3, 5e-3]
cbar_label = r"$\theta' \ / $ K"
dz = 3
Y_filter = [500., 1500.]
ylims = [0, 30]
xlims = [-150, 150]

# Things that differ with each subplot
results_dirs = ['skamarock_klemp_gw_adv', 'skamarock_klemp_gw_cons']
titles = ['Advective', 'Consistent']

# ---------------------------------------------------------------------------- #
# Things that are likely the same for all scripts
# ---------------------------------------------------------------------------- #

# This is declared BEFORE figure and ax are initialised
set_tomplot_style(18)

# ---------------------------------------------------------------------------- #
# Field plots
# ---------------------------------------------------------------------------- #

fig, axarray = plt.subplots(1,2,figsize=figsize,sharey='row')

for i, (ax, results_dir, title) in enumerate(zip(axarray, results_dirs, titles)):

    # ------------------------------------------------------------------------ #
    # Get coordinates
    # ------------------------------------------------------------------------ #

    results_filename = f'{results_stem}/{results_dir}/{mesh}/results/lfric_diag.nc'
    initial_filename = f'{results_stem}/{results_dir}/{mesh}/results/lfric_initial.nc'
    data_file = Dataset(results_filename, 'r')

    root_coords_name = data_file[field_name].dimensions[-1]
    coords_X_name = root_coords_name[1:]+'_x'
    coords_Y_name = root_coords_name[1:]+'_y'
    num_levels = np.shape(data_file[field_name])[1]

    levels = range(num_levels)
    coords_X_level = data_file[coords_X_name][:]
    coords_Y_level = data_file[coords_Y_name][:]

    # Scale coordinate fields to km
    coords_X_level *= 0.001
    coords_X_level = np.where(coords_X_level < -75, coords_X_level + 225, coords_X_level - 75)

    data_file.close()

    # ------------------------------------------------------------------------ #
    # Loop through levels to collate data
    # ------------------------------------------------------------------------ #

    for lev_idx, level in enumerate(levels):

        data_file = Dataset(results_filename, 'r')
        field_data_level = data_file[field_name][time_idx,level,:]
        data_file.close()

        initial_data_file = Dataset(initial_filename, 'r')
        initial_field_data_level = initial_data_file[field_name][level,:]
        initial_data_file.close()

        height_level = np.ones_like(field_data_level)*dz*lev_idx

        # Spatial filtering of data based on horizontal coordinate
        df = pd.DataFrame({'X': coords_X_level,
                           'Y': coords_Y_level,
                           'height': height_level,
                           'field_data': field_data_level,
                           'initial_field_data': initial_field_data_level})
        
        df = df[(df['Y'] > Y_filter[0]) & (df['Y'] < Y_filter[1])]
        df = df.sort_values('X')

        # Create final data arrays if we are in the bottom level
        if level == 0:
            len_data = len(df['field_data'].values)
            field_data = np.zeros((len_data, num_levels))
            coords_X = np.zeros((len_data, num_levels))
            coords_Y = np.zeros((len_data, num_levels))

        # Populate final data arrays
        coords_X[:,lev_idx] = df['X'].values
        coords_Y[:,lev_idx] = df['height'].values
        # Use original column
        field_data[:,lev_idx] = df['field_data'].values - df['initial_field_data'].values[0]

    # ------------------------------------------------------------------------ #
    # Make resulting subplot
    # ------------------------------------------------------------------------ #

    remove_contour = 0.0
    cmap, line_contours = tomplot_cmap(contours, color_scheme=colour_scheme,
                                       remove_contour=remove_contour)

    cf, cl = plot_contoured_field(ax, coords_X, coords_Y, field_data,
                                  'contour', contours, line_contours, cmap=cmap)

    tomplot_field_title(ax, title, minmax=False)

    ax.set_xlim(xlims)
    ax.set_xticks(xlims)
    ax.set_xticklabels(xlims)
    ax.set_ylim(ylims)
    ax.set_yticks(ylims)
    ax.set_yticklabels(ylims)

# Move the subplots to the left to make space for colorbar
fig.subplots_adjust(right=0.9, wspace=0.1)
# Add colorbar in its own axis
cbar_ax = fig.add_axes([0.92, 0.11, 0.02, 0.77])
cb = fig.colorbar(cf, cax=cbar_ax, ticks=cbar_ticks)
cb.set_label(cbar_label, labelpad=-40)

print(f'Plotting to {plotname}')
fig.savefig(plotname, bbox_inches='tight')
plt.close()

