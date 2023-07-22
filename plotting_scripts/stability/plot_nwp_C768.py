"""
Plots lon-lat slices of high resolution NWP data
"""
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np
import pandas as pd
from tomplot import plot_contoured_field, add_colorbar, plot_cubed_sphere_panels, \
                    tomplot_field_title, rounded_limits, set_tomplot_style

# ---------------------------------------------------------------------------- #
# Variables to alter based on the desired plot and test case
# ---------------------------------------------------------------------------- #

# Parameters for data extraction and figure creation that must be specified
results_dir = '/hpc/scratch_xcs/d03/frib/cylc-run/r43240_c768/work/1/run_lfric_atm_nwp_gal9_c768_C768_MG_dt-300p0_intel_64-bit_fast-debug_attempt2b'
results_filename = f'{results_dir}/lfric_gal_diagnostics.nc'
plotdir = '/data/users/tbendall/results/nwp_c768/'
plotstem = 'nwp_c768_oz'
figsize = (24,10)

# Things that will be the same for each subplot
time_idx = -1
levels = [10,20,30,40]
coord_filter_x = [100, 160]
coord_filter_y = [-60, 0]
title_stem = 'Attempt 2'

# Things that differ with each subplot
field_names = ['u_in_w3', 'v_in_w3', 'w_in_wth', 'theta', 'm_v', 'm_cl']

# ---------------------------------------------------------------------------- #
# Things that are likely the same for all scripts
# ---------------------------------------------------------------------------- #

# This is declared BEFORE figure and ax are initialised
set_tomplot_style(fontsize=16)

# ---------------------------------------------------------------------------- #
# Loop through levels
# ---------------------------------------------------------------------------- #

for level in levels:

    # ------------------------------------------------------------------------ #
    # Field plots
    # ------------------------------------------------------------------------ #

    fig = plt.figure(figsize=figsize)

    for i, field_name in enumerate(field_names):

        # To get cubed sphere panels, need to use Cartopy and create ax with PlateCarree
        ax = fig.add_subplot(2,3,1+i,projection=ccrs.PlateCarree(central_longitude=0))

        data_file = Dataset(results_filename, 'r')

        root_coords_name = data_file[field_name].dimensions[-1]
        coords_X_name = root_coords_name[1:]+'_x'
        coords_Y_name = root_coords_name[1:]+'_y'

        coords_X = data_file[coords_X_name][:]
        coords_Y = data_file[coords_Y_name][:]
        field_data = data_file[field_name][time_idx,level,:]

        data_file.close()


        nice_minmax = rounded_limits(field_data)
        if field_name == 'w_in_wth':
            nice_minmax = [-1, 1]
        cmap = 'RdYlBu_r'

        cf, _ = plot_contoured_field(ax, coords_X, coords_Y, field_data,
                                    "scatter", nice_minmax, cmap=cmap,
                                    plot_contour_lines=False,
                                    marker_scaling=0.01)

        add_colorbar(ax, cf, '')
        tomplot_field_title(ax, field_name.replace('_', ' '), minmax=True,
                             field_data=field_data)

        plot_cubed_sphere_panels(ax, linewidth=0.5)

        if coord_filter_x is not None:
            ax.set_xlim(coord_filter_x)
            ax.set_xticks(coord_filter_x)
            ax.set_xticklabels(coord_filter_x)
        else:
            ax.set_xticklabels()

        if coord_filter_y is not None:
            ax.set_ylim(coord_filter_y)
            ax.set_yticks(coord_filter_y)
            ax.set_yticklabels(coord_filter_y)
        else:
            ax.set_yticklabels()

    fig.suptitle(f'{title_stem} Time idx {time_idx:02d}, Level {level:02d}')
    plotname = f'{plotdir}/{plotstem}_t{time_idx:02d}_l{level:02d}.png'
    print(f'Plotting to {plotname}')
    fig.savefig(plotname, bbox_inches='tight')
    plt.close()

