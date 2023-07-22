"""
Plots vertical slices of high resolution NWP data
"""
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np
import pandas as pd
from tomplot import plot_contoured_field, add_colorbar, \
                    tomplot_field_title, rounded_limits, set_tomplot_style, \
                    lonlat_to_alphabeta, plot_cubed_sphere_panels, \
                    plot_cubed_sphere_slice

# ---------------------------------------------------------------------------- #
# Variables to alter based on the desired plot and test case
# ---------------------------------------------------------------------------- #

# Parameters for data extraction and figure creation that must be specified
results_dir = '/hpc/scratch_xcs/d03/frib/cylc-run/r43240_c768/work/1/run_lfric_atm_nwp_gal9_c768_C768_MG_dt-300p0_intel_64-bit_fast-debug_attempt2b'
results_filename = f'{results_dir}/lfric_gal_diagnostics.nc'
height_filename = '/data/users/tbendall/results/nwp_c768/lfric_initial.nc'
plotdir = '/data/users/tbendall/results/nwp_c768/'
plotstem = 'nwp_c768_oz_vert'
figsize = (24,10)

# Things that will be the same for each subplot
time_idx = -1
lon_point, lat_point = 126.03515625, -25.80685424
title_stem = 'Attempt 2'
alpha_or_beta = 'alpha'
alpha_filter = [np.pi/8, np.pi/4]  # List of lower/upper values for alpha, or None for whole panel
level_filter = [0, 30]  # List of lower/upper values for level, or None

# Things that differ with each subplot
field_names = ['u_in_w3', 'v_in_w3', 'w_in_wth', 'theta', 'm_v', 'm_cl']

# ---------------------------------------------------------------------------- #
# Coordinate filtering
# ---------------------------------------------------------------------------- #

eps = 1e-4
alpha, beta, panel_filter = lonlat_to_alphabeta(lon_point, lat_point)
if alpha_or_beta == 'beta':
    alpha_filter = [alpha-eps, alpha+eps]
    if beta_filter is None:
        beta_filter = [-np.pi/4, np.pi/4]
else:
    beta_filter = [beta-eps, beta+eps]  # List of lower/upper values for beta, or None
    if alpha_filter is None:
        alpha_filter = [-np.pi/4, np.pi/4]
not_alpha_or_beta = r'$\beta$' if alpha_or_beta == 'alpha' else r'$\alpha$'
alpha_or_beta_value = beta if alpha_or_beta else alpha

# ---------------------------------------------------------------------------- #
# Things that are likely the same for all scripts
# ---------------------------------------------------------------------------- #

# This is declared BEFORE figure and ax are initialised
set_tomplot_style(fontsize=16)

fig = plt.figure(figsize=figsize)

# ---------------------------------------------------------------------------- #
# Plot location of slice
# ---------------------------------------------------------------------------- #

projection = ccrs.PlateCarree(central_longitude=0)
ax = fig.add_subplot(3,3,1,projection=projection)
ax.set_position([0.14, 0.62, 0.15, 0.15])
ax.coastlines()
plot_cubed_sphere_panels(ax, color='grey', linewidth=0.5)
plot_cubed_sphere_slice(ax, alpha_filter, beta_filter, panel_filter,
                        color='red', linewidth=2.0)

# ---------------------------------------------------------------------------- #
# Field plots
# ---------------------------------------------------------------------------- #

for i, field_name in enumerate(field_names):

    ax = fig.add_subplot(3,3,4+i)

    data_file = Dataset(results_filename, 'r')

    root_coords_name = data_file[field_name].dimensions[-1]
    coords_X_name = root_coords_name[1:]+'_x'
    coords_Y_name = root_coords_name[1:]+'_y'

    if level_filter is not None:
        levels = range(level_filter[0], level_filter[1]+1)
        num_levels = len(levels)
    else:
        num_levels = np.shape(data_file[field_name])[1]
        levels = range(num_levels)

    # Work out which level this should be on
    if (data_file[field_name].dimensions[1] == 'half_levels'
        and data_file[field_name].dimensions[2] == 'nMesh2d_face'):
        height_name = 'height_w3'
    elif (data_file[field_name].dimensions[1] == 'full_levels'
          and data_file[field_name].dimensions[2] == 'nMesh2d_face'):
        height_name = 'height_wth'
    else:
        raise NotImplementedError(f'Dimensions for {field_name} are not '
            + 'implemented so cannot work out height field')

    longitude = data_file[coords_X_name][:]
    latitude = data_file[coords_Y_name][:]
    alpha_level, beta_level, panel_level = lonlat_to_alphabeta(longitude, latitude)
    data_file.close()

    # ------------------------------------------------------------------------ #
    # Loop through levels to collate data
    # ------------------------------------------------------------------------ #

    for lev_idx, level in enumerate(levels):

        data_file = Dataset(results_filename, 'r')
        field_data_level = data_file[field_name][time_idx,level,:]
        data_file.close()

        height_file = Dataset(height_filename, 'r')
        height_level = height_file[height_name][level,:]
        height_file.close()

        # Spatial filtering of data based on horizontal coordinate
        df = pd.DataFrame({'alpha': alpha_level,
                           'beta': beta_level,
                           'panel': panel_level,
                           'height': height_level,
                           'field_data': field_data_level})
        if alpha_filter is not None:
            df = df[(df['alpha'] > alpha_filter[0]) & (df['alpha'] < alpha_filter[1])]
        if beta_filter is not None:
            df = df[(df['beta'] > beta_filter[0]) & (df['beta'] < beta_filter[1])]
        if panel_filter is not None:
            df = df[df['panel'] == panel_filter]

        # Create final data arrays if we are in the bottom level
        if level == 0:
            len_data = len(df['field_data'].values)
            field_data = np.zeros((len_data, num_levels))
            coords_X = np.zeros((len_data, num_levels))
            coords_Y = np.zeros((len_data, num_levels))

        # Populate final data arrays
        coords_X[:,lev_idx] = df[alpha_or_beta].values
        coords_Y[:,lev_idx] = df['height'].values
        field_data[:,lev_idx] = df['field_data'].values

    # ------------------------------------------------------------------------ #
    # Make resulting subplot
    # ------------------------------------------------------------------------ #

    nice_minmax = rounded_limits(field_data)
    contours = np.linspace(nice_minmax[0], nice_minmax[1], 256)
    cmap = 'RdYlBu_r'

    cf, _ = plot_contoured_field(ax, coords_X, coords_Y, field_data,
                                "contour", contours, cmap=cmap,
                                plot_contour_lines=False)

    add_colorbar(ax, cf, '', cbar_ticks=nice_minmax)
    tomplot_field_title(ax, field_name.replace('_', ' '), minmax=True,
                          field_data=field_data)

    # ------------------------------------------------------------------------ #
    # Adjust limits
    # ------------------------------------------------------------------------ #

    if alpha_filter is not None and alpha_or_beta == 'alpha':
        ax.set_xlim(alpha_filter)
        ax.set_xticks(alpha_filter)
        ax.set_xticklabels(alpha_filter)
    if beta_filter is not None and alpha_or_beta == 'beta':
        ax.set_xlim(beta_filter)
        ax.set_xticks(beta_filter)
        ax.set_xticklabels(beta_filter)
    else:
        ax.set_xticklabels([])

    ax.set_ylim([np.min(coords_Y), np.max(coords_Y)])

# ---------------------------------------------------------------------------- #
# Add title and save whole figure
# ---------------------------------------------------------------------------- #

fig.suptitle(fr'{title_stem}, time idx {time_idx:02d}, {not_alpha_or_beta} = {alpha_or_beta_value:.2f} rad', y=0.7)
plotname = f'{plotdir}/{plotstem}_vert_t{time_idx:02d}.png'
print(f'Plotting to {plotname}')
fig.savefig(plotname, bbox_inches='tight')
plt.close()

