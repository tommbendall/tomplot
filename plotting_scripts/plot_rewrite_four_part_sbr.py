"""
Plots the 4-part SBR
"""

import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np
import pandas as pd
import cartopy.crs as ccrs
from tomplot import plot_contoured_field, plot_field_quivers, \
                    automake_field_title, automake_cmap, add_colorbar, \
                    automake_field_axis_labels, plot_cubed_sphere_panels

# ---------------------------------------------------------------------------- #
# Things that can be altered and parameters for the test case
# ---------------------------------------------------------------------------- #

results_dirname = 'four_part_sbr'
field_info = ('w2_vector1', 'w2_vector2')
cbar_label = r'$|F| \ /$ m s$^{-1}$'
colour_scheme = 'OrRd'
time_factor = 10.0
level = 0
contours = np.linspace(0,4,17)
cbar_ticks = [0, 4]
vector_magnitude_cutoff = 0.2

# ---------------------------------------------------------------------------- #
# Things that are likely the same
# ---------------------------------------------------------------------------- #

plotdir = '/data/users/tbendall/results/'+results_dirname+'/figures'
filename = '/data/users/tbendall/results/'+results_dirname+'/rewrite_lfric_diag.nc'
data_file = Dataset(filename, 'r')
time_idxs = range(len(data_file['time_instant'][:]))
data_file.close()

# This is declared BEFORE figure and ax are initialised
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
font = {'size':48}
plt.rc('font',**font)

# ---------------------------------------------------------------------------- #
# Loop through time points
# ---------------------------------------------------------------------------- #

for time_idx in time_idxs:

    fig = plt.figure(figsize=(30, 20))
    # To get cubed sphere panels, need to use Cartopy and create ax with PlateCarree
    ax = fig.add_subplot(1,1,1,projection=ccrs.PlateCarree(central_longitude=0))

    # ------------------------------------------------------------------------ #
    # Data extraction
    # ------------------------------------------------------------------------ #

    # This code is all adapted from plot_control
    data_file = Dataset(filename, 'r')

    # Extract data
    time = time_factor*data_file['time_instant'][time_idx]
    field_data_X = data_file['w2_vector1'][time_idx,level,:]
    field_data_Y = data_file['w2_vector2'][time_idx,level,:]
    field_data_mag = np.sqrt(field_data_X**2 + field_data_Y**2)

    # ------------------------------------------------------------------------ #
    # Make dataframe to hold this data (allowing us to filter velocities)
    # ------------------------------------------------------------------------ #

    df = pd.DataFrame({'coords_X': data_file['Mesh2d_edge_x'][:],
                       'coords_Y': data_file['Mesh2d_edge_y'][:],
                       'field_data_X': field_data_X,
                       'field_data_Y': field_data_Y,
                       'field_data_mag': field_data_mag
                       })

    filter_df = df[df['field_data_mag'] > vector_magnitude_cutoff]

    data_file.close()

    # ------------------------------------------------------------------------ #
    # Plot data
    # ----------------------------------------------------------------------- #

    # Quivers
    coords_X, coords_Y = filter_df['coords_X'].values, filter_df['coords_Y'].values
    field_data_X, field_data_Y = filter_df['field_data_X'].values, filter_df['field_data_Y'].values
    _ = plot_field_quivers(ax, coords_X, coords_Y, field_data_X, field_data_Y,
                           scale=1.5, minlength=vector_magnitude_cutoff)

    # Contours based on magnitude of vectors
    coords_X, coords_Y = df['coords_X'].values, df['coords_Y'].values
    field_data_mag = df['field_data_mag'].values
    cmap, _ = automake_cmap(contours, colour_scheme)
    cf, _ = plot_contoured_field(ax, coords_X, coords_Y, field_data_mag,
                                 'tricontour', contours, cmap=cmap,
                                 plot_contour_lines=False)

    plot_cubed_sphere_panels(ax, linewidth=0.5)

    # ------------------------------------------------------------------------ #
    # Add labels
    # ------------------------------------------------------------------------ #

    cax = fig.add_axes([0.95, 0.25, 0.02, 0.5])
    cb = fig.colorbar(cf, cax=cax, orientation='vertical',
                      label=cbar_label, ticks=cbar_ticks)

    title = f'Time {time:06.0f} s'
    automake_field_title(ax, title, titlepad=20, minmax=True, field_data=field_data_mag)

    data_metadata = {'xlims':[-180,180],
                     'xlabel':r'$\lambda \ /$ deg',
                     'ylims':[-70,70],
                     'ylabel':r'$\vartheta \ /$ deg',
                     'xticks':[-180,180],
                     'xticklabels':[-180,180],
                     'yticks':[-70,70],
                     'yticklabels':[-70,70],
                     'xlabelpad': -20,
                     'ylabelpad':-20}
    automake_field_axis_labels(ax, data_metadata)

    # ------------------------------------------------------------------------ #
    # Save figure
    # ------------------------------------------------------------------------ #

    plotname = f'{plotdir}/rewrite_four_part_sbr_time_{time_idx:02d}.png'
    print(f'Saving figure to {plotname}')
    fig.savefig(plotname, bbox_inches='tight')
    plt.close()
