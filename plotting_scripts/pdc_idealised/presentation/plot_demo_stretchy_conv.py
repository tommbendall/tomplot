"""
Makes plots for the spherical stretchy convergence test case
"""

from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
import cartopy.crs as ccrs
from tomplot import individual_field_contour_plot, extract_lfric_2D_data

# ---------------------------------------------------------------------------- #
# Things that can be altered and parameters for the test case
# ---------------------------------------------------------------------------- #

results_dirname = 'pdc_demo_data/stretchy_conv'
base_plotname = 'results/pdc_idealised_paper/figures/stretchy_rho_demo'
cbar_label = r'$m_X \ / $ kg kg$^{-1}$'
colour_scheme = 'OrRd'
extrusion_details = {'domain':'sphere', 'extrusion':'linear',
                     'zmin':0.0, 'zmax':2000, 'topological_dimension':3}

field_name = 'physics_rho'
field_min = 0.01
field_max = 0.07
step = 0.005
time_gap = 100
num_contours = int(np.floor((field_max - field_min) / step)) + 1
contours = np.linspace(field_min, field_max, num_contours)
titles = [r'Front', r'Back']

# ---------------------------------------------------------------------------- #
# Derived things from options
# ---------------------------------------------------------------------------- #

slice_name = 'xy'
slice_idx = 0

init_data_file = Dataset('/data/users/tbendall/'+results_dirname+'/lfric_initial.nc','r')
data_file = Dataset('/data/users/tbendall/'+results_dirname+'/lfric_diag.nc','r')

# Add 1 to include so that the initial data is time_idx = 0
time_idxs = range(len(data_file['time'][:])+1)

# This is declared BEFORE figure and ax are initialised
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
font = {'size':24}
plt.rc('font',**font)

# ---------------------------------------------------------------------------- #
# Field plots
# ---------------------------------------------------------------------------- #

import cartopy.crs as ccrs
crs = ccrs.PlateCarree()

for time_idx in time_idxs:

    fig = plt.figure(figsize=(12,6))

    plotname = f'{base_plotname}_{field_name}_{time_idx:02d}.jpg'
    print(f'Plotting to {plotname}')

    this_time_idx = 0 if time_idx == 0 else time_idx - 1
    this_data_file = init_data_file if time_idx == 0 else data_file

    for i, (spherical_centre, title) in enumerate(zip([(0.0,0.0),(np.pi,0.0)], titles)):

        # Hack to get colour showing on first plot (use last time point of other data!)
        if time_idx == 0 and i == 1:
            this_time_idx = -1
            this_data_file = data_file

        # Extract data
        coords_X, coords_Y, field_data, data_metadata = \
            extract_lfric_2D_data(this_data_file, field_name, this_time_idx,
                                  slice_name=slice_name, slice_idx=slice_idx,
                                  extrusion_details=extrusion_details,
                                  central_lon=spherical_centre[0])

        lon_centre, lat_centre = spherical_centre[0]*180.0/np.pi, spherical_centre[1]*180.0/np.pi
        ax = fig.add_subplot(1, 2, 1+i,
                             projection=ccrs.Orthographic(lon_centre, lat_centre))

        time = data_metadata['time']
        coord_labels = data_metadata['coord_labels']
        coord_lims = data_metadata['coord_lims']
        coord_ticks = data_metadata['coord_ticks']
        slice_label = data_metadata['slice_label']


        cf = individual_field_contour_plot(coords_X, coords_Y, field_data,
                                           slice_name=slice_name,
                                           plotname=plotname,
                                           title_method=None,
                                           slice_idx=slice_idx, slice_label=slice_label,
                                           contours=contours,
                                           ax=ax, projection='orthographic',
                                           spherical_centre=spherical_centre,
                                           colour_levels_scaling=1.4,
                                           no_cbar=True,
                                           restricted_cmap='top',
                                           colour_scheme=colour_scheme,
                                           extend_cmap=False,
                                           remove_contour=0.02,
                                           gridline_args={'xlocs':[0.,45.,135.,180.,-135.,-45.], 'linestyle':'--'})

        if i == 0:
            plt.text(-45.,60.0,r'45$^\mathrm{o}$W', ha='center', va='center', fontsize=10, transform=crs)
            plt.text(0.0,60.0,r'0$^\mathrm{o}$', ha='center', va='center', fontsize=10, transform=crs)
            plt.text(45.0,60.0,r'45$^\mathrm{o}$E', ha='center', va='center', fontsize=10, transform=crs)
        elif i == 1:
            plt.text(135.0,60.0,r'135$^\mathrm{o}$E', ha='center', va='center', fontsize=10, transform=crs)
            plt.text(180.0,60.0,r'180$^\mathrm{o}$', ha='center', va='center', fontsize=10, transform=crs)
            plt.text(-135.0,60.0,r'135$^\mathrm{o}$W', ha='center', va='center', fontsize=10, transform=crs)

        # I think titles don't work with orographic projection?
        ax.set_title(title)

    suptitle = r'$t=$ '+f'{(time_idx*time_gap):.0f} s'

    fig.suptitle(suptitle, y=0.88)

    # Move the subplots to the left to make space for colorbar
    fig.subplots_adjust(right=0.85)
    # Add colorbar in its own axis
    cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])
    fig.colorbar(cf, label=cbar_label, cax=cbar_ax)

    fig.savefig(plotname, bbox_inches='tight')
    plt.close()

