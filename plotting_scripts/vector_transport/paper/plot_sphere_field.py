"""
Plots the initial condition and state halfway through a run for the
spherical "hooks" convergence test case.
"""

import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np
import cartopy.crs as ccrs
from tomplot import individual_quiver_plot, extract_2D_data

# ---------------------------------------------------------------------------- #
# Things that can be altered and parameters for the test case
# ---------------------------------------------------------------------------- #

results_dirname = 'vector_transport_paper/demo_4_spherical_hooks'
plotname = 'fig_4_sphere_fields'
plot_times = 'all'
run_id = 0
field_info = ('F_0', 'F_0_zonal', 'F_0_meridional')
cbar_label = r'$|F| \ /$ m s$^{-1}$'
titles = [r'$t=0, \ t=2T$', r'$t=T$']
colour_scheme = 'OrRd'

# ---------------------------------------------------------------------------- #
# Things that are likely the same
# ---------------------------------------------------------------------------- #

plotdir = 'results/vector_transport_paper/figures'
field_name, field_X_name, field_Y_name = field_info
slice = 'xy'
slice_idx = 0
filename = 'results/'+results_dirname+'/nc_fields/field_output_'+str(run_id)+'.nc'
data_file = Dataset(filename, 'r')
time_idxs = [0,int(np.ceil((len(data_file['time'][:])-1)/2))]
data_file.close()

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

fig = plt.figure(figsize=(16,8))

plotpath = f'{plotdir}/{plotname}.jpg'

for i, (spherical_centre, time_idx, title) in enumerate(zip([(0.0,0.0),(np.pi,0.0)], time_idxs, titles)):
    # This code is all adapted from plot_control
    data_file = Dataset(filename, 'r')

    # Extract data
    coords_X, coords_Y, field_X_data, data_metadata = \
        extract_2D_data(data_file, field_X_name, time_idx,
                        slice_name=slice, slice_idx=slice_idx,
                        central_lon=spherical_centre[0])

    coords_X, coords_Y, field_Y_data, data_metadata = \
        extract_2D_data(data_file, field_Y_name, time_idx,
                        slice_name=slice, slice_idx=slice_idx,
                        central_lon=spherical_centre[0])

    time = data_metadata['time']
    coord_labels = data_metadata['coord_labels']
    coord_lims = data_metadata['coord_lims']
    coord_ticks = data_metadata['coord_ticks']
    slice_label = data_metadata['slice_label']

    data_file.close()

    lon_centre, lat_centre = spherical_centre[0]*180.0/np.pi, spherical_centre[1]*180.0/np.pi
    ax = fig.add_subplot(1, 2, 1+i,
                         projection=ccrs.Orthographic(lon_centre, lat_centre))

    cf = individual_quiver_plot(coords_X, coords_Y, field_X_data, field_Y_data,
                                time=time, slice_name=slice,
                                slice_idx=0, slice_label=slice_label,
                                contours=np.linspace(0.0, 3.0, 11),
                                contour_method='magnitude',
                                x_offset=1,
                                quiver_npts=2, scale=5e-6,
                                restrict_quivers=True,
                                restricted_cmap='top',
                                colour_levels_scaling=1.2,
                                ax=ax, projection='orthographic',
                                colour_scheme=colour_scheme,
                                spherical_centre=spherical_centre,
                                extend_cmap=False,
                                no_cbar=True,
                                gridline_args={'xlocs':[0.,45.,135.,180.,-135.,-45.], 'linestyle':'--'})

    if i == 0:
        plt.text(-45.,0.0,r'45$^\mathrm{o}$W', ha='center', va='center', fontsize=16, transform=crs)
        plt.text(0.0,0.0,r'0$^\mathrm{o}$', ha='center', va='center', fontsize=16, transform=crs)
        plt.text(45.0,0.0,r'45$^\mathrm{o}$E', ha='center', va='center', fontsize=16, transform=crs)
    elif i == 1:
        plt.text(135.0,0.0,r'135$^\mathrm{o}$E', ha='center', va='center', fontsize=16, transform=crs)
        plt.text(180.0,0.0,r'180$^\mathrm{o}$', ha='center', va='center', fontsize=16, transform=crs)
        plt.text(-135.0,0.0,r'135$^\mathrm{o}$W', ha='center', va='center', fontsize=16, transform=crs)

    ax.set_title(title)

# Move the subplots to the left to make space for colorbar
fig.subplots_adjust(right=0.9, wspace=0.05)
# Add colorbar in its own axis
cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])
cb = fig.colorbar(cf, cax=cbar_ax, ticks=[0, 3])
cb.set_label(cbar_label, labelpad=-10)

print(f'Plotting to {plotpath}')
fig.savefig(plotpath, bbox_inches='tight')
plt.close()
