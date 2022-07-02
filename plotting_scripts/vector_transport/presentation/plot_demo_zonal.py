"""
Demonstrates the orthographic projection
"""

import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np
import cartopy.crs as ccrs
from tomplot import individual_field_contour_plot, extract_2D_data

# ---------------------------------------------------------------------------- #
# Things that can be altered and parameters for the test case
# ---------------------------------------------------------------------------- #

results_dirname = 'demo_3_spherical_zonal'
plot_times = 'all'
run_id = 0
field_name  = 'F_0_zonal'
cbar_label = r'$F_\lambda \ /$ m s$^{-1}$'
titles = [r'Front', r'Back']
colour_scheme = 'RdBu_r'

# ---------------------------------------------------------------------------- #
# Things that are likely the same
# ---------------------------------------------------------------------------- #

plotdir = 'results/'+results_dirname+'/figures'
testname = 'zonal_demo'
slice = 'xy'
slice_idx = 0
filename = 'results/'+results_dirname+'/nc_fields/field_output_'+str(run_id)+'.nc'
data_file = Dataset(filename, 'r')
time_idxs = range(len(data_file['time'][:]))
data_file.close()

# This is declared BEFORE figure and ax are initialised
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
font = {'size':24}
plt.rc('font',**font)

# ---------------------------------------------------------------------------- #
# Field plots
# ---------------------------------------------------------------------------- #

for time_idx in time_idxs:
    print('Making contour plot %d' % time_idx)

    fig = plt.figure(figsize=(16,8))

    plotname = '%s/both_%s_%s_time_%s.png' % (plotdir, testname,
                                              field_name, str(time_idx))

    for i, (spherical_centre, title) in enumerate(zip([(0.0,-np.pi/4),(np.pi,-np.pi/4)], titles)):
        # This code is all adapted from plot_control
        data_file = Dataset(filename, 'r')

        # Extract data
        coords_X, coords_Y, field_data, data_metadata = \
            extract_2D_data(data_file, field_name, time_idx,
                            slice_name=slice, slice_idx=slice_idx,
                            central_lon=spherical_centre[0])

        data_file.close()

        time = data_metadata['time']
        coord_labels = data_metadata['coord_labels']
        coord_lims = data_metadata['coord_lims']
        coord_ticks = data_metadata['coord_ticks']
        slice_label = data_metadata['slice_label']

        lon_centre, lat_centre = spherical_centre[0]*180.0/np.pi, spherical_centre[1]*180.0/np.pi
        ax = fig.add_subplot(1, 2, 1+i,
                             projection=ccrs.Orthographic(lon_centre, lat_centre))

        cf = individual_field_contour_plot(coords_X, coords_Y, field_data,
                                           testname=testname, plotname=plotname,
                                           time=time, slice_name=slice,
                                           slice_idx=0, slice_label=slice_label,
                                           contours=np.linspace(-3.0, 3.0, 25),
                                           ax=ax, projection='orthographic',
                                           colour_scheme=colour_scheme,
                                           remove_contour='middle',
                                           spherical_centre=spherical_centre,
                                           extend_cmap=False,
                                           no_cbar=True)
        ax.set_title(title)

    fig.suptitle('Time: %06.2f s' % time, y=0.88)


    # Move the subplots to the left to make space for colorbar
    fig.subplots_adjust(right=0.85)
    # Add colorbar in its own axis
    cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])
    fig.colorbar(cf, label=cbar_label, cax=cbar_ax)


    fig.savefig(plotname, bbox_inches='tight')
    plt.close()
