"""
Demonstrates the orthographic projection
"""

import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np
import cartopy.crs as ccrs
from tomplot import individual_quiver_plot, extract_2D_data

# ---------------------------------------------------------------------------- #
# Things that can be altered and parameters for the test case
# ---------------------------------------------------------------------------- #

results_dirname = 'demo_3_spherical_zonal'
plot_times = 'all'
run_id = 0
field_info = ('F_0', 'F_0_zonal', 'F_0_meridional')
cbar_label = r'$|F| \ /$ m s$^{-1}$'
titles = [r'Front', r'Back']
colour_scheme = 'OrRd'

# ---------------------------------------------------------------------------- #
# Things that are likely the same
# ---------------------------------------------------------------------------- #

plotdir = 'results/'+results_dirname+'/figures'
field_name, field_X_name, field_Y_name = field_info
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
    print('Making quiver plot %d' % time_idx)

    fig = plt.figure(figsize=(16,8))

    plotname = '%s/quiver_%s_%s_time_%02d.png' % (plotdir, testname,
                                                  field_name, time_idx)

    for i, (spherical_centre, title) in enumerate(zip([(0.0,-np.pi/4),(np.pi,-np.pi/4)], titles)):
        # This code is all adapted from plot_control
        data_file = Dataset(filename, 'r')

        # Extract data
        coords_X, coords_Y, field_X_data, time, \
        coord_labels, coord_lims, coord_ticks,  \
        slice_label = extract_2D_data(data_file, field_X_name, time_idx,
                                      slice_name=slice, slice_idx=slice_idx,
                                      central_lon=spherical_centre[0])

        coords_X, coords_Y, field_Y_data, time, \
        coord_labels, coord_lims, coord_ticks,  \
        slice_label = extract_2D_data(data_file, field_Y_name, time_idx,
                                      slice_name=slice, slice_idx=slice_idx,
                                      central_lon=spherical_centre[0])

        data_file.close()

        lon_centre, lat_centre = spherical_centre[0]*180.0/np.pi, spherical_centre[1]*180.0/np.pi
        ax = fig.add_subplot(1, 2, 1+i,
                             projection=ccrs.Orthographic(lon_centre, lat_centre))

        cf = individual_quiver_plot(coords_X, coords_Y, field_X_data, field_Y_data,
                                    testname=testname, plotname=plotname,
                                    time=time, slice_name=slice,
                                    slice_idx=0, slice_label=slice_label,
                                    contours=np.linspace(0.0, 3.0, 16),
                                    contour_method='magnitude',
                                    quiver_npts=2, scale=5e-6,
                                    restrict_quivers=True,
                                    ax=ax, projection='orthographic',
                                    colour_scheme=colour_scheme,
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
