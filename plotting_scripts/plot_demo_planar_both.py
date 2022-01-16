"""
Makes quiver plots for the planar demo case
"""

from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
from tomplot import individual_quiver_plot, extract_2D_data

# ---------------------------------------------------------------------------- #
# Things that can be altered and parameters for the test case
# ---------------------------------------------------------------------------- #

results_dirnames = ['demo_1_RTCF1_planar', 'demo_2_RTCF2_planar']
plot_times = 'all'
run_id = 2
field_info  = ('F_0', 'F_0_x', 'F_0_y')
cbar_label = r'$|F| \ /$ m s$^{-1}$'
titles = [r'RTCF$_1$', r'RTCF$_2$']

# ---------------------------------------------------------------------------- #
# Things that are likely the same
# ---------------------------------------------------------------------------- #

plotdir = 'results/'+results_dirnames[0]+'/figures'
field_name, field_X_name, field_Y_name = field_info
testname = 'planar_demo'
slice = 'xy'
slice_idx = 0
filename = 'results/'+results_dirnames[0]+'/nc_fields/field_output_'+str(run_id)+'.nc'
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

    fig, axarray = plt.subplots(1, 2, sharex='col', sharey='row', figsize=(16,8))

    plotname = '%s/both_%s_%s_slice_%s_run_%s_time_%s.png' % (plotdir, testname, field_name, slice,
                                                              str(run_id), str(time_idx))


    for ax, dirname, quiver_npts, title in zip(axarray, results_dirnames, [2,1], titles):
        # This code is all adapted from plot_control
        filename = 'results/'+dirname+'/nc_fields/field_output_'+str(run_id)+'.nc'
        data_file = Dataset(filename, 'r')
        time_idxs = range(len(data_file['time'][:]))

        # Extract data
        coords_X, coords_Y, field_X_data, time, \
        coord_labels, coord_lims, coord_ticks,  \
        slice_label = extract_2D_data(data_file, field_X_name, time_idx,
                                      slice_name=slice, slice_idx=slice_idx)

        coords_X, coords_Y, field_Y_data, time, \
        coord_labels, coord_lims, coord_ticks,  \
        slice_label = extract_2D_data(data_file, field_Y_name, time_idx,
                                      slice_name=slice, slice_idx=slice_idx)

        data_file.close()

        cf = individual_quiver_plot(coords_X, coords_Y, field_X_data, field_Y_data,
                                    testname=testname, plotname=plotname, time=time,
                                    field_name=field_name, slice_name=slice,
                                    slice_idx=slice_idx, slice_label=slice_label, ax=ax,
                                    contours=np.arange(0.0, 2.4, 0.2), contour_method='magnitude',
                                    quiver_npts=quiver_npts, scale=0.5,
                                    title=title, no_cbar=True,
                                    xlabel=coord_labels[0], ylabel=coord_labels[1],
                                    xlims=coord_lims[0], ylims=coord_lims[1],
                                    xticklabels=coord_ticks[0], yticklabels=coord_ticks[1])

    fig.suptitle('Time: %.2f s' % time)

    # Move the subplots to the left to make space for colorbar
    fig.subplots_adjust(right=0.9)
    # Add colorbar in its own axis
    cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])
    fig.colorbar(cf, label=cbar_label, cax=cbar_ax)


    fig.savefig(plotname, bbox_inches='tight')
    plt.close()

