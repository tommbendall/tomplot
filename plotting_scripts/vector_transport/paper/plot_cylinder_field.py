"""
Plots initial condition and state halfway through a run for the cylindrical
deformational convergence test, showing quivers for the vector.
"""

import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np
from tomplot import individual_quiver_plot, extract_2D_data

# ---------------------------------------------------------------------------- #
# Things that can be altered and parameters for the test case
# ---------------------------------------------------------------------------- #

results_dirname = 'demo_5_cylindrical_curly'
plot_times = 'all'
field_labels = ['RTCF2']
plotname = 'fig_1_cylinder_fields'
cbar_label = r'$|F| \ / $ m s$^{-1}$'
titles = [r'$t=0$', r'$t=\tau/2$']
colour_scheme = 'OrRd'
run_id = 0
ylabels = [r'$z \ / $ m', None]
field_info = ('F_0', 'F_0_zonal', 'F_0_meridional')

# ---------------------------------------------------------------------------- #
# Things that are likely the same
# ---------------------------------------------------------------------------- #

plotdir = 'results/vector_transport_paper/figures'
field_name, field_X_name, field_Y_name = field_info
slice = 'xy'
slice_idx = 0
filename = 'results/'+results_dirname+'/nc_fields/field_output_'+str(run_id)+'.nc'
data_file = Dataset(filename, 'r')
# Times are initial time and time halfway through
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

fig, axarray = plt.subplots(1,2,figsize=(16,8),sharey='row')

plotpath = f'{plotdir}/{plotname}.jpg'

for i, (ax, time_idx, title, ylabel) in enumerate(zip(axarray, time_idxs, titles, ylabels)):
    # This code is all adapted from plot_control
    data_file = Dataset(filename, 'r')

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
                                ax=ax, time=time, slice_name=slice, title=title,
                                title_method=None,
                                slice_idx=0, slice_label=slice_label,
                                contours=np.arange(0.0, 5.0, 0.5),
                                contour_method='magnitude',
                                no_cbar=True, colour_levels_scaling=1.2,
                                restricted_cmap='top', colour_scheme=colour_scheme,
                                extend_cmap=False, quiver_npts=2, scale=0.8,
                                ylabelpad=-30, xlabel=r'$\varrho \phi \ / $ m',
                                ylabel=ylabel, titlepad=20,
                                xlims=coord_lims[0], ylims=coord_lims[1],
                                xticklabels=coord_ticks[0], yticklabels=coord_ticks[1])


# Move the subplots to the left to make space for colorbar
fig.subplots_adjust(right=0.9, wspace=0.1)
# Add colorbar in its own axis
cbar_ax = fig.add_axes([0.92, 0.11, 0.02, 0.77])
cb = fig.colorbar(cf, cax=cbar_ax, ticks=[0, 4.5])
cb.set_label(cbar_label, labelpad=-30)

print(f'Plotting to {plotpath}')
fig.savefig(plotpath, bbox_inches='tight')
plt.close()

