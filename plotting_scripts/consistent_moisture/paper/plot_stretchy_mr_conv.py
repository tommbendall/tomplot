"""
Plots initial state and state halfway through run from the convergence
configuration of the stretchy test.
"""

import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np
from tomplot import individual_field_contour_plot, extract_lfric_2D_data

# ---------------------------------------------------------------------------- #
# Things that can be altered and parameters for the test case
# ---------------------------------------------------------------------------- #

results_dirname = 'consistent_moisture_paper/stretchy-conv-cons-240'
plotname = 'fig_5_mr_stretchy_conv'
cbar_label = r'$m_v \ / $ kg kg$^{-1}$'
titles = [r'$t = 0$',r'$t=\tau/2$',r'$t=\tau$']
colour_scheme = 'OrRd'
extrusion_details = {'domain':'plane', 'extrusion':'linear',
                     'zmin':0.0, 'zmax':2000, 'topological_dimension':3}
field_name = 'm_v'

field_min = 0.01
field_max = 0.07
step = 0.005
xlabel = r'$x \ / $ km'
ylabels = [r'$z \ / $ km', None, None]
xlims = [-1,1]
ylims = [0,2]

num_contours = int(np.floor((field_max - field_min) / step)) + 1
contours = np.linspace(field_min, field_max, num_contours)

# ---------------------------------------------------------------------------- #
# Things that are likely the same
# ---------------------------------------------------------------------------- #

plotdir = 'results/consistent_moisture_paper/figures'
slice = 'xz'
slice_idx = 0

# This is declared BEFORE figure and ax are initialised
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
font = {'size':24}
plt.rc('font',**font)

# ---------------------------------------------------------------------------- #
# Field plots
# ---------------------------------------------------------------------------- #

fig, axarray = plt.subplots(1,3,figsize=(18,6),sharey='row')

plotpath = f'{plotdir}/{plotname}.jpg'

for i, (ax, title, ylabel) in enumerate(zip(axarray, titles, ylabels)):

    if i == 0: # Initial condition, open lfric_initial.nc
        filename = 'results/'+results_dirname+'/raw_data/lfric_initial.nc'
        data_file = Dataset(filename, 'r')
        time_idx = 0
    else:
        filename = 'results/'+results_dirname+'/raw_data/lfric_diag.nc'
        data_file = Dataset(filename, 'r')
        if i == 1:
            # Find halfway point in time
            time_idx = int(np.ceil((len(data_file['time_instant'][:])-2)/2))
        else:
            time_idx = -1


    # Extract data
    coords_X, coords_Y, field_data, time, \
    coord_labels, coord_lims, coord_ticks,  \
    slice_label = extract_lfric_2D_data(data_file, field_name, time_idx,
                                        slice_name=slice, slice_idx=slice_idx,
                                        extrusion_details=extrusion_details)

    data_file.close()

    # Scale coordinate fields to km
    coords_X *= 0.001
    coords_Y *= 0.001

    cf = individual_field_contour_plot(coords_X, coords_Y, field_data,
                                       ax=ax, time=time, slice_name=slice,
                                       title=title, title_method=None,
                                       slice_idx=0, slice_label=slice_label,
                                       contours=contours, no_cbar=True,
                                       colour_levels_scaling=1.4,
                                       restricted_cmap='top', colour_scheme=colour_scheme,
                                       extend_cmap=False, remove_contour=0.02,
                                       ylabelpad=-10, ylabel=ylabel, titlepad=20,
                                       xlims=xlims, ylims=ylims, xlabel=xlabel,
                                       xticklabels=xlims, yticklabels=ylims)


# Move the subplots to the left to make space for colorbar
fig.subplots_adjust(right=0.9, wspace=0.1)
# Add colorbar in its own axis
cbar_ax = fig.add_axes([0.92, 0.11, 0.02, 0.77])
cb = fig.colorbar(cf, cax=cbar_ax, ticks=[field_min, field_max])
cb.set_label(cbar_label, labelpad=-30)

print(f'Plotting to {plotpath}')
fig.savefig(plotpath, bbox_inches='tight')
plt.close()

