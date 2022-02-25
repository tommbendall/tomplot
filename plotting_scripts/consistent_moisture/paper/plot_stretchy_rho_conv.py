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
plotname = 'fig_2_rho_stretchy_conv'
cbar_label = r'$\rho_d \ / $ kg m$^{-3}$'
titles = [r'$t = 0$',r'$t=\tau/2$']
colour_scheme = 'RdBu_r'
extrusion_details = {'domain':'plane', 'extrusion':'linear',
                     'zmin':0.0, 'zmax':2000, 'topological_dimension':3}
field_name = 'rho'

field_mins = [0.5, 0.0]
field_maxes = [1.0, 4.0]
steps = [0.1, 0.25]
xlabel = r'$x \ / $ km'
ylabels = [r'$z \ / $ km', None]
xlims = [-1,1]
ylims = [0,2]

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

fig, axarray = plt.subplots(1,2,figsize=(16,8),sharey='row')

plotpath = f'{plotdir}/{plotname}.jpg'

for i, (ax, title, ylabel, field_min, field_max, step) in \
    enumerate(zip(axarray, titles, ylabels, field_mins, field_maxes, steps)):

    if i == 0: # Initial condition, open lfric_initial.nc
        filename = 'results/'+results_dirname+'/raw_data/lfric_initial.nc'
        data_file = Dataset(filename, 'r')
        time_idx = 0
    else:
        filename = 'results/'+results_dirname+'/raw_data/lfric_diag.nc'
        data_file = Dataset(filename, 'r')
        # Find halfway point in time
        time_idx = int(np.ceil((len(data_file['time_instant'][:])-2)/2))


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

    num_contours = int(np.floor((field_max - field_min) / step)) + 1
    contours = np.linspace(field_min, field_max, num_contours)

    cf = individual_field_contour_plot(coords_X, coords_Y, field_data,
                                       ax=ax, time=time, slice_name=slice,
                                       title=None, title_method=None,
                                       slice_idx=0, slice_label=slice_label,
                                       contours=contours, no_cbar=True,
                                       colour_levels_scaling=(1.25,1.3),
                                       restricted_cmap='both', colour_scheme=colour_scheme,
                                       extend_cmap=False,
                                       ylabelpad=-10, ylabel=ylabel,
                                       xlims=xlims, ylims=ylims, xlabel=xlabel,
                                       xticklabels=xlims, yticklabels=ylims)

    # The two plots use different colour bars. Put them above the plots
    #fig.subplots_adjust(top=0.9)
    if i == 0:
        cax = fig.add_axes([0.125, 0.0, 0.355, 0.04])
    else:
        cax = fig.add_axes([0.5475, 0.0, 0.355, 0.04])
    cb = fig.colorbar(cf, ticks=[field_min, field_max], cax=cax,
                      orientation='horizontal')
    cb.set_label(cbar_label, labelpad=-10)

    ax.set_title(title, pad=20)

print(f'Plotting to {plotpath}')
fig.savefig(plotpath, bbox_inches='tight')
plt.close()

