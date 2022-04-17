"""
Plots final states of two Bryan and Fritsch bubble tests.
"""

import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np
from tomplot import individual_field_contour_plot, extract_lfric_2D_data

# ---------------------------------------------------------------------------- #
# Things that can be altered and parameters for the test case
# ---------------------------------------------------------------------------- #

base_results_dirname = 'consistent_moisture_paper/bryan_fritsch-'
results_opts = ['adv','consist']
time_idx = -1
plotname = 'fig_9_bf_fields'
cbar_label = r"$\theta_e \ / $ K"
titles = ['Advective','Consistent']
colour_scheme = 'OrRd'
extrusion_details = {'domain':'plane', 'extrusion':'linear',
                     'zmin':0.0, 'zmax':10000, 'topological_dimension':3}
field_names = ['m_v','theta','exner_in_wth']
Lv = 2.501e6
cp = 1005.0
theta_e_min = 319.0
theta_e_max = 327.0
step = 1.0
xlabel = r'$x \ / $ km'
ylabels = [r'$z \ / $ km', None]
xlims = [-5,5]
ylims = [0,10]

contours = np.arange(theta_e_min, theta_e_max+step, step=step)

# Routine to compute diagnostic theta_e
def diagnostic_field(data_dict):
    m_v = data_dict['m_v']
    theta = data_dict['theta']
    exner = data_dict['exner_in_wth']

    T = theta * exner
    exp_arg = Lv * m_v / (cp * T)
    theta_e = theta * np.exp(exp_arg)

    return theta_e

# ---------------------------------------------------------------------------- #
# Things that are likely the same
# ---------------------------------------------------------------------------- #

results_dirnames = [base_results_dirname+opt for opt in results_opts]
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

for i, (ax, results_dirname, title, ylabel) in enumerate(zip(axarray, results_dirnames, titles, ylabels)):
    # This code is all adapted from plot_control
    filename = 'results/'+results_dirname+'/raw_data/lfric_diag.nc'
    data_file = Dataset(filename, 'r')

    # Set up dictionary to store multiple sets of data
    data_values = {}

    # Extract data
    for field_name in field_names:
        coords_X, coords_Y, field_data, time, \
        coord_labels, coord_lims, coord_ticks,  \
        slice_label = extract_lfric_2D_data(data_file, field_name, time_idx,
                                           slice_name=slice, slice_idx=slice_idx,
                                           extrusion_details=extrusion_details)

        data_values[field_name] = field_data

    data_file.close()

    field_data = diagnostic_field(data_values)

    # Scale coordinate fields to km
    coords_X *= 0.001
    coords_Y *= 0.001

    # Set title manually to change format of mins and maxes
    minmax_title = 'min: %3.2f K, max: %3.2f K' % (np.min(field_data), np.max(field_data))

    cf = individual_field_contour_plot(coords_X, coords_Y, field_data,
                                       ax=ax, time=time, slice_name=slice,
                                       title=minmax_title, title_method=None,
                                       text=title, text_pos=(-0, 11.4),
                                       slice_idx=0, slice_label=slice_label,
                                       contours=contours, no_cbar=True,
                                       colour_levels_scaling=1.2,
                                       restricted_cmap='top', colour_scheme=colour_scheme,
                                       extend_cmap=False, remove_contour=320.0,
                                       ylabelpad=-20, ylabel=ylabel, titlepad=20,
                                       xlims=xlims, ylims=ylims, xlabel=xlabel,
                                       xticklabels=xlims, yticklabels=ylims)


# Move the subplots to the left to make space for colorbar
fig.subplots_adjust(right=0.9, wspace=0.1)
# Add colorbar in its own axis
cbar_ax = fig.add_axes([0.92, 0.11, 0.02, 0.77])
cb = fig.colorbar(cf, cax=cbar_ax, ticks=[theta_e_min, theta_e_max])
cb.set_label(cbar_label, labelpad=-30)

print(f'Plotting to {plotpath}')
fig.savefig(plotpath, bbox_inches='tight')
plt.close()

