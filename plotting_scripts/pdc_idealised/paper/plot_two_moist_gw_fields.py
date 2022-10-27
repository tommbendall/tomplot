"""
Plots two moist gravity wave fields
"""

import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np
from tomplot import individual_field_contour_plot, extract_lfric_2D_data

# ---------------------------------------------------------------------------- #
# Things that can be altered and parameters for the test case
# ---------------------------------------------------------------------------- #

results_dirnames = ['moist_gw_data/D600_P600', 'moist_gw_data/D600_P150']
plotname = 'fig_X_two_moist_gw_fields'
cbar_label = r"$\theta_e' \ / $ K"
titles = [r'$\Delta x_{phys}=\Delta x_{dyn}$,',r'$\Delta x_{phys}=4\Delta x_{dyn}$,']
colour_scheme = 'RdBu_r'
extrusion_details = {'domain':'plane', 'extrusion':'linear',
                     'zmin':0.0, 'zmax':10000, 'topological_dimension':3}

field_min = -0.003
field_max = 0.003
xlabel = r'$x \ / $ km'
ylabels = [r'$z \ / $ km', None]
xlims = [-150,150]
ylims = [0,10]

contours = np.linspace(field_min, field_max, 13)

# ---------------------------------------------------------------------------- #
# Theta_e things
# ---------------------------------------------------------------------------- #

diag_field_names = ['m_v','exner','theta']

Lv = 2.501e6
cp = 1005.0

# Routine to compute diagnostic theta_e
def diagnostic_field(data_dict):
    m_v = data_dict['m_v']
    exner = data_dict['exner']
    theta = data_dict['theta']

    T = theta * exner
    exp_arg = Lv * m_v / (cp * T)
    theta_e = theta * np.exp(exp_arg)

    return theta_e

def exner_in_wth(exner):
    data_shape = np.shape(exner)
    exner_out = np.zeros((data_shape[0]+1, data_shape[1]))

    # Assume uniform extrusion
    exner_out[0,:] = 2.0*exner[0,:] - exner[1,:]
    exner_out[-1,:] = 2.0*exner[-1,:] - exner[-2,:]
    for j in range(1,data_shape[0]):
        exner_out[j,:] = 0.5*(exner[j,:]+exner[j-1,:])

    return exner_out

# ---------------------------------------------------------------------------- #
# Things that are likely the same
# ---------------------------------------------------------------------------- #

plotdir = '/data/users/tbendall/results/pdc_idealised_paper/figures'
slice_name = 'xz'
slice_idx = 0
time_idx = -1
init_time_idx = 0

# This is declared BEFORE figure and ax are initialised
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
font = {'size':24}
plt.rc('font',**font)

# ---------------------------------------------------------------------------- #
# Field plots
# ---------------------------------------------------------------------------- #

fig, axarray = plt.subplots(1,2,figsize=(13.5,5),sharey='row')

plotpath = f'{plotdir}/{plotname}.jpg'

for i, (ax, title, ylabel, results_dirname) \
    in enumerate(zip(axarray, titles, ylabels, results_dirnames)):

    init_data_file = Dataset('/data/users/tbendall/'+results_dirname+'/lfric_initial.nc','r')
    data_file = Dataset('/data/users/tbendall/'+results_dirname+'/lfric_diag.nc','r')

    # Set up dictionary to store multiple sets of data
    data_values = {}

    # Extract data
    for diag_field_name in diag_field_names:
        coords_X, coords_Y, field_data, time, \
        coord_labels, coord_lims, coord_ticks,  \
        slice_label = extract_lfric_2D_data(data_file, diag_field_name, time_idx,
                                        slice_name=slice_name, slice_idx=slice_idx,
                                        extrusion_details=extrusion_details)

        if diag_field_name == 'exner':
            data_values[diag_field_name] = exner_in_wth(field_data)
        else:
            data_values[diag_field_name] = field_data
    
    theta_data = diagnostic_field(data_values)

    # Get initial values 
    data_values = {}

    # Extract data
    for diag_field_name in diag_field_names:
        coords_X, coords_Y, field_data, data_metadata = \
            extract_lfric_2D_data(init_data_file, diag_field_name, init_time_idx,
                                  slice_name=slice_name, slice_idx=slice_idx,
                                  extrusion_details=extrusion_details)

        if diag_field_name == 'exner':
            data_values[diag_field_name] = exner_in_wth(field_data)
        else:
            data_values[diag_field_name] = field_data
    
    init_data = diagnostic_field(data_values)

    # Subtract base theta from everything
    field_data = np.zeros_like(theta_data)
    for i in range(len(theta_data[0,:])):
        field_data[:,i] = theta_data[:,i] - init_data[:,0]

    init_data_file.close()
    data_file.close()

    time = data_metadata['time']
    coord_labels = data_metadata['coord_labels']
    coord_lims = data_metadata['coord_lims']
    coord_ticks = data_metadata['coord_ticks']
    slice_label = data_metadata['slice_label']

    coords_X = coords_X*0.001
    coords_Y = coords_Y*0.001

    cf = individual_field_contour_plot(coords_X, coords_Y, field_data,
                                       ax=ax, time=time, slice_name=slice,
                                       title=title, title_method='minmax',
                                       slice_idx=0, slice_label=slice_label,
                                       contours=contours, no_cbar=True,
                                    #    colour_levels_scaling=(1.25,1.3),
                                    #    restricted_cmap='both',
                                       colour_scheme=colour_scheme,
                                       extend_cmap=False, titlepad=20,
                                       title_size=18,
                                       ylabelpad=-10, ylabel=ylabel,
                                       xlims=xlims, ylims=ylims, xlabel=xlabel,
                                       xticklabels=xlims, yticklabels=ylims)


# Move the subplots to the left to make space for colorbar
fig.subplots_adjust(right=0.9, wspace=0.125)
# Add colorbar in its own axis
cbar_ax = fig.add_axes([0.92, 0.11, 0.02, 0.77])
cb = fig.colorbar(cf, cax=cbar_ax, ticks=[field_min, field_max])
cb.set_label(cbar_label, labelpad=-50)

print(f'Plotting to {plotpath}')
fig.savefig(plotpath, bbox_inches='tight')
plt.close()

