"""
Makes field plots for the moist gravity wave test case
"""

from netCDF4 import Dataset
import numpy as np
from tomplot import individual_field_contour_plot, extract_lfric_2D_data


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
# Things that can be altered and parameters for the test case
# ---------------------------------------------------------------------------- #

results_dirname = 'pdc_demo_data/moist_gw_60'
plot_times = 'all'
base_plotname = '/data/users/tbendall/results/pdc_idealised_paper/figures/moist_gw_demo'
cbar_label = r"$\theta_e' \ / $ K"
colour_scheme = 'RdBu_r'
extrusion_details = {'domain':'plane', 'extrusion':'linear',
                     'zmin':0.0, 'zmax':10000, 'topological_dimension':3}

slice_name = 'xz'
slice_idx = 0
time_gap = 180.0
init_time_idx = 0
field_name = 'theta_e_pert'
diag_field_names = ['m_v','exner','theta']
xlabel = r'$x \ / $ km'
ylabel = r'$z \ / $ km'
xlims = [-150,150]
ylims = [0,10]
field_min = -0.01
field_max = 0.01
contour_levels = np.linspace(field_min, field_max, 41)
cbar_ticks = [field_min, field_max]


# ---------------------------------------------------------------------------- #
# Derived things from options
# ---------------------------------------------------------------------------- #

init_data_file = Dataset('/data/users/tbendall/results/'+results_dirname+'/lfric_initial.nc','r')
data_file = Dataset('/data/users/tbendall/results/'+results_dirname+'/lfric_diag.nc','r')

# Add 1 to include so that the initial data is time_idx = 0 and don't do every time!
time_idxs = range(0,len(data_file['time'][:])+1,2)

# ---------------------------------------------------------------------------- #
# Field plots
# ---------------------------------------------------------------------------- #

for time_idx in time_idxs:

    # ------------------------------------------------------------------------ #
    # Theta_e full field data
    # ------------------------------------------------------------------------ #
    # Set up dictionary to store multiple sets of data
    data_values = {}

    this_data_file = init_data_file if time_idx == 0 else data_file

    # Extract data
    for diag_field_name in diag_field_names:
        coords_X, coords_Y, field_data, data_metadata = \
            extract_lfric_2D_data(this_data_file, diag_field_name, time_idx-1,
                                  slice_name=slice_name, slice_idx=slice_idx,
                                  extrusion_details=extrusion_details)

        if diag_field_name == 'exner':
            data_values[diag_field_name] = exner_in_wth(field_data)
        else:
            data_values[diag_field_name] = field_data

    theta_data = diagnostic_field(data_values)

    # ------------------------------------------------------------------------ #
    # Theta_e background field data
    # ------------------------------------------------------------------------ #
    # Get initial values 
    data_values = {}

    # Extract background data
    for diag_field_name in diag_field_names:
        coords_X, coords_Y, field_data, data_metadata = \
            extract_lfric_2D_data(init_data_file, diag_field_name, init_time_idx,
                                  slice_name=slice_name, slice_idx=slice_idx,
                                  extrusion_details=extrusion_details)

        if diag_field_name == 'exner':
            data_values[diag_field_name] = exner_in_wth(field_data)
        else:
            data_values[diag_field_name] = field_data
    
    time = data_metadata['time']
    coord_labels = data_metadata['coord_labels']
    coord_lims = data_metadata['coord_lims']
    coord_ticks = data_metadata['coord_ticks']
    slice_label = data_metadata['slice_label']

    init_data = diagnostic_field(data_values)

    # Subtract base theta from everything
    field_data = np.zeros_like(theta_data)
    for i in range(len(theta_data[0,:])):
        field_data[:,i] = theta_data[:,i] - init_data[:,0]
    
    coords_X = coords_X*0.001
    coords_Y = coords_Y*0.001

    plotname = f'{base_plotname}_{field_name}_{int(time_idx/2):02d}.jpg'
    print(f'Plotting to {plotname}')

    individual_field_contour_plot(coords_X, coords_Y, field_data, figsize=(10,6),
                                  time=time, slice_name=slice_name,
                                  plotname=plotname,
                                  title=r'$t=$ '+f'{(time_idx*time_gap):.0f} s',
                                  title_method=None,
                                  slice_idx=slice_idx, slice_label=slice_label,
                                  contours=contour_levels,
                                  colour_scheme=colour_scheme,
                                  cbar_ticks=cbar_ticks, cbar_label=cbar_label,
                                  extend_cmap=False, titlepad=20, ylabel=ylabel,
                                  xlims=xlims, ylims=ylims, xlabel=xlabel,
                                  xticklabels=xlims, yticklabels=ylims,
                                  ylabelpad=-20, cbar_labelpad=-40)