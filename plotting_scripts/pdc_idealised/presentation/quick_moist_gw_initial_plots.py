"""
Makes plots for the spherical stretchy convergence test case
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

results_dirname = 'pdc_moist_gw'
plot_times = 'all'
base_plotname = 'results/pdc_idealised_paper/figures/quick_init'
cbar_label = r'$m_X \ / $ kg kg$^{-1}$'
colour_scheme = 'OrRd'
extrusion_details = {'domain':'plane', 'extrusion':'linear',
                     'zmin':0.0, 'zmax':10000, 'topological_dimension':3}

slice_name = 'xz'
slice_idx = 0
time_idx = 0
field_names = ['theta','theta_pert','theta_e','theta_e_pert','m_v','m_cl','m_v_pert','m_cl_pert','u1']
diag_field_names = ['m_v','exner','theta']


# ---------------------------------------------------------------------------- #
# Derived things from options
# ---------------------------------------------------------------------------- #

data_file = Dataset('results/'+results_dirname+'/lfric_initial.nc','r')

# ---------------------------------------------------------------------------- #
# Field plots
# ---------------------------------------------------------------------------- #

for field_name in field_names:

    if field_name == 'theta_pert':
        # Extract data
        coords_X, coords_Y, theta_data, time, \
        coord_labels, coord_lims, coord_ticks,  \
        slice_label = extract_lfric_2D_data(data_file, 'theta', time_idx,
                                            slice_name=slice_name, slice_idx=slice_idx,
                                            extrusion_details=extrusion_details)
        # Subtract base theta from everything
        field_data = np.zeros_like(theta_data)
        for i in range(len(theta_data[0,:])):
            field_data[:,i] = theta_data[:,i] - theta_data[:,0]
    

    elif field_name == 'm_v_pert':
        # Extract data
        coords_X, coords_Y, m_v_data, time, \
        coord_labels, coord_lims, coord_ticks,  \
        slice_label = extract_lfric_2D_data(data_file, 'm_v', time_idx,
                                            slice_name=slice_name, slice_idx=slice_idx,
                                            extrusion_details=extrusion_details)
        # Subtract base theta from everything
        field_data = np.zeros_like(m_v_data)
        for i in range(len(m_v_data[0,:])):
            field_data[:,i] = m_v_data[:,i] - m_v_data[:,0]

    elif field_name == 'm_cl_pert':
        # Extract data
        coords_X, coords_Y, m_cl_data, time, \
        coord_labels, coord_lims, coord_ticks,  \
        slice_label = extract_lfric_2D_data(data_file, 'm_cl', time_idx,
                                            slice_name=slice_name, slice_idx=slice_idx,
                                            extrusion_details=extrusion_details)
        # Subtract base theta from everything
        field_data = np.zeros_like(m_cl_data)
        for i in range(len(m_cl_data[0,:])):
            field_data[:,i] = m_cl_data[:,i] - m_cl_data[:,0]

    elif field_name in ['theta_e','theta_e_pert']:

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

        if field_name == 'theta_e':
            field_data = diagnostic_field(data_values)
        elif field_name == 'theta_e_pert':
            theta_data = diagnostic_field(data_values)
            # Subtract base theta from everything
            field_data = np.zeros_like(theta_data)
            for i in range(len(theta_data[0,:])):
                field_data[:,i] = theta_data[:,i] - theta_data[:,0]

    else:
        # Extract data
        coords_X, coords_Y, field_data, time, \
        coord_labels, coord_lims, coord_ticks,  \
        slice_label = extract_lfric_2D_data(data_file, field_name, time_idx,
                                            slice_name=slice_name, slice_idx=slice_idx,
                                            extrusion_details=extrusion_details)


    plotname = f'{base_plotname}_{field_name}.jpg'

    print(field_name, np.min(field_data), np.max(field_data))
    print(coord_lims)

    individual_field_contour_plot(coords_X, coords_Y, field_data,
                                  time=time, slice_name=slice_name,
                                  plotname=plotname,
                                  title=field_name.replace('_',''), title_method='minmax',
                                  slice_idx=slice_idx, slice_label=slice_label,
                                  colour_levels_scaling=1.4,
                                  restricted_cmap='top', colour_scheme=colour_scheme,
                                  extend_cmap=False, titlepad=20)