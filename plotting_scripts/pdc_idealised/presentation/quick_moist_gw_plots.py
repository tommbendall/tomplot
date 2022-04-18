"""
Makes plots for the spherical stretchy convergence test case
"""

from netCDF4 import Dataset
import numpy as np
from tomplot import individual_field_contour_plot, extract_lfric_2D_data

# ---------------------------------------------------------------------------- #
# Things that can be altered and parameters for the test case
# ---------------------------------------------------------------------------- #

results_dirname = 'pdc_moist_gw'
plot_times = 'all'
base_plotname = 'results/pdc_idealised_paper/figures/quick'
cbar_label = r'$m_X \ / $ kg kg$^{-1}$'
colour_scheme = 'OrRd'
extrusion_details = {'domain':'plane', 'extrusion':'linear',
                     'zmin':0.0, 'zmax':10000, 'topological_dimension':3}

slice_name = 'xz'
slice_idx = 0
time_idx = 0
field_names = ['rho','exner','theta','m_v','m_cl', 'u1', 'u2', 'u3', 'theta_pert']

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
                                  title=field_name.replace('_',''), title_method=None,
                                  slice_idx=slice_idx, slice_label=slice_label,
                                  colour_levels_scaling=1.4,
                                  restricted_cmap='top', colour_scheme=colour_scheme,
                                  extend_cmap=False, titlepad=20)