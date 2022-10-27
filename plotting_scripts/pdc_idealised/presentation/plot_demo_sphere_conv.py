"""
Makes plots for the spherical stretchy convergence test case
"""

from netCDF4 import Dataset
import numpy as np
from tomplot import individual_field_contour_plot, extract_lfric_2D_data

# ---------------------------------------------------------------------------- #
# Things that can be altered and parameters for the test case
# ---------------------------------------------------------------------------- #

results_dirname = 'pdc_idealised_paper/stretchy_sphere_conv_C96_C48/raw_data'
plot_times = 'all'
base_plotname = '/data/users/tbendall/results/pdc_idealised_paper/stretchy_sphere_conv_all/stretchy_rho'
cbar_label = r'$m_X \ / $ kg kg$^{-1}$'
colour_scheme = 'OrRd'
extrusion_details = {'domain':'sphere', 'extrusion':'linear',
                     'zmin':0.0, 'zmax':2000, 'topological_dimension':3}

field_name = 'physics_rho'
field_min = 0.01
field_max = 0.07
step = 0.005
num_contours = int(np.floor((field_max - field_min) / step)) + 1
contours = np.linspace(field_min, field_max, num_contours)

# ---------------------------------------------------------------------------- #
# Derived things from options
# ---------------------------------------------------------------------------- #

slice_name = 'xy'
slice_idx = 0
data_file = Dataset('/data/users/tbendall/results/'+results_dirname+'/lfric_diag.nc','r')
time_idxs = range(len(data_file['time'][:]))

# ---------------------------------------------------------------------------- #
# Field plots
# ---------------------------------------------------------------------------- #

for time_idx in time_idxs:
    # Extract data
    coords_X, coords_Y, field_data, data_metadata = \
        extract_lfric_2D_data(data_file, field_name, time_idx,
                              slice_name=slice_name, slice_idx=slice_idx,
                              extrusion_details=extrusion_details)

    title = f't = {time:.1f}'
    plotname = f'{base_plotname}_{time_idx}.jpg'

    time = data_metadata['time']
    coord_labels = data_metadata['coord_labels']
    coord_lims = data_metadata['coord_lims']
    coord_ticks = data_metadata['coord_ticks']
    slice_label = data_metadata['slice_label']

    individual_field_contour_plot(coords_X, coords_Y, field_data,
                                  time=time, slice_name=slice_name,
                                  plotname=plotname,
                                  title=title, title_method=None,
                                  slice_idx=slice_idx, slice_label=slice_label,
                                  contours=contours, cbar_label=cbar_label,
                                  colour_levels_scaling=1.4,
                                  restricted_cmap='top', colour_scheme=colour_scheme,
                                  extend_cmap=False, remove_contour=0.02, titlepad=20)