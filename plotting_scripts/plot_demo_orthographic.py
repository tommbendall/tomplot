"""
Demonstrates the orthographic projection
"""

from netCDF4 import Dataset
import numpy as np
from tomplot import make_field_plots, make_convergence_plots, make_time_series_plots

# ---------------------------------------------------------------------------- #
# Things that can be altered and parameters for the test case
# ---------------------------------------------------------------------------- #

results_dirname = 'demo_3_spherical_zonal'
plot_times = -1

# ---------------------------------------------------------------------------- #
# Derived things from options
# ---------------------------------------------------------------------------- #

data = Dataset('results/'+results_dirname+'/global_output.nc','r')
run_ids = data['run_id'][-2:]
data.close()

# ---------------------------------------------------------------------------- #
# Field plots
# ---------------------------------------------------------------------------- #

print('Making field plots')
for run_id in run_ids:
    print('Field plot %d' % run_id)
    field_names = ['F_0_zonal', 'F_0_meridional']

    make_field_plots(results_dirname, run_id, 'orthographic',
                     field_names, -1, 'xy', projection='orthographic',
                     spherical_centre=(0.0,0.0),
                     colour_scheme='RdBu_r', remove_contour='middle',
                     contours=np.linspace(-3.0, 3.0, 25))

    for lon in [-np.pi,-np.pi/2,0.0,np.pi/2,np.pi]:
        for lat in [-np.pi/2, -np.pi/4, 0.0,np.pi/4]:
            deg_lon, deg_lat = lon*180./np.pi, lat*180./np.pi
            testname = 'orthographic_%.0f_%.0f' % (deg_lon, deg_lat)
            make_field_plots(results_dirname, run_id, testname,
                             field_names, -1, 'xy', projection='orthographic',
                             spherical_centre=(lon,lat),
                             colour_scheme='RdBu_r', remove_contour='middle',
                             contours=np.linspace(-3.0, 3.0, 25))


    make_field_plots(results_dirname, run_id, 'plane',
                     field_names, -1, 'xy', projection=None,
                     remove_contour='middle',
                     colour_scheme='RdBu_r',
                     contours=np.linspace(-3.0, 3.0, 25))
