"""
Makes plots for the spherical convergence test case
"""

from netCDF4 import Dataset
import numpy as np
from tomplot import make_field_plots, make_convergence_plots, make_time_series_plots

# ---------------------------------------------------------------------------- #
# Things that can be altered and parameters for the test case
# ---------------------------------------------------------------------------- #

results_dirname = 'demo_6_convergence'
plot_times = -1
field_labels = ['RTCF2']
plot_testname = 'sphere_hooks'

# ---------------------------------------------------------------------------- #
# Derived things from options
# ---------------------------------------------------------------------------- #

data = Dataset('/data/users/tbendall/results/'+results_dirname+'/global_output.nc','r')
run_ids = data['run_id'][:]
num_setups = len(field_labels)
data.close()

# ---------------------------------------------------------------------------- #
# Convergence plots
# ---------------------------------------------------------------------------- #

if len(run_ids) > 1:
    print('Making convergence plots')
    field_names = ['F_'+str(i) for i in range(num_setups)]
    make_convergence_plots(results_dirname, 'rncells_per_dim', field_names, run_ids,
                           'L2_error', field_labels=field_labels,
                           testname=plot_testname)
