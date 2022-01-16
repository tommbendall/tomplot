"""
Makes plots for the cylindrical convergence test case
"""

from netCDF4 import Dataset
import numpy as np
from tomplot import make_field_plots, make_convergence_plots, make_time_series_plots

# ---------------------------------------------------------------------------- #
# Things that can be altered and parameters for the test case
# ---------------------------------------------------------------------------- #

results_dirname = 'conv_1_quads_cyl_curly'
plot_times = -1
field_labels = ['plain','vorticity','recovered']
plot_testname = 'cylinder_curly'

# ---------------------------------------------------------------------------- #
# Derived things from options
# ---------------------------------------------------------------------------- #

data = Dataset('results/'+results_dirname+'/global_output.nc','r')
run_ids = data['run_id'][:]
num_setups = len(field_labels)
data.close()

# ---------------------------------------------------------------------------- #
# Field plots
# ---------------------------------------------------------------------------- #

print('Making field plots')
for run_id in run_ids:
    for i in range(num_setups):
        print('Field plot %d %d' % (run_id, i))
        field_names = ['F_'+str(i)+'_zonal', 'F_'+str(i)+'_meridional']

        make_field_plots(results_dirname, run_id, plot_testname,
                         field_names, -1, 'xy', contours=np.linspace(-1.0, 4.0, 11))

# ---------------------------------------------------------------------------- #
# Convergence plots
# ---------------------------------------------------------------------------- #

if len(run_ids) > 1:
    print('Making convergence plots')
    field_names = ['F_'+str(i) for i in range(num_setups)]
    make_convergence_plots(results_dirname, 'rncells_per_dim', field_names, run_ids,
                           'L2_error', field_labels=field_labels,
                           testname=plot_testname)

    field_names = ['F_'+str(i)+'_zonal' for i in range(num_setups)]
    make_convergence_plots(results_dirname, 'rncells_per_dim', field_names, run_ids,
                           'L2_error', field_labels=field_labels,
                           testname=plot_testname+'_zonal')

    field_names = ['F_'+str(i)+'_meridional' for i in range(num_setups)]
    make_convergence_plots(results_dirname, 'rncells_per_dim', field_names, run_ids,
                           'L2_error', field_labels=field_labels,
                           testname=plot_testname+'_meridional')

# ---------------------------------------------------------------------------- #
# Time series plots
# ---------------------------------------------------------------------------- #

print('Making time series plots')
F_labels = ['F_'+str(i) for i in range(num_setups)]
make_time_series_plots(results_dirname, F_labels, run_ids[-1], 'L2',
                       field_labels=field_labels, testname=plot_testname)
