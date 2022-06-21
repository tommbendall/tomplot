"""
Makes plots for the spherical convergence test case
"""

from netCDF4 import Dataset
import numpy as np
from tomplot import make_field_plots, make_convergence_plots, make_time_series_plots

# ---------------------------------------------------------------------------- #
# Things that can be altered and parameters for the test case
# ---------------------------------------------------------------------------- #

base_results_dirname = 'conv_2_quads_dt'
field_labels = ['plain','rec','vort']
plot_testname = 'sphere_hooks'
dts = ['0pt1','0pt05','0pt02','0pt2','0pt5']

# ---------------------------------------------------------------------------- #
# Derived things from options
# ---------------------------------------------------------------------------- #

for dt in dts:

    results_dirname = base_results_dirname+dt

    data = Dataset('results/'+results_dirname+'/global_output.nc','r')
    run_ids = data['run_id'][:]
    num_setups = len(field_labels)
    data.close()

    # ------------------------------------------------------------------------ #
    # Convergence plots
    # ------------------------------------------------------------------------ #

    title = r'$\Delta t =$' + dt.replace('pt','.')

    if len(run_ids) > 1:
        print('Making convergence plots')
        field_names = ['F_'+str(i) for i in range(num_setups)]
        make_convergence_plots(results_dirname, 'rncells_per_dim', field_names, run_ids,
                               'L2_error', field_labels=field_labels, label_style='gradient_plain',
                               legend_bbox=(-0.2,1.2), legend_ncol=3,
                               testname=plot_testname, titles=title )
