"""
Makes plots for the cylindrical convergence test case
"""

from netCDF4 import Dataset
import numpy as np
from tomplot import make_field_plots, make_convergence_plots, make_time_series_plots

# ---------------------------------------------------------------------------- #
# Things that can be altered and parameters for the test case
# ---------------------------------------------------------------------------- #

shape = 'quads' # or 'tris'
results_dirname = 'conv_1_quads_cyl_curly' if shape == 'quads' else 'conv_1_tris_cyl_curly'
field_labels = ['plain','rec','vort'] if shape == 'quads' else ['plain','rec']
plot_testname = 'cylinder_curly'

# ---------------------------------------------------------------------------- #
# Derived things from options
# ---------------------------------------------------------------------------- #

data = Dataset('results/'+results_dirname+'/global_output.nc','r')
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
                           'L2_error', field_labels=field_labels, label_style='gradient_plain',
                           testname=plot_testname, legend_bbox=(-0.2,1.2), legend_ncol=3)
