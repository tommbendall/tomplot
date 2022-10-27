"""
Makes plots for the spherical convergence test case
"""

from netCDF4 import Dataset
import numpy as np
from tomplot import make_field_plots, make_convergence_plots, make_time_series_plots

# ---------------------------------------------------------------------------- #
# Things that can be altered and parameters for the test case
# ---------------------------------------------------------------------------- #

shape = 'quads' # 'quads' or 'tris'
results_dirname = 'conv_2_quads_hooks_trial' if shape == 'quads' else 'conv_2_tris_hooks_trial'
field_labels = ['Benchmark','Recovered','Vorticity'] if shape == 'quads' else ['Benchmark','Recovered']
leg_xcentre = 0.45 if shape == 'quads' else 0.5
plot_testname = 'sphere_hooks'
title = 'Quadrilateral cells' if shape == 'quads' else 'Triangular cells'

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
    make_convergence_plots(results_dirname, 'dx', field_names, run_ids,
                           'L2_error', field_labels=field_labels, label_style='gradient_plain',
                           testname=plot_testname, legend_bbox=(leg_xcentre,1.14), legend_ncol=3,
                           titlepad=60, titles=title, format='jpg', dpi=300,
                           ylabels=r'$\ln(||\mathbf{F}-\mathbf{F}_{true}||)$',
                           leg_col_spacing=0.1, leg_fontsize=20)
