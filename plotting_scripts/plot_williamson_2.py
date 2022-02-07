"""
Makes plots for the shallow water Williamson 2 test case,
"""

import numpy as np
from tomplot import make_field_plots, make_convergence_plots, make_time_series_plots

# ---------------------------------------------------------------------------- #
# Things that can be altered and parameters for the test case
# ---------------------------------------------------------------------------- #

shape = 'quads' # or 'tris'
base_results_dirname = 'conv_3_quads_sw2_RTCF1' if shape == 'quads' else 'conv_3_tris_sw2_RT1'
results_options = ['upwind','recovered','vorticity']
results_labels = ['plain','rec','vort']
base_angle_list = [0.0]
base_res_list = [20,24,30,36,40]

# ---------------------------------------------------------------------------- #
# Derived things from options
# ---------------------------------------------------------------------------- #

num_angles = len(base_angle_list)
num_resolutions = len(base_res_list)
res_list = [res for res in base_res_list for j in range(num_angles)]
angle_list = [base_angle_list[j] for i in range(num_resolutions) for j in range(num_angles)]
num_runs = len(res_list)

# ---------------------------------------------------------------------------- #
# Final convergence test
# ---------------------------------------------------------------------------- #

plot_dirs = [f'{base_results_dirname}_{result_option}' for result_option in results_options]

make_convergence_plots(plot_dirs, 'rncells_per_dim', 'D', range(num_runs),
                       'L2_error_normalised', testname='conv_comparison_D',
                       ylabels=r'$\ln\left[||D-D_{true}||_{L_2}/||D_{true}||_{L_2}\right]$',
                       field_labels=results_labels, legend_bbox=(-0.2,1.2), legend_ncol=3,
                       label_style='gradient_plain')
make_convergence_plots(plot_dirs, 'rncells_per_dim', 'u', range(num_runs),
                       'L2_error_normalised', testname='conv_comparison_u',
                       ylabels=r'$\ln\left[||\textbf{u}-\textbf{u}_{true}||_{L_2}/||\textbf{u}_{true}||_{L_2}\right]$',
                       field_labels=results_labels, legend_bbox=(-0.2,1.2), legend_ncol=3, label_style='gradient_plain')
