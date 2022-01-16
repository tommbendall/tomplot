"""
Makes plots for the shallow water Galewsky jet test case,
"""

import numpy as np
from tomplot import make_field_plots, make_convergence_plots, make_time_series_plots

# ---------------------------------------------------------------------------- #
# Things that can be altered and parameters for the test case
# ---------------------------------------------------------------------------- #

results_dirname = 'galewsky_RTCF1_vorticity'
plot_times = 'all'
base_angle_list = [np.pi/4]
base_res_list = [128]
perturb = True

# ---------------------------------------------------------------------------- #
# Derived things from options
# ---------------------------------------------------------------------------- #

num_angles = len(base_angle_list)
num_resolutions = len(base_res_list)
res_list = [res for res in base_res_list for j in range(num_angles)]
angle_list = [base_angle_list[j] for i in range(num_resolutions) for j in range(num_angles)]
num_runs = len(res_list)

# ---------------------------------------------------------------------------- #
# Field plots
# ---------------------------------------------------------------------------- #

print('Making field plots')
text_pos = (0.0, 0.53*np.pi)

for run_id in range(num_runs):
    text = r'$\alpha$: %1.2f rad, $n$: %d' % (angle_list[run_id], res_list[run_id])

    make_field_plots(results_dirname, run_id, 'galewsky',
                     'vorticity', plot_times, 'xy',
                     contours=np.arange(-1.8e-4,2.0e-4,step=2e-5),
                     ylims=[10*np.pi/180, 80*np.pi/180], yticklabels=[10,80],
                     ylabel=r'$\phi \ / $ deg', figsize=(15,6),
                     remove_contour=0.0,
                     time_method='days', text=text, text_pos=text_pos)

# ---------------------------------------------------------------------------- #
# Convergence plots
# ---------------------------------------------------------------------------- #
if not perturb and num_resolutions > 1:
    print('Making convergence plots')

    list_of_run_ids = [[i + j*num_angles for j in range(num_resolutions)] for i in range(num_angles)]

    if len(angle_list) > 1:
        field_labels = [r'$\alpha$: %1.2f,' % angle for angle in base_angle_list]
    else:
        field_labels = None

    list_of_variables = ['D' for i in angle_list]
    make_convergence_plots(results_dirname, 'dx', list_of_variables, list_of_run_ids,
                           'L2_error_normalised', testname='galewsky_D',
                           ylabels=r'$\ln\left[||D-D_{true}||_{L_2}/||D_{true}||_{L_2}\right]$',
                           field_labels=field_labels)
    list_of_variables = ['u' for i in angle_list]
    make_convergence_plots(results_dirname, 'dx', list_of_variables, list_of_run_ids,
                           'L2_error_normalised', testname='galewsky_u',
                           ylabels=r'$\ln\left[||\textbf{u}-\textbf{u}_{true}||_{L_2}/||\textbf{u}_{true}||_{L_2}\right]$',
                           field_labels=field_labels)

# -------------------------------------------------------------------- #
# Time series plots
# -------------------------------------------------------------------- #

if perturb:
    variables = ['D', 'vorticity']
    variable_labels = ['Normalised mass', r'Normalised $\int \zeta $d$x$']
    for variable, label in zip(variables, variable_labels):
        make_time_series_plots(results_dirname, variable, range(num_runs),
                               'total', testname='galewsky_'+variable,
                               ylabels=label, time_units='days', normalise=True)
