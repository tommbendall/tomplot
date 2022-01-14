"""
Makes plots for the shallow water Williamson 2 test case,
"""

import numpy as np
from tomplot import make_field_plots, make_convergence_plots, make_time_series_plots

# ---------------------------------------------------------------------------- #
# Things that can be altered and parameters for the test case
# ---------------------------------------------------------------------------- #

base_results_dirname = 'sw2_RTCF1'
plot_times = -1
base_angle_list = [0.0]
base_res_list = [20,24,30,36,40]
results_options = ['plain','recovered','vorticity']
results_labels = ['plain', 'recovered', 'vorticity']

# ---------------------------------------------------------------------------- #
# Derived things from options
# ---------------------------------------------------------------------------- #

num_angles = len(base_angle_list)
num_resolutions = len(base_res_list)
res_list = [res for res in base_res_list for j in range(num_angles)]
angle_list = [base_angle_list[j] for i in range(num_resolutions) for j in range(num_angles)]
num_runs = len(res_list)

# ---------------------------------------------------------------------------- #
# Error labels
# ---------------------------------------------------------------------------- #

error_label_dict = {('L1_error_normalised', 'D'): r'$||D-D_{true}||_{L_1}/||D_{true}||_{L_1}$',
                    ('L1_error_normalised', 'u'): r'$||\textbf{u}-\textbf{u}_{true}||_{L_1}/||\textbf{u}_{true}||_{L_1}$',
                    ('L2_error_normalised', 'D'): r'$||D-D_{true}||_{L_2}/||D_{true}||_{L_2}$',
                    ('L2_error_normalised', 'u'): r'$||\textbf{u}-\textbf{u}_{true}||_{L_2}/||\textbf{u}_{true}||_{L_2}$',
                    ('Linf_error_normalised', 'D'): r'$||D-D_{true}||_{L_\infty}/||D_{true}||_{L_\infty}$',
                    ('Linf_error_normalised', 'u'): r'$||\textbf{u}-\textbf{u}_{true}||_{L_\infty}/||\textbf{u}_{true}||_{L_\infty}$',
                    ('mean_normalised', 'D'): r'$ (\bar{D}-\bar{D}_{true}) / \bar{D}_0$',
                    ('variance_normalised', 'D'): r'$($Var$(D)-$Var$(D_{true}))/$Var$(D_0)$',
                    ('max_normalised', 'D'): r'$(\max D - \max D_{true})/(\max D_0 - \min D_0)$',
                    ('min_normalised', 'D'): r'$(\min D - \min D_{true})/(\max D_0 - \min D_0)$',
                    ('total', 'D'): 'Normalised mass',
                    ('total', 'energy'): 'Normalised energy',
                    ('total', 'potential_enstrophy'): 'Normalised potential enstrophy',
                    ('total', 'vorticity'): r'Normalised $\int \zeta $d$x$',
                    ('total', 'div_u'): r'Normalised $\int \nabla \cdot u $d$x$'}

# ---------------------------------------------------------------------------- #
# Field plots
# ---------------------------------------------------------------------------- #

text_pos = (0.0, 0.675*np.pi)

for result_option in results_options:
    results_dirname = f'{base_results_dirname}_{result_option}'

    print('Making field plots')

    for run_id in range(num_runs):
        text = r'$\alpha$: %1.2f rad, $n$: %d' % (angle_list[run_id], res_list[run_id])

        make_field_plots(results_dirname, run_id, 'williamson_2',
                         'D', plot_times, 'xy',
                         contours=np.arange(1000.,3250.,step=250.),
                         time_method='days', text=text, text_pos=text_pos)
        make_field_plots(results_dirname, run_id, 'williamson_2',
                         'D_diff', plot_times,
                         'xy', contours=np.arange(-500.,600.,step=100.),
                         time_method='days', text=text, text_pos=text_pos)

# ---------------------------------------------------------------------------- #
# Convergence plots
# ---------------------------------------------------------------------------- #
    if num_resolutions > 1:
        print('Making convergence plots')

        list_of_run_ids = [[i + j*num_angles for j in range(num_resolutions)] for i in range(num_angles)]

        if len(angle_list) > 1:
            field_labels = [r'$\alpha$: %1.2f,' % angle for angle in base_angle_list]
        else:
            field_labels = None

        list_of_variables = ['D' for i in range(num_angles)]
        make_convergence_plots(results_dirname, 'rncells_per_dim', list_of_variables, list_of_run_ids,
                               'L2_error_normalised', testname='conv_williamson_2_D',
                               ylabels=r'$\ln\left[||D-D_{true}||_{L_2}/||D_{true}||_{L_2}\right]$',
                               field_labels=field_labels)
        list_of_variables = ['u' for i in range(num_angles)]
        make_convergence_plots(results_dirname, 'rncells_per_dim', list_of_variables, list_of_run_ids,
                               'L2_error_normalised', testname='conv_williamson_2_u',
                               ylabels=r'$\ln\left[||\textbf{u}-\textbf{u}_{true}||_{L_2}/||\textbf{u}_{true}||_{L_2}\right]$',
                               field_labels=field_labels)

# ---------------------------------------------------------------------------- #
# Time series plots
# ---------------------------------------------------------------------------- #

    print('Making time series plots')

    field_labels = [r'$\alpha$: %1.2f, $n$: %d'
                    % (angle_list[run_id], res_list[run_id]) for run_id in range(num_runs)]

    errors = ['L2_error_normalised', 'L1_error_normalised', 'Linf_error_normalised',
                     'mean_normalised', 'variance_normalised', 'max_normalised', 'min_normalised']
    for error in errors:
        make_time_series_plots(results_dirname, 'D', range(num_runs),
                               error, field_labels=field_labels,
                               testname='williamson_2_D',
                               ylabels=error_label_dict[(error,'D')],
                               time_units='days')
    errors = ['L2_error_normalised', 'L1_error_normalised', 'Linf_error_normalised']
    for error in errors:
        make_time_series_plots(results_dirname, 'u', range(num_runs),
                               error, field_labels=field_labels,
                               testname='williamson_2_u',
                               ylabels=error_label_dict[(error,'u')],
                               time_units='days')

# ---------------------------------------------------------------------------- #
# Final convergence test
# ---------------------------------------------------------------------------- #

plot_dirs = [f'{base_results_dirname}_{result_option}' for result_option in results_options]

make_convergence_plots(plot_dirs, 'rncells_per_dim', 'D', range(num_runs),
                       'L2_error_normalised', testname='conv_comparison_D',
                       ylabels=r'$\ln\left[||D-D_{true}||_{L_2}/||D_{true}||_{L_2}\right]$',
                       field_labels=results_labels)
make_convergence_plots(plot_dirs, 'rncells_per_dim', 'u', range(num_runs),
                       'L2_error_normalised', testname='conv_comparison_u',
                       ylabels=r'$\ln\left[||\textbf{u}-\textbf{u}_{true}||_{L_2}/||\textbf{u}_{true}||_{L_2}\right]$',
                       field_labels=results_labels)
