"""
Makes plots for the shallow water Galewsky jet test case,
"""

import numpy as np
from tomplot import make_field_plots, make_convergence_plots, make_time_series_plots

# ---------------------------------------------------------------------------- #
# Things that can be altered and parameters for the test case
# ---------------------------------------------------------------------------- #

shape = 'quads' # or 'tris'
base_results_dirname = 'galewsky_'+shape
results_options = ['plain','recovered','vorticity']
plot_times = -1
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

for option in results_options:
    results_dirname = base_results_dirname + '_' + option

    print('Making field plots')
    text_pos = (-np.pi/2, 0.53*np.pi)

    for run_id in range(num_runs):
        text = option+r' $n$: %d' % (res_list[run_id])

        make_field_plots(results_dirname, run_id, 'galewsky',
                         'vorticity', plot_times, 'xy',
                         contours=np.arange(-1.6e-4,1.8e-4,step=2e-5),
                         spherical_centre=(-np.pi/2,np.pi/4),
                         plot_coords_1d=(np.linspace(-3*np.pi/2, np.pi/2, 201),
                                         np.linspace(1./9*np.pi/2, 8./9*np.pi/2, 101)),
                         ylims=[10*np.pi/180, 80*np.pi/180], yticklabels=[10,80],
                         xticklabels=[-270,90], xlabel=r'$\lambda \ / $ deg',
                         ylabel=r'$\vartheta \ / $ deg', figsize=(15,6),
                         cbar_label=r'$\zeta \ / $ s$^{-1}$',
                         cbar_labelpad=-60,
                         colour_scheme='RdBu_r',
                         remove_contour='middle',
                         extend_cmap=False, ylabelpad=-25,
                         cbar_ticks=[-1.6e-4, 1.6e-4],
                         cbar_format='%.1e',
                         time_method='days', text=text, text_pos=text_pos)
