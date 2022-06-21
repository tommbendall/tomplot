"""
Makes convergence and consistency plots for the transport-only test
"""

from netCDF4 import Dataset
import matplotlib.pyplot as plt
from tomplot import individual_convergence_plot, individual_time_series_plot

# ---------------------------------------------------------------------------- #
# Things that can be altered and parameters for the test case
# ---------------------------------------------------------------------------- #

plotname = 'fig_X_moist_gw_200_convergence'
ylabels = [r'$\ln(||\theta_e-\theta_e^{true}||)$','']

# Things for convergence plot
results_dirname = 'moist_gw_200_convergence'
all_field_labels = [[r'$\Delta x_D = \Delta x_P$', r'$\Delta x_P = 0.5$ km', r'$\Delta x_D = \Delta x_P/2$'],
                    [r'$\Delta x_P = \Delta x_D$', r'$\Delta x_D = 0.5$ km', r'$\Delta x_P = \Delta x_D/2$']]
field_name = ['theta_e']*3
conv_variables = ['dyn_dx', 'phys_dx']
xlabels = [r'$\Delta x_D$', r'$\Delta x_P$']
error = 'L2_error'
# Magic numbers from order that data was input
all_run_ids = [[[0,3,5,10,14,17,25,29], [2,7,12,15,24,27], [1,6,11,15,26,28]],
               [[0,3,5,10,14,17,25,29], [18,19,20,21,22,23], [4,8,9,13,16,21]]]

# ---------------------------------------------------------------------------- #
#
# ---------------------------------------------------------------------------- #

plotdir = 'results/pdc_idealised_paper/figures'

# This is declared BEFORE figure and ax are initialised
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
font = {'size':24}
plt.rc('font',**font)

fig, axarray = plt.subplots(1,2,figsize=(16,8),sharey='row')

plotpath = f'{plotdir}/{plotname}.jpg'

for i, (ax, xlabel, ylabel, run_ids, conv_variable, field_labels) \
    in enumerate(zip(axarray, xlabels, ylabels, all_run_ids, conv_variables, all_field_labels)):

# ---------------------------------------------------------------------------- #
# Convergence plots
# ---------------------------------------------------------------------------- #

    individual_convergence_plot(results_dirname, conv_variable, field_name,
                                run_ids, error, ax=ax,
                                field_labels=field_labels,
                                colours=['blue','red','black'],
                                label_style='gradient_plain',
                                legend_bbox=(0.5,1.12),
                                legend_ncol=3,
                                xlabel=xlabel, ylabel=ylabel,
                                leg_col_spacing=0.1, leg_fontsize=14)

# ---------------------------------------------------------------------------- #

print(f'Plotting to {plotpath}')
fig.subplots_adjust(wspace=0.28)
fig.savefig(plotpath, bbox_inches='tight')
plt.close()
