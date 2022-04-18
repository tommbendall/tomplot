"""
Makes convergence and consistency plots for the transport-only test
"""

from netCDF4 import Dataset
import matplotlib.pyplot as plt
from tomplot import individual_convergence_plot, individual_time_series_plot

# ---------------------------------------------------------------------------- #
# Things that can be altered and parameters for the test case
# ---------------------------------------------------------------------------- #

base_dirname = 'pdc_idealised_paper'
plotname = 'fig_X_convergence_consistency'
titles = ['Convergence','Consistency']
ylabels = [r'$\ln(||q-q^{true}||)$',r'$\left(||q(t)||-||q_0||\right)/||q_0||$']

# Things for convergence plot
conv_dirnames = ['stretchy_sphere_conv_adv_all','stretchy_sphere_conv_all']
conv_field_labels = ['Adv x 2', 'Adv max 96', 'Cons x 2', 'Cons max 96']
conv_field_names = ['phys-rho', 'phys-rho']
conv_variable = 'coarse_dx'
error = 'Rel-L2-error'
conv_run_ids = [[0,1,2,3,4], [4,5,6,7,8]]

# Things for consistency plot
cons_dirname = 'stretchy_sphere_cons'
cons_field_names = ['dyn-rho', 'phys-rho']
cons_diagnostic = 'L2-now'
cons_labels = [r'$\rho_d$', r'$q$']
cons_run_ids = [0]
cons_xticks = [0, 2000]
cons_yticks = [-0.1, 0.2]

# ---------------------------------------------------------------------------- #
#
# ---------------------------------------------------------------------------- #

plotdir = 'results/pdc_idealised_paper/figures'

# This is declared BEFORE figure and ax are initialised
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
font = {'size':24}
plt.rc('font',**font)

fig, axarray = plt.subplots(1,2,figsize=(16,8))

plotpath = f'{plotdir}/{plotname}.jpg'

for i, (ax, title, ylabel) in enumerate(zip(axarray, titles, ylabels)):

# ---------------------------------------------------------------------------- #
# Convergence plots
# ---------------------------------------------------------------------------- #

    if i == 0:

        results_dirnames = [f'{base_dirname}/{dirname}' for dirname in conv_dirnames]
        individual_convergence_plot(results_dirnames, conv_variable, conv_field_names,
                                    conv_run_ids, error, ax=ax,
                                    field_labels=conv_field_labels,
                                    colours=['blue','red','purple','black'],
                                    label_style='gradient_plain',
                                    legend_bbox=(0.5,1.12),
                                    legend_ncol=2, titlepad=55, title=title,
                                    xlabel=r'$\ln(\Delta x)$',
                                    ylabel=ylabel, leg_col_spacing=0.1, leg_fontsize=18)

# ---------------------------------------------------------------------------- #
# Consistency plot
# ---------------------------------------------------------------------------- #

    else:
        results_dirname = f'{base_dirname}/{cons_dirname}'
        individual_time_series_plot(results_dirname, cons_field_names, cons_run_ids,
                                    cons_diagnostic, ax=ax, field_labels=cons_labels,
                                    label_style='range_plain',  normalise=True,
                                    legend_bbox=(0.5,1.12), colours=['black','black'],
                                    linestyles=['--','-'],
                                    legend_ncol=2, titlepad=55, title=title,
                                    ylabel=ylabel, leg_col_spacing=0.5, leg_fontsize=18,
                                    xlims=cons_xticks, ylims=cons_yticks)

print(f'Plotting to {plotpath}')
fig.subplots_adjust(wspace=0.28)
fig.savefig(plotpath, bbox_inches='tight')
plt.close()
