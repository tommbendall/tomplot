"""
Makes convergence and consistency plots for the transport-only test
"""

from netCDF4 import Dataset
import matplotlib.pyplot as plt
from tomplot import individual_convergence_plot, individual_time_series_plot

# ---------------------------------------------------------------------------- #
# Things that can be altered and parameters for the test case
# ---------------------------------------------------------------------------- #

plotname = 'fig_X_moist_gw_new_convergence'
ylabel = r'$\ln(||\theta_e-\theta_e^{true}||)$'

# Things for convergence plot
results_dirname = 'moist_gw_new_convergence'
field_labels = [r'$\Delta x_P = \Delta x_D$', r'$\Delta x_P > \Delta x_D$', r'$\Delta x_P < \Delta x_D$']
field_name = ['theta_e']*3
conv_variable = 'phys_dx'
xlabel = r'$\Delta x_P$'
error = 'L2_error'
# Magic numbers from order that data was input
run_ids = [[0,1,2,3,5,7], [6,13,14], [4,10,11,12]]

# ---------------------------------------------------------------------------- #
#
# ---------------------------------------------------------------------------- #

plotdir = 'results/pdc_idealised_paper/figures'

# This is declared BEFORE figure and ax are initialised
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
font = {'size':24}
plt.rc('font',**font)

# ---------------------------------------------------------------------------- #
# Convergence plots
# ---------------------------------------------------------------------------- #

individual_convergence_plot(results_dirname, conv_variable, field_name,
                            run_ids, error,
                            field_labels=field_labels,
                            colours=['black','red','blue'],
                            label_style='plain', best_fit=True,
                            best_fit_deg=2, legend_bbox=(0.5,1.15),
                            linestyles=['--']*3,
                            legend_ncol=3, markers=['^','o','s'],
                            xlabel=xlabel, ylabel=ylabel, markersize=12,
                            format='jpg', plotdir=plotdir, testname=plotname)

# ---------------------------------------------------------------------------- #

