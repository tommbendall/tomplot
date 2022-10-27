"""
Makes convergence and consistency plots for the transport-only test
"""

from netCDF4 import Dataset
import matplotlib.pyplot as plt
from tomplot import individual_convergence_plot, individual_time_series_plot

# ---------------------------------------------------------------------------- #
# Things that can be altered and parameters for the test case
# ---------------------------------------------------------------------------- #

plotname = 'fig_X_moist_gw_convergence'
ylabel = r'$\ln(||\theta_e-\theta_e^{true}||)$'

# Things for convergence plot
results_dirname = '/data/users/tbendall/results/moist_gw_all_convergence/global_output.nc'
field_labels = [None, None, None, None, None,
                r'$\Delta x_P = \Delta x_D$',
                r'$\Delta x_P > \Delta x_D$',
                r'$\Delta x_P < \Delta x_D$']
field_name = ['theta_e']*len(field_labels)
conv_variable = 'phys_dx'
xlabel = r'$\ln(\Delta x_P \ / $ m$ )$'
error = 'L2_error'
# Magic numbers from order that data was input
run_ids = [ # All runs with a fixed nx_D (no markers, just lines)
           [9,2,10,11,12,13],         # nx_D = 100
           [14,15,5,16,17,18],        # nx_D = 200
           [19,20,21,22,6,23,24,25],  # nx_D = 300
           [26,27,28,29,7,30,31],     # nx_D = 400
           [32,33,34,35,8],           # nx_D = 600
           # Markers for all runs, no lines
           [0,1,2,3,4,5,6,7,8],       # nx_D == nx_P
           [9,14,15,19,20,21,22,      # nx_D > nx_P
            26,27,28,29,32,33,34,35],
           [10,11,12,13,16,17,18,23,  # nx_D < nx_P
            24,25,30,31]]

linestyles = ['--', '--', '--', '--', '--', '', '', '']
colours = ['black', 'black', 'black', 'black', 'black',
           'black', 'blue', 'red']
markers = ['', '', '', '', '',
           '^', 's', 'o']

# ---------------------------------------------------------------------------- #
#
# ---------------------------------------------------------------------------- #

plotdir = '/data/users/tbendall/results/pdc_idealised_paper/figures'

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
                            colours=colours, best_fit=False,
                            label_style='plain',
                            legend_bbox=(0.5,1.15),
                            linestyles=linestyles,
                            legend_ncol=3, markers=markers,
                            xlabel=xlabel, ylabel=ylabel, markersize=12,
                            format='jpg', plotdir=plotdir, testname=plotname,
                            override_dirname=True)

# ---------------------------------------------------------------------------- #

