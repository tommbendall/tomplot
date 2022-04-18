"""
Makes convergence and consistency plots for the transport-only test
"""

from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
from tomplot import individual_convergence_plot, individual_time_series_plot

# ---------------------------------------------------------------------------- #
# Things that can be altered and parameters for the test case
# ---------------------------------------------------------------------------- #

base_dirname = 'consistent_moisture_paper'
plotname = 'fig_8_conservation'
titles = ['Transport Test','Rising Bubble']
ylabels = [r'$\left(M(t) - M_0\right)/M_0$','']
all_dirnames = [['stretchy-conv-adv-120', 'stretchy-conv-cons-120'],
                ['bryan_fritsch-adv', 'bryan_fritsch-consist']]
field_labels = ['Advective', 'Consistent']
field_names = ['mr', 'total']
measures = ['Mass', 'After_timestep']
all_run_ids = [[4],[0]]

# ---------------------------------------------------------------------------- #
# Field plots
# ---------------------------------------------------------------------------- #

plotdir = 'results/consistent_moisture_paper/figures'

# This is declared BEFORE figure and ax are initialised
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
font = {'size':24}
plt.rc('font',**font)

fig, axarray = plt.subplots(1,2,figsize=(16,8))

plotpath = f'{plotdir}/{plotname}.jpg'

for i, (dirnames, ax, title, ylabel, field_name, measure, run_ids) \
  in enumerate(zip(all_dirnames, axarray, titles, ylabels, field_names, measures, all_run_ids)):

# ---------------------------------------------------------------------------- #
# Consistency plots
# ---------------------------------------------------------------------------- #

    results_dirnames = [f'{base_dirname}/{dirname}' for dirname in dirnames]
    individual_time_series_plot(results_dirnames, field_name, run_ids,
                                measure, ax=ax, field_labels=field_labels,
                                label_style='range_plain',
                                legend_bbox=(0.5,1.15), colours=['black','black'],
                                linestyles=['--','-'], normalise=True,
                                legend_ncol=2, titlepad=70, title=title,
                                ylabel=ylabel, leg_col_spacing=0.5, leg_fontsize=18)

print(f'Plotting to {plotpath}')
fig.savefig(plotpath, bbox_inches='tight')
plt.close()
