"""
Makes time series plots for the shallow water Galewsky jet test case,
"""

import numpy as np
import matplotlib.pyplot as plt
from tomplot import individual_time_series_plot

# ---------------------------------------------------------------------------- #
# Things that can be altered and parameters for the test case
# ---------------------------------------------------------------------------- #

plotdir = 'results/vector_transport_paper/figures'
results_options = ['plain','recovered','vorticity','vorticity_fancy_SUPG']
field_labels = ['Benchmark', 'Recovered', 'Vorticity w/o SUPG','Vorticity with SUPG']
diagnostic_fields = ['energy', 'potential_enstrophy']
diagnostic = 'total'
ylabels = ['(E(t)-E(0))/E(0)', '(Q(t)-Q(0))/Q(0)']
titles = ['Energy', 'Enstrophy']
linestyles = ['-.', '-', '--', ':']
no_legends = [False, True]
run_id = 0

base_results_dirname = 'galewsky_quads'
plotname = 'fig_7_galewsky_time_series'

# ---------------------------------------------------------------------------- #
# Time Series plots
# ---------------------------------------------------------------------------- #

results_dirnames = [f'vector_transport_paper/galewsky_quads_{option}' for option in results_options]

# This is declared BEFORE figure and ax are initialised
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
font = {'size':24}
plt.rc('font',**font)

fig, axarray = plt.subplots(1,2,figsize=(16,8))

plotpath = f'{plotdir}/{plotname}.jpg'

for i, (ax, title, ylabel, diagnostic_field, no_legend) \
  in enumerate(zip(axarray, titles, ylabels, diagnostic_fields, no_legends)):

    individual_time_series_plot(results_dirnames, diagnostic_field, run_id,
                                diagnostic, ax=ax, field_labels=field_labels,
                                format='jpg', ylabel=ylabel, title=title,
                                time_units='days', normalise=True,
                                linestyles=linestyles, linewidth=2.5,
                                legend_ncol=4, legend_bbox=(1.1, 1.15),
                                titlepad=70, no_legend=no_legend)

print(f'Plotting to {plotpath}')
fig.subplots_adjust(wspace=0.28)
fig.savefig(plotpath, bbox_inches='tight')
plt.close()
