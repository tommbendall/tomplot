"""
A tomplot example, for making a convergence plot.
"""

import matplotlib.pyplot as plt
from os.path import abspath, dirname
from tomplot import (set_tomplot_style, plot_convergence,
                     only_minmax_ticklabels, tomplot_legend_ax,
                     tomplot_legend_fig)
# ---------------------------------------------------------------------------- #
# Some dummy data
# ---------------------------------------------------------------------------- #
# This should be replaced with your own data if you start from this example!
# This is not actually meaningful data, just some made up numbers
dx_values = [0.01, 0.02, 0.04, 0.08]
data_1 = [1.40067551, 3.82953108, 14.88772346, 62.4848737]
data_2 = [0.809123, 2.4431, 3.96, 8.4]
data_3 = [25.132098, 27.24987, 32.60984, 50.00146]

# ---------------------------------------------------------------------------- #
# Directory for results and plots
# ---------------------------------------------------------------------------- #
# When copying this example these should not be relative to this file
plot_dir = f'{abspath(dirname(__file__))}/../tests/tmp_figures'
plot_name = f'{plot_dir}/example_convergence.png'
# There might be a results_dir here too

# ---------------------------------------------------------------------------- #
# Things that should be altered based on the plot
# ---------------------------------------------------------------------------- #
# Three data sets: three colours, markers and labels
all_error_data = [data_1, data_2, data_3]
colours = ['red', 'blue', 'purple']
markers = ['s', 'o', '^']
base_labels = ['Scheme 1', 'Scheme 2', 'Scheme 3']
labels = [base_label+', gradient =' for base_label in base_labels]
# Convergence plot options that we generally want to use
log_by = 'data'
xlabel = r"$\log(\Delta x)$"
ylabel = r"$\log(\Delta q)$"
legend_loc = 'upper center'
# ---------------------------------------------------------------------------- #
# Things that are likely the same for all plots
# ---------------------------------------------------------------------------- #
set_tomplot_style()
fig, ax = plt.subplots(1, 1, figsize=(5, 5))
# ---------------------------------------------------------------------------- #
# Data extraction
# ---------------------------------------------------------------------------- #
# You might put something here if we were actually extracting from a file

# ---------------------------------------------------------------------------- #
# Plot data
# ---------------------------------------------------------------------------- #
for error_data, colour, marker, label in \
        zip(all_error_data, colours, markers, labels):
    plot_convergence(ax, dx_values, error_data, label=label,
                     color=colour, marker=marker, log_by=log_by)

only_minmax_ticklabels(ax)
# tomplot_legend_fig(fig)
# plt.legend(loc=legend_loc)
ax.set_xlabel(xlabel)
ax.set_ylabel(ylabel)

tomplot_legend_ax(ax, location='bottom')

plt.grid()
# ---------------------------------------------------------------------------- #
# Save figure
# ---------------------------------------------------------------------------- #
print(f'Saving figure to {plot_name}')
fig.savefig(plot_name, bbox_inches='tight')
plt.close()


