"""
Makes convergence plots for the transport test cases
"""

from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
from tomplot import individual_convergence_plot

# ---------------------------------------------------------------------------- #
# Things that can be altered and parameters for the test case
# ---------------------------------------------------------------------------- #

error = 'L2_error'
variable = 'dx'
result_opts = ['conv_1_quads_cylinder','conv_2_quads_hooks']
field_labels = ['Benchmark','Recovered','Vorticity']
leg_xcentres = [0.45,0.5]
plotname = 'fig_2_convergence'
titles = ['Cylindrical domain','Spherical domain']
ylabels = [r'$\ln(||\mathbf{F}-\mathbf{F}_{true}||)$',None]

# ---------------------------------------------------------------------------- #
# Field plots
# ---------------------------------------------------------------------------- #

plotdir = '/data/users/tbendall/results/vector_transport_paper/figures'

# This is declared BEFORE figure and ax are initialised
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
font = {'size':24}
plt.rc('font',**font)

fig, axarray = plt.subplots(1,2,figsize=(16,8))

plotpath = f'{plotdir}/{plotname}.jpg'

for i, (ax, result_opt, leg_xcentre, title, ylabel) in \
    enumerate(zip(axarray, result_opts, leg_xcentres, titles, ylabels)):

# ---------------------------------------------------------------------------- #
# Get run ID and setups info
# ---------------------------------------------------------------------------- #

    results_dirname = f'vector_transport_paper/{result_opt}'
    data = Dataset(f'/data/users/tbendall/results/{results_dirname}/global_output.nc','r')
    run_ids = data['run_id'][:]
    num_setups = len(field_labels)
    data.close()

# ---------------------------------------------------------------------------- #
# Convergence plots
# ---------------------------------------------------------------------------- #

    field_names = ['F_'+str(i) for i in range(num_setups)]

    individual_convergence_plot(results_dirname, variable, field_names, run_ids,
                                error, ax=ax, field_labels=field_labels, label_style='gradient_plain',
                                legend_bbox=(leg_xcentre,1.12),
                                legend_ncol=3, titlepad=55, title=title,
                                ylabel=ylabel, leg_col_spacing=0.1, leg_fontsize=18)

print(f'Plotting to {plotpath}')
fig.savefig(plotpath, bbox_inches='tight')
plt.close()
