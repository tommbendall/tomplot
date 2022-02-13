"""
Makes convergence plots for the spherical test case
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
results_dirnames = ['conv_2_quads_hooks','conv_2_tris_hooks']
field_label_sets = [['Benchmark','Recovered','Vorticity'],['Benchmark','Recovered']]
leg_xcentres = [0.45,0.5]
plotname = 'fig_5_sphere_convergence'
titles = ['Quadrilateral cells','Triangular cells']
ylabels = [r'$\ln(||\mathbf{F}-\mathbf{F}_{true}||)$',None]

# ---------------------------------------------------------------------------- #
# Field plots
# ---------------------------------------------------------------------------- #

plotdir = 'results/vector_transport_paper'

# This is declared BEFORE figure and ax are initialised
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
font = {'size':24}
plt.rc('font',**font)

fig, axarray = plt.subplots(1,2,figsize=(16,8),sharey='row')

plotpath = f'{plotdir}/{plotname}.jpg'

for i, (ax, results_dirname, field_labels, leg_xcentre, title, ylabel) in \
    enumerate(zip(axarray, results_dirnames, field_label_sets, leg_xcentres, titles, ylabels)):

# ---------------------------------------------------------------------------- #
# Get run ID and setups info
# ---------------------------------------------------------------------------- #

    data = Dataset('results/'+results_dirname+'/global_output.nc','r')
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
