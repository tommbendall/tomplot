"""
Makes convergence plots for the velocity for Williamson 2 test case
"""

from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
from tomplot import individual_convergence_plot

# ---------------------------------------------------------------------------- #
# Things that can be altered and parameters for the test case
# ---------------------------------------------------------------------------- #

error = 'L2_error_normalised'
variable = 'dx'
field_name='u'
results_dirname_sets = [['conv_3_quads_sw2_upwind','conv_3_quads_sw2_recovered','conv_3_quads_sw2_vorticity'],
                        ['conv_3_tris_sw2_upwind','conv_3_tris_sw2_recovered','conv_3_tris_sw2_vorticity']]
field_label_sets = [['Benchmark','Recovered','Vorticity']]*2
plotname = 'fig_6_will2_u_convergence'
titles = ['Quadrilateral cells','Triangular cells']
ylabels = [r'$\ln\left[||\textbf{u}-\textbf{u}_{true}||_{L_2}/||\textbf{u}_{true}||_{L_2}\right]$',None]

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

for i, (ax, results_dirnames, field_labels, title, ylabel) in \
    enumerate(zip(axarray, results_dirname_sets, field_label_sets, titles, ylabels)):

# ---------------------------------------------------------------------------- #
# Get run ID and setups info
# ---------------------------------------------------------------------------- #

    data = Dataset('results/'+results_dirnames[0]+'/global_output.nc','r')
    run_ids = data['run_id'][:]
    num_setups = len(field_labels)
    data.close()

# ---------------------------------------------------------------------------- #
# Convergence plots
# ---------------------------------------------------------------------------- #

    field_names = ['F_'+str(i) for i in range(num_setups)]

    individual_convergence_plot(results_dirnames, variable, field_name, run_ids,
                                error, ax=ax, field_labels=field_labels, label_style='gradient_plain',
                                legend_bbox=(0.5,1.12),
                                legend_ncol=3, titlepad=55, title=title,
                                ylabel=ylabel, leg_col_spacing=0.1, leg_fontsize=18)

print(f'Plotting to {plotpath}')
fig.savefig(plotpath, bbox_inches='tight')
plt.close()
