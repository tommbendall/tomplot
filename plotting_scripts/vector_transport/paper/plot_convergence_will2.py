"""
Makes convergence plots for the Williamson 2 test case
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
field_names = ['u','D']
result_opts = ['conv_3_quads_sw2_up','conv_3_quads_sw2_rec','conv_3_quads_sw2_vort']
field_labels = ['Benchmark','Recovered','Vorticity']
plotname = 'fig_5_will2_convergence'
titles = ['Velocity','Height']
ylabels = [r'$\ln\left[||\textbf{u}-\textbf{u}_{true}||_{L_2}/||\textbf{u}_{true}||_{L_2}\right]$',
           r'$\ln\left[||h-h_{true}||_{L_2}/||h_{true}||_{L_2}\right]$']
leg_xcentres = [0.45,0.5]

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

for i, (ax, field_name, title, ylabel, leg_xcentre) in \
    enumerate(zip(axarray, field_names, titles, ylabels, leg_xcentres)):

# ---------------------------------------------------------------------------- #
# Get run ID and setups info
# ---------------------------------------------------------------------------- #

    results_dirnames = [f'vector_transport_paper/{result_opt}' for result_opt in result_opts]
    data = Dataset(f'/data/users/tbendall/results/{results_dirnames[0]}/global_output.nc','r')
    run_ids = data['run_id'][:]
    num_setups = len(field_labels)
    data.close()

# ---------------------------------------------------------------------------- #
# Convergence plots
# ---------------------------------------------------------------------------- #

    individual_convergence_plot(results_dirnames, variable, field_name, run_ids,
                                error, ax=ax, field_labels=field_labels, label_style='gradient_plain',
                                legend_bbox=(leg_xcentre,1.12),
                                legend_ncol=3, titlepad=55, title=title,
                                ylabel=ylabel, leg_col_spacing=0.1, leg_fontsize=18)

print(f'Plotting to {plotpath}')
fig.savefig(plotpath, bbox_inches='tight')
plt.close()
