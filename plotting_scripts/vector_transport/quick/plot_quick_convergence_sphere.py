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
base_results_dir = 'conv_2_quads'
result_opts = ['dt0pt05','dt0pt2','dt0pt5','update_vort']
field_labels = ['Benchmark','Recovered','Vorticity','SUPG','Fancy']
base_plotname = 'quick_conv_sph'
ylabel = r'$\ln(||\mathbf{F}-\mathbf{F}_{true}||)$'

# ---------------------------------------------------------------------------- #
# Field plots
# ---------------------------------------------------------------------------- #

plotdir = 'results/vector_transport_paper/quick_figures'

# This is declared BEFORE figure and ax are initialised
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
font = {'size':24}
plt.rc('font',**font)


for i, result_opt in enumerate(result_opts):

    fig, ax = plt.subplots(1,1,figsize=(8,8))

    plotpath = f'{plotdir}/{base_plotname}_{result_opt}.jpg'

# ---------------------------------------------------------------------------- #
# Get run ID and setups info
# ---------------------------------------------------------------------------- #

    results_dirname = f'vector_transport_paper/{base_results_dir}_{result_opt}'
    data = Dataset(f'results/{results_dirname}/global_output.nc','r')
    run_ids = data['run_id'][3:]
    num_setups = len(field_labels)
    data.close()

# ---------------------------------------------------------------------------- #
# Convergence plots
# ---------------------------------------------------------------------------- #

    field_names = ['F_'+str(i) for i in range(num_setups)]

    individual_convergence_plot(results_dirname, variable, field_names, run_ids,
                                error, field_labels=field_labels, ax=ax,
                                label_style='gradient_plain',
                                legend_bbox=(0.5,1.12),
                                legend_ncol=5, ylabel=ylabel, leg_col_spacing=0.1,
                                leg_fontsize=16)

    print(f'Plotting to {plotpath}')
    fig.savefig(plotpath, bbox_inches='tight')
    plt.close()
