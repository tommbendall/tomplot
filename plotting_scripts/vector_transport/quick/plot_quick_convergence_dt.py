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
field_labels = ['0.005','0.01','0.02','0.05', '0.002']
base_plotname = 'quick_conv_cyl_dt'
ylabel = r'$\ln(||\mathbf{F}-\mathbf{F}_{true}||)$'

# ---------------------------------------------------------------------------- #
# Field plots
# ---------------------------------------------------------------------------- #

plotdir = '/data/users/tbendall/results/vector_transport_paper/quick_figures'

# This is declared BEFORE figure and ax are initialised
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
font = {'size':24}
plt.rc('font',**font)

all_run_ids = [[0, 1, 2, 3, 4, 5, 6, 26, 27, 28],
               [7, 8, 9, 10, 11, 12, 13],
               [14, 15, 16, 17, 18, 19, 20],
               [21, 22, 23, 24, 25],
               [29, 30, 31, 32, 33, 34, 35, 36, 37, 38]]

for field in ['F_rec', 'F_vort']:
    
    fig, ax = plt.subplots(1,1,figsize=(8,8))

    plotpath = f'{plotdir}/{base_plotname}_{field}.jpg'
    results_dirname = '/data/users/tbendall/results/vector_transport_paper/conv_1_dt_all'

    field_names = [field]*len(all_run_ids)

    individual_convergence_plot(results_dirname, variable, field_names, all_run_ids,
                                error, field_labels=field_labels,
                                label_style='gradient_plain', ax=ax,
                                legend_bbox=(0.5,1.12),
                                legend_ncol=len(field_labels)+1, ylabel=ylabel, leg_col_spacing=0.1,
                                leg_fontsize=16, best_fit=True,
                                comparison_lines=[2.0])

    print(f'Plotting to {plotpath}')
    fig.savefig(plotpath, bbox_inches='tight')
    plt.close()
