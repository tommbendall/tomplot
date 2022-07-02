"""
Plots three Straka fields
"""

import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np
from tomplot import individual_field_contour_plot, extract_lfric_2D_data

# ---------------------------------------------------------------------------- #
# Things that can be altered and parameters for the test case
# ---------------------------------------------------------------------------- #

results_dirnames = ['pdc_demo_data/straka_D200_P400', 
                    'pdc_demo_data/straka_D200_P200',
                    'pdc_demo_data/straka_D200_P100']
plotname = 'straka_demo'
cbar_label = r"$\theta \ / $ K"
titles = [r'$\Delta x_{phys}=\frac{1}{4}\Delta x_{dyn}$',
          r'$\Delta x_{phys}=\Delta x_{dyn}$',
          r'$\Delta x_{phys}=4\Delta x_{dyn}$,']
colour_scheme = 'Blues_r'
extrusion_details = {'domain':'plane', 'extrusion':'linear',
                     'zmin':0.0, 'zmax':6400, 'topological_dimension':3}

field_name = 'theta'
field_min = 290.0
field_max = 301.0
xlabel = r'$x \ / $ km'
ylabels = [r'$z \ / $ km', None, None]
xlims = [0,16]
ylims = [0,5]

contours = np.linspace(field_min, field_max, 12)

# ---------------------------------------------------------------------------- #
# Things that are likely the same
# ---------------------------------------------------------------------------- #

plotdir = 'results/pdc_idealised_paper/figures'
slice_name = 'xz'
slice_idx = 0
time_idx = -1
init_time_idx = 0

# This is declared BEFORE figure and ax are initialised
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
font = {'size':24}
plt.rc('font',**font)

# ---------------------------------------------------------------------------- #
# Field plots
# ---------------------------------------------------------------------------- #

fig, axarray = plt.subplots(1,3,figsize=(18,5),sharey='row')

plotpath = f'{plotdir}/{plotname}.jpg'

for i, (ax, title, ylabel, results_dirname) \
    in enumerate(zip(axarray, titles, ylabels, results_dirnames)):

    data_file = Dataset('/data/users/tbendall/'+results_dirname+'/lfric_diag.nc','r')

    coords_X, coords_Y, field_data, data_metadata = \
        extract_lfric_2D_data(data_file, field_name, time_idx,
                              slice_name=slice_name, slice_idx=slice_idx,
                              extrusion_details=extrusion_details)

    data_file.close()

    time = data_metadata['time']
    coord_labels = data_metadata['coord_labels']
    coord_lims = data_metadata['coord_lims']
    coord_ticks = data_metadata['coord_ticks']
    slice_label = data_metadata['slice_label']

    coords_X = coords_X*0.001
    coords_Y = coords_Y*0.001

    cf = individual_field_contour_plot(coords_X, coords_Y, field_data,
                                       ax=ax, time=time, slice_name=slice,
                                       title=title, title_method=None,
                                       slice_idx=0, slice_label=slice_label,
                                       contours=contours, no_cbar=True,
                                       colour_scheme=colour_scheme,
                                       extend_cmap=False, titlepad=20,
                                       remove_contour=300.0,
                                       ylabelpad=-10, ylabel=ylabel,
                                       xlims=xlims, ylims=ylims, xlabel=xlabel,
                                       xticklabels=xlims, yticklabels=ylims)


# Move the subplots to the left to make space for colorbar
fig.subplots_adjust(right=0.9, wspace=0.125)
# Add colorbar in its own axis
cbar_ax = fig.add_axes([0.92, 0.11, 0.02, 0.77])
cb = fig.colorbar(cf, cax=cbar_ax, ticks=[field_min, field_max])
cb.set_label(cbar_label, labelpad=-30)

print(f'Plotting to {plotpath}')
fig.savefig(plotpath, bbox_inches='tight')
plt.close()

