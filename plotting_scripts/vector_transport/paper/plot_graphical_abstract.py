"""
Makes field plots for the shallow water Galewsky jet test case for the
graphic abstract.
"""

import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from tomplot import extract_2D_data, individual_field_contour_plot

# ---------------------------------------------------------------------------- #
# Things that can be altered and parameters for the test case
# ---------------------------------------------------------------------------- #

plotdir = '/data/users/tbendall/results/vector_transport_paper/figures'
results_options = ['plain','recovered','vorticity','vorticity_fancy_SUPG']
titles = ['Standard upwind scheme', 'Recovery into higher-order function space',
          'Mixed velocity-vorticity form without stabilisation',
          'Mixed velocity-vorticity form with SUPG stabilisation']
time_idx = -1
colour_scheme = 'RdBu_r'
field_name = 'vorticity'
spherical_centre = (-np.pi/2,np.pi/4)
slice = 'xy'
slice_idx = 0
run_id = 0

base_results_dirname = 'galewsky_quads'
plotname = 'fig_0_graphical_abstract'

# ---------------------------------------------------------------------------- #
# Field plots
# ---------------------------------------------------------------------------- #

results_dirnames = [f'vector_transport_paper/galewsky_quads_{option}' for option in results_options]

# This is declared BEFORE figure and ax are initialised
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
font = {'size':32}
plt.rc('font',**font)

fig, axarray = plt.subplots(4,1,figsize=(16,16),sharex='col')

plotpath = f'{plotdir}/{plotname}.jpg'

for i, (ax, dirname, title, xticklabels, xlabel, xticks) in \
    enumerate(zip(axarray, results_dirnames, titles, [None,None,None,[-270,90]],
                    [None,None,None,r'$\lambda \ / $ deg'], [None,None,None,[-3*np.pi/2,np.pi/2]])):
    # This code is all adapted from plot_control
    filename = f'/data/users/tbendall/results/{dirname}/nc_fields/field_output_{str(run_id)}.nc'
    data_file = Dataset(filename, 'r')

    # Extract data
    coords_X, coords_Y, field_data, data_metadata = \
        extract_2D_data(data_file, field_name, time_idx,
                        slice_name=slice, slice_idx=slice_idx,
                        central_lon=spherical_centre[0],
                        plot_coords_1d=(np.linspace(-3*np.pi/2, np.pi/2, 201),
                                        np.linspace(1./9*np.pi/2, 8./9*np.pi/2, 101)))

    data_file.close()

    time = data_metadata['time']
    coord_labels = data_metadata['coord_labels']
    coord_lims = data_metadata['coord_lims']
    coord_ticks = data_metadata['coord_ticks']
    slice_label = data_metadata['slice_label']


    yticks = [10*np.pi/180, 80*np.pi/180] if i == 0 else None
    yticklabels = [10,80] if i == 0 else ['','']
    ylabel = r'$\vartheta \ / $ deg' if i == 0 else None

    cf = individual_field_contour_plot(coords_X, coords_Y, field_data, ax=ax,
                                        title=title, title_method=None,
                                        titlepad=10, title_size=32,
                                        contours=np.arange(-1.8e-4,2.0e-4,step=2e-5),
                                        spherical_centre=spherical_centre,
                                        ylims=[10*np.pi/180, 80*np.pi/180], yticklabels=yticklabels,
                                        ylabel=ylabel, yticks=yticks,
                                        xticklabels=xticklabels, xlabel=xlabel,
                                        xticks=xticks,
                                        restricted_cmap='both',
                                        colour_levels_scaling=(1.35,1.25),
                                        no_cbar=True,colour_scheme=colour_scheme,
                                        remove_contour='middle',
                                        extend_cmap=False, ylabelpad=-25)

# Move the subplots to the left to make space for colorbar
# and make them slightly closer together
fig.subplots_adjust(right=0.9, hspace=0.2)
# Make collective colour bar
cbar_label=r'$\zeta \ / $ s$^{-1}$'
cbar_labelpad=-80
cbar_ticks=[-1.8e-4, 1.8e-4]
cbar_format='%.1e'
# Add colorbar in its own axis
cbar_ax = fig.add_axes([0.92, 0.11, 0.03, 0.77])
cb = fig.colorbar(cf, orientation='vertical', cax=cbar_ax,
                    format=cbar_format, ticks=cbar_ticks)
cb.set_label(cbar_label, labelpad=cbar_labelpad)

print(f'Plotting to {plotpath}')
fig.savefig(plotpath, bbox_inches='tight')
plt.close()

