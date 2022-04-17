"""
Makes plots for the shallow water Galewsky jet test case,
"""

import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from tomplot import extract_2D_data, individual_field_contour_plot, individual_time_series_plot

# ---------------------------------------------------------------------------- #
# Things that can be altered and parameters for the test case
# ---------------------------------------------------------------------------- #

plotdir = 'results/vector_transport_paper/figures'
results_options = ['plain','recovered','vorticity','vorticity_fancy_SUPG']
titles = ['Benchmark,', 'Recovered,', 'Vorticity,','SUPG']
diagnostic_fields = ['energy', 'potential_enstrophy']
time_idx = -1
colour_scheme = 'RdBu_r'
field_name = 'vorticity'
spherical_centre = (-np.pi/2,np.pi/4)
slice = 'xy'
slice_idx = 0
run_id = 0
diagnostic = 'total'
sw_labels = {('total', 'energy'): 'Normalised energy',
             ('total', 'potential_enstrophy'): 'Normalised potential enstrophy'}

base_results_dirname = 'galewsky_quads'
plotname = 'fig_6_galewsky'

# ---------------------------------------------------------------------------- #
# Field plots
# ---------------------------------------------------------------------------- #

results_dirnames = [f'results/vector_transport_paper/galewsky_quads_{option}' for option in results_options]

# This is declared BEFORE figure and ax are initialised
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
font = {'size':24}
plt.rc('font',**font)

fig, axarray = plt.subplots(4,1,figsize=(16,20),sharex='col')

plotpath = f'{plotdir}/{plotname}.jpg'

for i, (ax, dirname, title, xticklabels, xlabel, xticks) in \
    enumerate(zip(axarray, results_dirnames, titles, [None,None,None,[-270,90]],
                    [None,None,None,r'$\lambda \ / $ deg'], [None,None,None,[-3*np.pi/2,np.pi/2]])):
    # This code is all adapted from plot_control
    filename = f'{dirname}/nc_fields/field_output_{str(run_id)}.nc'
    data_file = Dataset(filename, 'r')

    # Extract data
    coords_X, coords_Y, field_data, time, \
    coord_labels, coord_lims, coord_ticks,  \
    slice_label = extract_2D_data(data_file, field_name, time_idx,
                                    slice_name=slice, slice_idx=slice_idx,
                                    central_lon=spherical_centre[0],
                                    plot_coords_1d=(np.linspace(-3*np.pi/2, np.pi/2, 201),
                                                    np.linspace(1./9*np.pi/2, 8./9*np.pi/2, 101)))

    data_file.close()

    cf = individual_field_contour_plot(coords_X, coords_Y, field_data, ax=ax,
                                        title=title, title_method='minmax',
                                        titlepad=10,
                                        contours=np.arange(-1.8e-4,2.0e-4,step=2e-5),
                                        spherical_centre=spherical_centre,
                                        ylims=[10*np.pi/180, 80*np.pi/180], yticklabels=[10,80],
                                        ylabel=r'$\vartheta \ / $ deg',
                                        xticklabels=xticklabels, xlabel=xlabel,
                                        xticks=xticks,
                                        restricted_cmap='both',
                                        colour_levels_scaling=(1.35,1.25),
                                        no_cbar=True,colour_scheme=colour_scheme,
                                        remove_contour='middle',
                                        extend_cmap=False, ylabelpad=-25)

# Move the subplots to the left to make space for colorbar
# and make them slightly closer together
fig.subplots_adjust(bottom=0.1, hspace=0.14)
# Make collective colour bar
cbar_label=r'$\zeta \ / $ s$^{-1}$'
cbar_labelpad=-25
cbar_ticks=[-1.8e-4, 1.8e-4]
cbar_format='%.1e'
# Add colorbar in its own axis
cbar_ax = fig.add_axes([0.16, 0.04, 0.7, 0.02])
cb = fig.colorbar(cf, orientation='horizontal', cax=cbar_ax,
                    format=cbar_format, ticks=cbar_ticks)
cb.set_label(cbar_label, labelpad=cbar_labelpad)

print(f'Plotting to {plotpath}')
fig.savefig(plotpath, bbox_inches='tight')
plt.close()

for diagnostic_field in diagnostic_fields:

    individual_time_series_plot(results_dirnames, diagnostic_field, run_id,
                                diagnostic, testname=f'fig_7_gal_{diagnostic_field}',
                                plotdir=plotdir, field_labels=results_options,
                                format='jpg', ylabel=sw_labels[(diagnostic,diagnostic_field)],
                                time_units='days', normalise=True)


