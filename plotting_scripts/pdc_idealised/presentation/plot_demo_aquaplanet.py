"""
Plots an animation of aquaplanets evolving
"""

import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np
from tomplot import individual_field_contour_plot, extract_lfric_2D_data

# ---------------------------------------------------------------------------- #
# Things that can be altered and parameters for the test case
# ---------------------------------------------------------------------------- #

init_case = False
nice_radiation = True
results_dirnames = ['pdc_demo_data/aquaplanet_D48_P96', 
                    'pdc_demo_data/aquaplanet_D48_P48',
                    'pdc_demo_data/aquaplanet_D48_P24']
plotname = 'aquaplanet_demo'
cbar_label = r"$\theta \ / $ K"
titles = [r'$\Delta x_{phys}=\frac{1}{2}\Delta x_{dyn}$',
          r'$\Delta x_{phys}=\Delta x_{dyn}$',
          r'$\Delta x_{phys}=2\Delta x_{dyn}$,']
colour_scheme = 'Blues_r'
extrusion_details = {'domain':'sphere', 'extrusion':'linear',
                     'zmin':0.0, 'zmax':40000, 'topological_dimension':3}

field_names = ['m_v', 'm_cl','sw_down_surf']
colour_schemes = ['Blues_r', 'Greys', 'Greys_r']
transparencies = [1.0, 0.6, 0.3]
slice_idxs = [1,8,0]
all_contours = [np.linspace(0.0, 0.025, 41),
                np.linspace(5e-6, 5e-4, 40),
                np.linspace(0.0, 1200.0, 13)]

# ---------------------------------------------------------------------------- #
# Things that are likely the same
# ---------------------------------------------------------------------------- #

plotdir = '/data/users/tbendall/results/pdc_idealised_paper/figures'
slice_name = 'xy'
time_gap = 1.5
spherical_centre = (np.pi,0.0)

# This is declared BEFORE figure and ax are initialised
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
font = {'size':24}
plt.rc('font',**font)

import cartopy.crs as ccrs
crs = ccrs.PlateCarree()

# Get time indices
data_file = Dataset(f'/data/users/tbendall/{results_dirnames[0]}/lfric_diag.nc','r')
if init_case:
    time_idxs = range(0,len(data_file['time'][:])+1,1)
else:
    time_idxs = range(0,len(data_file['time'][:]),1)
data_file.close()

# ---------------------------------------------------------------------------- #
# Field plots
# ---------------------------------------------------------------------------- #

for time_idx in time_idxs:

    plotpath = f'{plotdir}/{plotname}_{time_idx:02d}.jpg'
    fig = plt.figure(figsize=(15,5))

    if init_case and time_idx > 0:
        this_time_idx = time_idx - 1
    else:
        this_time_idx = time_idx

    for i, (title, results_dirname) \
        in enumerate(zip(titles, results_dirnames)):

        if init_case and time_idx == 0:
            data_name = 'lfric_initial.nc'
        else:
            data_name = 'lfric_diag.nc'
        data_file = Dataset(f'/data/users/tbendall/{results_dirname}/{data_name}','r')

        lon_centre, lat_centre = spherical_centre[0]*180.0/np.pi, spherical_centre[1]*180.0/np.pi
        ax = fig.add_subplot(1, 3, 1+i,
                             projection=ccrs.Orthographic(lon_centre, lat_centre))


        for j, (field_name, contours, slice_idx, transparency, colour_scheme) in \
            enumerate(zip(field_names, all_contours, slice_idxs, transparencies, colour_schemes)):

            # No radiation data in initial file
            if time_idx == 0 and j == 2 and init_case:
                pass
            else:
                if nice_radiation and time_idx in [5,13] and j == 2:
                    coords_X, coords_Y, prev_field_data, data_metadata = \
                        extract_lfric_2D_data(data_file, field_name, (this_time_idx-1),
                                              slice_name=slice_name, slice_idx=slice_idx,
                                              extrusion_details=extrusion_details,
                                              central_lon=spherical_centre[0],
                                              num_points=100)
                    coords_X, coords_Y, next_field_data, data_metadata = \
                        extract_lfric_2D_data(data_file, field_name, (this_time_idx+1),
                                              slice_name=slice_name, slice_idx=slice_idx,
                                              extrusion_details=extrusion_details,
                                              central_lon=spherical_centre[0],
                                              num_points=100)

                    field_data = 0.5*(prev_field_data + next_field_data)

                else:
                    coords_X, coords_Y, field_data, data_metadata = \
                        extract_lfric_2D_data(data_file, field_name, this_time_idx,
                                              slice_name=slice_name, slice_idx=slice_idx,
                                              extrusion_details=extrusion_details,
                                              central_lon=spherical_centre[0],
                                              num_points=100)

                time = data_metadata['time']
                coord_labels = data_metadata['coord_labels']
                coord_lims = data_metadata['coord_lims']
                coord_ticks = data_metadata['coord_ticks']
                slice_label = data_metadata['slice_label']


                cf = individual_field_contour_plot(coords_X, coords_Y, field_data,
                                                   ax=ax, title_method=None,
                                                   contours=contours, no_cbar=True,
                                                   transparency=transparency,
                                                   colour_scheme=colour_scheme,
                                                   contour_lines=False,
                                                   extend_cmap=False,
                                                   projection='orthographic',
                                                   spherical_centre=spherical_centre)

        data_file.close()

        # I think titles don't work with orographic projection?
        ax.set_title(title, pad=20)

    suptitle = r'$t=$ '+f'{(time_idx*time_gap):.1f} hr'
    fig.suptitle(suptitle, y=1.02)

    print(f'Plotting to {plotpath}')
    fig.savefig(plotpath, bbox_inches='tight')
    plt.close()

