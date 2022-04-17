"""
Demonstrates the orthographic projection
"""

import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np
import cartopy.crs as ccrs
from tomplot import individual_quiver_plot, extract_2D_data

# ---------------------------------------------------------------------------- #
# Things that can be altered and parameters for the test case
# ---------------------------------------------------------------------------- #

results_dirname = 'demo_4_spherical_hooks'
plotname = 'fig_3_sphere_paths'
time_idx = 0
run_id = 0
field_info = ('F_0', 'F_0_zonal', 'F_0_meridional')
cbar_label = r'$|F| \ /$ m s$^{-1}$'
titles = [r'Front', r'Back']
colour_scheme = 'OrRd'

# ---------------------------------------------------------------------------- #
# Things that are likely the same
# ---------------------------------------------------------------------------- #

plotdir = 'results/vector_transport_paper/figures'
field_name, field_X_name, field_Y_name = field_info
slice = 'xy'
slice_idx = 0
filename = 'results/'+results_dirname+'/nc_fields/field_output_'+str(run_id)+'.nc'

# This is declared BEFORE figure and ax are initialised
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
font = {'size':24}
plt.rc('font',**font)

# ---------------------------------------------------------------------------- #
# Field plots
# ---------------------------------------------------------------------------- #

# Add paths
import cartopy.crs as ccrs
crs = ccrs.PlateCarree()

fig = plt.figure(figsize=(16,8))

plotpath = f'{plotdir}/{plotname}.jpg'

for i, (spherical_centre, title) in enumerate(zip([(0.0,-np.pi/6),(np.pi,-np.pi/6)], titles)):
    # This code is all adapted from plot_control
    data_file = Dataset(filename, 'r')

    # Extract data
    coords_X, coords_Y, field_X_data, time, \
    coord_labels, coord_lims, coord_ticks,  \
    slice_label = extract_2D_data(data_file, field_X_name, time_idx,
                                  slice_name=slice, slice_idx=slice_idx,
                                  central_lon=spherical_centre[0])

    coords_X, coords_Y, field_Y_data, time, \
    coord_labels, coord_lims, coord_ticks,  \
    slice_label = extract_2D_data(data_file, field_Y_name, time_idx,
                                  slice_name=slice, slice_idx=slice_idx,
                                  central_lon=spherical_centre[0])

    data_file.close()

    lon_centre, lat_centre = spherical_centre[0]*180.0/np.pi, spherical_centre[1]*180.0/np.pi
    ax = fig.add_subplot(1, 2, 1+i,
                         projection=ccrs.Orthographic(lon_centre, lat_centre))

    coords_X *= 180.0/np.pi
    coords_Y *= 180.0/np.pi
    zeros = np.zeros_like(coords_X)
    ax.scatter(coords_X, coords_Y, color='white', transform=crs)
    gl = ax.gridlines(xlocs=[0.,45.,135.,180.,-135.,-45.])

    if i == 0:
        plt.text(0.2,0.94,r'45$^\mathrm{o}$W', ha='center', va='bottom', fontsize=20, transform=ax.transAxes)
        plt.text(0.5,1.0,r'0$^\mathrm{o}$', ha='center', va='bottom', fontsize=20, transform=ax.transAxes)
        plt.text(0.8,0.94,r'45$^\mathrm{o}$E', ha='center', va='bottom', fontsize=20, transform=ax.transAxes)
    elif i == 1:
        plt.text(0.2,0.94,r'135$^\mathrm{o}$E', ha='center', va='bottom', fontsize=20, transform=ax.transAxes)
        plt.text(0.5,1.0,r'180$^\mathrm{o}$', ha='center', va='bottom', fontsize=20, transform=ax.transAxes)
        plt.text(0.8,0.94,r'135$^\mathrm{o}$W', ha='center', va='bottom', fontsize=20, transform=ax.transAxes)


    path_1_lon = np.linspace(0.0, 180.0, 51)
    path_1_lat = np.linspace(-30.0, -30.0, 51)
    ax.plot(path_1_lon, path_1_lat, marker='', linestyle='-', lw=15,
            color='grey', transform=crs)


    new_lon = np.pi/180.0*np.linspace(-90.0, 90.0, 19)
    new_lat = np.pi/180.0*np.linspace(-60.0, -60.0, 19)
    path_2_lon = np.arctan2(-np.cos(new_lon)*np.cos(new_lat), np.sin(new_lat))
    path_2_lat = np.arctan2(np.sin(new_lon)*np.cos(new_lat),
        (np.sqrt(np.sin(new_lat)**2 + np.cos(new_lon)**2*np.cos(new_lat)**2)))
    path_2_lon *= 180.0/np.pi
    path_2_lat *= 180.0/np.pi
    ax.plot(path_2_lon, path_2_lat, marker='', linestyle='-', lw=15,
            color='grey', transform=crs)
    if i == 1:
        ax.plot(path_2_lon[0], path_2_lat[0], marker='o', linestyle='',
                color='black', ms=30, transform=crs)
        ax.text(path_2_lon[0], path_2_lat[0], r'$\textbf{1}$', color='white',
                ha='center', va='center', transform=crs)

    path_3_lon = np.linspace(180.0, 360.0, 51)
    path_3_lat = np.linspace(30.0, 30.0, 51)
    ax.plot(path_3_lon, path_3_lat, marker='', linestyle='-', lw=15,
            color='grey', transform=crs)
    if i == 1:
        ax.plot(path_3_lon[0], path_3_lat[0], marker='o', linestyle='',
                color='black', ms=30, transform=crs)
        ax.text(path_3_lon[0], path_3_lat[0], r'$\textbf{2}$', color='white',
                ha='center', va='center', transform=crs)

    new_lon = np.pi/180.0*np.linspace(90.0, -90.0, 19)
    new_lat = np.pi/180.0*np.linspace(60.0, 60.0, 19)
    path_4_lon = np.arctan2(np.cos(new_lon)*np.cos(new_lat), np.sin(new_lat))
    path_4_lat = np.arctan2(np.sin(new_lon)*np.cos(new_lat),
        (np.sqrt(np.sin(new_lat)**2 + np.cos(new_lon)**2*np.cos(new_lat)**2)))
    path_4_lon *= 180.0/np.pi
    path_4_lat *= 180.0/np.pi
    ax.plot(path_4_lon, path_4_lat, marker='', linestyle='-', lw=15,
            color='grey', transform=crs)
    if i == 0:
        ax.plot(path_4_lon[0], path_4_lat[0], marker='o', linestyle='',
                color='black', ms=30, transform=crs)
        ax.text(path_4_lon[0], path_4_lat[0], r'$\textbf{3}$', color='white',
                ha='center', va='center', transform=crs)
    if i == 0:
        ax.plot(path_4_lon[-1], path_4_lat[-1], marker='o', linestyle='',
                color='black', ms=30, transform=crs)
        ax.text(path_4_lon[-1], path_4_lat[-1], r'$\textbf{0}$', color='white',
                ha='center', va='center', transform=crs)

print(f'Plotting to {plotpath}')
fig.savefig(plotpath, bbox_inches='tight')
plt.close()