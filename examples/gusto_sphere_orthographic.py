"""
A tomplot example, using an orthographic spherical projection with Gusto data.
"""

import cartopy.crs as ccrs
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from os.path import abspath, dirname
from tomplot import (set_tomplot_style, tomplot_contours, tomplot_cmap,
                     plot_contoured_field, add_colorbar_fig,
                     tomplot_field_title, extract_gusto_coords,
                     extract_gusto_field, regrid_horizontal_slice)

# ---------------------------------------------------------------------------- #
# Directory for results and plots
# ---------------------------------------------------------------------------- #
# When copying this example these should not be relative to this file
results_dir = f'{abspath(dirname(__file__))}/../tests/data'
plot_dir = f'{abspath(dirname(__file__))}/../tests/tmp_figures'
results_file_name = f'{results_dir}/gusto_sphere_2d_field_output.nc'
plot_name = f'{plot_dir}/example_gusto_orthographic.png'

# ---------------------------------------------------------------------------- #
# Things that should be altered based on the plot
# ---------------------------------------------------------------------------- #
num_subplots = 2
# Specify lists for variables that are different between subplots
spherical_centres = [(90, 60), (-90, -15)]  # Centres of spherical projections
subtitles = [f'Centered on {spherical_centre}' for spherical_centre in spherical_centres]
# Things that are the same for both subplots
field_name = 'D'  # Plot depth on front and back of sphere
colour_scheme = 'RdYlBu_r'
field_label = r'$D \ / $ m'
time_idx = -1
contour_method = 'contour'  # Best method for orthographic projection
remove_contour = 10000.
# We need to regrid onto lon-lat grid -- specify that here
lon_1d = np.linspace(-180.0, 180.0, 80)
lat_1d = np.linspace(-90, 90, 80)
coords_lon, coords_lat = np.meshgrid(lon_1d, lat_1d, indexing='ij')

# ---------------------------------------------------------------------------- #
# Things that are likely the same for all plots
# ---------------------------------------------------------------------------- #
set_tomplot_style()
data_file = Dataset(results_file_name, 'r')
# Don't create fig with subplots
# When using cartopy we need to specify projections for each axes separately
fig = plt.figure(figsize=(10, 5))

# ---------------------------------------------------------------------------- #
# Data extraction
# ---------------------------------------------------------------------------- #
# Since we are looping through front/back of sphere, we can extact data first
orig_field_data = extract_gusto_field(data_file, field_name, time_idx=time_idx)
orig_coords_lon, orig_coords_lat = extract_gusto_coords(data_file, field_name)
time = data_file['time'][time_idx]

# For orthographic projection to look good, need to regrid onto lon-lat grid
field_data = regrid_horizontal_slice(coords_lon, coords_lat, orig_coords_lon,
                                     orig_coords_lat, orig_field_data,
                                     periodic_fix='sphere')

contours = tomplot_contours(field_data)
cmap, lines = tomplot_cmap(contours, colour_scheme, remove_contour=remove_contour)

# Loop through subplots
for i, (spherical_centre, subtitle) in \
    enumerate(zip(spherical_centres, subtitles)):

    # Specify projection and set up axes here
    projection = ccrs.Orthographic(spherical_centre[0], spherical_centre[1])
    ax = fig.add_subplot(1, num_subplots, 1+i, projection=projection)

    # ------------------------------------------------------------------------ #
    # Plot data
    # ------------------------------------------------------------------------ #
    cf, _ = plot_contoured_field(ax, coords_lon, coords_lat, field_data,
                                 contour_method, contours, cmap=cmap,
                                 line_contours=lines,
                                 projection=projection)
    ax.set_title(subtitle)

# Share colorbar and title
add_colorbar_fig(fig, cf, None)
# Generate title with return_title=True
suptitle = tomplot_field_title(None, f't = {time:.1f}', return_title=True,
                               minmax=True, field_data=field_data)
fig.suptitle(suptitle)
# ---------------------------------------------------------------------------- #
# Save figure
# ---------------------------------------------------------------------------- #
print(f'Saving figure to {plot_name}')
fig.savefig(plot_name, bbox_inches='tight')
plt.close()


