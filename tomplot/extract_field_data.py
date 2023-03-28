"""
This file provides some plotting tools
"""
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np
from .plot_coordinates import *


def extract_1D_data(data_file, field_name, time_idx,
                    slice_name=None, slice_idx=None, num_points=None,
                    plot_coords_1d=None):
    """
    Extract 1D field data from a netCDF field output file, ready to plot.

    First the coordinate grid for the plotting is set-up, and then the field is
    interpolated onto this.

    :arg data_file:  the netCDF data file
    :arg field_name: a string giving the name of the field to plot
    :arg time_idx:   the index in time at which to extract the field
    :arg slice_name: (Optional). The dimension along which to slice the data.
                     Required if the domain has dimension greater than 1.
    :arg slice_idx:  (Optional). The index of the other dimensions on the
                     interpolation grid at which to slice
    :arg num_points: (Optional). The number of coordinate points (in 1D) to use
                     in the plotting grid. If not specified then this will be
                     determined from the data, and if it can't then set to be 100.
    :arg plot_coords_1d: (Optional) An array of 1D arrays of the points at which
                         to plot. If not specified then this will be determined
                         from the data
    """

    #--------------------------------------------------------------------------#
    # Checks
    #--------------------------------------------------------------------------#

    # TODO: fill in checks here

    if slice_name not in ['x','y','z']:
        raise ValueError('For 1D plots slice must be x, y or z')

    time = data_file.variables['time'][time_idx]

    #--------------------------------------------------------------------------#
    # Construct coordinate points
    #--------------------------------------------------------------------------#

    # Get space name
    coords_name = data_file.groups[field_name].variables['field_values'].dimensions[0]
    # This should be a string 'coords_XXX'
    space_name = coords_name[7:]

    plot_coords, data_coords, interp_coords, \
    coord_label, coord_lims, coord_ticks, \
    slice_label = get_coords_1d(data_file, space_name,
                                slice_name, slice_idx, num_points,
                                plot_coords_1d=plot_coords_1d)

    #--------------------------------------------------------------------------#
    # Interpolate field data
    #--------------------------------------------------------------------------#

    from scipy.interpolate import griddata

    field_full = data_file.groups[field_name].variables['field_values'][:,time_idx]
    field_data = griddata(data_coords, field_full, interp_coords, method='linear')
    field_near = griddata(data_coords, field_full, interp_coords, method='nearest')
    field_data[np.isnan(field_data)] = field_near[np.isnan(field_data)]

    # Guarantee that plot_coords and field_data will be the same shape
    plot_coords = np.reshape(plot_coords, np.shape(field_data))

    data_metadata = {'time': time, 'slice_label': slice_label,
                     'coord_labels': coord_label, 'coord_lims': coord_lims,
                     'coord_ticks': coord_ticks}


    return plot_coords, field_data, data_metadata


def extract_2D_data(data_file, field_name, time_idx,
                    slice_name=None, slice_idx=None, num_points=None,
                    central_lon=0.0, plot_coords_1d=None,
                    interpolate=True):
    """
    Extract 2D field data from a netCDF field output file, ready to plot.

    First the coordinate grid for the plotting is set-up, and then the field is
    interpolated onto this.

    :arg data_file:  the netCDF data file
    :arg field_name: a string giving the name of the field to plot
    :arg time_idx:   the index in time at which to extract the field
    :arg slice_name: (Optional). The dimension along which to slice the data.
                     Required if the domain has dimension greater than 1.
    :arg slice_idx:  (Optional). The index of the other dimensions on the
                     interpolation grid at which to slice
    :arg num_points: (Optional). The number of coordinate points (in 1D) to use
                     in the plotting grid. If not specified then this will be
                     determined from the data, and if it can't then set to be 100.
    :arg central_lon: (Optional) longitude at centre of plot. Only implemented
                      for spherical domains.
    :arg plot_coords_1d: (Optional) An array of 1D arrays of the points at which
                         to plot. If not specified then this will be determined
                         from the data
    :arg interpolate: Boolean, whether to interpolate to a meshgrid or not.
    """

    #--------------------------------------------------------------------------#
    # Checks
    #--------------------------------------------------------------------------#

    # TODO: fill in checks here

    if slice_name not in ['xy','yz','xz', None]:
        raise ValueError('For 2D plots slice must be xy, yz or xz')

    time = data_file.variables['time'][time_idx]

    #--------------------------------------------------------------------------#
    # Construct coordinate points
    #--------------------------------------------------------------------------#

    # Get space name
    coords_name = data_file.groups[field_name].variables['field_values'].dimensions[0]
    # This should be a string 'coords_XXX'
    space_name = coords_name[7:]

    if interpolate:
        plot_coords, data_coords, interp_coords, \
        coord_labels, coord_lims, coord_ticks, \
        slice_label = get_coords_2d(data_file, space_name, slice_name, slice_idx,
                                    num_points, central_lon=central_lon,
                                    plot_coords_1d=plot_coords_1d)

        #--------------------------------------------------------------------------#
        # Interpolate field data
        #--------------------------------------------------------------------------#

        from scipy.interpolate import griddata

        field_full = data_file.groups[field_name].variables['field_values'][:,time_idx]
        field_data = griddata(data_coords, field_full, interp_coords, method='linear')
        field_near = griddata(data_coords, field_full, interp_coords, method='nearest')
        field_data[np.isnan(field_data)] = field_near[np.isnan(field_data)]

        data_metadata = {'time': time, 'slice_label': slice_label,
                        'coord_labels': coord_labels, 'coord_lims': coord_lims,
                        'coord_ticks': coord_ticks}

    else:
        field_data = data_file.groups[field_name].variables['field_values'][:,time_idx]
        data_metadata = None

        if f'lon_{space_name}' in data_file.variables.keys():
            data_x = data_file.variables[f'lon_{space_name}'][:]
            # Adjust based on central longitude
            data_x[data_x < (central_lon - np.pi)] += 2*np.pi
            data_x[data_x > (central_lon + np.pi)] -= 2*np.pi

        elif f'x_{space_name}' in data_file.variables.keys():
            data_x = data_file.variables[f'x_{space_name}'][:]
        else:
            raise KeyError('Could not work out x coordinates')

        # Try to figure out what the domain is
        if slice_name == 'xy':
            if f'lat_{space_name}' in data_file.variables.keys():
                data_y = data_file.variables[f'lat_{space_name}'][:]
            elif f'y_{space_name}' in data_file.variables.keys():
                data_y = data_file.variables[f'y_{space_name}'][:]
            else:
                raise KeyError('Could not work out y coordinates')

            plot_coords = (data_x, data_y)
        elif slice_name == 'xz':
            if f'z_{space_name}' in data_file.variables.keys():
                data_z = data_file.variables[f'z_{space_name}'][:]
            else:
                raise KeyError('Could not work out z coordinates')

            plot_coords = (data_x, data_z)
        else:
            raise NotImplementedError(f'Slice {slice_name} not implemented')

    return plot_coords[0], plot_coords[1], field_data, data_metadata
