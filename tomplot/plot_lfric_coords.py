"""
This provides routines for setting up coordinates for plotting fields on, with
data extracted from LFRic
"""
import numpy as np
from .plot_decorations import get_domain_label
from .data_coords import get_lfric_data_coords
from .domain_properties import get_lfric_domain_properties


def get_lfric_coords_2d(data, hori_placement, vert_placement, slice_name,
                        slice_idx=1, num_points=None, extrusion_details=None,
                        length_units='m', angular_units='rad', central_lon=0.0,
                        plot_coords_1d=None):

    if extrusion_details is None:
        # TODO: could find a way of working this out from data file
        extrusion_details = {'domain':'plane',
                             'topological_dimension':3,
                             'zmin':0.0,
                             'zmax':1000.0,
                             'extrusion':'linear'}

    domain = extrusion_details['domain']
    dim = extrusion_details['topological_dimension']

    if vert_placement is None and slice_name != 'xy':
        raise ValueError(f'Cannot take slice {slice_name} of 2D data')

    if dim == 3:
        dim_names = ['x','y','z']

        if slice_name is None:
            raise ValueError('slice_name must be specified for 3D data')

    else:
        raise ValueError(f'Setting up LFRic coordinates not implemented '+
                         f'for domains of dimensions {dim}')

    #--------------------------------------------------------------------------#
    # Work out which coordinates we are using
    #--------------------------------------------------------------------------#

    # TODO: can all this code be combined with transportdrake code?
    names, domain_extents, ticklabels = get_lfric_domain_properties(data,
                                                                    extrusion_details,
                                                                    central_lon)

    coord_lims = [domain_extents[slice_name[0]], domain_extents[slice_name[1]]]
    coord_ticks = [ticklabels[slice_name[0]], ticklabels[slice_name[1]]]
    coord_labels = []
    for i in range(2):
        domain_label, units = get_domain_label(names[slice_name[i]], length_units, angular_units)
        coord_labels.append(domain_label+r'$ \ / \ $'+units)

    #--------------------------------------------------------------------------#
    # Determine number of points in each direction
    #--------------------------------------------------------------------------#

    if plot_coords_1d is not None:
        num = tuple([len(i) for i in plot_coords_1d])
    elif num_points is not None:
        # num points can be an integer or a 2D tuple
        if type(num_points) == int:
            num = tuple([num_points for i in range(dim)])
        elif len(num_points) != dim:
            raise TypeError('num_points must be an integer or list/tuple/array '+
                             'of length matching topological dimension')
        elif len(num_points) == 2:
            num = (num_points[0], num_points[1])
        else:
            num = (num_points[0], num_points[1], num_points[2])

    else:
        num_list = []

        for dim_name in dim_names:
            if dim_name == 'z':
                if vert_placement is None:
                    num_part = 0
                else:
                    num_part = len(np.unique(data[vert_placement][:]))
            elif domain == 'sphere':
                num_part = 101
            elif (names[dim_name] == 'x'):
                num_part = len(np.unique(data[hori_placement+'_x'][:]))
            elif (names[dim_name] == 'y'):
                num_part = len(data[hori_placement+'_y'][:])
            else:
                raise ValueError('Something has gone wrong')
            num_list.append(num_part)

        # Turn list into tuple
        num = tuple(num_list)

    # Get number of plotting points for each direction
    if slice_name == 'xy':
        num_plot_points_x = num[0]
        num_plot_points_y = num[1]
    elif slice_name == 'xz':
        num_plot_points_x = num[0]
        if dim == 2:
            num_plot_points_y = num[1]
        else:
            num_plot_points_y = num[2]
    elif slice_name == 'yz':
        num_plot_points_x = num[1]
        num_plot_points_y = num[2]
    else:
        raise ValueError('Slice name %s is not recognised' % slice_name)

    num_plot_points = num_plot_points_x * num_plot_points_y

    #--------------------------------------------------------------------------#
    # Make coordinates for plotting
    #--------------------------------------------------------------------------#

    if dim == 2:
        raise NotImplementedError('LFRic coords not implemented for dim 2')

    else:
        if plot_coords_1d is None:
            plot_coords_1d = (np.linspace(domain_extents['x'][0],
                                          domain_extents['x'][1], num[0]),
                              np.linspace(domain_extents['y'][0],
                                          domain_extents['y'][1], num[1]),
                              np.linspace(domain_extents['z'][0],
                                          domain_extents['z'][1], num[2]))

        if slice_name == 'xy':
            plot_coords = np.meshgrid(plot_coords_1d[0], plot_coords_1d[1])

            # Turn this into a list of points for interpolating to
            interp_coords_2d = np.array([[[coord_x, coord_y]
                                          for coord_x, coord_y in zip(plot_coords[0][:][j], plot_coords[1][:][j])]
                                         for j in range(num_plot_points_y)])

            # Must be at constant z
            domain_label, units = get_domain_label(names['z'], length_units, angular_units)
            if vert_placement is None:
                slice_label = ''
            else:
                slice_label = domain_label+r'$=$ %.1f' % plot_coords_1d[2][slice_idx]
                slice_label += ' '+units
            coord_lims = [domain_extents['x'], domain_extents['y']]
            coord_ticks = [ticklabels['x'], ticklabels['y']]

        elif slice_name == 'xz':
            plot_coords = np.meshgrid(plot_coords_1d[0], plot_coords_1d[2])

            # Turn this into a list of points for interpolating to
            interp_coords_2d = (np.array(plot_coords_1d[0]),
                                np.ones_like(plot_coords_1d[0])*plot_coords_1d[1][slice_idx])

            # Must be at constant y
            domain_label, units = get_domain_label(names['y'], length_units, angular_units)
            slice_label = domain_label+r'$=$ %.1f' % plot_coords_1d[1][slice_idx]
            slice_label += ' '+units
            coord_lims = [domain_extents['x'], domain_extents['z']]
            coord_ticks = [ticklabels['x'], ticklabels['z']]

        elif slice_name == 'yz':
            plot_coords = np.meshgrid(plot_coords_1d[1], plot_coords_1d[2])

            # Turn this into a list of points for interpolating to
            interp_coords_2d = (plot_coords_1d[0][slice_idx]*np.ones_like(plot_coords_1d[1]),
                                np.array(plot_coords_1d[1]))

            # Must be at constant x
            domain_label, units = get_domain_label(names['x'], length_units, angular_units)
            slice_label = domain_label+r'$=$ %.1f' % plot_coords_1d[0][slice_idx]
            slice_label += ' '+units
            coord_lims = [domain_extents['y'], domain_extents['z']]
            coord_ticks = [ticklabels['y'], ticklabels['z']]

        else:
            raise ValueError('Slice %s not recognised')

    #--------------------------------------------------------------------------#
    # Extract data coordinates
    #--------------------------------------------------------------------------#

    data_coords_2d = get_lfric_data_coords(data, hori_placement,
                                           extrusion_details, central_lon)

    return plot_coords, data_coords_2d, interp_coords_2d, \
           coord_labels, coord_lims, coord_ticks, slice_label

