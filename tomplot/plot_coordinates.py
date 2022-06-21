"""
This provides routines for setting up coordinates for plotting fields on
"""
import numpy as np
from .plot_decorations import get_domain_label
from .domain_properties import get_domain_properties
from .data_coords import get_td_data_coords

def get_coords_1d(data, space_name, slice_name, slice_idx=0, num_points=None,
                  length_units='m', angular_units='rad', plot_coords_1d=None,
                  central_lon=0.0):
    """
    Returns coordinates for the mesh and the coordinates of the data points
    for use in creating 1D slice plots.

    :arg data:       netCDF data file
    :arg space_name: the name of the space to get the coordinates for
    :arg slice_name: the name of the slice
    :arg slice_idx:  the index to slice along
    :arg num_points: the number of points for the plotting mesh
    :arg plot_coords_1d: (Optional) An array of 1D arrays of the points at which
                         to plot. If not specified then this will be determined
                         from the data
    """

    # TODO: allow different slice indices for 1D slices of 3D data

    #--------------------------------------------------------------------------#
    # Work out which coordinates we are using
    #--------------------------------------------------------------------------#

    domain = data.variables['domain'][0]
    names, domain_extents, ticklabels = get_domain_properties(data, central_lon)

    if slice_name is None:
        slice_name = 'x'

    #--------------------------------------------------------------------------#
    # Set up coordinates for plot
    #--------------------------------------------------------------------------#

    dim = data.variables['topological_dimension'][0]

    if plot_coords_1d is not None:
        num_plot_points = tuple([len(i) for i in plot_coords_1d])
    elif num_points is not None:
        num_plot_points = num_points
    elif slice_name == 'z':
        num_plot_points = int(data.variables['nz'][0])
    elif ((names[slice_name] == 'x') or
          (slice_name == 'x' and domain in ['circle', 'cylinder'])):
        num_plot_points = int(data.variables['nx'][0])
    elif ((names[slice_name] == 'y') or
          (slice_name == 'y' and domain == 'cylinder')):
        num_plot_points = int(data.variables['ny'][0])
    else:
        # We are a sphere or a torus
        num_plot_points = 101

    if slice_idx == 'midpoint':
        slice_idx = int(num_plot_points / 2)


    if dim == 1:
        if plot_coords_1d is None:
            plot_coords_1d = (np.linspace(domain_extents['x'][0],
                                          domain_extents['x'][1], num_plot_points),)
        plot_coords = plot_coords_1d
        interp_coords = plot_coords_1d

        # We have not sliced
        slice_label = None

    elif dim == 2:
        if 'y' in names.keys():
            if plot_coords_1d is None:
                plot_coords_1d = (np.linspace(domain_extents['x'][0],
                                              domain_extents['x'][1], num_plot_points),
                                  np.linspace(domain_extents['y'][0],
                                              domain_extents['y'][1], num_plot_points))

        else:
            if plot_coords_1d is None:
                plot_coords_1d = (np.linspace(domain_extents['x'][0],
                                              domain_extents['x'][1], num_plot_points),
                                  np.linspace(domain_extents['z'][0],
                                              domain_extents['z'][1], num_plot_points))

        if slice_name == 'x':
            plot_coords = plot_coords_1d[0]
            interp_coords = np.array([[coord, plot_coords_1d[1][slice_idx]] for coord in plot_coords])

            if 'y' in names.keys():
                domain_label, units = get_domain_label(names['y'], length_units, angular_units)
            else:
                domain_label, units = get_domain_label(names['z'], length_units, angular_units)
            slice_label = domain_label+r'$=$ %.1f' % plot_coords_1d[1][slice_idx]
            slice_label += ' '+units

        else:
            plot_coords = plot_coords_1d[1]
            interp_coords = np.array([[plot_coords_1d[0][slice_idx], coord] for coord in plot_coords])

            # Must be at constant x
            domain_label, units = get_domain_label(names['x'], length_units, angular_units)
            slice_label = domain_label+r'$=$ %.1f' % plot_coords_1d[0][slice_idx]
            slice_label += ' '+units

    else:
        plot_coords_1d = (np.linspace(domain_extents['x'][0],
                                      domain_extents['x'][1], num_plot_points),
                          np.linspace(domain_extents['y'][0],
                                      domain_extents['y'][1], num_plot_points),
                          np.linspace(domain_extents['z'][0],
                                      domain_extents['z'][1], num_plot_points))

        if slice_name == 'x':
            plot_coords = plot_coords_1d[0]

            # Turn this into a list of points for interpolating to
            interp_coords = np.array([[coord_x, plot_coords_1d[1][slice_idx], plot_coords_1d[2][slice_idx]]
                                       for coord_x in plot_coords])

            # Must be at constant y and z
            domain_label, units = get_domain_label(names['y'], length_units, angular_units)
            slice_label = domain_label+r'$=$ %.1f' % plot_coords_1d[1][slice_idx]
            slice_label += ' '+units
            domain_label, units = get_domain_label(names['z'], length_units, angular_units)
            slice_label += ', '+domain_label+r'$=$ %.1f' % plot_coords_1d[2][slice_idx]
            slice_label += ' '+units

        elif slice_name == 'y':
            plot_coords = plot_coords_1d[1]

            # Turn this into a list of points for interpolating to
            interp_coords = np.array([[plot_coords_1d[0][slice_idx], coord_y, plot_coords_1d[2][slice_idx]]
                                       for coord_y in plot_coords])

            # Must be at constant x and z
            domain_label, units = get_domain_label(names['x'], length_units, angular_units)
            slice_label = domain_label+r'$=$ %.1f' % plot_coords_1d[0][slice_idx]
            slice_label += ' '+units
            domain_label, units = get_domain_label(names['z'], length_units, angular_units)
            slice_label += ', '+domain_label+r'$=$ %.1f' % plot_coords_1d[2][slice_idx]
            slice_label += ' '+units

        else:
            plot_coords = plot_coords_1d[2]

            # Turn this into a list of points for interpolating to
            interp_coords = np.array([[plot_coords_1d[0][slice_idx], plot_coords_1d[1][slice_idx], coord_z]
                                       for coord_z in plot_coords])

            # Must be at constant x and y
            domain_label, units = get_domain_label(names['x'], length_units, angular_units)
            slice_label = domain_label+r'$=$ %.1f' % plot_coords_1d[0][slice_idx]
            slice_label += ' '+units
            domain_label, units = get_domain_label(names['y'], length_units, angular_units)
            slice_label += ', '+domain_label+r'$=$ %.1f' % plot_coords_1d[1][slice_idx]
            slice_label += ' '+units

    # Set up coordinate limits and labels
    coord_lims = domain_extents[slice_name]
    coord_ticks = ticklabels[slice_name]
    domain_label, units = get_domain_label(names[slice_name], length_units, angular_units)
    coord_label = domain_label+r'$ \ / \ $'+units

    #--------------------------------------------------------------------------#
    # Extract data coordinates
    #--------------------------------------------------------------------------#

    data_coords = get_td_data_coords(data, space_name, central_lon)

    return plot_coords, data_coords, interp_coords, coord_label, coord_lims, coord_ticks, slice_label


def get_coords_2d(data, space_name, slice_name, slice_idx=0, num_points=None,
                  length_units='m', angular_units='rad', central_lon=0.0,
                  plot_coords_1d=None):
    """
    Returns coordinates for the mesh and the coordinates of the data points
    for use in creating 2D slice plots.

    :arg data:        netCDF data file
    :arg space_name:  the name of the space to get the coordinates for
    :arg slice_name:  the name of the slice
    :arg slice_idx:   the index to slice along
    :arg num_points:  the number of points for the plotting mesh
    :arg central_lon: (Optional) longitude at centre of plot. Only implemented
                      for spherical domains.
    :arg plot_coords_1d: (Optional) An array of 1D arrays of the points at which
                         to plot. If not specified then this will be determined
                         from the data
    """

    domain = data.variables['domain'][0]
    extruded = True if data.variables['extrusion'][0] == 'True' else False
    dim = data.variables['topological_dimension'][0]

    if dim < 2:
        raise ValueError('The domain must be two dimensional or higher')

    if dim == 2:
        if extruded:
            dim_names = ['x', 'z']
            if slice_name is None:
                slice_name = 'xz'
        else:
            dim_names = ['x', 'y']
            if slice_name is None:
                slice_name = 'xy'
    else:
        dim_names = ['x', 'y', 'z']

        if slice_name is None:
            raise ValueError('slice_name must be specified for 3D data')

    #--------------------------------------------------------------------------#
    # Work out which coordinates we are using
    #--------------------------------------------------------------------------#

    domain = data.variables['domain'][0]
    names, domain_extents, ticklabels = get_domain_properties(data, central_lon)

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
                num_part = int(data.variables['nz'][0])
            elif ((names[dim_name] == 'x') or
                  (dim_name == 'x' and domain in ['circle', 'cylinder'])):
                num_part = int(data.variables['nx'][0])
            elif ((names[dim_name] == 'y') or
                  (dim_name == 'y' and domain == 'cylinder')):
                num_part = int(data.variables['ny'][0])
            else:
                # We are a sphere or a torus
                num_part = 101
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

        # We are not slicing
        slice_label = None

        if 'y' in names.keys():
            if plot_coords_1d is None:
                plot_coords_1d = (np.linspace(domain_extents['x'][0],
                                              domain_extents['x'][1], num[0]),
                                  np.linspace(domain_extents['y'][0],
                                              domain_extents['y'][1], num[1]))
            coord_lims = [domain_extents['x'], domain_extents['y']]
            coord_ticks = [ticklabels['x'], ticklabels['y']]

        else:
            if plot_coords_1d is None:
                plot_coords_1d = (np.linspace(domain_extents['x'][0],
                                              domain_extents['x'][1], num[0]),
                                  np.linspace(domain_extents['z'][0],
                                              domain_extents['z'][1], num[1]))
            coord_lims = [domain_extents['x'], domain_extents['z']]
            coord_ticks = [ticklabels['x'], ticklabels['z']]

        plot_coords = np.meshgrid(plot_coords_1d[0], plot_coords_1d[1])

        # Turn this into a list of points for interpolating to
        interp_coords = np.array([[[coord_x, coord_y] for coord_x, coord_y in
                                    zip(plot_coords[0][:][j], plot_coords[1][:][j])]
                                    for j in range(num_plot_points_y)])

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
            interp_coords = np.array([[[coord_x, coord_y, plot_coords_1d[2][slice_idx]] for coord_x, coord_y in
                                        zip(plot_coords[0][:][j], plot_coords[1][:][j])]
                                        for j in range(num_plot_points_y)])

            # Must be at constant z
            domain_label, units = get_domain_label(names['z'], length_units, angular_units)
            slice_label = domain_label+r'$=$ %.1f' % plot_coords_1d[2][slice_idx]
            slice_label += ' '+units
            coord_lims = [domain_extents['x'], domain_extents['y']]
            coord_ticks = [ticklabels['x'], ticklabels['y']]

        elif slice_name == 'xz':
            plot_coords = np.meshgrid(plot_coords_1d[0], plot_coords_1d[2])

            # Turn this into a list of points for interpolating to
            interp_coords = np.array([[[coord_x, plot_coords_1d[1][slice_idx], coord_z] for coord_x, coord_z in
                                       zip(plot_coords[0][:][j], plot_coords[2][:][j])]
                                       for j in range(num_plot_points_y)])

            # Must be at constant y
            domain_label, units = get_domain_label(names['y'], length_units, angular_units)
            slice_label = domain_label+r'$=$ %.1f' % plot_coords_1d[1][slice_idx]
            slice_label += ' '+units
            coord_lims = [domain_extents['x'], domain_extents['z']]
            coord_ticks = [ticklabels['x'], ticklabels['z']]

        elif slice_name == 'yz':
            plot_coords = np.meshgrid([plot_coords_1d[0][slice_idx]], plot_coords_1d[1], plot_coords_1d[2])

            # Turn this into a list of points for interpolating to
            interp_coords = np.array([[[plot_coords_1d[0][slice_idx], coord_y, coords_z]  for coord_y, coord_z in
                                       zip(plot_coords[1][:][j], plot_coords[2][:][j])]
                                       for j in range(num_plot_points_y)])

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

    data_coords = get_td_data_coords(data, space_name, central_lon)

    # FIXME: Need a routine for adjusting coordinates

    # Adjust cylinder phi coordinates to have dimensions of length
    # FIXME: Only implemented for cylinder
    if domain == 'cylinder' and names['x'] == 'phi':
        # Assume phi data is first in the array
        # Adjust coordinates for all points outside of branch point
        radius = data.variables['base_radius'][0]
        plot_coords[0] *= radius
        interp_coords.T[0] *= radius
        coord_lims[0] = domain_extents['x']*radius
        coord_labels[0] = r'$r\phi \ / $ m'
        coord_ticks[0] = None
        # data_coords already adjusted!

    return plot_coords, data_coords, interp_coords, coord_labels, coord_lims, coord_ticks, slice_label
