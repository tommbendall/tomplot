"""
This returns the coordinates of the data points for a specified field from
some output data file.
"""
import numpy as np
from .domain_properties import get_domain_properties, get_lfric_domain_properties

def get_td_data_coords(data, space_name, central_lon=0.0):

    #--------------------------------------------------------------------------#
    # Work out which coordinates we are using
    #--------------------------------------------------------------------------#

    domain = data.variables['domain'][0]
    names, _, _ = get_domain_properties(data, central_lon)

    #--------------------------------------------------------------------------#
    # Extract data coordinates
    #--------------------------------------------------------------------------#

    if len(names.keys()) == 1:
        data_coords = (data.variables[names['x']+'_'+space_name][:],)

    elif len(names.keys()) == 2:
        if 'y' in names.keys():
            data_coords = np.array([[data_x, data_y] for data_x, data_y in
                                    zip(data.variables[names['x']+'_'+space_name][:],
                                        data.variables[names['y']+'_'+space_name][:])])
        else:
            data_coords = np.array([[data_x, data_z] for data_x, data_z in
                                   zip(data.variables[names['x']+'_'+space_name][:],
                                       data.variables[names['z']+'_'+space_name][:])])

    else:
        data_coords = np.array([[data_x, data_y, data_z] for data_x, data_y, data_z in
                                zip(data.variables[names['x']+'_'+space_name][:],
                                    data.variables[names['y']+'_'+space_name][:],
                                    data.variables[names['z']+'_'+space_name][:])])

    # FIXME: Need a routine for adjusting coordinates

    # Adjust longitudinal coordinates if necessary
    # FIXME: Only implemented for sphere
    if domain == 'sphere' and names['x'] == 'lon':
        # Assume longitude data is first in the array
        # Adjust coordinates for all points outside of branch point
        data_coords.T[0][data_coords.T[0] < (central_lon - np.pi)] += 2*np.pi
        data_coords.T[0][data_coords.T[0] > (central_lon + np.pi)] -= 2*np.pi

    # Adjust cylinder phi coordinates to have dimensions of length
    # FIXME: Only implemented for cylinder
    if domain == 'cylinder' and names['x'] == 'phi':
        # Assume phi data is first in the array
        # Adjust coordinates for all points outside of branch point
        radius = data.variables['base_radius'][0]
        data_coords.T[0] *= radius

    return data_coords


def get_lfric_data_coords(data, hori_placement, extrusion_details=None,
                          central_lon=0.0):

    if extrusion_details is None:
        # TODO: could find a way of working this out from data file
        extrusion_details = {'domain':'plane',
                             'topological_dimension':3,
                             'zmin':0.0,
                             'zmax':1000.0,
                             'extrusion':'linear'}

    domain = extrusion_details['domain']

    #--------------------------------------------------------------------------#
    # Work out which coordinates we are using
    #--------------------------------------------------------------------------#

    names, _, _ = get_lfric_domain_properties(data,
                                              extrusion_details,
                                              central_lon)

    #--------------------------------------------------------------------------#
    # Extract data coordinates
    #--------------------------------------------------------------------------#

    if len(names.keys()) == 2:
        raise NotImplementedError('Set up of data coords not implemented')
    else:
        if domain == 'sphere':
            # Convert to radians
            data_coords_2d = (data.variables[hori_placement+'_'+names['x']][:]*np.pi/180.0,
                              data.variables[hori_placement+'_'+names['y']][:]*np.pi/180.0)
        else:
            data_coords_2d = (data.variables[hori_placement+'_'+names['x']][:],
                              data.variables[hori_placement+'_'+names['y']][:])

    if (domain == 'sphere' and abs(central_lon) > 1e-12):
        # Adjust longitudinal coordinates if necessary
        # Assume longitude data is first in the array
        # Adjust coordinates for all points outside of branch point
        data_coords_2d[0][data_coords_2d[0] < (central_lon - np.pi)] += 2*np.pi
        data_coords_2d[0][data_coords_2d[0] > (central_lon + np.pi)] -= 2*np.pi


    return data_coords_2d
