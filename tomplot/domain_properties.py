"""
This returns details about the domain from a data file.
"""
import numpy as np

def work_out_gusto_domain(data):

    # Guess space names
    if 'coords_DG' in data.dimensions.keys():
        space_name = 'DG'
    elif 'coords_DG0' in data.dimensions.keys():
        space_name = 'DG0'
    elif 'coords_CG2' in data.dimensions.keys():
        space_name = 'CG2'
    elif 'coords_theta' in data.dimensions.keys():
        space_name = 'theta'
    else:
        raise NotImplementedError('Need a space name to work out domain details')

    lon_bool = f'lon_{space_name}' in data.variables.keys()
    lat_bool = f'lat_{space_name}' in data.variables.keys()
    x_bool = f'x_{space_name}' in data.variables.keys()
    y_bool = f'y_{space_name}' in data.variables.keys()
    z_bool = f'z_{space_name}' in data.variables.keys()

    if lon_bool and lat_bool:
        domain = 'sphere'
        dim = 3 if z_bool else 2
    else:
        if not x_bool:
            raise KeyError('Could not work out coordinates')
        if not y_bool and not z_bool:
            raise ValueError('Cannot work in 1D')
        if y_bool:
            domain = 'plane'
            dim = 3 if z_bool else 2
        else:
            domain = 'interval'
            dim = 2
    extruded = z_bool

    return domain, extruded, dim

def get_domain_properties(data, domain, extruded, central_lon=0.0):
    """
    Get dictionaries of information about the domain from the data metadata.

    :arg central_lon: (Optional) Longitude of centre of plot. Currently only
                      implemented for spherical domains. Default is 0.0.
                      Should be given in radians.
    """

    if 'domain' not in data.variables.keys() and domain != 'sphere':
        # Guess a space name
        if 'coords_DG' in data.dimensions.keys():
            space_name = 'DG'
        elif 'coords_DG0' in data.dimensions.keys():
            space_name = 'DG0'
        elif 'coords_CG2' in data.dimensions.keys():
            space_name = 'CG2'
        elif 'coords_theta' in data.dimensions.keys():
            space_name = 'theta'
        else:
            raise NotImplementedError('Need a space name to work out domain details')
        
        if domain == 'plane':
            names = {'x':'x', 'y':'y'}
            domain_extents = {'x': (0, np.max(data.variables[f'x_{space_name}'])),
                              'y': (0, np.max(data.variables[f'y_{space_name}']))}
            ticklabels = {'x':None, 'y':None, 'z':None}

            if extruded:
                raise NotImplementedError
        else:
            raise NotImplementedError
    
    else:

        #----------------------------------------------------------------------#
        # Work out which coordinates we are using
        #----------------------------------------------------------------------#

        if domain == 'interval':
            names = {'x':'x'}
            domain_extents = {'x': (0, data.variables['Lx'][0])}
            ticklabels = {'x':None, 'z':None}
            if extruded:
                names['z'] = 'z'
                domain_extents['z'] = (0, data.variables['H'][0])

        elif domain == 'circle':
            names = {'x':'phi'}
            domain_extents = {'x': (-np.pi, np.pi)}
            ticklabels = {'x':(r'$-\pi$', r'$\pi$'),
                          'z':None}

            if extruded:
                names['z'] = 'r'
                domain_extents['z'] = (0, data.variables['H'][0])

        elif domain == 'plane':
            names = {'x':'x', 'y':'y'}
            domain_extents = {'x': (0, data.variables['Lx'][0]),
                              'y': (0, data.variables['Ly'][0])}
            ticklabels = {'x':None, 'y':None, 'z':None}

            if extruded:
                names['z'] = 'z'
                domain_extents['z'] = (0, data.variables['H'][0])

        elif domain == 'cylinder':
            names = {'x':'phi', 'y':'z'}
            domain_extents = {'x': (-np.pi, np.pi),
                              'y': (0, data.variables['L'][0])}
            ticklabels = {'x':(r'$-\pi$', r'$\pi$'),
                          'y':None,
                          'z':None}

            if extruded:
                names['z'] = 'h'
                domain_extents['z'] = (0, data.variables['H'][0])

        elif domain == 'torus':
            names = {'x':'lambda','y':'sigma'}
            domain_extents = {'x': (-np.pi, np.pi),
                              'y': (-np.pi, np.pi)}
            ticklabels = {'x':(r'$-\pi$', r'$\pi$'),
                          'y':(r'$-\pi$', r'$\pi$'),
                          'z':None}

            if extruded:
                names['z'] = 'h'
                domain_extents['z'] = (0, data.variables['H'][0])

        elif domain == 'sphere':
            names = {'x':'lon', 'y':'lat'}
            domain_extents = {'x': (-np.pi+central_lon, np.pi+central_lon),
                              'y': (-np.pi/2, np.pi/2)}
            if abs(central_lon) < 1e-6:
                xticklabels = (r'$-\pi$', r'$\pi$')
            elif abs(central_lon - np.pi) < 1e-6:
                xticklabels = (r'0', r'$2\pi$')
            elif abs(central_lon + np.pi) < 1e-6:
                xticklabels = (r'$-2\pi$', r'$0$')
            elif abs(central_lon + np.pi/2) < 1e-6:
                xticklabels = (r'$-3\pi/2$', r'$\pi/2$')
            elif abs(central_lon - np.pi/2) < 1e-6:
                xticklabels = (r'$-\pi/2$', r'$3\pi/2$')
            else:
                xticklabels = None

            ticklabels = {'x':xticklabels,
                          'y':(r'$-\pi/2$', r'$\pi/2$'),
                          'z':None}
            if extruded:
                names['z'] = 'h'
                domain_extents['z'] = (0, data.variables['H'][0])

        elif domain == 'box':
            names = {'x':'x', 'y':'y', 'z':'z'}
            domain_extents = {'x': (0, data.variables['Lx'][0]),
                              'y': (0, data.variables['Ly'][0]),
                              'z': (0, data.variables['H'][0])}
            ticklabels = {'x': None, 'y':None, 'z':None}

        else:
            raise ValueError('Domain %s not recognised' % domain)

    return names, domain_extents, ticklabels


def get_lfric_domain_properties(data, extrusion_details, central_lon=0.0):

    domain = extrusion_details['domain']

    #--------------------------------------------------------------------------#
    # Work out which coordinates we are using
    #--------------------------------------------------------------------------#

    if domain == 'plane':
        names = {'x':'x', 'y':'y'}

        # Find appropriate names for mesh entities
        # Needing to support different versions of LFRic data
        mesh_entity_names = [('edge_edge','edge_node'),('edge','node')]
        mesh_entity_found = False
        for mesh_entity_pair in mesh_entity_names:
            if (f'Mesh2d_{mesh_entity_pair[0]}_x' in data.variables.keys() and
                f'Mesh2d_{mesh_entity_pair[1]}_x' in data.variables.keys()):
                entity = mesh_entity_pair
                mesh_entity_found = True
                break
        if not mesh_entity_found:
            raise KeyError('Mesh entities for determining details of LFRic plane not found. '+
                           'Likely another mesh entity needs adding to get_lfric_domain_properties')

        domain_extents = {'x': (np.minimum(np.min(data[f'Mesh2d_{entity[0]}_x'][:]),
                                           np.min(data[f'Mesh2d_{entity[1]}_x'][:])),
                                np.maximum(np.max(data[f'Mesh2d_{entity[0]}_x'][:]),
                                           np.max(data[f'Mesh2d_{entity[1]}_x'][:]))),
                          'y': (np.minimum(np.min(data[f'Mesh2d_{entity[0]}_y'][:]),
                                           np.min(data[f'Mesh2d_{entity[1]}_y'][:])),
                                np.maximum(np.max(data[f'Mesh2d_{entity[0]}_y'][:]),
                                           np.max(data[f'Mesh2d_{entity[1]}_y'][:])))}
        ticklabels = {'x':None, 'y':None, 'z':None}

        if extrusion_details['extrusion'] is not None:
            names['z'] = 'z'
            domain_extents['z'] = (extrusion_details['zmin'],
                                   extrusion_details['zmax'])
        
    elif domain == 'sphere':
        names = {'x':'x', 'y':'y'}
        domain_extents = {'x': (-np.pi+central_lon, np.pi+central_lon),
                          'y': (-np.pi/2, np.pi/2)}
        if abs(central_lon) < 1e-6:
            xticklabels = (r'$-\pi$', r'$\pi$')
        elif abs(central_lon - np.pi) < 1e-6:
            xticklabels = (r'0', r'$2\pi$')
        elif abs(central_lon + np.pi) < 1e-6:
            xticklabels = (r'$-2\pi$', r'$0$')
        elif abs(central_lon + np.pi/2) < 1e-6:
            xticklabels = (r'$-3\pi/2$', r'$\pi/2$')
        elif abs(central_lon - np.pi/2) < 1e-6:
            xticklabels = (r'$-\pi/2$', r'$3\pi/2$')
        else:
            xticklabels = None

        ticklabels = {'x':xticklabels,
                      'y':(r'$-\pi/2$', r'$\pi/2$'),
                      'z':None}
        if extrusion_details['extrusion'] is not None:
            names['z'] = 'h'
            domain_extents['z'] = (extrusion_details['zmin'],
                                   extrusion_details['zmax'])

    else:
        raise ValueError('Domain %s not recognised' % domain)

    return names, domain_extents, ticklabels
