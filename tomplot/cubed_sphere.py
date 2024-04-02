"""Routines relating to the cubed-sphere grid."""

import numpy as np

__all__ = ["lonlat_to_alphabeta", "alphabeta_to_lonlat",
           "plot_cubed_sphere_panels", "plot_cubed_sphere_slice"]


def lonlat_to_alphabeta(lon, lat):
    """
    Converts longitude, latitude values into alpha, beta values.

    Args:
        lon (float, `np.ndarray`): longitude value or array of longitude values.
        lat (float, `np.ndarray`): latitude value or array of latitude values.

    Returns:
        alpha, beta, panel: the equiangular cubed sphere coordinates (including
            panel_id) corresponding to the input coordinates.
    """

    rad_lon = np.pi*lon/180.
    rad_lat = np.pi*lat/180.

    # Compute Cartesian coordinates
    X = np.cos(rad_lon)*np.cos(rad_lat)
    Y = np.sin(rad_lon)*np.cos(rad_lat)
    Z = np.sin(rad_lat)

    # Work out which panel we are on
    panel_conditions = [np.logical_and(np.cos(rad_lon) >= 1.0 / np.sqrt(2),
                                       np.abs(np.sin(rad_lat)) <= np.cos(rad_lon)*np.cos(rad_lat)),
                        np.logical_and(np.sin(rad_lon) >= 1.0 / np.sqrt(2),
                                       np.abs(np.sin(rad_lat)) <= np.cos(rad_lon-np.pi/2)*np.cos(rad_lat)),
                        np.logical_and(np.cos(rad_lon) <= -1.0 / np.sqrt(2),
                                       np.abs(np.sin(rad_lat)) <= np.cos(rad_lon-np.pi)*np.cos(rad_lat)),
                        np.logical_and(np.sin(rad_lon) <= -1.0 / np.sqrt(2),
                                       np.abs(np.sin(rad_lat)) <= np.cos(rad_lon+np.pi/2)*np.cos(rad_lat)),
                        (np.sin(rad_lat) > 0.0),
                        (np.sin(rad_lat) < 0.0)]

    panel_ids = [1, 2, 3, 4, 5, 6]
    panel = np.select(panel_conditions, panel_ids)

    # Rotate -- these have to match the order of the panel_ids
    panel_rotation_conditions = [(panel == panel_id) for panel_id in panel_ids]
    panel_rotations = [(X, Y, Z), (Y, -X, Z), (-X, Z, Y),
                       (-Y, Z, -X), (Z, Y, -X), (-Z, -X, Y)]
    rot_X, rot_Y, rot_Z = np.select(panel_rotation_conditions, panel_rotations)

    alpha = np.arctan2(rot_Y, rot_X)
    beta = np.arctan2(rot_Z, rot_X)

    return alpha, beta, panel


def alphabeta_to_lonlat(alpha, beta, panel, rotated_pole=None):
    """
    Converts equiangular cubed sphere coordinates (alpha, beta and panel ID)
    into longitude, latitude values.

    Args:
        alpha (float, `np.ndarray`): alpha coord value or array of alpha values.
        beta (float, `np.ndarray`): beta coord value or array of beta values.
        panel (int, `np.ndarray`): ID of cubed sphere panel.
        rotated_pole (tuple, optional): the coordinates of the North pole
            relative to a standard reference frame. These are a tuple of floats,
            (longitude, latitude).

    Returns:
        lon, lat: the longitude, latitude coordinates corresponding to the input
            coordinates.
    """

    # Compute Cartesian coordinates as if we are on Panel 1
    varrho = np.sqrt(1 + np.tan(alpha)**2 + np.tan(beta)**2)
    X = 1 / varrho
    Y = np.tan(alpha) / varrho
    Z = np.tan(beta) / varrho

    if type(panel) is int:
        panel = np.ones(np.shape(alpha), dtype=int)*panel

    # Rotate -- these have to match the order of the panel_ids
    panel_ids = [1, 2, 3, 4, 5, 6]
    panel_rotation_conditions = [(panel == panel_id) for panel_id in panel_ids]
    panel_rotations = [(X, Y, Z), (-Y, X, Z), (-X, Z, Y),
                       (-Z, -X, Y), (-Z, Y, X), (-Y, Z, -X)]
    rot_X, rot_Y, rot_Z = np.select(panel_rotation_conditions, panel_rotations)

    if rotated_pole is not None:
        lon, lat, _ = rotated_lonlatr_coords((rot_X, rot_Y, rot_Z), rotated_pole)
    else:
        # Compute longitude-latitude from global Cartesian coordinates
        lon = np.arctan2(rot_Y, rot_X)
        lat = np.arctan2(rot_Z, np.sqrt(rot_X**2 + rot_Y**2))

    deg_lon = 180.*lon/np.pi
    deg_lat = 180.*lat/np.pi

    return deg_lon, deg_lat


def plot_cubed_sphere_panels(ax, units='deg', color='black', linewidth=None,
                             projection=None, rotated_pole=None, lon_shift=None):
    """
    For plotting the cubed sphere panel edges on a horizontal plot.

    Args:
        ax (:class:`AxesSubplot`): the pyplot axes object used for the plot.
        units (str, optional): which units the underlying coordinates are in.
            Must be 'deg' or 'rad'. Defaults to 'deg'.
        color (str, optional): the colour of the lines used to plot the panel
            edges. Defaults to 'black'.
        linewidth (float, optional): the width of the lines used to plot the
            panel edges. Defaults to None.
        projection (:class:`Projection`, optional): a cartopy projection object.
            Defaults to None.
        rotated_pole (tuple, optional): the coordinates of the North pole
            relative to a standard reference frame. These are a tuple of floats,
            (longitude, latitude).
    """

    import cartopy.crs as ccrs

    if units not in ['deg', 'rad']:
        raise ValueError('Units for plotting cubed sphere panel edges should '
                         + f'be "deg" or "rad", not {units}')

    transform = ccrs.Geodetic()

    if units == 'deg':
        unit_factor = 1.0
    elif units == 'rad':
        unit_factor = np.pi/180.0

    vertices_alphabetapanel = [[(np.pi/4, np.pi/4, 1), (-np.pi/4, np.pi/4, 1)],
                               [(-np.pi/4, np.pi/4, 1), (-np.pi/4, -np.pi/4, 1)],
                               [(-np.pi/4, -np.pi/4, 1), (np.pi/4, -np.pi/4, 1)],
                               [(np.pi/4, -np.pi/4, 1), (np.pi/4, np.pi/4, 1)],
                               [(np.pi/4, np.pi/4, 2), (-np.pi/4, np.pi/4, 2)],
                               [(np.pi/4, -np.pi/4, 2), (np.pi/4, np.pi/4, 2)],
                               [(-np.pi/4, -np.pi/4, 2), (np.pi/4, -np.pi/4, 2)],
                               [(-np.pi/4, np.pi/4, 3), (-np.pi/4, -np.pi/4, 3)],
                               [(-np.pi/4, -np.pi/4, 3), (np.pi/4, -np.pi/4, 3)],
                               [(np.pi/4, -np.pi/4, 3), (np.pi/4, np.pi/4, 3)],
                               [(-np.pi/4, np.pi/4, 4), (-np.pi/4, -np.pi/4, 4)],
                               [(np.pi/4, -np.pi/4, 4), (np.pi/4, np.pi/4, 4)]]

    edges_alphabetapanel = [[np.linspace(edge[0][0], edge[1][0], 101),  # alpha values
                             np.linspace(edge[0][1], edge[1][1], 101),  # beta values
                             np.ones(101, dtype=int)*edge[0][2]]        # panels
                            for edge in vertices_alphabetapanel]

    edges_lonlat = [alphabeta_to_lonlat(edge[0], edge[1], edge[2],
                                        rotated_pole=rotated_pole)
                    for edge in edges_alphabetapanel]

    for edge in edges_lonlat:
        x_coords, y_coords = edge
        if lon_shift is not None:
            x_coords += lon_shift
        ax.plot(x_coords*unit_factor, y_coords*unit_factor, color=color,
                linewidth=linewidth, transform=transform)


def plot_cubed_sphere_slice(ax, alpha_coords, beta_coords, panel, units='deg',
                            color='black', linewidth=None, projection=None):
    """
    For plotting a line in the equiangular (alpha, beta) cubed sphere
    coordinates on a horizontal plane.

    Args:
        ax (:class:`AxesSubplot`): the pyplot axes object used for the plot.
        alpha_coords (`np.array`): array of alpha coordinates of the points
            along the edge.
        beta_coords (`np.array`): array
        panel (`np.array`): _description_
        units (str, optional): which units the underlying coordinates are in.
            Must be 'deg' or 'rad'. Defaults to 'deg'.
        color (str, optional): the colour of the line to plot. Defaults to
            'black'.
        linewidth (float, optional): the width of the line to plot. Defaults to
            None.
        projection (:class:`Projection`, optional): a cartopy projection object.
            Defaults to None.
    """

    import cartopy.crs as ccrs

    if units not in ['deg', 'rad']:
        raise ValueError('Units for plotting cubed sphere panel edges should '
                         + f'be "deg" or "rad", not {units}')

    transform = projection if projection is not None else ccrs.Geodetic()

    lon_coords, lat_coords = alphabeta_to_lonlat(alpha_coords, beta_coords, panel)

    ax.plot(lon_coords, lat_coords, color=color, linewidth=linewidth, transform=transform)


def rotated_lonlatr_coords(xyz, new_pole):
    """
    Returns the rotated (lon,lat,r) coordinates, given a rotation axis and the
    (X,Y,Z) coordinates.

    Args:
        xyz (tuple of :class:`np.ndarray`): Original geocentric Cartesian
            coordinates.
        new_pole (tuple): a tuple of floats (lon, lat) of the new pole, in the
            original coordinate system. The longitude and latitude must be
            expressed in radians.

        tuple of :class:`np.ndarray`: the rotated (lon,lat,r) coordinates.
    """

    rot_axis, rot_angle = pole_rotation(new_pole)

    # If numpy, shape (x,y,z) array as a vector
    old_xyz_vector = np.transpose(xyz)

    # Do rotations to get new coordinates
    new_xyz_vector = rodrigues_rotation(old_xyz_vector, rot_axis, rot_angle)

    new_xyz = np.transpose(new_xyz_vector)

    new_lonlatr = lonlatr_from_xyz(new_xyz[0], new_xyz[1], new_xyz[2])

    return new_lonlatr


def rodrigues_rotation(old_vector, rot_axis, rot_angle):
    u"""
    Performs the rotation of a vector v about some axis k by some angle ϕ, as
    given by Rodrigues' rotation formula:                                     \n

    v_rot = v * cos(ϕ) + (k cross v) sin(ϕ) + k (k . v)*(1 - cos(ϕ))          \n

    Returns a new vector. All components must be (x,y,z) components.

    Args:
        old_vector (:class:`np.ndarray`): the original vector or vector-valued
            field to be rotated, to be expressed in geocentric Cartesian (x,y,z)
            components in the original coordinate system.
        rot_axis (tuple): the vector representing the axis to rotate around,
            expressed in geocentric Cartesian (x,y,z) components (in the frame
            before the rotation).
        rot_angle (float): the angle to rotate by.

    Returns:
        :class:`np.ndarray`: the rotated vector or vector-valued field.
    """

    cos = np.cos
    sin = np.sin
    cross = np.cross

    # Numpy vector may need reshaping
    if np.shape(rot_axis) != np.shape(old_vector):
        # Construct shape for tiling vector
        tile_shape = [dim for dim in np.shape(old_vector)[:-1]]+[1]
        # Tile rot_axis vector to create an ndarray
        rot_axis = np.tile(rot_axis, tile_shape)

        # Replace dot product routine with something that does elementwise dot
        def dot(a, b):
            dotted_vectors = np.einsum('ij,ij->i', a, b)
            # Add new axis to allow further multiplication by a vector
            return dotted_vectors[:, np.newaxis]
    else:
        dot = np.dot

    new_vector = (old_vector * cos(rot_angle)
                  + cross(rot_axis, old_vector) * sin(rot_angle)
                  + rot_axis * dot(rot_axis, old_vector) * (1 - cos(rot_angle)))

    return new_vector


def pole_rotation(new_pole):
    """
    Computes the rotation axis and angle associated with rotating the pole from
    lon = 0 and lat = pi / 2 to a new longitude and latitude. Returns the
    rotation axis and angle for use in the Rodrigues rotation.

    Args:
        new_pole (tuple): a tuple of floats (lon, lat) of the new pole. The
            longitude and latitude must be expressed in radians.

    Returns:
        tuple: (rot_axis, rot_angle). This describes the rotation axis (a tuple
            or :class:`as_vector` of (x, y, z) components of the rotation axis,
            and a float describing the rotation angle.
    """

    # We assume that the old pole is old_lon_p = 0 and old_lat_p = pi / 2.
    old_lat_p = np.pi / 2

    # Then moving the pole to new_lon_p, new_lat_p is akin to rotating the pole
    # about lon_rot = new_lon + pi / 2, lat_rot = 0
    new_lon_p, new_lat_p = new_pole
    lon_rot = new_lon_p + np.pi / 2
    lat_rot = 0.0

    # The rotation angle is only in the latitudinal direction
    rot_angle = old_lat_p - new_lat_p

    # Turn rotation axis into a vector
    # it points in the radial direction and has a length of one
    rot_axis = xyz_vector_from_lonlatr(0, 0, 1, (lon_rot, lat_rot, 1),
                                       position_units='lonlatr_rad')

    return rot_axis, rot_angle


def lonlatr_from_xyz(x, y, z, angle_units='rad'):
    """
    Returns the spherical lon, lat and r coordinates from the global geocentric
    Cartesian x, y, z coordinates.

    Args:
        x (:class:`np.ndarray`): x-coordinate.
        y (:class:`np.ndarray`): y-coordinate.
        z (:class:`np.ndarray`): z-coordinate.
        angle_units (str, optional): the units to use for the angle. Valid
            options are 'rad' (radians) or 'deg' (degrees). Defaults to 'rad'.

    Returns:
        tuple of :class`np.ndarray`: the tuple of (lon, lat, r) coordinates in
            the appropriate form given the provided arguments.
    """

    if angle_units not in ['rad', 'deg']:
        raise ValueError(f'angle_units arg {angle_units} not valid')

    # Determine whether to use firedrake or numpy functions
    atan2 = np.arctan2
    sqrt = np.sqrt
    pi = np.pi

    if angle_units == 'deg':
        unit_factor = 180./pi
    if angle_units == 'rad':
        unit_factor = 1.0

    lon = atan2(y, x)
    r = sqrt(x**2 + y**2 + z**2)
    l = sqrt(x**2 + y**2)
    lat = atan2(z, l)

    return lon*unit_factor, lat*unit_factor, r

def xyz_vector_from_lonlatr(lon_component, lat_component, r_component,
                            position_vector, position_units="xyz"):
    """
    Returns the Cartesian geocentric x, y and z components of a vector from a
    vector whose components are in lon, lat and r spherical coordinates. If
    dealing with Firedrake, a vector expression is returned.

    Args:
        lon_component (:class:`np.ndarray`): the zonal component of the input
            vector.
        lat_component (:class:`np.ndarray`): the meridional component of the
            input vector.
        r_component (:class:`np.ndarray`): the radial component of the input
            vector.
        position_vector (:class:`np.ndarray`): the position vector, either as
            (x, y, z) or (lon, lat, radius) coordinates, subject to the
            `position_units` argument. Should match the shape of the input
            vector.
        position_units (str, optional): in which units the provided position
            vector is. Valid options are ["xyz", "lonlatr_rad", "lonlatr_deg"].
            Defaults to "xyz".

    Returns:
        :class:`np.ndarray`: (x, y, z) components of the input vector.
    """

    # Check position units argument is valid
    if position_units not in ["xyz", "lonlatr_rad", "lonlatr_deg"]:
        raise ValueError('xyz_vector_from_lonlatr: the `position_units` arg '
                         + 'must be one of "xyz", "lonlatr_rad", "lonlatr_deg "'
                         + f'but {position_units} was provided')

    # Determine whether to use firedrake or numpy functions
    pi = np.pi
    cos = np.cos
    sin = np.sin

    # Convert position to lonlatr_rad
    if position_units == 'xyz':
        lon, lat, _ = lonlatr_from_xyz(position_vector[0], position_vector[1],
                                       position_vector[2])
    elif position_units == 'lonlatr_rad':
        lon, lat, _ = position_vector
    elif position_units == 'lonlatr_deg':
        lon, lat = position_vector[0]*pi/180., position_vector[1]*pi/180.

    # f is our vector, e_i is the ith unit basis vector
    # f = f_r * e_r + f_lon * e_lon + f_lat * e_lat
    # We want f = f_x * e_x + f_y * e_y + f_z * e_z

    # f_x = dot(f, e_x)
    # e_x = cos(lon)*cos(lat) * e_r - sin(lon) * e_lon - cos(lon)*sin(lat) * e_lat
    x_component = (cos(lon)*cos(lat) * r_component
                   - sin(lon) * lon_component
                   - cos(lon)*sin(lat) * lat_component)

    # f_y = dot(f, e_y)
    # e_y = sin(lon)*cos(lat) * e_r + cos(lon) * e_lon - sin(lon)*sin(lat) * e_lat
    y_component = (sin(lon)*cos(lat) * r_component
                   + cos(lon) * lon_component
                   - sin(lon)*sin(lat) * lat_component)

    # f_z = dot(f, e_z)
    # e_z = sin(lat) * e_r + cos(lat) * e_lat
    z_component = (sin(lat) * r_component
                   + cos(lat) * lat_component)

    return (x_component, y_component, z_component)