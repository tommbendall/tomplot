"""Routines relating to the icosahedral sphere grid."""

import numpy as np

__all__ = ["plot_icosahedral_sphere_panels"]


def plot_icosahedral_sphere_panels(ax, units='deg', color='black', linewidth=None):
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
    """

    import cartopy.crs as ccrs

    if units not in ['deg', 'rad']:
        raise ValueError('Units for plotting icosahedral sphere panel edges '
                         + f'should be "deg" or "rad", not {units}')

    transform = ccrs.Geodetic()

    if units == 'deg':
        unit_factor = 180.0/np.pi
    elif units == 'rad':
        unit_factor = 1.0

    # (x, y, z) coordinates of the vertices of the icosahedral sphere
    phi = (1 + np.sqrt(5)) / 2
    xyz_vertices = [
        (-1, phi, 0),
        (1, phi, 0),
        (-1, -phi, 0),
        (1, -phi, 0),
        (0, -1, phi),
        (0, 1, phi),
        (0, -1, -phi),
        (0, 1, -phi),
        (phi, 0, -1),
        (phi, 0, 1),
        (-phi, 0, -1),
        (-phi, 0, 1)
    ]

    # (lon, lat) coordinates of the vertices of the icosahedral sphere
    lonlat_vertices = [
        (np.arctan2(y, x), np.arctan2(z, np.sqrt(x**2 + y**2)))
        for x, y, z in xyz_vertices
    ]

    # Connectivity of the vertices
    edge_connections = [
        (0, 11),
        (5, 11),
        (0, 5),
        (0, 1),
        (1, 5),
        (1, 7),
        (0, 7),
        (0, 10),
        (10, 11),
        (1, 9),
        (5, 9),
        (4, 5),
        (4, 11),
        (2, 10),
        (2, 11),
        (6, 7),
        (6, 10),
        (7, 10),
        (1, 8),
        (7, 8),
        (3, 4),
        (3, 9),
        (4, 9),
        (2, 3),
        (2, 4),
        (2, 6),
        (3, 6),
        (6, 8),
        (3, 8),
        (8, 9)
    ]

    # Create points lying along a great arc between the vertices
    edges_lonlat = [
        great_circle_linspace(
            lonlat_vertices[i][0], lonlat_vertices[j][0],
            lonlat_vertices[i][1], lonlat_vertices[j][1], 21
        ) for i, j in edge_connections
    ]

    for edge in edges_lonlat:
        x_coords, y_coords = edge
        ax.plot(x_coords*unit_factor, y_coords*unit_factor, color=color,
                linewidth=linewidth, transform=transform)


def great_circle_linspace(lon1, lon2, lat1, lat2, num_points):
    """
    Generate points along a great circle path between two longitudes.

    Args:
        lon1 (float): Longitude of the first point in radians.
        lon2 (float): Longitude of the second point in radians.
        lat1 (float): Latitude of the first point in radians.
        lat2 (float): Latitude of the second point in radians.
        num_points (int): Number of points to generate.

    Returns:
        np.ndarray: Array of longitudes and latitudes along the great circle path.
    """

    # Convert to Cartesian coordinates
    x1 = np.cos(lat1) * np.cos(lon1)
    y1 = np.cos(lat1) * np.sin(lon1)
    z1 = np.sin(lat1)
    x2 = np.cos(lat2) * np.cos(lon2)
    y2 = np.cos(lat2) * np.sin(lon2)
    z2 = np.sin(lat2)

    # Compute points lying on vector directly between points
    x_points = np.linspace(x1, x2, num_points)
    y_points = np.linspace(y1, y2, num_points)
    z_points = np.linspace(z1, z2, num_points)

    # Convert back to spherical coordinates
    lon_points = np.arctan2(y_points, x_points)
    lat_points = np.arctan2(z_points, np.sqrt(x_points**2 + y_points**2))

    return lon_points, lat_points