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


def alphabeta_to_lonlat(alpha, beta, panel):
    """
    Converts equiangular cubed sphere coordinates (alpha, beta and panel ID)
    into longitude, latitude values.

    Args:
        alpha (float, `np.ndarray`): alpha coord value or array of alpha values.
        beta (float, `np.ndarray`): beta coord value or array of beta values.
        panel (int, `np.ndarray`): ID of cubed sphere panel.

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

    # Compute longitude-latitude from global Cartesian coordinates
    lon = np.arctan2(rot_Y, rot_X)
    lat = np.arctan2(rot_Z, np.sqrt(rot_X**2 + rot_Y**2))

    deg_lon = 180.*lon/np.pi
    deg_lat = 180.*lat/np.pi

    return deg_lon, deg_lat


def plot_cubed_sphere_panels(ax, units='deg', color='black', linewidth=None,
                             projection=None):
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

    edges_lonlat = [alphabeta_to_lonlat(edge[0], edge[1], edge[2]) for edge in edges_alphabetapanel]

    for edge in edges_lonlat:
        x_coords, y_coords = edge
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
