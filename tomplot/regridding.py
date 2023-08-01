"""
Routines for regridding data.
"""

__all__ = ['regrid_horizontal_slice', 'regrid_regular_horizontal_slice',
           'regrid_vertical_slice']

import numpy as np
import pandas as pd
from scipy.interpolate import RectBivariateSpline, griddata


def regrid_horizontal_slice(new_coords_X, new_coords_Y,
                            old_coords_X, old_coords_Y, field_data,
                            method='combined_linear', periodic_fix=None):
    """
    Regrids a horizontal slice of field data to new horizontal coordinates.
    Uses the scipy `griddata` routine.

    Args:
        new_coords_X (`numpy.ndarray`): array of X coordinates of the new grid.
        new_coords_Y (`numpy.ndarray`): array of Y coordniates of the new grid.
        old_coords_X (`numpy.ndarray`): array of X coordinates of the old grid.
        old_coords_Y (`numpy.ndarray`): array of Y coordinates of the old grid.
        field_data (`numpy.ndarray`): the field data array on the old grid.
        method (str, optional): method to use for interpolating. Default is
            "combined_linear", which also performs a "nearest" interpolation to
            try to prevent NaN values. Valid options are "linear", "cubic",
            "nearest", "combined_linear" and "combined_cubic".
        periodic_fix (str, optional): whether to extend the original data to
            get good interpolation at periodic boundaries. For non-periodic
            meshes this may give strange results so should be avoided. Allowed
            values are currently "sphere" and None. Default is None.

        Returns:
        `numpy.ndarray`: the regridded field data.
        """

    methods = ['combined_linear', 'combined_cubic', 'linear', 'cubic', 'nearest']
    assert method in methods, f'regrid_horizontal_slice: method {method} not valid'

    # ------------------------------------------------------------------------ #
    # Periodic fix
    # ------------------------------------------------------------------------ #
    assert periodic_fix in ['sphere', None], 'regrid_horizontal_slice: ' \
        + f'periodic_fix can be "sphere" or None, not {periodic_fix}'

    if periodic_fix == 'sphere':
        # Extend mesh over longitude coordinate
        min_x, max_x = np.min(old_coords_X), np.max(old_coords_X)
        x_diff = max_x - min_x

        # x_diff will likely not be exactly 360 or 2*pi. Determine which to use
        # based on actual value
        if x_diff > 100:
            x_diff = 360.0
        else:
            x_diff = 2*np.pi

        # Put existing flattened data into a pandas Dataframe
        data_dict = {'X': old_coords_X.flatten(),
                     'Y': old_coords_Y.flatten(),
                     'field': field_data.flatten()}
        old_data = pd.DataFrame(data_dict)

        # Extend at minimum x
        min_x_data = old_data[old_data['X'] < min_x + 0.05*x_diff]
        min_x_data['X'].values[:] = min_x_data['X'].values[:] + x_diff

        # Extend at maximum x
        max_x_data = old_data[old_data['X'] > max_x - 0.05*x_diff]
        max_x_data['X'].values[:] = max_x_data['X'].values[:] - x_diff

        # Update original data frame
        old_data = pd.concat([old_data, min_x_data, max_x_data], ignore_index=True)

        # Pull out data values
        field_data = old_data['field'].values[:]
        old_coords_X = old_data['X'].values[:]
        old_coords_Y = old_data['Y'].values[:]

    # Bundle coordinates into a tuple
    new_coords = (new_coords_X, new_coords_Y)
    old_coords = (old_coords_X, old_coords_Y)

    # ------------------------------------------------------------------------ #
    # Interpolation
    # ------------------------------------------------------------------------ #
    # For a single method, just pass method to griddata
    if method in ['nearest', 'linear', 'cubic']:
        new_field_data = griddata(old_coords, field_data, new_coords, method=method)

    # For other methods, try to prevent NaNs by combining with "nearest" method
    elif method == 'combined_linear':
        nearest_field_data = griddata(old_coords, field_data, new_coords, method='nearest')
        new_field_data = griddata(old_coords, field_data, new_coords, method='linear')
        new_field_data[np.isnan(new_field_data)] = nearest_field_data[np.isnan(new_field_data)]

    elif method == 'combined_cubic':
        nearest_field_data = griddata(old_coords, field_data, new_coords, method='cubic')
        new_field_data = griddata(old_coords, field_data, new_coords, method='linear')
        new_field_data[np.isnan(new_field_data)] = nearest_field_data[np.isnan(new_field_data)]

    return new_field_data


def regrid_regular_horizontal_slice(new_coords_X, new_coords_Y,
                                    old_coords_X_1d, old_coords_Y_1d,
                                    field_data, spline_degree=3):
    """
    Regrids a horizontal slice of field data to new horizontal coordinates.
    Assumes that the old coordinates are 1D arrays, and uses the scipy
    RectBivariateSpline method. This is more efficient than other regridding
    so should be used where this is possible.

    Args:
        new_coords_X (`numpy.ndarray`): array of X coordinates of the new grid
        new_coords_Y (`numpy.ndarray`): array of Y coordniates of the new grid
        old_coords_X_1d (`numpy.ndarray`): 1D array of X coords for the old grid
        old_coords_Y_1d (`numpy.ndarray`): 1D array of Y coords for the old grid
        field_data (`numpy.ndarray`): the field data array on the old grid
        spline_degree (int, optional): the spline degree to pass to the
            interpolation object. Defaults to 3. The same degree is used for
            both directions.

    Returns:
        `numpy.ndarray`: the regridded field data.
    """

    interp_func = RectBivariateSpline(old_coords_X_1d, old_coords_Y_1d,
                                      field_data, kx=spline_degree,
                                      ky=spline_degree)
    new_field_data = interp_func(new_coords_X, new_coords_Y)

    return new_field_data


def regrid_vertical_slice(new_coords_1d, slice_along, slice_at,
                          old_coords_X, old_coords_Y, old_coords_Z,
                          field_data, method='combined_linear'):
    """
    Regrids field data onto a single vertical slice. This is done taking a
    series of 1D horizontal slices at each vertical level.

    Args:
        new_coords_1d (`numpy.ndarray`): 1D array of horizontal coordinates to
            use for the new grid.
        slice_along (str): which coordinate axis to take the slice along. Valid
            values are 'x', 'y', 'lon' or 'lat'.
        slice_at (float): the coordinate value at which to slice the data.
        old_coords_X (`numpy.ndarray`): array of X coordinates of the old grid.
        old_coords_Y (`numpy.ndarray`): array of Y coordinates of the old grid.
        old_coords_Z (`numpy.ndarray`): array of Y coordinates of the old grid.
        field_data (`numpy.ndarray`): the field data array on the old grid.
        method (str, optional): method to use for interpolating. Default is
            "combined_linear", which also performs a "nearest" interpolation to
            try to prevent NaN values. Valid options are "linear", "cubic",
            "nearest", "combined_linear" and "combined_cubic".

    Returns:
        tuple of `numpy.ndarray`: (new_field, new_coords_hori, new_coords_Z).
            The regridded data and its coordinates.
    """

    methods = ['combined_linear', 'combined_cubic', 'linear', 'cubic', 'nearest']
    assert method in methods, f'regrid_vertical_slice: method {method} not valid'

    assert slice_along in ['x', 'y', 'lon', 'lat'], \
        f'slice_along variable must correspond to a coordinate. {slice_along}' \
        + ' is not a valid value'

    # ------------------------------------------------------------------------ #
    # Make empty data arrays to be filled
    # ------------------------------------------------------------------------ #
    # Second dimension of data array is number of levels
    num_levels = np.shape(field_data)[1]
    coords_z_1d = old_coords_Z[0, :]
    new_coords_hori, new_coords_Z = np.meshgrid(new_coords_1d, coords_z_1d, indexing='ij')
    new_field_data = np.zeros_like(new_coords_hori)

    # ------------------------------------------------------------------------ #
    # Loop through levels, doing 1D interpolation
    # ------------------------------------------------------------------------ #
    for lev_idx in range(num_levels):
        # Make 1D mesh to interpolate to

        hori_coords_x, hori_coords_y = \
            np.meshgrid(new_coords_1d, np.array([slice_at]), indexing='ij')

        # Perform horizontal interpolation
        new_field_data_level = regrid_horizontal_slice(hori_coords_x,
                                                       hori_coords_y,
                                                       old_coords_X[:, 0],
                                                       old_coords_Y[:, 0],
                                                       field_data[:, lev_idx],
                                                       method=method)
        new_field_data[:, lev_idx] = new_field_data_level[:, 0]

    return new_field_data, new_coords_hori, new_coords_Z
