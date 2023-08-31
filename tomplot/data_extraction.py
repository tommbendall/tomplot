"""
Routines for extracting data from Gusto and LFRic field data files.
"""

import numpy as np
import pandas as pd
import warnings
from .cubed_sphere import lonlat_to_alphabeta

__all__ = ["extract_gusto_field", "extract_gusto_coords",
           "extract_lfric_field", "extract_lfric_coords",
           "extract_lfric_heights", "extract_lfric_vertical_slice",
           "extract_gusto_vertical_slice", "reshape_gusto_data"]


def extract_gusto_field(dataset, field_name, time_idx=None):
    """
    Extracts the data corresponding to a Gusto field.

    Args:
        dataset (:class:`Dataset`): the netCDF dataset containing the data.
        field_name (str): the name of the field to be extracted.
        time_idx (int, optional): index of the time point to extract at.
            Defaults to None, in which case all time points are extracted.
    """

    # If time_idx or level are None, we would index array as :
    # Make equivalent slice objects to these to simplify code below
    if time_idx is None:
        time_idx = slice(None, None)

    field_data = dataset[field_name]['field_values'][:, time_idx]

    return field_data


def extract_gusto_coords(dataset, field_name, units=None):
    """
    Extracts the arrays of coordinate data for a specified field from a Gusto
    field data file.

    Args:
        dataset (:class:`Dataset`): the netCDF dataset containing the data.
        field_name (str): the name of the field to be extracted.
        units (str, optional): a string specifying which units to return the
            coordinates in. Must match the specified domain. Defaults to None,
            in which case longitude/latitude coordinates are returned in degrees
            while Cartesian or radial coordinates are returned in metres. Valid
            options are 'deg', 'rad', 'm' or 'km'.

    Returns:
        tuple of `numpy.ndarray`s: the coordinate data.
    """
    # ------------------------------------------------------------------------ #
    # Checks on units argument and set default for this domain
    # ------------------------------------------------------------------------ #
    # Work out which units to return the coords in, based on domain metadata
    domain = dataset.variables['domain_type'][:]
    if domain == 'spherical_shell':
        assert units in [None, 'deg', 'rad'], 'extract_gusto_coords: for a ' \
            + f'spherical shell domain units must be "deg" or "rad" not {units}'
        if units is None:
            units = 'deg'
    elif domain in ['vertical_slice', 'plane', 'extruded_plane', 'interval']:
        assert units in [None, 'm', 'km'], 'extract_gusto_coords: for a ' \
            + f'{domain} domain units must be "m" or "km" not {units}'
        if units is None:
            units = 'km'
    elif domain == 'extruded_spherical_shell':
        assert units in [None, 'deg', 'rad'], 'extract_gusto_coords: for an ' \
            + f'extruded sphere, units must be "deg" or "rad" not {units}'
        if units is None:
            units = 'deg'
    else:
        raise NotImplementedError(f'extract_gusto_coords: domain {domain} '
                                  + ' either not implemented or recognised')

    # ------------------------------------------------------------------------ #
    # Set unit factors
    # ------------------------------------------------------------------------ #
    if units in ['m', 'rad']:
        unit_factor = 1
    elif units == 'km':
        unit_factor = 0.001
    elif units == 'deg':
        unit_factor = 180.0 / np.pi
    else:
        raise NotImplementedError(
            f'extract_gusto_coords: invalid units {units} arg')

    # ------------------------------------------------------------------------ #
    # Determine which coordinates to extract
    # ------------------------------------------------------------------------ #
    coord_variable = dataset[field_name]['field_values'].dimensions[0]
    # Variable name is "coords_*". We want to find the suffix
    coord_space = coord_variable[7:]

    # ------------------------------------------------------------------------ #
    # Work out which coordinates to expect based on the domain metadata
    # ------------------------------------------------------------------------ #
    if domain == 'spherical_shell':
        coords_X = dataset.variables[f'lon_{coord_space}'][:]*unit_factor
        coords_Y = dataset.variables[f'lat_{coord_space}'][:]*unit_factor
        return coords_X, coords_Y
    elif domain == 'vertical_slice':
        coords_X = dataset.variables[f'x_{coord_space}'][:]*unit_factor
        coords_Z = dataset.variables[f'z_{coord_space}'][:]*unit_factor
        return coords_X, coords_Z
    elif domain == 'interval':
        coords_X = dataset.variables[f'x_{coord_space}'][:]*unit_factor
        return coords_X
    elif domain == 'plane':
        coords_X = dataset.variables[f'x_{coord_space}'][:]*unit_factor
        coords_Y = dataset.variables[f'y_{coord_space}'][:]*unit_factor
        return coords_X, coords_Y
    elif domain == 'extruded_plane':
        coords_X = dataset.variables[f'x_{coord_space}'][:]*unit_factor
        coords_Y = dataset.variables[f'y_{coord_space}'][:]*unit_factor
        coords_Z = dataset.variables[f'z_{coord_space}'][:]*unit_factor
        return coords_X, coords_Y, coords_Z
    elif domain == 'extruded_spherical_shell':
        coords_X = dataset.variables[f'lon_{coord_space}'][:]*unit_factor
        coords_Y = dataset.variables[f'lat_{coord_space}'][:]*unit_factor
        try:
            # Support old data which outputted radius
            coords_Z = dataset.variables[f'r_{coord_space}'][:]  # No unit factor
            coords_Z -= np.min(coords_Z)  # Subtract radius
        except KeyError:
            coords_Z = dataset.variables[f'h_{coord_space}'][:]  # No unit factor
        return coords_X, coords_Y, coords_Z
    else:
        raise NotImplementedError(f'extract_gusto_coords: domain {domain} '
                                  + ' either not implemented or recognised')


def extract_lfric_field(dataset, field_name, time_idx=None, level=None):
    """
    Extracts the data corresponding to a LFRic field.

    Args:
        dataset (:class:`Dataset`): the netCDF dataset containing the data.
        field_name (str): the name of the field to be extracted.
        time_idx (int, optional): index of the time point to extract at.
            Defaults to None, in which case all time points are extracted.
        level (int, optional): index of the vertical level to extract at (if
            there is any). Defaults to None, in which case all of the data is
            extracted.
    """

    # If time_idx or level are None, we would index array as : and extract the
    # whole field - make equivalent slice objects to this to simplify code below
    if time_idx is None:
        time_idx = slice(None, None)
    if level is None:
        level = slice(None, None)

    # 3D data
    if len(dataset[field_name].dimensions) == 3:
        field_data = dataset[field_name][time_idx, level, :]

    # 2D time-varying field
    elif (len(dataset[field_name].dimensions) == 2
          and dataset[field_name].dimensions[0] == 'time'):
        field_data = dataset[field_name][time_idx, :]

    # 3D non-time-varying field
    elif len(dataset[field_name].dimensions) == 2:
        field_data = dataset[field_name][level, :]

    # 2D non-time-varying field
    elif len(dataset[field_name].dimensions) == 1:
        field_data = field_data = dataset[field_name][:]

    else:
        raise RuntimeError(
            'extract_lfric_field: unable to work out how to handle data with '
            + f'{len(dataset[field_name].dimensions)} dimensions')

    return field_data


def extract_lfric_coords(dataset, field_name, units=None):
    """
    Extracts the arrays of coordinate data for a specified field from an LFRic
    field data file.

    Args:
        dataset (:class:`Dataset`): the netCDF dataset containing the data.
        field_name (str): the name of the field to be plotted.
        units (str, optional): a string specifying which units to return the
            coordinates in. Must match the specified domain. Defaults to None,
            which will not alter the coordinates from the data file. Valid
            options are 'deg', 'rad', 'm' or 'km'.

    Returns:
        tuple of `numpy.ndarray`s: the coordinate data.
    """

    # Determine unit factor -- we don't know what the domain is so assume that
    # a valid option has been passed
    if units in [None, 'm', 'deg']:
        unit_factor = 1
    elif units == 'km':
        unit_factor = 0.001
    elif units == 'rad':
        unit_factor = np.pi / 180.0
    else:
        raise ValueError(f'extract_lfric_coords: invalid units {units} arg')

    # Get name of coordinate data, e.g. "nMesh2d_edge"
    root_coords_name = dataset[field_name].dimensions[-1]
    # Corresponding variable name is "Mesh2d_edge_x"
    coords_X_name = root_coords_name[1:]+'_x'
    coords_Y_name = root_coords_name[1:]+'_y'

    coords_X = dataset[coords_X_name][:]*unit_factor
    coords_Y = dataset[coords_Y_name][:]*unit_factor

    return coords_X, coords_Y


def extract_lfric_heights(height_dataset, field_dataset, field_name, level=None):
    """
    Extracts the arrays of height data for a specified field from an LFRic
    data file.

    Args:
        height_dataset (:class:`Dataset`): the netCDF dataset containing the
            height data to be extracted.
        field_dataset (:class:`Dataset`): the netCDF dataset containing the
            data of the field to be plotted.
        field_name (str): the name of the field to be plotted.
        level (int, optional): the vertical level to extract the heights at.
            Defaults to None, in which case the whole height field is extracted.

    Returns:
        `numpy.ndarray`: the height data.
    """

    # Work out which level this should be on
    if (field_dataset[field_name].dimensions[-2] == 'half_levels'
            and field_dataset[field_name].dimensions[-1] == 'nMesh2d_face'):
        height_name = 'height_w3'
    elif (field_dataset[field_name].dimensions[-2] == 'full_levels'
          and field_dataset[field_name].dimensions[-1] == 'nMesh2d_face'):
        height_name = 'height_wth'
    else:
        raise NotImplementedError(f'Dimensions for {field_name} are not '
                                  + 'implemented so cannot work out height field')

    # If time_idx or level are None, we would index array as : and extract the
    # whole field - make equivalent slice objects to this to simplify code below
    if level is None:
        level = slice(None, None)

    # Data may be time-varying or not, so these cases need handling separately
    if len(height_dataset[height_name].dimensions) == 2:
        # Initial data with no time value
        heights = height_dataset[height_name][level, :]
    elif len(height_dataset[height_name].dimensions) == 3:
        # 3D data -- take the first time value
        heights = height_dataset[height_name][0, level, :]
    else:
        raise RuntimeError(f'Trying to extract {height_name} field, but '
                           + 'the array is not 2D or 3D, so cannot proceed')

    return heights


def extract_lfric_vertical_slice(field_dataset, field_name, time_idx,
                                 slice_along, slice_at=0.0, height_dataset=None,
                                 levels=None, panel=None, tolerance=1e-4):
    """
    Extracts the field and coordinates for a vertical slice of LFRic data.

    Args:
        field_dataset (:class:`Dataset`): the netCDF dataset containing the
            data of the field to be extracted.
        field_name (str): the name of the field to be extracted.
        time_idx (int): the time index at which to extract the data. Can be None
            if not relevant.
        slice_along (str): which coordinate axis to take the slice along. Valid
            values are 'x', 'y', 'lon', 'lat', 'alpha' or 'beta'.
        slice_at (float, optional): the coordinate value at which to slice the
            data. Defaults to 0.0.
        height_dataset (:class:`Dataset`, optional): the netCDF dataset
            containing the height data to be extracted. If not specified, model
            level is used instead as a vertical coordinate. Defaults to None.
        levels (iter, optional): iterable of integers, indicating which model
            levels to extract the slice at. Defaults to None, in which case data
            is extracted at all levels.
        panel (int, optional): if "alpha" or "beta" are specified as the coord
            to slice along, then this argument must also be provided, in order
            to indicate which panel to slice along. Defaults to None.
        tolerance (float, optional): tolerance to use in determining how close
            data points are to the points to slice along. Defaults to 1e-4.
    """

    assert slice_along in ['x', 'y', 'lon', 'lat', 'alpha', 'beta'], \
        f'slice_along variable must correspond to a coordinate. {slice_along}' \
        + ' is not a valid value'

    if slice_along in ['alpha', 'beta']:
        assert panel is not None, 'extract_lfric_vertical_slice: if slicing ' \
            + 'along alpha or beta coordinates, panel argument must be provided'

    # Combine coord systems into general X/Y system to simplify later code
    if slice_along in ['x', 'lon', 'alpha']:
        local_slice_along = 'X'
    if slice_along in ['y', 'lat', 'beta']:
        local_slice_along = 'Y'
    local_X_coord = 'Y' if local_slice_along == 'X' else 'X'

    # ------------------------------------------------------------------------ #
    # Determine number of vertical levels for this field
    # ------------------------------------------------------------------------ #
    if levels is None:
        if (len(field_dataset[field_name].dimensions) == 2):
            num_levels = np.shape(field_dataset[field_name])[0]
        elif (len(field_dataset[field_name].dimensions) == 3):
            num_levels = np.shape(field_dataset[field_name])[1]
        else:
            raise RuntimeError('extract_lfric_vertical_slice: cannot work with '
                               + 'data of this shape')

        levels = range(num_levels)
    else:
        num_levels = len(levels)

    # ------------------------------------------------------------------------ #
    # Extract horizontal coordinates for this field
    # ------------------------------------------------------------------------ #
    coords_X, coords_Y = extract_lfric_coords(field_dataset, field_name)
    # Convert to alpha/beta if that is specified
    if slice_along in ['alpha', 'beta']:
        # Assume coords are already lon/lat
        coords_X, coords_Y, panel_ids = lonlat_to_alphabeta(coords_X, coords_Y)

    # ------------------------------------------------------------------------ #
    # Loop through levels and build up data
    # ------------------------------------------------------------------------ #
    for lev_idx, level in enumerate(levels):
        # Extract field data for this level
        field_data_level = extract_lfric_field(field_dataset, field_name,
                                               time_idx=time_idx, level=level)

        # Make dictionary of data, to pass to pandas DataFrame
        data_dict = {'field_data': field_data_level,
                     'X': coords_X, 'Y': coords_Y}

        # Add height data, if specified
        if height_dataset is not None:
            height_data_level = extract_lfric_heights(
                height_dataset, field_dataset, field_name, level)
            data_dict['height'] = height_data_level

        # Add panel IDs, if necessary
        if slice_along in ['alpha', 'beta']:
            data_dict['panel'] = panel_ids

        # -------------------------------------------------------------------- #
        # Make data frame and filter it
        # -------------------------------------------------------------------- #
        df = pd.DataFrame(data_dict)

        # Are there values close to the specified "slice_at" value?
        slice_df = df[np.abs(df[local_slice_along] - slice_at) < tolerance]
        # Additionally filter based on panel
        if slice_along in ['alpha', 'beta']:
            slice_df = slice_df[slice_df['panel'] == panel]

        if len(slice_df) == 0:
            # If there aren't points, we shall take nearest two sets of points
            warnings.warn('extract_lfric_vertical_slice: No data points found '
                          + f'at {slice_along} = {slice_at}. Finding nearest '
                          + 'points instead. This implies that you will need to'
                          + ' regrid')
            # Try to get 5% of values
            max_coord = np.max(df[local_slice_along].values)
            min_coord = np.min(df[local_slice_along].values)
            tolerance = 0.05 * (max_coord - min_coord)
            slice_df = df[np.abs(df[local_slice_along] - slice_at) < tolerance]

            # Additionally filter based on panel
            if slice_along in ['alpha', 'beta']:
                slice_df = slice_df[slice_df['panel'] == panel]

        # -------------------------------------------------------------------- #
        # Add data to a larger array
        # -------------------------------------------------------------------- #
        # Create data arrays if we're doing the first level
        if lev_idx == 0:
            len_data = len(slice_df['field_data'].values)
            field_data = np.zeros((len_data, num_levels))
            coords_X_final = np.zeros((len_data, num_levels))
            coords_Y_final = np.zeros((len_data, num_levels))
            coords_Z_final = np.zeros((len_data, num_levels))

        # Populate data arrays
        field_data[:, lev_idx] = slice_df['field_data'].values
        coords_X_final[:, lev_idx] = slice_df[local_X_coord].values
        coords_Y_final[:, lev_idx] = slice_df[local_slice_along].values
        coords_Z_final[:, lev_idx] = level if height_dataset is None else slice_df['height'].values

    return field_data, coords_X_final, coords_Y_final, coords_Z_final


def extract_gusto_vertical_slice(field_dataset, field_name, time_idx,
                                 slice_along, slice_at=0.0,
                                 panel=None, tolerance=1e-4):
    """
    Extracts the field and coordinates for a vertical slice of Gusto data.

    Args:
        field_dataset (:class:`Dataset`): the netCDF dataset containing the
            data of the field to be extracted.
        field_name (str): the name of the field to be extracted.
        time_idx (int): the time index at which to extract the data. Can be None
            if not relevant.
        slice_along (str): which coordinate axis to take the slice along. Valid
            values are 'x', 'y', 'lon', 'lat', 'alpha' or 'beta'.
        slice_at (float, optional): the coordinate value at which to slice the
            data. Defaults to 0.0.
        panel (int, optional): if "alpha" or "beta" are specified as the coord
            to slice along, then this argument must also be provided, in order
            to indicate which panel to slice along. Defaults to None.
        tolerance (float, optional): tolerance to use in determining how close
            data points are to the points to slice along. Defaults to 1e-4.
    """

    assert slice_along in ['x', 'y', 'lon', 'lat', 'alpha', 'beta'], \
        f'slice_along variable must correspond to a coordinate. {slice_along}' \
        + ' is not a valid value'

    if slice_along in ['alpha', 'beta']:
        assert panel is not None, 'extract_gusto_vertical_slice: if slicing ' \
            + 'along alpha or beta coordinates, panel argument must be provided'

    # Combine coord systems into general X/Y system to simplify later code
    if slice_along in ['x', 'lon', 'alpha']:
        local_slice_along = 'X'
    if slice_along in ['y', 'lat', 'beta']:
        local_slice_along = 'Y'
    local_X_coord = 'Y' if local_slice_along == 'X' else 'X'

    # ------------------------------------------------------------------------ #
    # Extract horizontal coordinates for this field
    # ------------------------------------------------------------------------ #
    all_coords = extract_gusto_coords(field_dataset, field_name)
    # Is data 2D or 3D?
    if len(all_coords) == 2:
        raise NotImplementedError('Taking a vertical slice of Gusto data is '
                                  + 'currently not supported. However 2D data '
                                  + 'can just be plotted directly')
    elif len(all_coords) == 3:
        coords_X, coords_Y, coords_Z = all_coords

        # Convert to alpha/beta if that is specified
        if slice_along in ['alpha', 'beta']:
            # Assume coords are already lon/lat
            coords_X, coords_Y, panel_ids = lonlat_to_alphabeta(coords_X, coords_Y)
        else:
            panel_ids = None

    # ------------------------------------------------------------------------ #
    # Extract whole field and reshape it to
    # ------------------------------------------------------------------------ #

    field_data = extract_gusto_field(field_dataset, field_name, time_idx)

    if panel_ids is not None:
        field_data, coords_X, coords_Y, coords_Z, panel_ids = \
            reshape_gusto_data(field_data, coords_X, coords_Y, coords_Z,
                               other_arrays=[panel_ids])
    else:
        field_data, coords_X, coords_Y, coords_Z = \
            reshape_gusto_data(field_data, coords_X, coords_Y, coords_Z)

    # We can now determine number of levels in data
    num_levels = np.shape(field_data)[1]
    levels = range(num_levels)

    # ------------------------------------------------------------------------ #
    # Loop through levels and build up data
    # ------------------------------------------------------------------------ #
    for lev_idx, _ in enumerate(levels):
        # Make dictionary of data, to pass to pandas DataFrame
        data_dict = {'field_data': field_data[:, lev_idx],
                     'X': coords_X[:, lev_idx],
                     'Y': coords_Y[:, lev_idx],
                     'Z': coords_Z[:, lev_idx]}

        # Add panel IDs, if necessary
        if slice_along in ['alpha', 'beta']:
            data_dict['panel'] = panel_ids[:, lev_idx]

        # -------------------------------------------------------------------- #
        # Make data frame and filter it
        # -------------------------------------------------------------------- #

        df = pd.DataFrame(data_dict)

        # Are there values close to the specified "slice_at" value?
        slice_df = df[np.abs(df[local_slice_along] - slice_at) < tolerance]
        # Additionally filter based on panel
        if slice_along in ['alpha', 'beta']:
            slice_df = slice_df[slice_df['panel'] == panel]

        if len(slice_df) == 0:
            # If there aren't points, we shall take nearest two sets of points
            warnings.warn('extract_gusto_vertical_slice: No data points found '
                          + f'at {slice_along} = {slice_at}. Finding nearest '
                          + 'points instead. This implies that you will need to'
                          + ' regrid')
            # Try to get 5% of values
            max_coord = np.max(df[local_slice_along].values)
            min_coord = np.min(df[local_slice_along].values)
            tolerance = 0.05 * (max_coord - min_coord)
            slice_df = df[np.abs(df[local_slice_along] - slice_at) < tolerance]

            # Additionally filter based on panel
            if slice_along in ['alpha', 'beta']:
                slice_df = slice_df[slice_df['panel'] == panel]

        # -------------------------------------------------------------------- #
        # Add data to a larger array
        # -------------------------------------------------------------------- #
        # Create data arrays if we're doing the first level
        if lev_idx == 0:
            len_data = len(slice_df['field_data'].values)
            field_data_final = np.zeros((len_data, num_levels))
            coords_X_final = np.zeros((len_data, num_levels))
            coords_Y_final = np.zeros((len_data, num_levels))
            coords_Z_final = np.zeros((len_data, num_levels))

        # Populate data arrays
        field_data_final[:, lev_idx] = slice_df['field_data'].values
        coords_X_final[:, lev_idx] = slice_df[local_X_coord].values
        coords_Y_final[:, lev_idx] = slice_df[local_slice_along].values
        coords_Z_final[:, lev_idx] = slice_df['Z'].values

    return field_data_final, coords_X_final, coords_Y_final, coords_Z_final


def reshape_gusto_data(field_data, coords_X, coords_Y_or_Z, coords_Z_3d=None,
                       other_arrays=None):
    """
    Reshapes an unstructured data array to introduce vertical structure.

    Gusto data is entirely unstructured, including in the vertical direction. To
    take advantage of some other plotting routines, it can be helpful to reshape
    this to be unstructured in the horizontal but structured in the vertical.

    Args:
        field_data (`numpy.ndarray`): unstructured field data array.
        coords_X (`numpy.ndarray`): unstructured array of X coordinates.
        coords_Y_or_Z (`numpy.ndarray`): unstructured array of Y or Z coords.
            When the data is 2D (e.g. a vertical slice) this should be Z coords.
        coords_Z_3d (`numpy.ndarray`, optional): unstructured array of Z coords.
            Defaults to None, as this should not be specified for 2D data. For
            3D data, this should be specified.
        other_arrays (list, optional): iterable of other data arrays to be
            reshaped in the same way. Must be of the same shape as the field
            data. Defaults to None.

    Returns:
        tuple of `numpy.ndarrays`: (field_data, coords_X, coords_Y, coords_Z)
            that are now structured to have horizontal data in a separate
            dimension from vertical data. If the `other_arrays` argument is
            provided, this is also returned. If the data is 2D, this will be
            (field_data, coords_X, coords_Z)
    """

    if coords_Z_3d is None:
        data_is_3d = False
        coords_Y = None
        coords_Z = coords_Y_or_Z
    else:
        data_is_3d = True
        coords_Y = coords_Y_or_Z
        coords_Z = coords_Z_3d

    if other_arrays is None:
        other_arrays = []

    # ------------------------------------------------------------------------ #
    # Round data to ensure sorting in dataframe is OK
    # ------------------------------------------------------------------------ #

    # Work out digits to round to, based on number of points and range of coords
    num_points = np.size(coords_X)
    data_range = np.max(coords_X) - np.min(coords_X)
    digits = int(np.floor(-np.log10(data_range / num_points)) + 3)
    coords_X = coords_X.round(digits)

    if data_is_3d:
        data_range = np.max(coords_Y) - np.min(coords_Y)
        digits = int(np.floor(-np.log10(data_range / num_points)) + 3)
        coords_Y = coords_Y.round(digits)

    # ------------------------------------------------------------------------ #
    # Make data frame
    # ------------------------------------------------------------------------ #

    data_dict = {'field': field_data, 'X': coords_X, 'Z': coords_Z}
    if coords_Y is not None:
        data_dict['Y'] = coords_Y

    for idx, other_array in enumerate(other_arrays):
        assert np.shape(other_array) == np.shape(field_data), \
            'reshape_gusto_data: shape of other array must match the field'
        data_dict[f'other_{idx}'] = other_array

    # Put everything into a pandas dataframe
    data = pd.DataFrame(data_dict)

    # Sort array by X and Y coordinates
    if data_is_3d:
        data = data.sort_values(by=['X', 'Y', 'Z'])
        first_X, first_Y = data['X'].values[0], data['Y'].values[0]
        first_point = data[(np.isclose(data['X'], first_X))
                           & (np.isclose(data['Y'], first_Y))]

    else:
        data = data.sort_values(by=['X', 'Z'])
        first_X = data['X'].values[0]
        first_point = data[np.isclose(data['X'], first_X)]

    # Number of levels should correspond to the number of points with the first
    # coordinate values
    num_levels = len(first_point)
    assert len(data) % num_levels == 0, 'Unable to nicely divide data into levels'

    # ------------------------------------------------------------------------ #
    # Create new arrays to store structured data
    # ------------------------------------------------------------------------ #

    num_hori_points = int(len(data) / num_levels)
    field_data = np.zeros((num_hori_points, num_levels))
    coords_X = np.zeros((num_hori_points, num_levels))
    if data_is_3d:
        coords_Y = np.zeros((num_hori_points, num_levels))
    coords_Z = np.zeros((num_hori_points, num_levels))
    new_other_arrays = [np.zeros((num_hori_points, num_levels)) for _ in other_arrays]

    # ------------------------------------------------------------------------ #
    # Fill arrays, on the basis of the dataframe already being sorted
    # ------------------------------------------------------------------------ #

    for lev_idx in range(num_levels):
        data_slice = slice(lev_idx, num_hori_points*num_levels+lev_idx, num_levels)
        field_data[:, lev_idx] = data['field'].values[data_slice]
        coords_X[:, lev_idx] = data['X'].values[data_slice]
        if data_is_3d:
            coords_Y[:, lev_idx] = data['Y'].values[data_slice]
        coords_Z[:, lev_idx] = data['Z'].values[data_slice]
        for idx, _ in enumerate(other_arrays):
            new_other_arrays[idx][:, lev_idx] = data[f'other_{idx}'].values[data_slice]

    if len(new_other_arrays) > 0 and data_is_3d:
        return field_data, coords_X, coords_Y, coords_Z, new_other_arrays
    elif data_is_3d:
        return field_data, coords_X, coords_Y, coords_Z
    elif len(new_other_arrays) > 0:
        return field_data, coords_X, coords_Z, new_other_arrays
    else:
        return field_data, coords_X, coords_Z
