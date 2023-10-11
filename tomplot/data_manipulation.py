"""
Common routines for manipulating field data and the corresponding coordinates.
"""
import pandas as pd
import numpy as np


def area_restriction(field_data, coords_X, coords_Y, coord_lims):
    """
    A function that restricts a data field and returns the new restricted field
    as well as its corresponding co-ordinates.

    Args:
        field_data (:class:`numpy.ndarray`): The field to be filtered.
        coords_X (':class:`numpy.ndarray`): The field along the X axis.
        coords_Y (':class:`numpy.ndarray`): The field along the Y axis.
        coord_lims (dict): A dicitionary containing the
            coordinates as keys and the values are tuples containing the limits.

    Returns:
        new_field_data (':class:`numpy.ndarray`): The new restricted data.
        new_X_coords (':class:`numpy.ndarray`): The X coordinates for the restricted data
        new_Y_coords (':class:`numpy.ndarray`): The Y coordinates for the restricted data
    """

    if len(np.shape(field_data)) != 1:
        raise ValueError('area_restriction: input data must be 1D to be filtered by pandas data frame.')
    if len(np.shape(coords_X)) != 1:
        raise ValueError('area_restriction: input data must be 1D to be filtered by pandas data frame.')
    if len(np.shape(coords_Y)) != 1:
        raise ValueError('area_restriction: input data must be 1D to be filtered by pandas data frame.')
    if len(coord_lims) == 0:
        raise ValueError('area_restriction: limit dictionary is empty, please provide an axis and limits')

    data_dict = {'field': field_data, 'X': coords_X, 'Y': coords_Y}
    df = pd.DataFrame(data_dict)
    for key in coord_lims:
        if key not in ('X', 'Y'):
            raise KeyError('Key error, Please choose a valid axis: X, Y')
        min, max = coord_lims[key]
        df = df[(df[key] >= min) & (df[key] <= max)]

    new_field_data = df['field']
    new_coords_X = df['X']
    new_coords_Y = df['Y']

    return (new_field_data, new_coords_X, new_coords_Y)
