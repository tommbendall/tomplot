import pandas as pd
import numpy as np


def area_restriction(field_data, coords_X, coords_Y, X_lim=None, Y_lim=None):
    """
    A function that restricts a data field and returns the new restricted field 
    as well as its corresponding co-ordinates

    Args:
        field_data (':class:`numpy.ndarray`): The field to be filtered.
        coords_X (':class:`numpy.ndarray`): The field along the X axis.
        coords_Y (':class:`numpy.ndarray`): The field along the Y axis.
        X_lim: tuple containing the (lower, upper) limits for the X axis
        Y_lim: tuple containing the (lower, upper) limits for the Y axis

    Returns:
        new_field_data (':class:`numpy.ndarray`):  
        new_X_coords (':class:`numpy.ndarray`):
        new_Y_coords (':class:`numpy.ndarray`):
    """

    if len(np.shape(field_data)) != 1:
        raise ValueError('area_restriction: input data must be 1D to be filtered by pandas data frame')
    if len(np.shape(coords_X)) != 1:
        raise ValueError('area_restriction: input data must be 1D to be filtered by pandas data frame')
    if len(np.shape(coords_Y)) != 1:
        raise ValueError('area_restriction: input data must be 1D to be filtered by pandas data frame')

    data_dict = {'field': field_data, 'X': coords_X, 'Y': coords_Y}
    df = pd.DataFrame(data_dict)

    if X_lim is not None:
        X_min, X_max = X_lim
        df = df[(df["X"] >= X_min) & (df["X"] <= X_max)]

    if Y_lim is not None:
        Y_min, Y_max = Y_lim
        df = df[(df["Y"] >= Y_min) & (df["Y"] <= Y_max)]

    new_field_data = df['field']
    new_coords_X = df['X']
    new_coords_Y = df['Y']

    return (new_field_data, new_coords_X, new_coords_Y)
