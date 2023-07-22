"""
Routines for extracting data from Gusto and LFRic field data files.
"""

__all__ = ["extract_gusto_coords", "extract_lfric_coords"]

def extract_gusto_coords(dataset, field_name, slice_along=None):
    """
    Extracts the arrays of coordinate data for a specified field from a Gusto
    field data file.

    Args:
        dataset (:class:`Dataset`): the netCDF dataset containing the data.
        field_name (str): the name of the field to be plotted.
        slice_along (str, optional): a string specifying which direction to
            slice 3D data, and hence which coordinates are to be returned.
            Defaults to None.

    Returns:
        tuple of `numpy.ndarray`s: the coordinate data.
    """

    coord_variable = dataset[field_name]['field_values'].dimensions[0]
    # Variable name is "coords_*". We want to find the suffix
    coord_space = coord_variable[7:]

    # Work out which coordinates to expect based on the domain metadata
    domain = dataset.variables['domain_type'][:]
    if domain == 'spherical_shell':
        coords_X = dataset.variables[f'lon_{coord_space}'][:]
        coords_Y = dataset.variables[f'lat_{coord_space}'][:]
        return coords_X, coords_Y
    else:
        raise NotImplementedError(f'extract_gusto_coords: domain {domain} '
                                  +' either not implemented or recognised')


def extract_lfric_coords(dataset, field_name):
    """
    Extracts the arrays of coordinate data for a specified field from an LFRic
    field data file.

    Args:
        dataset (:class:`Dataset`): the netCDF dataset containing the data.
        field_name (str): the name of the field to be plotted.

    Returns:
        tuple of `numpy.ndarray`s: the coordinate data.
    """

