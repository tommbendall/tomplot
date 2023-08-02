"""
This tests the routine to reshape 3D Gusto data.
"""

from tomplot import reshape_gusto_data, extract_gusto_field, extract_gusto_coords
import numpy as np
import pytest
from os.path import abspath, dirname
from netCDF4 import Dataset


@pytest.mark.parametrize("domain", ["2D", "3D"])
def test_reshape_gusto_data(domain):

    time_idx = 0
    results_dir = f'{abspath(dirname(__file__))}/data'

    if domain == "2D":
        field_name = 'rho'
        results_file_name = f'{results_dir}/gusto_slice_2d_field_output.nc'
    else:
        field_name = 'theta'
        results_file_name = f'{results_dir}/gusto_sphere_3d_field_output.nc'

    data_file = Dataset(results_file_name, 'r')

    # ------------------------------------------------------------------------ #
    # Extract data
    # ------------------------------------------------------------------------ #

    field_full = extract_gusto_field(data_file, field_name, time_idx)

    if domain == "2D":
        coords_X_full, coords_Z_full = \
            extract_gusto_coords(data_file, field_name)

        # Reshape
        _, _, coords_Z_full = \
            reshape_gusto_data(field_full, coords_X_full, coords_Z_full)

    else:
        coords_X_full, coords_Y_full, coords_Z_full = \
            extract_gusto_coords(data_file, field_name)

        # Reshape
        _, _, _, coords_Z_full = \
            reshape_gusto_data(field_full, coords_X_full,
                               coords_Y_full, coords_Z_full)

    data_file.close()

    # ------------------------------------------------------------------------ #
    # Check that the data is now ordered
    # ------------------------------------------------------------------------ #

    first_heights = coords_Z_full[0, :]
    sorted_heights = np.sort(first_heights)

    assert np.all(np.isclose(first_heights, sorted_heights)), \
        f'reshaped {domain} gusto data does not appear to be reshaped correctly'
