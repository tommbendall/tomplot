"""
This tests the routine to restrict Gusto data.
"""

from tomplot import (area_restriction, extract_gusto_field, 
                     extract_gusto_coords, reshape_gusto_data)
import numpy as np
import pytest
from os.path import abspath, dirname
from netCDF4 import Dataset

@pytest.mark.parametrize('limits', [{'X': (0, 90)}, 
                                    {'Y': (0, 90)}, 
                                    {'X': (0, 90), 'Y': (0, 90)}])
def test_area_restriction(limits):

    time_idx = 0
    level = 0
    results_dir = f'{abspath(dirname(__file__))}/data'
    field_name = 'theta'
    results_file_name = f'{results_dir}/gusto_sphere_3d_field_output.nc'
    data_file = Dataset(results_file_name, 'r')


    # ------------------------------------------------------------------------ #
    # Extract data
    # ------------------------------------------------------------------------ #
    coords_X_full, coords_Y_full, coords_Z_full = \
        extract_gusto_coords(data_file, field_name)
    field_full = extract_gusto_field(data_file, field_name, time_idx)
    # Reshape
    _, _, _, coords_Z_full = \
        reshape_gusto_data(field_full, coords_X_full,
                            coords_Y_full, coords_Z_full)
    _, coords_X ,coords_Y = area_restriction(field_full[:,level], 
                                                      coords_X_full[:,level],
                                                      coords_Y_full[:,level], 
                                                      limits) 
    data_file.close() 
    # ------------------------------------------------------------------------ #
    # Check the limits have been applied
    # ------------------------------------------------------------------------ #
    for key in limits:
        min, max = limits[key]
        if key == 'X':
            assert np.min(coords_X) == min
            assert np.max(coords_X) == max

        if key == 'Y':
            assert np.min(coords_Y) == min
            assert np.max(coords_Y) == max
