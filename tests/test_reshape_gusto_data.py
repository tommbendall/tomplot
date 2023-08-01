"""
This tests the routine to reshape 3D Gusto data.
"""

from tomplot import reshape_gusto_data
import numpy as np
import pytest

@pytest.mark.parametrize("domain", ["2D", "3D"])
def test_reshape_gusto_data(domain):

    # Test a series of raw data and check if we get nice rounded values
    tol = 1e-12

    # ------------------------------------------------------------------------ #
    # Specify raw values and answers
    # ------------------------------------------------------------------------ #

    # ------------------------------------------------------------------------ #
    # Specify raw values and answers
    # ------------------------------------------------------------------------ #

    assert np.all(diff) < tol, 'Nice contours for situation ' \
        + f'{situation} are incorrect'
