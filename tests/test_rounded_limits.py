"""
This tests the routine to find nice rounded limits from a data array.
"""

from tomplot import rounded_limits
import numpy as np
import pytest

situations = ["single_digits", "negative_symmetric", "negative_asymmetric",
              "theta_type", "single_digits_to_teens", "below_and_above_one",
              "single_digits_to_twenties"]

@pytest.mark.parametrize("situation", situations)
def test_rounded_limits(situation):

    # Test a series of raw data and check if we get nice rounded values
    tol = 1e-12

    # ------------------------------------------------------------------------ #
    # Specify raw values and answers
    # ------------------------------------------------------------------------ #

    if situation == "single_digits":
        raw_min, raw_max = 1.1, 8.6
        answer_min, answer_max = 1.0, 9.0
    elif situation == "negative_symmetric":
        raw_min, raw_max =  -4.3, 3.5
        answer_min, answer_max = -5.0, 5.0
    elif situation == "negative_asymmetric":
        raw_min, raw_max = -0.7, 6.5
        answer_min, answer_max = -1.0, 7.0
    elif situation == "theta_type":
        raw_min, raw_max =  316.4, 389.3
        answer_min, answer_max = 310, 390
    elif situation == "single_digits_to_teens":
        raw_min, raw_max =  1.1, 18.6
        answer_min, answer_max = 0.0, 20.0
    elif situation == "single_digits_to_twenties":
        raw_min, raw_max =  1.1, 27.4
        answer_min, answer_max = 0.0, 30.0
    elif situation == "below_and_above_one":
        raw_min, raw_max = 0.42, 1.15
        answer_min, answer_max = 0.4, 1.2
    else:
        raise ValueError(f'test_rounded_limits: situation {situation} not implemented')

    # ------------------------------------------------------------------------ #
    # Specify raw values and answers
    # ------------------------------------------------------------------------ #
    data_array = np.array([raw_max, raw_min])
    rounded_min, rounded_max = rounded_limits(data_array)

    assert abs(rounded_min - answer_min) < tol, 'Rounded min for situation ' \
        + f'{situation} is incorrect, got {rounded_min} but expected {answer_min}'
    assert abs(rounded_max - answer_max) < tol, 'Rounded max for situation ' \
        + f'{situation} is incorrect, got {rounded_max} but expected {answer_max}'