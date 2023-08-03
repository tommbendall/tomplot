"""
This tests the routine to find nice contours for a data array.
"""

from tomplot import tomplot_contours
import numpy as np
import pytest

situations = ["single_digits", "negative_symmetric", "negative_asymmetric",
              "theta_type", "single_digits_to_teens", "below_and_above_one",
              "single_digits_to_twenties", "tricky_depth", "constant", "zero"]


@pytest.mark.parametrize("situation", situations)
def test_tomplot_contours(situation):

    # Test a series of raw data and check if we get nice rounded values
    tol = 1e-12

    # ------------------------------------------------------------------------ #
    # Specify raw values and answers
    # ------------------------------------------------------------------------ #

    if situation == "single_digits":
        raw_min, raw_max = 1.1, 8.6
        answer_min, answer_max = 1.0, 9.0
        answer_step = 0.5
        min_num_bins = 10
    elif situation == "negative_symmetric":
        raw_min, raw_max = -4.3, 3.5
        answer_min, answer_max = -5.0, 5.0
        answer_step = 1.0
        min_num_bins = 10
    elif situation == "negative_asymmetric":
        raw_min, raw_max = -0.7, 6.5
        answer_min, answer_max = -1.0, 7.0
        answer_step = 1.0
        min_num_bins = 5
    elif situation == "theta_type":
        raw_min, raw_max = 316.4, 389.3
        answer_min, answer_max = 310, 390
        answer_step = 5
        min_num_bins = 15
    elif situation == "single_digits_to_teens":
        raw_min, raw_max = 1.1, 18.6
        answer_min, answer_max = 0.0, 20.0
        answer_step = 2.5
        min_num_bins = 7
    elif situation == "single_digits_to_twenties":
        raw_min, raw_max = 1.1, 27.4
        answer_min, answer_max = 0.0, 30.0
        answer_step = 2.5
        min_num_bins = 10
    elif situation == "below_and_above_one":
        raw_min, raw_max = 0.42, 1.15
        answer_min, answer_max = 0.4, 1.2
        answer_step = 0.05
        min_num_bins = 10
    elif situation == 'tricky_depth':
        # Corresponds to depths in Galewsky jet
        raw_min, raw_max = 8850, 10170
        answer_min, answer_max = 8800, 10200
        answer_step = 100.
        min_num_bins = 12
    elif situation == 'constant':
        raw_min, raw_max = 0.58735, 0.58735
        answer_min, answer_max = 0.586, 0.588
        answer_step = 0.002 / 3
        min_num_bins = 3
    elif situation == 'zero':
        raw_min, raw_max = 0.0, 0.0
        answer_min, answer_max = -0.001, 0.001
        answer_step = 0.002 / 3
        min_num_bins = 3
    else:
        raise ValueError(f'test_rounded_limits: situation {situation} not implemented')

    # ------------------------------------------------------------------------ #
    # Specify raw values and answers
    # ------------------------------------------------------------------------ #
    data_array = np.array([raw_min, raw_max])
    answer_contours = np.arange(answer_min, answer_max+answer_step/2.0, step=answer_step)
    nice_contours = tomplot_contours(data_array, min_num_bins=min_num_bins)
    diff = nice_contours - answer_contours

    assert np.all(diff) < tol, 'Nice contours for situation ' \
        + f'{situation} are incorrect'
