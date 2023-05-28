"""
This tests the function for making a 2D field plot with a removed contour
"""

import matplotlib.pyplot as plt
from tomplot import (automake_cmap, plot_contoured_field,
                     add_colorbar, automake_field_title)
import numpy as np
import pytest


def remove_contour_settings(situation, setup_func):

    tracer_max = 305.0
    tracer_background = 300.0

    if situation == "middle":
        num_bins = 10
        num_bins_over_2 = int(num_bins/2)
        ic_name = 'dipole'
        remove_contour = "middle"
        contour_max = tracer_max - tracer_background
        contour_min = -contour_max
        contours = np.linspace(contour_min, contour_max, num_bins+1)
        contour_step = (contour_max - contour_min) / num_bins
        answer_contours = [contour_min + contour_step*i for i in range(num_bins_over_2)] \
            + [contour_step*(1+i) for i in range(num_bins_over_2)]
        title = 'Remove middle contour'

    elif situation == "odd":
        num_bins = 7
        ic_name = 'two_gaussian'
        remove_contour = tracer_background
        contour_step = (tracer_max - tracer_background) / (num_bins - 2)
        contour_min = tracer_background - 2*contour_step
        contour_max = tracer_max
        contours = np.arange(contour_min, contour_max+contour_step, step=contour_step)
        answer_contours = [contour_min, contour_min+contour_step] + \
                          [tracer_background + (i+1)*contour_step for i in range(num_bins-2)]
        title = 'Remove contour with odd number of contours'

    elif situation == "even":
        num_bins = 6
        ic_name = 'two_gaussian'
        remove_contour = tracer_background
        contour_step = (tracer_max - tracer_background) / (num_bins - 1)
        contour_min = tracer_background - contour_step
        contour_max = tracer_max
        contours = np.arange(contour_min, contour_max+contour_step, step=contour_step)
        answer_contours = [contour_min] + \
                          [tracer_background + (i+1)*contour_step for i in range(num_bins-1)]
        title = 'Remove contour with even number of contours'

    elif situation == "odd_zero":
        num_bins = 20
        num_bins_over_2 = int(num_bins/2)
        ic_name = 'dipole'
        remove_contour = 0.0
        contour_max = tracer_max - tracer_background
        contour_min = -contour_max
        contour_step = (contour_max - contour_min) / num_bins
        contours = np.linspace(contour_min, contour_max, num_bins+1)
        answer_contours = [contour_min + contour_step*i for i in range(num_bins_over_2)] \
            + [contour_step*(1+i) for i in range(num_bins_over_2)]
        title = 'Remove zero contour with odd number of contours'

    elif situation == "even_zero":
        num_bins = 12
        tracer_max -= tracer_background
        tracer_background = 0.0
        ic_name = 'two_gaussian'
        remove_contour = 0.0
        contour_step = (tracer_max - tracer_background) / (num_bins - 2)
        contour_min = tracer_background - 2*contour_step
        contour_max = tracer_max
        contours = np.arange(contour_min, contour_max+contour_step, step=contour_step)
        answer_contours = [contour_min, contour_min+contour_step] + \
                          [tracer_background + (i+1)*contour_step for i in range(num_bins-2)]
        title = 'Remove zero contour with even number of contours'

    else:
        raise ValueError(f'Situation {situation} for removing contours not implemented')

    setup = setup_func(ic_name, tracer_background, tracer_max)

    return setup, contours, np.array(answer_contours), remove_contour, title


@pytest.mark.parametrize("situation", ["middle", "odd", "even", "odd_zero", "even_zero"])
def test_remove_contour(situation, plot_setup):

    tol = 1e-14

    # ------------------------------------------------------------------------ #
    # Settings for test
    # ------------------------------------------------------------------------ #

    setup, contours, answer_contours, remove_contour, title = \
        remove_contour_settings(situation, plot_setup)

    coords_X, coords_Y = setup.coords_X, setup.coords_Y
    field_data = setup.field_data
    colour_scheme = setup.colour_scheme

    # ------------------------------------------------------------------------ #
    # Remove contour
    # ------------------------------------------------------------------------ #

    cmap, line_contours = automake_cmap(contours, colour_scheme,
                                        remove_contour=remove_contour)

    assert (len(line_contours) == len(answer_contours)), 'Length of line ' + \
        f'contours is not correct for situation {situation}'
    diff_contours = line_contours - answer_contours
    assert np.all(diff_contours < tol), \
        f'Removal of line contour for situation {situation} is not correct'

    # ------------------------------------------------------------------------ #
    # Make plot
    # ------------------------------------------------------------------------ #

    plt.close()

    _, ax = plt.subplots(1, 1, figsize=(5, 5))

    cf, _ = plot_contoured_field(ax, coords_X, coords_Y, field_data,
                                 'contour', contours, cmap=cmap,
                                 line_contours=line_contours)

    add_colorbar(ax, cf, '')
    automake_field_title(ax, title)

    plot_name = f'remove_contour_{situation}.png'
    setup.make_plots(plot_name)
