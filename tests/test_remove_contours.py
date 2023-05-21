"""
This tests the function for making a 2D field plot with a removed contour
"""

import matplotlib.pyplot as plt
from tomplot import (automake_cmap, plot_contoured_field,
                     add_colorbar, automake_field_title)
import numpy as np
import pytest
import sys


def dipole_initial_condition(x, y, Lx, Ly, tracer_background, tracer_max):
    r = np.sqrt(((x - Lx/2) / (Lx/5))**2 + ((y - Ly/2) / (Ly/5))**2)
    return 2*((x - Lx/2) / (Lx/5))*(tracer_max - tracer_background) * np.exp(-r**2)

def two_gaussian_initial_condition(x, y, Lx, Ly, tracer_background, tracer_max):
    r1 = np.sqrt(((x - Lx/2) / (Lx/5))**2 + ((y - Ly/3) / (Ly/5))**2)
    r2 = np.sqrt(((x - 3*Lx/4) / (Lx/10))**2 + ((y - 3*Ly/4) / (Ly/10))**2)
    return (tracer_background + (tracer_max - tracer_background) * np.exp(-r1**2) +
            (tracer_max - tracer_background) / 20.0 * np.exp(-r2**2))

def remove_contour_settings(situation):

    tracer_max = 305.0
    tracer_background = 300.0

    if situation == "middle":
        num_bins = 10
        num_bins_over_2 = int(num_bins/2)
        initial_condition = dipole_initial_condition
        remove_contour = "middle"
        contour_max = tracer_max - tracer_background
        contour_min = -contour_max
        contours = np.linspace(contour_min, contour_max, num_bins+1)
        contour_step = (contour_max - contour_min) / num_bins
        answer_contours = [contour_min + contour_step*i for i in range(num_bins_over_2)] \
                           + [contour_step*(1+i) for i in range(num_bins_over_2)]
        colour_scheme = 'RdBu_r'
        title = 'Remove middle contour'

    elif situation == "odd":
        num_bins = 7
        initial_condition = two_gaussian_initial_condition
        remove_contour = tracer_background
        contour_step = (tracer_max - tracer_background) / (num_bins - 2)
        contour_min = tracer_background - 2*contour_step
        contour_max = tracer_max
        contours = np.arange(contour_min, contour_max+contour_step, step=contour_step)
        answer_contours = [contour_min, contour_min+contour_step] + \
                          [tracer_background + (i+1)*contour_step for i in range(num_bins-2)]
        colour_scheme = 'Blues'
        title = 'Remove contour with odd number of contours'

    elif situation == "even":
        num_bins = 6
        initial_condition = two_gaussian_initial_condition
        remove_contour = tracer_background
        contour_step = (tracer_max - tracer_background) / (num_bins - 1)
        contour_min = tracer_background - contour_step
        contour_max = tracer_max
        contours = np.arange(contour_min, contour_max+contour_step, step=contour_step)
        answer_contours = [contour_min] + \
                          [tracer_background + (i+1)*contour_step for i in range(num_bins-1)]
        colour_scheme = 'Blues'
        title = 'Remove contour with even number of contours'

    elif situation == "odd_zero":
        num_bins = 20
        num_bins_over_2 = int(num_bins/2)
        initial_condition = dipole_initial_condition
        remove_contour = 0.0
        contour_max = tracer_max - tracer_background
        contour_min = -contour_max
        contour_step = (contour_max - contour_min) / num_bins
        contours = np.linspace(contour_min, contour_max, num_bins+1)
        answer_contours = [contour_min + contour_step*i for i in range(num_bins_over_2)] \
                           + [contour_step*(1+i) for i in range(num_bins_over_2)]
        colour_scheme = 'RdBu_r'
        title = 'Remove zero contour with odd number of contours'

    elif situation == "even_zero":
        num_bins = 12
        tracer_max -= tracer_background
        tracer_background = 0.0
        initial_condition = two_gaussian_initial_condition
        remove_contour = 0.0
        contour_step = (tracer_max - tracer_background) / (num_bins - 2)
        contour_min = tracer_background - 2*contour_step
        contour_max = tracer_max
        contours = np.arange(contour_min, contour_max+contour_step, step=contour_step)
        answer_contours = [contour_min, contour_min+contour_step] + \
                          [tracer_background + (i+1)*contour_step for i in range(num_bins-2)]
        colour_scheme = 'Blues'
        title = 'Remove zero contour with even number of contours'

    else:
        raise ValueError(f'Situation {situation} for removing contours not implemented')

    return contours, np.array(answer_contours), colour_scheme, initial_condition, \
        remove_contour, tracer_background, tracer_max, title


@pytest.mark.parametrize("situation", ["middle", "odd", "even", "odd_zero", "even_zero"])
def test_remove_contour(situation):

    tol = 1e-14

    # ------------------------------------------------------------------------ #
    # Settings for test
    # ------------------------------------------------------------------------ #

    contours, answer_contours, colour_scheme, initial_condition, remove_contour, \
        tracer_background, tracer_max, title = remove_contour_settings(situation)

    Lx = Ly = 10.0
    x_1d = np.linspace(0, Lx, 15)
    y_1d = np.linspace(0, Ly, 15)
    coords_X, coords_Y = np.meshgrid(x_1d, y_1d, indexing='ij')

    field_data = initial_condition(coords_X, coords_Y, Lx, Ly, tracer_background, tracer_max)

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

    # if '--show-plots' in sys.argv:
    plt.show()
