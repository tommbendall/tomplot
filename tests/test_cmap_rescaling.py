"""
This tests the rescaling of a colour map
"""

import matplotlib.pyplot as plt
from tomplot import (tomplot_cmap, plot_contoured_field,
                     add_colorbar_ax, tomplot_field_title)
import numpy as np
import pytest


def rescaling_settings(cmap_arrangement, cmap_rescale_type, setup_func):

    tracer_max = 305.0
    tracer_background = 300.0

    if cmap_arrangement in ['divergent', 'divergent_removed']:
        if cmap_arrangement == 'divergent':
            num_bins = 20
            remove_contour = None
            title = 'Rescaled divergent cmap'
            lower_true_cmap = -10
            upper_true_cmap = 30
        else:
            num_bins = 10
            remove_contour = 'middle'
            title = 'Rescaled divergent cmap with contour removed'
            lower_true_cmap = 5
            upper_true_cmap = 15

        equiv_remove_contour = remove_contour
        contour_max = tracer_max - tracer_background
        contour_min = -contour_max
        contours = np.linspace(contour_min, contour_max, num_bins+1)
        ic_name = 'dipole'
        cmap_rescaling = 0.5
        equiv_contours = np.linspace(2.0*contour_min, 2.0*contour_max, 2*num_bins+1)
        assert cmap_rescale_type == 'both', 'cmap_rescale_type ' + \
            f'{cmap_rescale_type} does not make sense for divergent cmap'

    elif cmap_arrangement == "linear_difficult":
        contours = np.linspace(-0.2, 1.2, 15)
        cmap_rescaling = 0.75
        remove_contour = 0.0
        equiv_contours = np.linspace(-0.2, -0.2 + cmap_rescaling*1.4, 15)
        equiv_remove_contour = equiv_contours[2]
        ic_name = 'two_gaussian'
        lower_true_cmap = 0
        upper_true_cmap = 14
        title = 'Difficult linear with removed contour'
        tracer_background = 0.0
        tracer_max = 1.0

    else:
        contour_max = tracer_max
        if cmap_arrangement == 'linear_odd':
            num_bins = 11
            contour_step = (tracer_max - tracer_background) / (num_bins - 1)
            contour_min = tracer_background - contour_step
            cmap_rescaling = 11.0/13.0
            if cmap_rescale_type == 'both':
                equiv_contours = np.arange(contour_min - contour_step,
                                           contour_max+2*contour_step, step=contour_step)
                title = 'Linear cmap with odd bins, scaled from top and bottom'
                lower_true_cmap = 1
                upper_true_cmap = num_bins
            elif cmap_rescale_type == 'bottom':
                equiv_contours = np.arange(contour_min - 2*contour_step,
                                           contour_max+contour_step, step=contour_step)
                title = 'Linear cmap with odd bins, scaled from bottom'
                lower_true_cmap = 2
                upper_true_cmap = num_bins+1
            elif cmap_rescale_type == 'top':
                equiv_contours = np.arange(contour_min, contour_max+3*contour_step,
                                           step=contour_step)
                title = 'Linear cmap with odd bins, scaled from top'
                lower_true_cmap = 0
                upper_true_cmap = num_bins-1

        elif cmap_arrangement == 'linear_even':
            num_bins = 12
            contour_step = (tracer_max - tracer_background) / (num_bins - 2)
            contour_min = tracer_background - 2.0*contour_step
            cmap_rescaling = 8.0/12.0
            if cmap_rescale_type == 'both':
                equiv_contours = np.arange(contour_min - 2*contour_step,
                                           contour_max+3*contour_step, step=contour_step)
                title = 'Linear cmap with even bins, scaled from top and bottom'
                lower_true_cmap = 2
                upper_true_cmap = num_bins+1
            elif cmap_rescale_type == 'bottom':
                equiv_contours = np.arange(contour_min - 4*contour_step,
                                           contour_max+contour_step, step=contour_step)
                title = 'Linear cmap with even bins, scaled from bottom'
                lower_true_cmap = 4
                upper_true_cmap = num_bins+3
            elif cmap_rescale_type == 'top':
                equiv_contours = np.arange(contour_min, contour_max+5*contour_step,
                                           step=contour_step)
                title = 'Linear cmap with even bins, scaled from top'
                lower_true_cmap = 0
                upper_true_cmap = num_bins-1

        # Other settings for linear cmaps
        contours = np.arange(contour_min, contour_max+contour_step, step=contour_step)
        ic_name = 'two_gaussian'
        remove_contour = None
        equiv_remove_contour = None

    setup = setup_func(ic_name, tracer_background, tracer_max)

    true_cmap, _ = tomplot_cmap(equiv_contours, setup.colour_scheme,
                                remove_contour=equiv_remove_contour)
    true_cmap_values = [true_cmap(i) for i in range(lower_true_cmap,
                                                    upper_true_cmap+1)]

    return setup, contours, title, true_cmap_values, cmap_rescaling, remove_contour


# Generate all combos of (linear_odd, linear_even) with (both, bottom, top)
def cmap_combos(cmap_arrangements, cmap_rescale_types):
    import itertools
    return list(itertools.product(cmap_arrangements, cmap_rescale_types))


@pytest.mark.parametrize("cmap_arrangement, cmap_rescale_type",
                         # All combinations of odd/even with both/bottom/top
                         cmap_combos(["linear_odd", "linear_even"],
                                     ["both", "bottom", "top"])
                         # Also add a diverging cmap (only rescaling "both" makes sense here)
                         + [("divergent", "both"), ("divergent_removed", "both"),
                            ("linear_difficult", "top")])
# @pytest.mark.parametrize("cmap_arrangement, cmap_rescale_type",
#                         [("linear_difficult", "top")])
def test_cmap_rescaling(cmap_arrangement, cmap_rescale_type, plot_setup):

    # Not currently bothered about getting exact colours, just general scaling
    tol = 0.12

    # ------------------------------------------------------------------------ #
    # Settings for test
    # ------------------------------------------------------------------------ #

    setup, contours, title, true_cmap_values, cmap_rescaling, remove_contour \
        = rescaling_settings(cmap_arrangement, cmap_rescale_type, plot_setup)

    coords_X, coords_Y = setup.coords_X, setup.coords_Y
    field_data = setup.field_data
    colour_scheme = setup.colour_scheme

    # ------------------------------------------------------------------------ #
    # Generate rescaled cmap
    # ------------------------------------------------------------------------ #

    cmap, line_contours = tomplot_cmap(contours, colour_scheme,
                                       cmap_rescale_type=cmap_rescale_type,
                                       cmap_rescaling=cmap_rescaling,
                                       remove_contour=remove_contour)

    # ------------------------------------------------------------------------ #
    # Compare with true cmap
    # ------------------------------------------------------------------------ #

    cmap_values = [cmap(i) for i in range(cmap.N)]
    colour_diff = np.array(true_cmap_values[0]) - np.array(cmap_values[0])
    assert np.all(colour_diff < tol), 'Colours at initial position differ, ' + \
        f'for {cmap_rescale_type} with {cmap_arrangement} cmap'
    colour_diff = np.array(true_cmap_values[-1]) - np.array(cmap_values[-1])
    assert np.all(colour_diff < tol), 'Colours at final position differ, ' + \
        f'for {cmap_rescale_type} with {cmap_arrangement} cmap'

    # ------------------------------------------------------------------------ #
    # Make plot
    # ------------------------------------------------------------------------ #

    plt.close()

    _, ax = plt.subplots(1, 1, figsize=(5, 5))

    cf, _ = plot_contoured_field(ax, coords_X, coords_Y, field_data,
                                 'contour', contours, cmap=cmap,
                                 line_contours=line_contours)

    add_colorbar_ax(ax, cf, '')
    tomplot_field_title(ax, title)

    plot_name = f'cmap_rescaling_{cmap_arrangement}_{cmap_rescale_type}.png'
    setup.make_plots(plot_name)
