"""Routines relating to nicely placing colorbars."""

import matplotlib.pyplot as plt
from .tomplot_tools import work_out_cmap_extension
import numpy as np
import warnings

__all__ = ["add_colorbar_ax", "add_colorbar_fig"]


# This is a separate routine to "plot_contoured_field" so that:
# (a) colorbars can be easily applied to some (but not all) subplot axes
# (b) the number of arguments to "plot_contoured_field" can be minimised
# This routine motivates "plot_contoured_field" returning the `cf` object.
def add_colorbar_ax(ax, cf, cbar_label=None, cbar_format=None, cbar_ticks=None,
                    cbar_labelpad=None, extra_labelpad=0.0, location=None,
                    **colorbar_kwargs):
    """
    Adds a colour bar to a filled contour plot, for a particular axes.

    Args:
        ax (:class:`AxesSubplot`): the pyplot axes object to plot the field on.
            Note that this is currently unused.
        cf (`matplotlib.contour`): the filled contour object for the plot. This
            is used to obtain the colours to be displayed in the colorbar.
        cbar_label (str, optional): the label to be applied to the colorbar.
            Defaults to None.
        cbar_format (float, optional): the formatting to used for the ticklabels
            attached to the colorbar. Defaults to None.
        cbar_ticks (iter, optional): which ticks to be attached to the colorbar.
            Defaults to None.
        cbar_labelpad (float, optional): the padding to apply to the colorbar
            label. Defaults to None.
        extra_labelpad (float, optional): any extra label padding to be applied.
            If cbar_labelpad is not None, this is ignored. Defaults to 0.0.
        location (str, optional): the location, relative to the parent axes,
            where the colorbar axes is created. Defaults to None.
        **colorbar_kwargs: keyword arguments to be passed to plt.colorbar.
    """

    # ------------------------------------------------------------------------ #
    # Work out whether colorbar needs extending
    # ------------------------------------------------------------------------ #
    if 'extend' in colorbar_kwargs.keys():
        extend = colorbar_kwargs['extend']
        del colorbar_kwargs['extend']
    elif hasattr(cf, "levels"):
        extend = work_out_cmap_extension(cf.cmap, cf.levels)
    else:
        extend = 'neither'

    # ------------------------------------------------------------------------ #
    # Make colorbar
    # ------------------------------------------------------------------------ #
    cbar_ticks, cbar_format, cbar_format_str = tomplot_cbar_format(cf, cbar_ticks, cbar_format)
    cb = plt.colorbar(cf, format=cbar_format, ticks=cbar_ticks,
                      location=location, extend=extend, **colorbar_kwargs)
    if cbar_label is not None:
        cbar_labelpad = tomplot_cbar_labelpad(cbar_labelpad, extra_labelpad,
                                              cbar_format_str, location, cb)
        cb.set_label(cbar_label, labelpad=cbar_labelpad)


def add_colorbar_fig(fig, cf, cbar_label=None, location='right', ax_idxs=None,
                     cbar_format=None, cbar_ticks=None, cbar_labelpad=None,
                     extra_labelpad=0.0, cbar_padding=0.0, cbar_thickness=0.02,
                     **colorbar_kwargs):
    """
    Adds a colour bar to a filled contour plot, to the whole figure.

    Args:
        fig (:class:`AxesSubplot`): the pyplot axes object to plot the field on.
            Note that this is currently unused.
        cf (`matplotlib.contour`): the filled contour object for the plot. This
            is used to obtain the colours to be displayed in the colorbar.
        cbar_label (str, optional): the label to be applied to the colorbar.
            Defaults to None.
        location (str, optional): where the colorbar should be placed. Should
            be one of "right", "left", "top" or "bottom". Defaults to "right".
        ax_idxs (list of int, optional): a list of indices of the figure's axes
            to use for determining the position of the colorbar. Defaults to
            None, in which case all axes are used.
        cbar_format (float, optional): the formatting to used for the ticklabels
            attached to the colorbar. Defaults to None.
        cbar_ticks (iter, optional): which ticks to be attached to the colorbar.
            Defaults to None.
        cbar_labelpad (float, optional): the padding to apply to the colorbar
            label. Defaults to None.
        extra_labelpad (float, optional): any extra label padding to be applied.
            If cbar_labelpad is not None, this is ignored. Defaults to 0.0.
        cbar_padding (float, optional): the padding to apply to the colorbar
            relative to the array of axes. Defaults to 0.0.
        cbar_thickness (float, optional): how thick to make the colorbar.
            Defaults to 0.02.
        **colorbar_kwargs: keyword arguments to be passed to fig.colorbar.
        """

    # ------------------------------------------------------------------------ #
    # Work out axes for colorbar
    # ------------------------------------------------------------------------ #

    # Check position variable is valid
    assert location in ['right', 'left', 'top', 'bottom'], \
        f'add_colorbar_fig: invalid location variable {location}'

    # Move the subplots to make space for colorbar
    if location == 'right':
        fig.subplots_adjust(right=0.9, wspace=0.1)
    elif location == 'left':
        fig.subplots_adjust(left=0.2, wspace=0.1)
    elif location == 'top':
        fig.subplots_adjust(top=0.8, hspace=0.1)
    elif location == 'bottom':
        fig.subplots_adjust(bottom=0.2, hspace=0.1)

    if ax_idxs is not None:
        ax_to_use = [fig.get_axes()[idx] for idx in ax_idxs]
    else:
        ax_to_use = fig.get_axes()

    # Find corners of axes for whole plot
    min_x = min_y = 1.0  # Initial big values for minimum
    max_x = max_y = 0.0  # Initial small values for maximum
    for ax in ax_to_use:
        ax_coords = ax.get_position()
        min_x = np.minimum(min_x, ax_coords.xmin)
        max_x = np.maximum(max_x, ax_coords.xmax)
        min_y = np.minimum(min_y, ax_coords.ymin)
        max_y = np.maximum(max_y, ax_coords.ymax)

    # Set cbar a bit away from main axes
    if location in ['right', 'top']:
        cbar_delta = 0.025 + cbar_padding
    else:
        # Need to create more spaces for numbers
        cbar_delta = 0.06 + cbar_padding

    if location == 'right':
        orientation = 'vertical'
        cbar_coords = \
            [max_x + cbar_delta,                   # x-coord of bottom left corner
             min_y,                                # y-coord of bottom left corner
             cbar_thickness,                       # x-extent of cbar
             max_y - min_y]                        # y-extent of cbar
    elif location == 'left':
        orientation = 'vertical'
        cbar_coords = \
            [min_x - cbar_delta - cbar_thickness,  # x-coord of bottom left corner
             min_y,                                # y-coord of bottom left corner
             cbar_thickness,                       # x-extent of cbar
             max_y - min_y]                        # y-extent of cbar
    elif location == 'top':
        orientation = 'horizontal'
        cbar_coords = \
            [min_x,                                # x-coord of bottom left corner
             max_y + cbar_delta,                   # y-coord of bottom left corner
             max_x - min_x,                        # x-extent of cbar
             cbar_thickness]                       # y-extent of cbar
    elif location == 'bottom':
        orientation = 'horizontal'
        cbar_coords = \
            [min_x,                                # x-coord of bottom left corner
             min_y - cbar_delta - cbar_thickness,  # y-coord of bottom left corner
             max_x - min_x,                        # x-extent of cbar
             cbar_thickness]                       # y-extent of cbar

    # Add colorbar in its own axis
    cbar_ax = fig.add_axes(cbar_coords)

    # ------------------------------------------------------------------------ #
    # Work out whether colorbar needs extending
    # ------------------------------------------------------------------------ #
    if 'extend' in colorbar_kwargs.keys():
        extend = colorbar_kwargs['extend']
        del colorbar_kwargs['extend']
    elif hasattr(cf, "levels"):
        extend = work_out_cmap_extension(cf.cmap, cf.levels)
    else:
        extend = 'neither'

    # ------------------------------------------------------------------------ #
    # Make colorbar
    # ------------------------------------------------------------------------ #
    cbar_ticks, cbar_format, cbar_format_str = tomplot_cbar_format(cf, cbar_ticks, cbar_format)
    cb = fig.colorbar(cf, cax=cbar_ax, format=cbar_format, ticks=cbar_ticks,
                      orientation=orientation, ticklocation=location,
                      extend=extend, **colorbar_kwargs)

    if cbar_label is not None:
        cbar_labelpad = tomplot_cbar_labelpad(cbar_labelpad, extra_labelpad,
                                              cbar_format_str, location, cb)
        cb.set_label(cbar_label, labelpad=cbar_labelpad)


def tomplot_cbar_format(cf, cbar_ticks=None, cbar_format=None):
    """
    Determines the ticks and format to use for the colorbar. If these aren't
    specified by the user, the cbar_ticks will correspond to the highest and
    lowest contour values. The default formatting is to 3 significant figures.

    This currently does nothing when used with a scatter plot, as it is unclear
    how to determine the minimum and maximum data values.

    Args:
        cf (`matplotlib.contour`): the filled contour object for the plot. This
            is used to obtain the colours to be displayed in the colorbar.
        cbar_ticks (iter, optional): which ticks to be attached to the colorbar.
            Defaults to None.
        cbar_format (float, optional): the formatting to used for the ticklabels
            attached to the colorbar. Defaults to None.

    Returns:
        (list, str): the new cbar_ticks and the new cbar_format.
    """

    cbar_format_str = None

    if cbar_ticks is None:
        if hasattr(cf, "levels"):
            # Take the minimum and maximum contours
            cbar_ticks = [cf.levels[0], cf.levels[-1]]

        # Unclear how to find min and max values when using scatter plot,
        # so do nothing here

    if cbar_format is None and cbar_ticks is not None:
        min_val, max_val = np.min(cbar_ticks), np.max(cbar_ticks)
        # Use scientific notation if the minimum is less than -999 or
        # between -0.01 and 0.01
        min_sci = (min_val < -9999 or (max_val < 0.1 and min_val < 0.1 and min_val > -0.1))

        # Use scientific notation if the maximum is less than -999 or
        # between -0.01 and 0.01
        max_sci = (max_val > 9999 or (min_val > -0.1 and max_val < 0.1 and max_val > -0.1))

        if min_sci or max_sci:
            cbar_format_str = '.2e'
            cbar_format = "{x:.2e}"
        else:
            cbar_format_str = '.2f'
            cbar_format = "{x:.2f}"

    elif cbar_format is not None:
        cbar_format_str = cbar_format
        cbar_format = "{x:"+cbar_format+"}"

    return cbar_ticks, cbar_format, cbar_format_str


def tomplot_cbar_labelpad(cbar_labelpad, extra_labelpad, cbar_format_str,
                          location, colorbar):
    """
    Determines a reasonable padding for the colorbar label. This assumes that
    only the upper and lower limits of the colorbar have tick labels.

    Args:
        cbar_labelpad (float): manually specified labelpad to use.
        extra_labelpad (float): any extra label padding that can be added to the
            calculated padding.
        cbar_format_str (str): the formatting string for the colorbar ticklabels
        location (str): location of the colorbar.
        colorbar (`matplotlib.colorbar.Colorbar`): the colorbar to be labelled.
    """

    assert location in ['top', 'bottom', 'right', 'left', None], \
        f'location {location} not valid.'

    if cbar_format_str is None:
        cbar_format_str = '.2f'

    if cbar_labelpad is not None:
        return cbar_labelpad

    # Do nothing for horizontal colorbars
    if location in ['top', 'bottom']:
        return None

    # Do nothing if we have multiple ticks
    if len(colorbar.get_ticks()) > 2:
        return None

    # Translate cbar_format into number of digits
    digits = 0  # Initial value
    cbar_format_str_split = cbar_format_str.split('.')
    if len(cbar_format_str_split) != 2:
        warnings.warn('Unable to work out digits for cbar labelpadding. '
                      + 'Falling back to default.')
        return None

    # Don't take part before decimal point
    cbar_format_str = cbar_format_str_split[1]

    for character in ['b', 'c', 'd', 'f', 'F', 's', '%']:
        if character in cbar_format_str:
            cbar_format_str = cbar_format_str.replace(character, '')
            # Probably a float, work out number of digits from cbar ticks
            max_tick_value = np.max(np.abs(colorbar.get_ticks()))
            if max_tick_value >= 10:
                digits = int(np.floor(np.log10(max_tick_value))) + 1
            else:
                digits = 1

    # Replace expontential identifiers with more digits
    for character in ['e', 'E', 'g', 'G']:
        if character in cbar_format_str:
            cbar_format_str = cbar_format_str.replace(character, '4')
            digits = 1

    for character in cbar_format_str:
        digits += int(character)

    # Assume we are using the default font size for the cbar ticks
    fontsize = plt.rcParams['font.size']

    # Multiply num digits by fontsize
    factor = 0.45  # empirically determined
    cbar_labelpad = -factor*digits*fontsize + extra_labelpad

    return cbar_labelpad
