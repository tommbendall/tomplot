"""
Tools provided by tomplot for making nice figures, that may be common for
multiple types of figure.
"""

import numpy as np
import matplotlib.cm as cm
from matplotlib.colors import ListedColormap

__all__ = ['automake_field_axis_labels', 'automake_field_title', 'automake_cmap',
           'automake_minmax']


def automake_field_axis_labels(ax, data_metadata):
    """
    Sets labels, ticks, ticklabels and limits for the axes for a field plot,
    based on the metadata that tomplot uses to describe the plotted field data.

    Args:
        ax (:class:`AxesSubplot`): the pyplot axes object to plot the field on.
        data_metadata (dict): a dictionary, whose keys are strings describing
            aspects of the axes to be set, and whose values are the values to
            set them to. Active keys are: "xlims", "ylims", "xlabel", "ylabel",
            "xticks", "yticks", "xticklabels", "yticklabels", "xlabelpad" and
            "ylabelpad".

    Returns:
        :class:`AxesSubplot`: the pyplot axes object which has been adjusted.
    """

    if 'xlims' in data_metadata.keys():
        ax.set_xlim(data_metadata['xlims'])
    if 'ylims' in data_metadata.keys():
        ax.set_ylim(data_metadata['ylims'])

    if 'xlabel' in data_metadata.keys():
        if 'xlabelpad' in data_metadata.keys():
            xlabelpad = data_metadata['xlabelpad']
        else:
            xlabelpad = None
        ax.set_xlabel(data_metadata['xlabel'], labelpad=xlabelpad)
    if 'ylabel' in data_metadata.keys():
        if 'ylabelpad' in data_metadata.keys():
            ylabelpad = data_metadata['ylabelpad']
        else:
            ylabelpad = None
        ax.set_ylabel(data_metadata['ylabel'], labelpad=ylabelpad)

    if 'xticks' in data_metadata.keys():
        ax.set_xticks(data_metadata['xticks'])
    if 'yticks' in data_metadata.keys():
        ax.set_yticks(data_metadata['yticks'])

    if 'xticklabels' in data_metadata.keys():
        ax.set_xticklabels(data_metadata['xticklabels'])
    if 'yticklabels' in data_metadata.keys():
        ax.set_yticklabels(data_metadata['yticklabels'])

    return ax


def automake_field_title(ax, title, titlepad=None, fontsize=None,
                         minmax=False, minmax_format='.2f', field_data=None):
    """
    Adds a title to a subplot.

    If the `minmax` argument is set to True, the field's minimum and maximum
    values are included in the title.

    Args:
        ax (:class:`AxesSubplot`): the pyplot axes object to plot the field on.
        title (str): the title to add to the subplot, or the base of the title
            to add if other elements (such as mins/maxes) are being added.
        titlepad (float, optional): amount of padd. Defaults to None.
        fontsize(float, optional): size of the text to use. Defaults to None.
        minmax(bool, optional): whether to append the field's min and max values
            to the title. Defaults to False.
        minmax_format (str, optional): format to use for the field's min and max
            values in the title. Defaults to ".2f". Examples would be ".2f" or
            "3.1e".
        field_data (numpy.ndarray, optional): the field data used to give the
            min and max values. Defaults to None, but must be specified if the
            `minmax` argument is set to True.

    Returns:
        :class:`AxesSubplot`: the pyplot axes object which has been adjusted.
    """

    if minmax:
        if field_data is None:
            raise ValueError('If generating title using "minmax", '
                             + 'field data must be provided')
        field_min = np.min(field_data)
        field_max = np.max(field_data)

        format_str = '{:'+minmax_format+'}'

    if title in [None, '']:
        full_title = ''
    elif not minmax:
        full_title = title
    else:
        full_title = f'{title}, min: {format_str.format(field_min)}, ' \
            + f'max: {format_str.format(field_max)}'

    ax.set_title(full_title, pad=titlepad, fontsize=fontsize)

    return ax


def automake_cmap(contours, color_scheme='Blues',
                  cmap_rescale_type=None, cmap_rescaling=0.8,
                  remove_contour=None, extend_cmap=False):
    """
    Generates a color map based on contours and a named color scheme. This
    provides options to linearly rescale the colors in the color map, and to
    blend neighbouring colors to "remove" a contour.

    Args:
        contours (iter): list or array of the edges of color bins to be used
            in generating the color map.
        color_scheme (str, optional): specifies which color scheme to use in
            making the color map. Defaults to `Blues`. For named maps, see
            https://matplotlib.org/stable/gallery/color/colormap_reference.html
        cmap_rescale_type (str, optional): whether to rescale the colors in the
            color map, and if so from which end of the color map. Allowed args
            are None, 'top', 'bottom' and 'both'. Defaults to None.
        cmap_rescaling (float, optional): the rescaling to apply to the color
            map. If `cmap_rescale_type` is 'both', then a tuple of two values
            can be provided here. Defaults to 0.8, which "shrinks" the color map
            towards the central colors.
        remove_contour (float/str, optional): whether to blend two color
            bins, thus "removing" a contour. Can be a float, corresponding to
            the value of contour to be removed, or the string 'middle', which
            removes the central contour. Defaults to None.
        extend_cmap (bool, optional): whether to allow the color map to be
            extended beyond the values specified in `contours`. Defaults to
            False.

    Raises:
        ValueError: if `remove_contour` is 'middle' but the number of contours
            is not an odd integer.

    Returns:
        `matplotlib.Colormap`: the color map to be generated.
        iter: a list or numpy array of contour lines, taking into account any
            removed contour.
    """

    if cmap_rescale_type not in [None, 'both', 'top', 'bottom']:
        raise ValueError(f'cmap_rescale_type {cmap_rescale_type} not recognised')

    # First reproduce the contours for the lines from the existing set of contours
    line_contours = contours.copy()

    # Scale the colour map to remove the most extreme colours
    if cmap_rescale_type is not None:
        if cmap_rescale_type == 'both':
            if type(cmap_rescaling) not in [tuple, list]:
                cmap_rescaling = [cmap_rescaling]*2
            lower_colour = 0.5 * (1.0 - cmap_rescaling[0])
            upper_colour = 0.5 * (1.0 + cmap_rescaling[1])
            avg_colour_scaling = 0.5*(cmap_rescaling[0]+cmap_rescaling[1])
        elif cmap_rescale_type == 'top':
            lower_colour = 0.0
            upper_colour = cmap_rescaling
            avg_colour_scaling = cmap_rescaling
        elif cmap_rescale_type == 'bottom':
            lower_colour = 1.0 - cmap_rescaling
            upper_colour = 1.0
            avg_colour_scaling = cmap_rescaling

        # Make a colour map for more contours than we'll use
        actual_num_colour_levels = len(contours)
        pure_num_colour_levels = np.ceil(actual_num_colour_levels/avg_colour_scaling)

        pure_cmap = cm.get_cmap(color_scheme, pure_num_colour_levels)
        new_colours = pure_cmap(np.linspace(lower_colour, upper_colour),
                                actual_num_colour_levels)
        cmap = ListedColormap(new_colours)

    else:
        cmap = cm.get_cmap(color_scheme, len(contours))

    # Remove a particular contour
    if remove_contour is not None:
        if remove_contour == 'middle':
            # Find middle contour
            # Only works for an odd number of lines!
            if len(contours) % 2 == 1:
                level_to_remove = int(np.floor((len(contours) - 1) / 2))
            else:
                raise ValueError('Can only use remove_contour method "middle" '
                                 + 'when there are an odd number of contour lines')
        elif isinstance(remove_contour, float) or isinstance(remove_contour, int):
            remove_contour = float(remove_contour)
            # Search through contours to find this specific contour
            contour_found = False
            for i, contour in enumerate(contours):
                if np.isclose(contour, remove_contour):
                    level_to_remove = i
                    contour_found = True
                    break

            if not contour_found:
                # If we get here then we have not found this contour
                raise ValueError('contour %.3f was not found' % remove_contour)

        else:
            raise ValueError('remove_contour %s not recognised' % remove_contour)

        # Remove contour line
        try:
            # For lists
            line_contours.pop(level_to_remove)
        except AttributeError:
            # For numpy arrays
            line_contours = np.delete(line_contours, level_to_remove)

        # Remove colour from colour map
        cmap = remove_colour(cmap, level_to_remove, len(contours))

    return cmap, line_contours


def remove_colour(old_cmap, level_to_remove, num_levels):
    """
    Replaces two colours in a colour map with a shared colour, in effect
    removing a contour from the colour map.

    Args:
        old_cmap (`matplotlib.Colormap`): a colour map object
        level_to_remove (int): the index of the contour to be removed
        num_levels(int): the total number of contour levels

    Returns:
        `matplotlib.Colormap`: a new colour map.
    """

    newcolours = old_cmap(np.linspace(0, 1.0, num_levels-1))
    merged_colour = old_cmap(level_to_remove/num_levels)
    newcolours[level_to_remove-1] = merged_colour
    newcolours[level_to_remove] = merged_colour
    new_cmap = ListedColormap(newcolours)

    return new_cmap


def roundup(number, digits=0):
    """
    Rounds a number up, to a specified precision.
    
    Args:
        number (float): the number to be rounded.
        digits (int): the digits up to which to perform the rounding.

    Returns:
        float: the rounded number.
    """
    import math

    n = 10**-digits
    return round(math.ceil(number / n) * n, digits)


def rounddown(number, digits=0):
    """
    Rounds a number down, to a specified precision.
    
    Args:
        number (float): the number to be rounded.
        digits (int): the digits up to which to perform the rounding.

    Returns:
        float: the rounded number.
    """
    import math

    n = 10**-digits
    return round(math.floor(number / n) * n, digits)


def automake_limits(data):
    """
    Finds rounded min and max values to give nice limits. Can be used to create
    nice contours for previously uninspected data.

    Args:
        data (`numpy.array`): 
    """
    import math

    raw_max = np.amax(data)
    raw_min = np.amin(data)
    if (raw_max - raw_min > 0):
        digits = math.floor(-np.log10(raw_max - raw_min))
    else:
        digits = 1
    if (raw_min < 0 and raw_max > 0):
        # Assume diverging profile, so make symmetrical around 0
        max_abs = np.maximum(np.abs(raw_min), np.abs(raw_max))
        data_max = roundup(max_abs, digits=digits)
        data_min = - data_max
    else:
        data_min = rounddown(raw_min, digits=digits)
        data_max = roundup(raw_max, digits=digits)

    # We may need to adjust this if we've made a colour bar that it too wide
    raw_diff = raw_max - raw_min
    col_diff = data_max - data_min
    if raw_diff < 0.5*col_diff:
        stored_data_min = data_min
        if (raw_min < 0 and raw_max > 0):
            bins = 8 if raw_diff < 0.25*col_diff else 4
        else:
            bins = 5 if raw_diff < 0.2*col_diff else 4
        for bin_edge in np.linspace(data_min, data_max, bins+1):
            if raw_min < bin_edge:
                break
            else:
                data_min = bin_edge
        for bin_edge in np.linspace(data_max, stored_data_min, bins+1):
            if raw_max > bin_edge:
                break
            else:
                data_max = bin_edge

    return data_min, data_max
