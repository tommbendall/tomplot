"""
Tools provided by tomplot for making nice figures, that may be common for
multiple types of figure.
"""
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
from matplotlib.colors import ListedColormap

__all__ = ['set_tomplot_style', 'tomplot_field_title', 'tomplot_cmap',
           'tomplot_field_markersize', 'tomplot_contours', 'rounded_limits',
           'work_out_cmap_extension']


def set_tomplot_style(fontsize=16, family='serif'):
    """
    Set some initial matplotlib variables to use a nice tomplot style. By
    default this will use latex formatting.

    Args:
        fontsize (int, optional): default fontsize to use. Defaults to 48.
        family (str, optional): family of font to use. Defaults to "serif".
    """

    plt.rc('text', usetex=True)
    font_opts = {'size': fontsize, 'family': family}
    plt.rc('font', **font_opts)


def tomplot_field_title(ax, title, titlepad=None, return_title=False,
                        fontsize=None, minmax=False, minmax_format='.2f',
                        field_data=None):
    """
    Adds a title to a subplot, or generates a title to be returned depending on
    the `return_title` argument.

    If the `minmax` argument is set to True, the field's minimum and maximum
    values are included in the title.

    Args:
        ax (:class:`AxesSubplot`): the pyplot axes object to plot the field on.
        title (str): the title to add to the subplot, or the base of the title
            to add if other elements (such as mins/maxes) are being added.
        titlepad (float, optional): amount of padd. Defaults to None.
        return_title (bool, optional): whether to return the generated title, or
            the axes object. If True, does not add the title to the axes.
            Defaults to False.
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
        str or :class:`AxesSubplot`: either the generated field title, or the
            pyplot axes object which has been adjusted.
    """

    if minmax:
        if field_data is None:
            raise ValueError('If generating title using "minmax", '
                             + 'field data must be provided')
        field_min = np.min(field_data)
        field_max = np.max(field_data)

        if minmax_format == '.2f' and (-abs(field_min) > 0.01 or abs(field_max) < 0.01):
            minmax_format = '.2e'

        format_str = '{:'+minmax_format+'}'

    if title in [None, '']:
        full_title = ''
    elif not minmax:
        full_title = title
    else:
        full_title = f'{title}, min: {format_str.format(field_min)}, ' \
            + f'max: {format_str.format(field_max)}'

    if return_title:
        return full_title
    else:
        ax.set_title(full_title, pad=titlepad, fontsize=fontsize)

        return ax


def tomplot_cmap(contours, color_scheme='Blues',
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
            False. If True, sets "over" to yellow and "under" to magenta.

    Raises:
        ValueError: if `remove_contour` is 'middle' but the number of contours
            is not an odd integer.
        NotImplementedError: if `cmap_rescale_type` is not None and extend_cmap
            is True. These options are not currently compatible, as rescaling
            the cmap makes it difficult to work out if the cmap has been
            extended.

    Returns:
        `matplotlib.Colormap`: the color map to be generated.
        iter: a list or numpy array of contour lines, taking into account any
            removed contour.
    """

    if cmap_rescale_type not in [None, 'both', 'top', 'bottom']:
        raise ValueError(f'cmap_rescale_type {cmap_rescale_type} not recognised')

    if cmap_rescale_type is not None and extend_cmap:
        raise NotImplementedError('tomplot_cmap: cmap_rescale_type cannot be '
                                  + 'used with extend_cmap')

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

        pure_cmap = mpl.colormaps[color_scheme].resampled(pure_num_colour_levels)
        new_colours = pure_cmap(np.linspace(lower_colour, upper_colour),
                                actual_num_colour_levels)
        cmap = ListedColormap(new_colours)

        # Allow us to work out if this has been rescaled
        setattr(cmap, '_tomplot_rescaling', True)

    else:
        cmap = mpl.colormaps[color_scheme].resampled(len(contours)-1)
        setattr(cmap, '_tomplot_rescaling', False)

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

    if extend_cmap:
        cmap.set_under('magenta')
        cmap.set_over('yellow')

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


def work_out_cmap_extension(cmap, contours):
    """
    Determines, from a colormap and set of contours, whether the colormap has
    been extended.

    Args:
        cmap (`matplotlib.Colormap`): the color map to be inspected.
        contours (iter): a list or numpy array of the contour levels that
            correspond to the colormap.

    Returns:
        str: description of how the cmap has been extended. If it has not been
            extended, this returns 'neither'.
    """

    if cmap is None or type(cmap) is str:
        # This will not be an extended cmap
        return 'neither'

    num_colour_bins = len(contours) - 1
    extended_under = False
    extended_over = False
    # If the "below" color is very different from the bottom colour, it has
    # been extended
    tol = 0.01
    if np.any(np.array(cmap(-1)) - np.array(cmap(0))) > tol:
        extended_under = True
    if np.any(np.array(cmap(100*num_colour_bins)) - np.array(cmap(num_colour_bins-1))) > tol:
        extended_over = True

    if hasattr(cmap, '_tomplot_rescaling'):
        # Using extend_cmap is not currently implemented with rescaling the cmap
        if cmap._tomplot_rescaling:
            return 'neither'

    if extended_under and extended_over:
        return 'both'
    elif extended_under:
        return 'min'
    elif extended_over:
        return 'max'
    else:
        return 'neither'


def tomplot_field_markersize(data, marker_scaling=1.0, ax=None):
    """
    Generates a markersize to use when using the "scatter" method for plotting
    2D fields.

    Args:
        data (`numpy.array`): the data to be plotted.
        marker_scaling (float, optional): scaling to be applied to marker size,
            to allow manual adjustment. Defaults to 1.0.
        ax (`matplotlib.Axes`, optional): the axes on which the data will be
            plotted. Defaults to None.

    Returns:
        float: the markersize.
    """

    if ax is not None:
        # Scale points based on number of points and figure size
        fig = ax.figure
        figsize = fig.get_size_inches()
        # Determine the number of subplots
        if hasattr(ax, "get_subplotspec"):
            subplot_shape = (ax.get_subplotspec().rowspan.stop,
                             ax.get_subplotspec().colspan.stop)
        else:
            subplot_shape = (1, 1)  # Default, assume 1 axes

        if subplot_shape == (1, 1):
            subplot_size = figsize
        else:
            subplot_size = (figsize[0]/subplot_shape[0], figsize[1]/subplot_shape[1])

        # Scaling dependent on tightest direction
        if len(np.shape(data)) == 2:
            # 2D data, assume points in X and Y direction
            point_density = np.maximum(np.shape(data)[0] / subplot_size[0],
                                       np.shape(data)[1] / subplot_size[1])
        else:
            # 1D data, so take square root
            point_density = np.maximum(int(len(data)**0.5) / subplot_size[0],
                                       int(len(data)**0.5) / subplot_size[1])
    else:
        # Just do scaling based on number of points
        n_points = np.max(np.shape(data))
        point_density = int(n_points**0.5)

    # Now derive the value (don't return a value less than 1)
    return max((72 / point_density)**1.6, 1)*marker_scaling


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


def tomplot_contours(data, min_num_bins=10, divergent_flag=False):
    """
    Generates some nice rounded contour values from a dataset. This can be used
    to make a quick plot from uninspected data. Contours are always picked to
    be nice values: e.g. multiples of 1, 2, 2.5 or 5.

    Args:
        data (`numpy.array`): the data to find nice contours of.
        min_num_bins (int, optional): the minimum number of bins to have (two
            contours corresponds to one bin). Defaults to 10.
        divergent_flag (boolean, optional): whether to enforce that a divergent
            profile is expected, to give results that are symmetric around 0.
            Defaults to False.

    Returns:
        `numpy.array`: a 1D array of points to use as contours.
    """

    # ------------------------------------------------------------------------ #
    # Round mins and maxes to nice numbers
    # ------------------------------------------------------------------------ #

    import math

    raw_max = np.amax(data)
    raw_min = np.amin(data)

    data_min = raw_min
    data_max = raw_max

    # Make initial discrete range definitely bigger than the real range
    discrete_diff = 3.0 * (raw_max - raw_min) + 1.0
    digit_inc = 1

    # ------------------------------------------------------------------------ #
    # Loop through rounding digits until we have nice rounded mins/maxes
    # ------------------------------------------------------------------------ #
    # If the real range if less than half of the discrete range,
    # round to the next digit
    while discrete_diff > 2.0*(raw_max - raw_min):

        if (raw_max - raw_min > 0):
            digits = math.floor(-np.log10(raw_max - raw_min)) + digit_inc
        else:
            digits = 1

        if (raw_min < 0 and raw_max > 0):
            # How symmetric are these around zero?
            max_abs = np.maximum(np.abs(raw_min), np.abs(raw_max))
            min_abs = np.minimum(np.abs(raw_min), np.abs(raw_max))
            if divergent_flag or min_abs > 0.5*max_abs:
                data_max = roundup(max_abs, digits=digits)
                data_min = - data_max
            else:
                data_max = roundup(raw_max, digits=digits)
                data_min = rounddown(raw_min, digits=digits)
        else:
            data_min = rounddown(raw_min, digits=digits)
            data_max = roundup(raw_max, digits=digits)

        discrete_diff = data_max - data_min
        digit_inc += 1  # Loop again to next digit if necessary

    # ------------------------------------------------------------------------ #
    # Find a nice step for contours
    # ------------------------------------------------------------------------ #
    # data_min and data_max are now decided, so find contours between them
    max_step = (data_max - data_min) / min_num_bins

    # To be safe, only make contours when we have non-zero data
    if max_step > 1e-32:
        step_digits = math.floor(-np.log10(max_step))
        # Find an even bigger maximum step by rounding max_step up
        max_max_step = roundup(max_step, digits=step_digits)

        # Our step should be a nice fraction of this
        # Loop through nice fractions to find first number with enough bins
        step_found = False
        for nice_fraction in [1.0, 0.5, 0.25, 0.2, 0.1]:
            step = nice_fraction*max_max_step
            # Num bins should be a nice integer number, but "round" to nearest
            # integer just to make sure
            num_bins = round((data_max - data_min) / step)
            if num_bins >= min_num_bins:
                step_found = True
                break

        if not step_found:
            raise RuntimeError('Unable to find nice contours')

        contours = np.linspace(data_min, data_max, num_bins+1)

    # Data is entirely zero, so return only two contours
    else:
        contours = np.array([data_min, data_max])

    return contours


def rounded_limits(data, divergent_flag=False):
    """
    Finds rounded min and max values to give nice limits. Can be used to create
    nice contours for previously uninspected data.

    Args:
        data (`numpy.array`): the data to find rounded mins and maxes of.
        divergent_flag (boolean, optional): whether to enforce that a divergent
            profile is expected, to give results that are symmetric around 0.
            Defaults to False.

    Returns:
        tuple: rounded (min and max) values of the data
    """

    nice_contours = tomplot_contours(data, divergent_flag=divergent_flag)
    data_min = nice_contours[0]
    data_max = nice_contours[-1]

    return data_min, data_max
