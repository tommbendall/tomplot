"""
This file provides some plotting tools
"""
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import ListedColormap
from netCDF4 import Dataset
import numpy as np


def specified_contours():
    """
    Returns a dictionary with details of contours for specific fields
    and test cases.
    """

    raise NotImplementedError('specified_contours is not yet implements')


def get_colour(testname, field, i):
    """
    Returns a string denoting the colour to be used in the plot. This can be
    a pre-specified testname/field combination, or instead will be indicated
    by the integer i.

    :arg testname: a string giving the test name
    :arg field:    the string giving the name of the field
    :arg i:        an integer
    """

    # Can store specific test cases and fields in this dictionary
    specific_colours = {'bob':{'rho':'black',
                               'theta':'red',
                               'v':'blue'}}

    preset_colours = ['black', 'red', 'blue', 'purple', 'green',
                      'orange', 'cyan', 'pink', 'brown', 'gray',
                      'brown', 'turquoise']

    if ((testname in specific_colours.keys()) and
        (field in specific_colours[testname].keys())):

        return specific_colours[testname][field]

    else:
        return preset_colours[i % len(preset_colours)]

def get_marker(testname, field, i):
    """
    Returns a string denoting the marker to be used in the plot. This can be
    a pre-specified testname/field combination, or instead will be indicated
    by the integer i.

    :arg testname: a string giving the test name
    :arg field:    the string giving the name of the field
    :arg i:        an integer
    """

    # Can store specific test cases and fields in this dictionary
    specific_markers = {'bob':{'rho':'+',
                               'theta':'o',
                               'v':'x'}}

    preset_markers = ['+', 'o', 'x', '^', 's',
                      'v', '*', 'D', '<', '>',
                      '1', 'X']

    if ((testname in specific_markers.keys()) and
        (field in specific_markers[testname].keys())):

        return specific_markers[testname][field]

    else:
        return preset_markers[i % len(preset_markers)]


def get_label(field):
    """
    Returns a label based on the field name.

    :arg field:    the string giving the name of the field
    """

    specific_labels = {'rho':    r'$\rho \ / $ kg m$^{-3}$',
                       'theta':  r'$\theta \ / $ K',
                       'mr_v':   r'$m_v \ / $ kg kg$^{-1}$',
                       'mr_cl':  r'$m_{cl} \ / $ kg kg$^{-1}$',
                       'mr_r':   r'$m_r \ / $ kg kg$^{-1}$',
                       'mr_sat': r'$m_{sat} \ / $ kg kg$^{-1}$',
                       'u_x':    r'$u \ / $ m s$^{-1}$',
                       'u_y':    r'$v \ / $ m s$^{-1}$',
                       'u_z':    r'$w \ / $ m s$^{-1}$',
                       'q':      r'$q \ / $ kg kg$^{-1}$',
                       'T':      r'$T \ / $ K',
                       'D':      r'$D \ / $ m',
                       'D_diff': r'$(D-D_{true}) \ / $ m'}

    if field in specific_labels.keys():
        return specific_labels[field]
    else:
        return field.replace('_', '\_')


def get_xlabel(variable, plot_type):
    """
    Returns an x-label for the axes

    :arg variable:  the x-variable for the plot
    :arg plot_type: 'convergence' or 'time_series' or 'slice'
    """

    if plot_type == 'convergence':
        if variable == 'dx':
            return r'$\ln(\Delta x)$'
        elif variable == 'dz':
            return r'$\ln(\Delta z)$'
        elif variable == 'dt':
            return r'$\ln(\Delta t)$'
        elif variable == 'rncells_per_dim':
            return r'$\ln\sqrt{1/n}$'
        else:
            return variable
    elif plot_type == 'grid parameters':
        if variable == 'dx':
            return r'$\Delta x$'
        elif variable == 'dz':
            return r'$\Delta z$'
        elif variable == 'dt':
            return r'$\Delta t$'
        else:
            return variable

    else:
        raise NotImplementedError('Other plot types have not been implemented yet')


def get_ylabel(variable, plot_type):
    """
    Returns an y-label for the axes

    :arg variable:  the y-variable for the plot
    :arg plot_type: 'convergence' or 'time_series' or 'slice'
    """

    if plot_type == 'convergence':
        if variable == 'L2_error':
            return r'$\ln(||q_{true}-q||)$'
        elif variable == 'L2_error_normalised':
            return r'$\ln(||q_{true}-q||/||q_{true}||)$'
        elif variable == 'dissipation_error':
            return r'$\ln($dissipation error$)$'
        elif variable == 'dispersion_error':
            return r'$\ln($dispersion error$)$'
        else:
            return variable
    elif plot_type == 'time series':
        if variable == 'L2':
            return r'$||q||_{L^2}$'
        elif variable == 'max':
            return r'$\max{q}$'
        elif variable == 'min':
            return r'$\min{q}$'
        elif variable == 'L2_error_normalised':
            return r'$||q-q_{true}||_{L_2}/||q_{true}||_{L_2}$'
        elif variable == 'L1_error_normalised':
            return r'$||q-q_{true}||_{L_1}/||q_{true}||_{L_1}$'
        elif variable == 'Linf_error_normalised':
            return r'$||q-q_{true}||_{L_\infty}/||q_{true}||_{L_\infty}$'
        elif variable == 'mean_normalised':
            return r'$ (\bar{q}-\bar{q}_{true}) / \bar{q}_0$'
        elif variable == 'variance_normalised':
            return r'$($Var$(q)-$Var$(q_{true}))/$Var$(q_0)$'
        elif variable == 'max_normalised':
            return r'$(\max q - \max q_{true})/(\max q_0 - \min q_0)$'
        elif variable == 'min_normalised':
            return r'$(\min q - \min q_{true})/(\max q_0 - \min q_0)$'

        else:
            return variable
    else:
        raise NotImplementedError('Other plot types have not been implemented yet')


def get_domain_label(direction, length_units='m', angular_units='rad'):
    """
    Returns label and units for axes that describe a spatial direction
    for plots of fields.

    :arg direction:     a string giving the direction.
    :arg length_units:  (Optional) which units to use for length dimensions.
                        Default is 'm' for metres.
    :arg angular_units: (Optional) which units to use for angular dimensions.
                        Default is 'rad' for radians.
    """

    if direction == 'x':
        return r'$x$', length_units
    elif direction == 'y':
        return r'$y$', length_units
    elif direction == 'z':
        return r'$z$', length_units
    elif direction == 'r':
        return r'$r$', length_units
    elif direction == 'h':
        return r'$h$', length_units
    elif direction == 'tau':
        return r'$\tau$', length_units
    elif direction == 'lon':
        return r'$\lambda$', angular_units
    elif direction == 'lat':
        return r'$\varphi$', angular_units
    elif direction == 'phi':
        return r'$\phi$', angular_units
    elif direction == 'sigma':
        return r'$\sigma$', angular_units
    else:
        raise ValueError('direction %s not recognised' % direction)


def colourmap_and_contours(contours, colour_scheme='Blues',
                           restricted_cmap=True,
                           remove_contour=False):
    """
    Makes the colour map and contour levels for a 2D field plot.

    :arg contours:        a list of the contour levels to use for the colour field.
    :arg colour_scheme:   (Optional) the string code to pass to pyplot for setting
                          up the colour scheme. Default is "Blues".
    :arg restricted_cmap: (Optional) Boolean describing whether to restrict the
                          natural pyplot colour cmap. If so, the colour map will
                          be rescaled to remove the most extreme colours.
    :arg remove_contour:  (Optional) argument for whether to remove a specific
                          contour. The colours touching this contour will then
                          be merged. This can be a float to remove a specific
                          contour value, an integer to remove a particular index
                          in the contour list, or "middle" (which only works when
                          there are an odd number of contours) to remove the
                          middle contour. Default is False.
    """

    # First reproduce the contours for the lines from the existing set of contours
    line_contours = contours.copy()

    # Scale the colour map to remove the most extreme colours
    if restricted_cmap:
        colour_levels_scaling = 1.2  # We scale colour map to avoid extreme colours

        # Make a colour map for more contours than we'll use
        actual_num_colour_levels = len(contours)
        pure_num_colour_levels = np.ceil(actual_num_colour_levels*colour_levels_scaling)

        pure_cmap = cm.get_cmap(colour_scheme, pure_num_colour_levels)
        new_colours = pure_cmap(np.linspace(0, 1/colour_levels_scaling),
                                actual_num_colour_levels)
        cmap = ListedColormap(new_colours)

    else:
        cmap = cm.get_cmap(colour_scheme, len(contours))

    # Remove a particular contour
    if remove_contour != False:
        if remove_contour == 'middle':
            # Find middle contour
            # Only works for an odd number of lines!
            if len(contours) % 2 == 1:
                level_to_remove = int(np.floor((len(contours) - 1) / 2))
            else:
                raise ValueError('Can only use remove_contour method "middle" '+
                                 'when there are an odd number of contour lines')
        elif isinstance(remove_contour, float):
            # Search through contours to find this specific contour
            contour_found = False
            for i, contour in enumerate(contours):
                if np.isclose(contour, remove_contour):
                    level_to_remove = i
                    contour_found = True
                    break

            if not contour_found:
                # If we get here then we have not found this contour
                raise ValueError('contour %.1f was not found' % remove_contour)

        elif isinstance(remove_contour, int):
            level_to_remove = remove_contour
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

    :arg cmap:            a colour map object
    :arg level_to_remove: the index of the contour to be removed
    :arg num_levels:      the total number of contour levels
    """

    newcolours = old_cmap(np.linspace(0, 1.0, num_levels-1))
    merged_colour = old_cmap(level_to_remove/num_levels)
    newcolours[level_to_remove-1] = merged_colour
    newcolours[level_to_remove] = merged_colour
    new_cmap = ListedColormap(newcolours)

    return new_cmap


def axes_limits_labels_and_titles(ax, xlabel=None, xlabelpad=None, xlims=None,
                                  xticks=None, xticklabels=None,
                                  ylabel=None, ylabelpad=None, ylims=None,
                                  yticks=None, yticklabels=None,
                                  title=None, title_method='full', titlepad=None,
                                  slice_label=None,
                                  time=None, time_method='seconds',
                                  field_min=None, field_max=None):
    """
    Decorates axes with supplied labels and limits, and adds a title.
    """
    if xlims is not None:
        ax.set_xlim(xlims)
    if ylims is not None:
        ax.set_ylim(ylims)
    if xticks is not None:
        ax.set_xticks(xticks)
    elif xlims is not None:
        ax.set_xticks(xlims)
    if xticklabels is not None:
        ax.set_xticklabels(xticklabels)
    if yticks is not None:
        ax.set_yticks(yticks)
    elif ylims is not None:
        ax.set_yticks(ylims)
    if yticklabels is not None:
        ax.set_yticklabels(yticklabels)
    if xlabel is not None:
        ax.set_xlabel(xlabel, labelpad=xlabelpad)
    if ylabel is not None:
        if ylabelpad is None and yticklabels is not None:
            if yticklabels[0] == '$-\\pi/2$':
                ylabelpad = -40
        ax.set_ylabel(ylabel, labelpad=ylabelpad)
    if title is not None:
        ax.set_title(title, pad=titlepad)
    elif (title_method is not None and title_method != 'none'):
        # We have a method for generating a title if one hasn't been provided
        if title_method == 'slice':
            title = slice_label
        elif title_method == 'time':
            title = 'time: '+get_time_string(time, time_method)
        elif (title_method == 'minmax' or (title_method == 'full' and slice_label is None and time is None)):
            title = 'min: %2.2e, max: %2.2e' % (field_min, field_max)
        elif title_method == 'full':
            if slice_label is None and time is None:
                if field_min is None:
                    title = ''
                else:
                    title = 'min: %2.2e, max: %2.2e' % (field_min, field_max)
            elif slice_label is None:
                if field_min is None:
                    title = 'time: '+get_time_string(time, time_method)
                else:
                    title = 'min: %2.2e, max: %2.2e, ' % (field_min, field_max)
                    title += 'time: '+get_time_string(time, time_method)
            elif time is None:
                if field_min is None:
                    title = slice_label
                else:
                    title = slice_label+', min: %2.2e, max: %2.2e' % (field_min, field_max)
            else:
                if field_min is None:
                    title = slice_label
                else:
                    title = slice_label+', min: %2.2e, max: %2.2e, ' % (field_min, field_max)
                title += 'time: '+get_time_string(time, time_method)
        else:
            raise ValueError('title_method %s not recognised' % title_method)
        ax.set_title(title, pad=titlepad)


def get_time_string(time, time_method):

    if time_method == 'seconds':
        return '%.2f s' % time
    elif time_method == 'days':
        return '%1.1f days' % (time / (24.*60.*60.))
    else:
        raise ValueError('time_method %s not recognised' % time_method)
