"""
This file provides a generic routine for plotting a field on a 2D surface,
generally as (filled) contours or with coloured points as a scatter plot.

Some auxiliary routines are also provided.
"""
import matplotlib.pyplot as plt

__all__ = ["plot_contoured_field", "add_colorbar", "label_contour_lines"]

def plot_contoured_field(ax, coords_X, coords_Y, field_data, contours, method,
                         plot_filled_contours=True, plot_contour_lines=True,
                         projection=None,
                         # Options relating to filled contours
                         cmap=None, remove_lines=False, transparency=1.0,
                         # Options relating to line contours
                         contour_linestyles=None, contour_linewidths=1,
                         contour_linecolors='black'):
    """
    Plots a 2D field using (filled) contours or coloured points.

    This routine is for plotting fields on a 2D surface, acting as a wrapper to
    the matplotlib `contour`, `contourf`, `tricontour` and `tricontourf`
    routines. It also allows plotting of fields as points via a `scatter` plot.
    Which underlying routine is called depends upon the `method` argument, which
    can be "contour", "tricontour" or "scatter".

    Filled (coloured) contours and contour lines will be plotted by default.
    Different Cartopy projections can also be passed through, for instance when
    wanting to plot on the surface of a sphere.

    Args:
        ax (:class:`AxesSubplot`): the pyplot axes object to plot the field on.
        coords_X (`numpy.ndarray`): an array containing the coordinates of the
            field data, corresponding to the plot's X axis. If the method is
            "contour", then this must form a `meshgrid` with the Y coordinates.
            The shape of this array must correspond to that of the field data.
        coords_Y (`numpy.ndarray`): an array containing the coordinates of the
            field data, corresponding to the plot's Y axis. If the method is
            "contour", then this must form a `meshgrid` with the X coordinates.
            The shape of this array must correspond to that of the field data.
        field_data (`numpy.ndarray`): the field data to be plotted. The shape of
            this array must match that of the coordinates.
        contours (iter): an iterable object containing the values to use to
            contour the field. In the case of a scatter plot, these values
            determine the bins used for colouring the points.
        method (str): determines the method to use to plot the field. Valid
            options are "contour", "tricontour" or "scatter" (which cannot be
            used to plot contour lines).
        plot_filled_contours (bool, optional): whether to plot (coloured) filled
            contours. Defaults to True.
        plot_contour_lines (bool, optional): whether to plot the contour lines.
            Defaults to True.
        projection (:class:`Projection`, optional): a cartopy projection object.
            Defaults to None.
        cmap (`matplotlib.cmap`, optional): the colour map to be used for the
            the coloured field. Defaults to None, in which case the default
            matplotlib option is called (usually viridis).
        remove_lines (bool, optional): when using a filled contour method, gaps
            may be left between the filled colours (corresponding to the
            contour lines). This option fills in those lines. Defaults to False.
        transparency (float, optional): a number between 0 and 1 representing
            the transparency of the filling/colouring to be applied, with 0
            corresponding to transparent. Defaults to 1.0 (i.e. opaque).
        contour_linestyles (_type_, optional): the line style to be used for
            contour lines. Defaults to None, in which case positive contours
            have a solid line and negative contours are dashed.
        contour_linewidths (int, optional): the line width to use for contours.
            Can be a number or an array matching the contours, Defaults to 1.
        contour_linecolors (str, optional): the colour to be used for the
            contour lines. Can be a string or an array of size matching the
            contours. Defaults to 'black'.

    Raises:
        ValueError: if the `remove_lines` option and `plot_contour_lines` are
            both set to True.
        ValueError: if the `remove_lines` option is used with the "scatter"
            method.
        ValueError: if the `plot_contour_lines` option is set to True but the
            "scatter" method is used.

    Returns:
        tuple: the generated filled contour and line contour set objects. These
            can be used later, for instance to generate a colorbar. If one of
            these objects is not created, None will be returned.
    """

    #--------------------------------------------------------------------------#
    # Checks
    #--------------------------------------------------------------------------#

    if method not in ['contour', 'tricontour', 'scatter']:
        raise ValueError(
            f'Method {method} to field contour plot not recognised. ' +
            'Valid options are "contour", "tricontour" or "scatter"')

    if not (plot_filled_contours or plot_contour_lines):
        raise ValueError('field_contour_plot invalid option: '+
                         'both plot_filled_contours and plot_contour_lines '+
                         'are set to False. One must be True!')

    if remove_lines and plot_contour_lines:
        raise ValueError('field_contour_plot invalid option: '+
                         'the remove_lines argument can only be True if the '+
                         'plot_contour_lines argument is False')

    if remove_lines and method == 'scatter':
        raise ValueError('field_contour_plot invalid option: '+
                         'the remove_lines argument cannot be used with '+
                         'the "scatter" method')

    #--------------------------------------------------------------------------#
    # Things related to a Cartopy projection
    #--------------------------------------------------------------------------#

    # By default, these variables are None and do nothing
    if projection is not None:
        ax.set_global()
        transform_extent = (0, 360, -90, 90)
        transform_crs = projection.PlateCarree()
    else:
        transform_extent = None
        transform_crs = None

    #--------------------------------------------------------------------------#
    # Plot coloured fillings between contours
    #--------------------------------------------------------------------------#

    if plot_filled_contours:
        # Plot field using cmap, depending on the method
        if method == 'contour':
            cf = ax.contourf(coords_X, coords_Y, field_data, contours,
                             cmap=cmap, alpha=transparency, origin='lower',
                             extent=transform_extent, transform=transform_crs)

        elif method == 'tricontour':
            cf = ax.tricontourf(coords_X, coords_Y, field_data, contours,
                                cmap=cmap, alpha=transparency, origin='lower',
                                extent=transform_extent, transform=transform_crs)

        elif method == 'scatter':
            raise NotImplementedError('Scatter method not yet implemented')

        # Contour lines may appear as gaps in plot. These can be filled here
        if remove_lines:
            for c in cf.collections:
                c.set_edgecolor("face")

    else:
        cf = None

    #--------------------------------------------------------------------------#
    # Plot contour lines
    #--------------------------------------------------------------------------#

    if plot_contour_lines:
        if method == 'contour':
            cl = ax.contour(coords_X, coords_Y, field_data, contours,
                            linewidths=contour_linewidths,
                            colors=contour_linecolors,
                            linestyles=contour_linestyles, origin='lower',
                            extent=transform_extent, transform=transform_crs)

        elif method == 'tricontour':
            cl = ax.tricontour(coords_X, coords_Y, field_data, contours,
                               linewidths=contour_linewidths,
                               colors=contour_linecolors,
                               linestyles=contour_linestyles, origin='lower',
                               extent=transform_extent, transform=transform_crs)

        elif method == 'scatter':
            raise ValueError(
                'Scatter method not compatible with plotting line contours. '+
                'If you want to use the scatter method then set the '+
                'plot_contour_lines argument to False.')

    else:
        cl = None

    #--------------------------------------------------------------------------#
    # Finish
    #--------------------------------------------------------------------------#

    return cf, cl

# This is a separate routine to "plot_contoured_field" so that:
# (a) colorbars can be easily applied to some (but not all) subplot axes
# (b) the number of arguments to "plot_contoured_field" can be minimised
# This routine motivates "plot_contoured_field" returning the `cf` object.
# TODO: allow the position of the colorbar to be easily specified and adjusted
# TODO: make a global version of this routine (in which the position can be
#       easily determined)
def add_colorbar(ax, cf, cbar_label, cbar_format=None, cbar_ticks=None,
                 cbar_labelpad=None):
    """
    Adds a colour bar to a filled contour plot.

    Args:
        ax (:class:`AxesSubplot`): the pyplot axes object to plot the field on.
            Note that this is currently unused.
        cf (`matplotlib.contour`): the filled contour object for the plot. This
            is used to obtain the colours to be displayed in the colorbar.
        cbar_label (str): the label to be applied to the colorbar.
        cbar_format (float, optional): the formatting to used for the ticklabels
            attached to the colorbar. Defaults to None.
        cbar_ticks (iter, optional): which ticks to be attached to the colorbar.
            Defaults to None.
        cbar_labelpad (float, optional): the padding to apply to the colorbar
            label. Defaults to None.
    """

    cb = plt.colorbar(cf, format=cbar_format, ticks=cbar_ticks)
    if cbar_label is not None:
        cb.set_label(cbar_label, labelpad=cbar_labelpad)

    # TODO: can we add this to the ax instead?

# This is a separate routine to "plot_contoured_field" so that the number of
# arguments to "plot_contoured_field" can be minimised, because I find that I
# don't usually want to label contours, and doing so requires specifying a lot
# of variables. Since this simply calls `clabel`, there is an argument that we
# don't need this routine, but maybe it has an easier interface than `clabel`.
# This routine motivates "plot_contoured_field" returning the `cl` object.
def label_contour_lines(ax, cl, clabel_levels=None, clabel_fontsize=None,
                        clabel_locations=None, clabel_format=None):
    """
    Adds labels to contour lines on a 2D contour plot.

    Args:
        ax (:class:`AxesSubplot`): the pyplot axes object containing the
            contoured field plot.
        cl (`matplotlib.contour`): the contoured lines object for the plot.
        clabel_levels (iter, optional): a list of floats corresponding to the
            contour levels to label. Defaults to None.
        clabel_fontsize (float, optional): the size of contour labels. Defaults
            to None.
        clabel_locations (iter, optional): specific locations in the axes
            coordinates for the contour labels. Defaults to None.
        clabel_format (str, optional): the format to use for the labels.
            Defaults to None.
    """

    ax.clabel(cl, levels=clabel_levels, fontsize=clabel_fontsize,
              manual=clabel_locations, fmt=clabel_format)
