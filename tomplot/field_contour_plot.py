"""
This file provides a generic routine for plotting a field on a 2D surface,
generally as (filled) contours or with coloured points as a scatter plot.

Some auxiliary routines are also provided.
"""
import numpy as np
import warnings
from .tomplot_tools import tomplot_field_markersize, work_out_cmap_extension

__all__ = ["plot_contoured_field", "label_contour_lines"]


def plot_contoured_field(ax, coords_X, coords_Y, field_data, method, contours,
                         line_contours=None, projection=None,
                         plot_filled_contours=True, plot_contour_lines=True,
                         # Options relating to filled contours
                         cmap=None, remove_lines=False, transparency=1.0,
                         # Options relating to line contours
                         contour_linestyles=None, contour_linewidths=1,
                         contour_linecolors='black',
                         # Options related to scatter method
                         markersize=None, marker_scaling=1.0,
                         # Keyword arguments to be passed to underlying matplotlib routines
                         **kwargs):
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
        method (str): determines the method to use to plot the field. Valid
            options are "contour", "tricontour" or "scatter" (which cannot be
            used to plot contour lines).
        contours (iter): an iterable object containing the values to use to
            contour the field. In the case of a scatter plot, these values
            determine the bins used for colouring the points. If the
            `line_contours` argument is specified, then this is only the set of
            contours used for the filled contours.
        line_contours (iter, optional): an iterable object containing the values
            to use to plot contour lines. Defaults to None, in which case
            contour lines are plotted with the `contours` argument.
        projection (:class:`Projection`, optional): a cartopy projection object.
            Defaults to None.
        plot_filled_contours (bool, optional): whether to plot (coloured) filled
            contours. Defaults to True.
        plot_contour_lines (bool, optional): whether to plot the contour lines.
            Defaults to True.
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
        markersize (float, optional): the size of markers to use for the
            "scatter" method.
        marker_scaling (float, optional): an optional scaling to apply to the
            markersize when using the "scatter" method, and auto-generating the
            markersizes. Defaults to 1.0.
        **kwargs: keyword arguments to be passed to the underlying plotting
            routines (contourf, tricontourf, contour, tricontour, scatter)

    Raises:
        ValueError: if the `remove_lines` option and `plot_contour_lines` are
            both set to True.
        ValueError: if the `remove_lines` option is used with the "scatter"
            method.
        ValueError: if the `plot_contour_lines` option is set to True but the
            "scatter" method is used.
        ValueError: if the `markersize` option is set but the "scatter" method
            is not used.

    Returns:
        tuple: the generated filled contour and line contour set objects. These
            can be used later, for instance to generate a colorbar. If one of
            these objects is not created, None will be returned.
    """

    # ------------------------------------------------------------------------ #
    # Checks
    # ------------------------------------------------------------------------ #

    if method not in ['contour', 'tricontour', 'scatter']:
        raise ValueError(
            f'Method {method} to field contour plot not recognised. '
            + 'Valid options are "contour", "tricontour" or "scatter"')

    if not (plot_filled_contours or plot_contour_lines):
        raise ValueError('field_contour_plot invalid option: '
                         + 'both plot_filled_contours and plot_contour_lines '
                         + 'are set to False. One must be True!')

    if remove_lines and plot_contour_lines:
        raise ValueError('field_contour_plot invalid option: '
                         + 'the remove_lines argument can only be True if the '
                         + 'plot_contour_lines argument is False')

    if remove_lines and method == 'scatter':
        raise ValueError('field_contour_plot invalid option: '
                         + 'the remove_lines argument cannot be used with '
                         + 'the "scatter" method')

    if markersize is not None and method != 'scatter':
        raise ValueError('field_contour_plot invalid option: '
                         + 'markersize is set but can only be used with the '
                         + '"scatter" method')

    if method == 'tricontour':
        if len(np.shape(coords_X)) > 1:
            coords_X = np.array(coords_X).flatten()
            warnings.warn('WARNING: field_contour_plot with "tricontour" '
                          + 'method requires flattened data. Flattening '
                          + 'coords_X for you.')
        if len(np.shape(coords_Y)) > 1:
            coords_Y = np.array(coords_Y).flatten()
            warnings.warn('WARNING: field_contour_plot with "tricontour" '
                          + 'method requires flattened data. Flattening '
                          + 'coords_Y for you.')
        if len(np.shape(field_data)) > 1:
            field_data = np.array(field_data).flatten()
            warnings.warn('WARNING: field_contour_plot with "tricontour" '
                          + 'method requires flattened data. Flattening '
                          + 'field_data for you.')

    # ------------------------------------------------------------------------ #
    # Things related to a Cartopy projection
    # ------------------------------------------------------------------------ #

    # By default, these variables are None and do nothing
    if projection is not None:
        import cartopy.crs as ccrs
        ax.set_global()
        transform_extent = (-180, 180, -90, 90)
        if method == "scatter":
            transform_crs = ccrs.Geodetic()
        else:
            transform_crs = ccrs.PlateCarree()
    else:
        transform_extent = None
        transform_crs = None

    # ------------------------------------------------------------------------ #
    # Determine how to handle values beyond the contour range
    # ------------------------------------------------------------------------ #

    if method in ['contour', 'tricontour']:
        if 'extend' in kwargs.keys():
            cmap_extension = kwargs['extend']
            del kwargs['extend']
        else:
            cmap_extension = work_out_cmap_extension(cmap, contours)

    # ------------------------------------------------------------------------ #
    # Plot coloured fillings between contours
    # ------------------------------------------------------------------------ #

    if plot_filled_contours:
        # Plot field using cmap, depending on the method
        if method == 'contour':
            cf = ax.contourf(coords_X, coords_Y, field_data, contours,
                             cmap=cmap, alpha=transparency, origin='lower',
                             extent=transform_extent, transform=transform_crs,
                             extend=cmap_extension, **kwargs)

        elif method == 'tricontour':
            cf = ax.tricontourf(coords_X, coords_Y, field_data, contours,
                                cmap=cmap, alpha=transparency, origin='lower',
                                extent=transform_extent, transform=transform_crs,
                                extend=cmap_extension, **kwargs)

        elif method == 'scatter':
            if markersize is None:
                markersize = tomplot_field_markersize(
                    field_data, marker_scaling=marker_scaling, ax=ax)

            if transform_crs is None:
                # Need to do this separately to ensure that points are plotted
                cf = ax.scatter(coords_X, coords_Y, c=field_data, s=markersize,
                                vmin=contours[0], vmax=contours[-1], cmap=cmap,
                                alpha=transparency, **kwargs)
            else:
                cf = ax.scatter(coords_X, coords_Y, c=field_data, s=markersize,
                                vmin=contours[0], vmax=contours[-1], cmap=cmap,
                                alpha=transparency, transform=transform_crs,
                                **kwargs)

        # Contour lines may appear as gaps in plot. These can be filled here
        if remove_lines:
            for c in cf.collections:
                c.set_edgecolor("face")

    else:
        cf = None

    # ----------------------------------------------------------------------- #
    # Plot contour lines
    # ----------------------------------------------------------------------- #

    if plot_contour_lines:
        # If line contours aren't specified, get these from filled contours
        if line_contours is None:
            line_contours = contours

        if method == 'contour':
            cl = ax.contour(coords_X, coords_Y, field_data, line_contours,
                            linewidths=contour_linewidths,
                            colors=contour_linecolors,
                            linestyles=contour_linestyles, origin='lower',
                            extent=transform_extent, transform=transform_crs,
                            **kwargs)

        elif method == 'tricontour':
            cl = ax.tricontour(coords_X, coords_Y, field_data, line_contours,
                               linewidths=contour_linewidths,
                               colors=contour_linecolors,
                               linestyles=contour_linestyles, origin='lower',
                               extent=transform_extent, transform=transform_crs,
                               **kwargs)

        elif method == 'scatter':
            raise ValueError(
                'Scatter method not compatible with plotting line contours. '
                + 'If you want to use the scatter method then set the '
                + 'plot_contour_lines argument to False.')

    else:
        cl = None

    # ------------------------------------------------------------------------ #
    # Finish
    # ------------------------------------------------------------------------ #

    return cf, cl


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
