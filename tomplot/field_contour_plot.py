"""
This file provides some plotting tools
"""
import matplotlib.pyplot as plt
import numpy as np
from .plot_decorations import *


def individual_field_contour_plot(coords_X, coords_Y, field_data,
                                  testname=None, plotname=None,
                                  figsize=(8,8), field_name=None,
                                  extra_field_data=None, extra_field_name=None,
                                  slice_name=None, slice_idx=None,
                                  projection=None,
                                  spherical_centre=(0.0, 0.0),
                                  cbar_label=None,
                                  cbar_labelpad=None,
                                  cbar_ticks=None,
                                  cbar_format=None,
                                  no_cbar=False,
                                  contours=None,
                                  tricontour=False,
                                  extra_contours=None,
                                  contour_lines=True,
                                  cmap=None,
                                  line_contours=None,
                                  clabel=False,
                                  clabel_levels=None, clabel_fmt=None,
                                  clabel_locations=None, clabel_fontsize=None,
                                  colour_scheme='Blues', restricted_cmap=None,
                                  colour_levels_scaling=1.2,
                                  extend_cmap=True, remove_contour=False,
                                  transparency=1.0,
                                  linestyle=None, linewidth=1,
                                  linecolours='black',
                                  fontsize=24, title=None, title_method='full',
                                  titlepad=30, ax=None, title_size=24,
                                  slice_label=None, time=None, time_method='seconds',
                                  text=None, text_pos=None,
                                  xlims=None, ylims=None, xticks=None, yticks=None,
                                  xticklabels=None, yticklabels=None,
                                  xlabel=None, ylabel=None, xlabelpad=-20, ylabelpad=None,
                                  dpi=None, gridline_args=None,
                                  remove_lines=False):
    """
    Makes an individual coloured 2D contour plot of a field from a netCDF
    field output file.

    :arg coords_X:   a 2D array of the X-coordinates at which to plot the field.
    :arg coords_Y:   a 2D array of the Y-coordinates at which to plot the field.
    :arg field_data: a 2D array giving the values of the field to be plotted.
    :arg testname:   (Optional) a string describing the test. Usually used for
                     accessing some pre-set specific plot settings.
    :arg field_name: (Optional). The name of the field. Used for providing text
                     information on the plot.
    :arg time:       (Optional). The simulation time at corresponding to the
                     data. Used for providing text information on the plot.
    :arg slice_name: (Optional). The dimension along which to slice the data.
                     Required if the domain has dimension greater than 1.
                     Used for providing text information on the plot.
    :arg slice_idx:  (Optional). The index of the other dimensions on the
                     interpolation grid at which to slice. Used for providing
                     text information on the plot.
    """

    #--------------------------------------------------------------------------#
    # Checks
    #--------------------------------------------------------------------------#

    # TODO: fill in checks here


    ax_provided = (ax is not None)
    if not ax_provided and plotname is None:
        raise ValueError('If ax is not provided to make 1D plot then a plotname must be specifed')

    #--------------------------------------------------------------------------#
    # Make figure
    #--------------------------------------------------------------------------#

    if ax is None:
        # This is declared BEFORE figure and ax are initialised
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif')
        font = {'size':fontsize}
        plt.rc('font',**font)

        fig = plt.figure(1, figsize=figsize)

        if projection is None:
            ax = fig.add_subplot(111)
        elif projection == 'orthographic':
            import cartopy.crs as ccrs
            ax = fig.add_subplot(1, 1, 1,
                                 projection=ccrs.Orthographic(spherical_centre[0]*180./np.pi,
                                                              spherical_centre[1]*180./np.pi))
        elif projection == 'robsinson':
            import cartopy.crs as ccrs
            ax = fig.add_subplot(1, 1, 1,
                                 projection=ccrs.Robinson(spherical_centre[0]*180./np.pi))
        else:
            raise ValueError('Projection %s not implemented' % projection)

    if projection is None:
        crs = None
        extent = None
    elif projection in ['orthographic','robinson']:
        import cartopy.crs as ccrs
        ax.set_global()
        coords_X *= 360.0/(2*np.pi)
        coords_Y *= 360.0/(2*np.pi)
        crs = ccrs.PlateCarree()
        extent = (0, 360, -90, 90)
    else:
        raise ValueError('Projection %s not implemented' % projection)


    #--------------------------------------------------------------------------#
    # Determine contours
    #--------------------------------------------------------------------------#

    field_min = np.amin(field_data)
    field_max = np.amax(field_data)

    if contours is None:
        # Deal with case of no actual field values
        if abs(field_max - field_min) < 1e-15:
            field_min -= 1e-14
            field_max += 1e-14
            num_steps = 1
        else:
            num_steps = 10

        field_step = (field_max - field_min) / num_steps
        colour_contours = np.arange(field_min, field_max+field_step, step=field_step)
    else:
        colour_contours = contours

    if extra_field_data is not None:
        if extra_contours is None:
            extra_field_min = np.amin(extra_field_data)
            extra_field_max = np.amax(extra_field_data)

            # Deal with case of no actual field values
            if abs(extra_field_max - extra_field_min) < 1e-15:
                extra_field_min -= 1e-14
                extra_field_max += 1e-14
                num_steps = 1
            else:
                num_steps = 10

            extra_field_step = (extra_field_max - extra_field_min) / num_steps
            extra_contours = np.arange(extra_field_min, extra_field_max+extra_field_step, step=extra_field_step)

    #--------------------------------------------------------------------------#
    # Set up colour bar details
    #--------------------------------------------------------------------------#

    if cmap is None:
        cmap, line_contours_new = colourmap_and_contours(colour_contours,
                                                         colour_scheme=colour_scheme,
                                                         restricted_cmap=restricted_cmap,
                                                         colour_levels_scaling=colour_levels_scaling,
                                                         remove_contour=remove_contour)

    if line_contours is None:
        line_contours = line_contours_new

    extend = 'both' if extend_cmap else 'neither'

    #--------------------------------------------------------------------------#
    # Plot
    #--------------------------------------------------------------------------#

    if tricontour:
        cf = ax.tricontourf(coords_X, coords_Y, field_data, colour_contours,
                            cmap=cmap, extent=extent, transform=crs,
                            origin='lower', alpha=transparency)
    else:
        cf = ax.contourf(coords_X, coords_Y, field_data, colour_contours,
                         cmap=cmap, extent=extent, transform=crs,
                         origin='lower', alpha=transparency)

    if remove_lines:
        for c in cf.collections:
            c.set_edgecolor("face")

    if extend_cmap:
        # Set colours for over and under shoots
        cf.cmap.set_under('magenta')
        cf.cmap.set_over('yellow')

    if contour_lines:
        if linestyle is None and extra_field_data is not None:
            linestyles = 'solid'
        if tricontour:
            cl = ax.tricontour(coords_X, coords_Y, field_data, line_contours,
                               linewidths=linewidth, colors=linecolours,
                               linestyles=linestyle, transform=crs, extent=extent,
                               origin='lower')
        else:
            cl = ax.contour(coords_X, coords_Y, field_data, line_contours,
                            linewidths=linewidth, colors=linecolours,
                            linestyles=linestyle, transform=crs, extent=extent,
                            origin='lower')

        if clabel == True:
            ax.clabel(cl, levels=clabel_levels, fontsize=clabel_fontsize, manual=clabel_locations, fmt=clabel_fmt)

    if extra_field_data is not None:
        if tricontour:
            raise NotImplementedError
        cle = ax.contour(coords_X, coords_Y, extra_field_data,
                         extra_contours, linewidths=linewidth,
                         colors=linecolours, linestyles='dashed',
                         transform=crs, extent=extent, origin='lower')

    # FIXME: should this not be associated with the axis?
    if cbar_label is None and field_name is not None:
        cbar_label = get_label(field_name)
    if not no_cbar:
        cb = plt.colorbar(cf, format=cbar_format, ticks=cbar_ticks)
        if cbar_label is not None:
            cb.set_label(cbar_label, labelpad=cbar_labelpad)

    # Add axes labels, set limits and add a title
    if projection is None:
        axes_limits_labels_and_titles(ax, xlabel=xlabel, xlabelpad=xlabelpad, xlims=xlims,
                                      xticks=xticks, xticklabels=xticklabels,
                                      ylabel=ylabel, ylabelpad=ylabelpad, ylims=ylims,
                                      yticks=yticks, yticklabels=yticklabels,
                                      title=title, title_method=title_method, titlepad=titlepad,
                                      slice_label=slice_label, time=time, time_method=time_method,
                                      field_min=field_min, field_max=field_max, title_size=title_size)
    elif gridline_args is not None:
        ax.gridlines(**gridline_args)

    if text is not None:
        if text_pos is None or not isinstance(text_pos, tuple):
            raise ValueError('If text has been provided, a text_pos tuple must also be provided')
        else:
            ax.text(text_pos[0], text_pos[1], text, fontsize=fontsize*1.25,
                    ha='center', va='center')

    #--------------------------------------------------------------------------#
    # Finish
    #--------------------------------------------------------------------------#

    if ax_provided:
        return cf
    else:

        print(f'Saving figure to {plotname}')
        fig.savefig(plotname, bbox_inches='tight', dpi=dpi)

        plt.close()
