"""
This file makes a 2D quiver plot
"""
import matplotlib.pyplot as plt
import numpy as np
from .plot_decorations import *


def individual_quiver_plot(coords_X, coords_Y, field_data_X, field_data_Y,
                           testname=None, plotname=None,
                           figsize=(8,8), field_name=None,
                           extra_field_data=None, extra_field_name=None,
                           slice_name=None, slice_idx=None,
                           quiver_npts=1, x_offset=None, y_offset=None,
                           units='xy', scale=None, angles='xy', scale_units='xy',
                           restrict_quivers=False,
                           projection=None,
                           spherical_centre=(0.0, 0.0),
                           cbar_label=None,
                           cbar_labelpad=None,
                           cbar_ticks=None,
                           cbar_format=None,
                           no_cbar=False,
                           contour_method=None,
                           contours=None,
                           contour_lines=True,
                           colour_scheme='Blues', restricted_cmap=None,
                           colour_levels_scaling=1.2,
                           extend_cmap=True, remove_contour=False,
                           linestyle=None, linewidth=1,
                           linecolours='black',
                           fontsize=24, title=None, title_method='full',
                           titlepad=30, ax=None,
                           slice_label=None, time=None, time_method='seconds',
                           text=None, text_pos=None,
                           xlims=None, ylims=None, xticks=None, yticks=None,
                           xticklabels=None, yticklabels=None,
                           xlabel=None, ylabel=None, xlabelpad=-20, ylabelpad=None,
                           dpi=None, gridline_args=None):
    """
    Makes an individual quiver plot of a field from a netCDF
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
        ax = fig.add_subplot(111)

        if projection is None:
            ax = fig.add_subplot(111)
        elif projection == 'orthographic':
            import cartopy.crs as ccrs
            ax = fig.add_subplot(1, 1, 1,
                                 projection=ccrs.Orthographic(spherical_centre[0]*180./np.pi,
                                                              spherical_centre[1]*180./np.pi))
        else:
            raise ValueError('Projection %s not implemented' % projection)

    if projection is None:
        crs = None
        extent = None
    elif projection == 'orthographic':
        import cartopy.crs as ccrs
        ax.set_global()
        coords_X *= 360.0/(2*np.pi)
        coords_Y *= 360.0/(2*np.pi)
        crs = ccrs.PlateCarree()
        extent = (0, 360, -90, 90)
    else:
        raise ValueError('Projection %s not implemented' % projection)

    #--------------------------------------------------------------------------#
    # Determine contour data
    #--------------------------------------------------------------------------#

    if contour_method is None:
        pass

    elif contour_method == 'magnitude':
        contour_data = np.sqrt(field_data_X**2 + field_data_Y**2)

    elif contour_method == 'x':
        contour_data = field_data_X

    elif contour_method == 'y':
        contour_data = field_data_Y

    elif contour_method == 'extra':
        contour_data = extra_field_data

    else:
        raise ValueError('Contour method %s not valid' % contour_method)

    #--------------------------------------------------------------------------#
    # Determine contour specifics
    #--------------------------------------------------------------------------#

    if contour_method is None:
        field_min = None
        field_max = None
    else:

        field_min = np.amin(contour_data)
        field_max = np.amax(contour_data)

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

        #----------------------------------------------------------------------#
        # Set up colour bar details
        #----------------------------------------------------------------------#

        cmap, line_contours = colourmap_and_contours(colour_contours,
                                                     colour_scheme=colour_scheme,
                                                     restricted_cmap=restricted_cmap,
                                                     colour_levels_scaling=colour_levels_scaling,
                                                     remove_contour=remove_contour)

        extend = 'both' if extend_cmap else 'neither'

    #--------------------------------------------------------------------------#
    # Plot contours
    #--------------------------------------------------------------------------#

    if contour_method is not None:
        cf = ax.contourf(coords_X, coords_Y, contour_data, colour_contours,
                         cmap=cmap, extend=extend, transform=crs, extent=extent,
                         origin='lower')

        if extend_cmap:
            # Set colours for over and under shoots
            cf.cmap.set_under('magenta')
            cf.cmap.set_over('yellow')

        if contour_lines:
            if linestyle is not None:
                cl = ax.contour(coords_X, coords_Y, contour_data, line_contours,
                                linewidths=linewidth, colors=linecolours,
                                linestyles=linestyle, transform=crs, extent=extent,
                                origin='lower')

        # FIXME: should this not be associated with the axis?
        if cbar_label is None and field_name is not None:
            cbar_label = get_label(field_name)
        if not no_cbar:
            cb = plt.colorbar(cf, format=cbar_format, ticks=cbar_ticks)
            if cbar_label is not None:
                cb.set_label(cbar_label, labelpad=cbar_labelpad)

    if gridline_args is not None:
        ax.gridlines(**gridline_args)

    #--------------------------------------------------------------------------#
    # Restrict extent of quivers if required
    #--------------------------------------------------------------------------#

    if restrict_quivers:
        # Assume that we are chopping off the top and bottom 10% of values
        # TODO: add option for this
        # TODO: do this in its own function?
        top_cutoff = np.min(coords_Y) + 0.9*(np.max(coords_Y) - np.min(coords_Y))
        bot_cutoff = np.min(coords_Y) + 0.1*(np.max(coords_Y) - np.min(coords_Y))

        # Use numpy function to do elementwise masking
        filter_array = np.logical_and(bot_cutoff <= coords_Y, coords_Y <= top_cutoff)
        # X coordinate shouldn't be restricted
        filter_shape = (np.shape(coords_Y)[0],
                        # Results in a 1D array
                        int(np.shape(coords_Y[filter_array])[0] / np.shape(coords_Y)[0]))

        # Need to reshape arrays after filtering
        field_data_X = np.reshape(field_data_X[filter_array], filter_shape)
        field_data_Y = np.reshape(field_data_Y[filter_array], filter_shape)
        coords_X = np.reshape(coords_X[filter_array], filter_shape)
        coords_Y = np.reshape(coords_Y[filter_array], filter_shape)

    #--------------------------------------------------------------------------#
    # Plot quivers
    #--------------------------------------------------------------------------#

    if type(quiver_npts) in (tuple,list):
        quiver_npts_x = quiver_npts[0]
        quiver_npts_y = quiver_npts[1]
    elif type(quiver_npts) in [int, float]:
        quiver_npts_x = quiver_npts
        quiver_npts_y = quiver_npts
    else:
        raise TypeError(f'Type {type(quiver_npts)} is not supported')
    

    x_slice = slice(x_offset, None, quiver_npts_x)
    y_slice = slice(y_offset, None, quiver_npts_y)

    if crs is None:
        # separately handle this case as None transform results in no arrows
        ax.quiver(coords_X[x_slice,y_slice], coords_Y[x_slice,y_slice],
                  field_data_X[x_slice,y_slice], field_data_Y[x_slice,y_slice],
                  units=units, scale=scale, scale_units=scale_units, angles=angles,
                  zorder=3)
    else:
        ax.quiver(coords_X[x_slice,y_slice], coords_Y[x_slice,y_slice],
                  field_data_X[x_slice,y_slice], field_data_Y[x_slice,y_slice],
                  units=units, scale=scale, scale_units=scale_units, angles=angles, transform=crs,
                  zorder=3)


    # Add axes labels, set limits and add a title
    axes_limits_labels_and_titles(ax, xlabel=xlabel, xlabelpad=xlabelpad, xlims=xlims,
                                  xticks=xticks, xticklabels=xticklabels,
                                  ylabel=ylabel, ylabelpad=ylabelpad, ylims=ylims,
                                  yticks=yticks, yticklabels=yticklabels,
                                  title=title, title_method=title_method, titlepad=titlepad,
                                  slice_label=slice_label, time=time, time_method=time_method,
                                  field_min=field_min, field_max=field_max)

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
        if contour_method is None:
            return None
        else:
            # Return cf to give access to colorbars
            return cf
    else:

        fig.savefig(plotname, bbox_inches='tight', dpi=dpi)

        plt.close()
