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
                           quiver_npts=1,
                           units='xy', scale=None, angles='xy', scale_units='xy',
                           cbar_label=None,
                           no_cbar=False,
                           contour_method=None,
                           contours=None,
                           contour_lines=True,
                           colour_scheme='Blues', restricted_cmap=True,
                           extend_cmap=True, remove_contour=False,
                           linestyle=None, linewidth=1,
                           linecolours='black',
                           fontsize=24, title=None, title_method='full',
                           titlepad=30, ax=None,
                           slice_label=None, time=None, time_method='seconds',
                           text=None, text_pos=None,
                           xlims=None, ylims=None, xticks=None, yticks=None,
                           xticklabels=None, yticklabels=None,
                           xlabel=None, ylabel=None, xlabelpad=-20, ylabelpad=None):
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
                                                     remove_contour=remove_contour)

        extend = 'both' if extend_cmap else 'neither'

    #--------------------------------------------------------------------------#
    # Plot contours
    #--------------------------------------------------------------------------#

    if contour_method is not None:
        cf = ax.contourf(coords_X, coords_Y, contour_data, colour_contours,
                         cmap=cmap, extend=extend)

        if extend_cmap:
            # Set colours for over and under shoots
            cf.cmap.set_under('magenta')
            cf.cmap.set_over('yellow')

        if contour_lines:
            if linestyle is not None:
                cl = ax.contour(coords_X, coords_Y, contour_data, line_contours,
                                linewidths=linewidth, colors=linecolours,
                                linestyles=linestyle)

        # FIXME: should this not be associated with the axis?
        if cbar_label is None:
            cbar_label = get_label(field_name)
        if not no_cbar:
            plt.colorbar(cf, label=cbar_label)

    #--------------------------------------------------------------------------#
    # Plot quivers
    #--------------------------------------------------------------------------#

    ax.quiver(coords_X[::quiver_npts,::quiver_npts], coords_Y[::quiver_npts,::quiver_npts],
              field_data_X[::quiver_npts,::quiver_npts], field_data_Y[::quiver_npts,::quiver_npts],
              units=units, scale=scale, scale_units=scale_units, angles=angles)

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

        fig.savefig(plotname, bbox_inches='tight')

        plt.close()
