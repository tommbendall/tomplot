"""
This file provides some plotting tools
"""
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np
from .plot_decorations import *

def individual_convergence_plot(dirnames, variable, fields, run_ids, error,
                                best_fit=True, testname=None, plotdir=None,
                                override_dirname=False, field_labels=None,
                                figsize=(8,8), colours=None, markers=None,
                                linestyles=None, linewidth=2, markersize=8,
                                fontsize=24, title=None, comparison_lines=None,
                                ax=None, grid=True, xlabel='default',
                                ylabel='default', best_fit_deg=1,
                                xlim=None, ylim=None, legend_bbox=(1.0,1.0),
                                legend_ncol=1, label_style='gradient_full',
                                format='png', dpi=None, titlepad=None,
                                leg_col_spacing=None, leg_fontsize=None):
    """
    Makes an individual convergence plot for errors from a global netCDF
    diagnostics file.

    :arg dirnames:          list of names of the directories with the data
    :arg variable:          a string giving the name of the variable that we
                            measure convergence with respect to. e.g. 'dx'.
    :arg fields:            a single string, or a list of all string giving the
                            fields to be plotted.
    :arg run_ids:           a list of IDs for the runs to be plotted. This
                            should be of the same length as the list of fields.
                            Can be replaced by a single integer.
    :arg error:             a string giving the error diagnostic to be plotted
    :arg best_fit:          (Optional) Boolean to plot best fit lines through
                            the error points. If True then the gradients of
                            these lines will be added to the legend labels.
                            Default is True.
    :arg testname:          (Optional) a string used for identifying particular
                            plot parameters and labelling the figures.
    :arg plotdir:           (Optional) the name of the directory to output the
                            plots to. If None (the default value) then plotdir
                            will be the same as dirname.
    :arg override_dirname:  (Optional) Boolean for determining whether to
                            override the default path to the dirname.
    :arg field_labels:      (Optional) a list of labels corresponding to the
                            fields to be used in the plot legends.
    :arg figsize:           (Optional) the size of the figure to use. If not
                            specified then the default is (8,8).
    :arg colours:           (Optional) a list of colours to use for the markers
                            and best fit lines. If not specified then will be
                            auto-picked.
    :arg markers:           (Optional) the markers to use in the convergence
                            plot. If not specified then it will be auto-picked.
    :arg linestyles:        (Optional) the line styles to use for the best fit
                            lines on the convergence plot. Default is '-'
    :arg linewidth:         (Optional) the widths to use for the best fit lines.
                            Default is 2.
    :arg markersize:        (Optional) the size of the markers to use in the
                            convergence plot. Default is 8.
    :arg fontsize:          (Optional) the fontsize to use for axis labels, tick
                            labels, plot titles and the legend labels.
    :arg title:             (Optional) a title of the plot. Default is None.
    :arg comparison_lines:  (Optional) a list of floats that can be used for
                            giving comparison convergence lines, corresponding
                            to specific convergence rates. e.g. [1.0, 2.0].
    :arg ax:                (Optional) a matplotlib.pyplot ax object. If this is
                            provided, then the plotting will be performed on
                            this ax and this will then be returned. This could
                            be useful for making subplots. Otherwise a new fig
                            and ax will be set up and no plot will be returned.
    :arg grid:              (Optional) Boolean determining whether to show the
                            background gridlines or not. Default is True.
    :arg xlabel:            (Optional) the label for the x-axis. If not provided
                            then this will be determined from 'variable'.
    :arg ylabel:            (Optional) the label for the y-axis. If not provided
                            then this will be determined from 'error'.
    :arg xlim:              (Optional) the limits for the x-axis. If not given
                            then just uses the pyplot defaults.
    :arg ylim:              (Optional) the limits for the y-axis. If not given
                            then just uses the pyplot defaults.
    """

    # Each field will represent a single line on the plot
    # Each run ID will be at a different resolution -- each is a point on the plot

    #--------------------------------------------------------------------------#
    # Checks
    #--------------------------------------------------------------------------#

    # Make results dirnames into an array if only a single string is provided
    if isinstance(dirnames, str):
        dirnames = [dirnames]

    # Make fields into an array if only a single string is provided
    if isinstance(fields, str):
        fields = [fields]

    # Make run_ids into an array if only a single integer is provided
    if isinstance(run_ids, int):
        run_ids = [[run_ids]]*len(fields)
    else:
        # Assume that run_ids is already a list
        # If it's a list of lists then we're already good
        if not isinstance(run_ids[0], list):
            # Then we assume that the same list is used for every field
            run_ids = [run_ids]*len(fields)

    if len(run_ids) != len(fields):
        raise ValueError('The lengths of run_ids list and fields list are not equal')

    if colours is not None:
        if len(colours) != len(fields)*len(dirnames):
            raise ValueError('The list of colours should have the same '+
                             'length as the list of fields. Found %d but expected %d' %
                             (len(colours), len(fields)*len(dirnames)))

    if markers is not None:
        if len(markers) != len(fields)*len(dirnames):
            raise ValueError('The list of markers should have the same '+
                             'length as the list of fields. Found %d but expected %d' %
                             (len(markers), len(fields)*len(dirnames)))

    if linestyles is not None:
        if len(linestyles) != len(fields)*len(dirnames):
            raise ValueError('The list of linestyles should have the same '+
                             'length as the list of fields. Found %d but expected %d' %
                             (len(linestyles), len(fields)*len(dirnames)))
    
    if type(best_fit) == list:
        if len(best_fit) != len(best_fit)*len(dirnames):
            raise ValueError('The list of best fits should have the same '+
                             'length as the list of fields. Found %d but expected %d' %
                             (len(best_fit), len(fields)*len(dirnames)))
    else:
        best_fit = [best_fit]*(len(fields)*len(dirnames))

    ax_provided = (ax is not None)

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
    # Loop through fields adding lines to plot
    #--------------------------------------------------------------------------#

    for j, dirname in enumerate(dirnames):

        if override_dirname:
            filename = dirname
        else:
            filename = 'results/'+dirname+'/global_output.nc'

        data_file = Dataset(filename, 'r')

        for i, (field, run_id_list) in enumerate(zip(fields, run_ids)):

            k = i + j*len(fields)

            #------------------------------------------------------------------#
            # Extract data
            #------------------------------------------------------------------#

            # Should be an array indexed by run_id
            variable_data = np.log(data_file.variables[variable][run_id_list])

            # Use only the errors from the final time step
            error_data = np.log(data_file.groups[field]['errors'].variables[error][run_id_list,-1])

            #------------------------------------------------------------------#
            # Determine marker colours, shapes and labels
            #------------------------------------------------------------------#

            colour = colours[k] if colours is not None else get_colour(testname, field, k)
            marker = markers[k] if markers is not None else get_marker(testname, field, k)
            linestyle = linestyles[k] if linestyles is not None else '-'

            if field_labels is not None:
                # Label is just read in
                label = field_labels[k]
            else:
                label = get_label(field)


            #------------------------------------------------------------------#
            # Plot errors
            #------------------------------------------------------------------#

            if best_fit[k]:
                # Get line of best fit first to amend the label
                best_fit_line = np.poly1d(np.polyfit(variable_data, error_data, deg=best_fit_deg))
                ax.plot(variable_data, best_fit_line(variable_data),
                        linestyle=linestyle, color=colour, lw=linewidth)

                if label_style == 'gradient_full':
                    label = label+' gradient: %1.3f' % best_fit_line[1]
                elif label_style == 'gradient_plain':
                    label = label+': %1.3f' % best_fit_line[1]
                elif label_style != 'plain':
                    raise ValueError('label_style not recognised')

            # Put a line through data points if a linestyle is specified
            if not best_fit[k] and linestyles is not None:
                data_ls = linestyles[k]
            else:
                data_ls = ''

            # Plot error points
            ax.plot(variable_data, error_data, color=colour,
                    marker=marker, label=label, linestyle=data_ls, ms=markersize)

        data_file.close()

    #--------------------------------------------------------------------------#
    # Decorations
    #--------------------------------------------------------------------------#

    if xlabel == 'default':
        xlabel = get_xlabel(variable, 'convergence')
    if xlabel is not None:
        ax.set_xlabel(xlabel)
    if ylabel == 'default':
        ylabel = get_ylabel(error, 'convergence')
        ax.set_ylabel(ylabel)
    if ylabel is not None:
        ax.set_ylabel(ylabel)


    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)
    if title is not None:
        ax.set_title(title, pad=titlepad)

    if grid:
        ax.grid('on')

    #--------------------------------------------------------------------------#
    # Comparison lines
    #--------------------------------------------------------------------------#

    # This is done after plot decorations to ensure
    # the lines pass through the centre of the figure

    if comparison_lines is not None:
        # Get centre of figure
        x_lower, x_upper = ax.get_xlim()[0], ax.get_xlim()[1]
        y_lower, y_upper = ax.get_ylim()[0], ax.get_ylim()[1]
        x_center = 0.5 * (x_lower + x_upper)
        y_center = 0.5 * (y_lower + y_upper)
        x = np.array([x_lower, x_upper])

        comparison_linestyles = ['dashed', 'dotted', 'dashdot']

        for j, gradient in enumerate(comparison_lines):

            label = get_xlabel(variable, 'grid parameters')
            if isinstance(gradient, int):
                if gradient != 1:
                    label += r'$^{%d}$' % gradient
            elif isinstance(gradient, float):
                label += r'$^{%.1f}$' % gradient
            else:
                raise ValueError('Gradients provided in comparison_lines must be '+
                                 'floats or integers, not %s' % type(gradient))

            # y = m*x + c
            y = gradient * (x - x_center) + y_center
            ax.plot(x, y, color='black', linestyle=comparison_linestyles[j % 3],
                    label=label, linewidth=linewidth, marker='')

    #--------------------------------------------------------------------------#
    # Legend
    #--------------------------------------------------------------------------#

    handles, labels = ax.get_legend_handles_labels()
    lgd = ax.legend(handles, labels, loc='upper center', ncol=legend_ncol,
                    bbox_to_anchor=legend_bbox, edgecolor='black',
                    fontsize=leg_fontsize, handletextpad=0.0,
                    columnspacing=leg_col_spacing)

    #--------------------------------------------------------------------------#
    # Save and finish plot
    #--------------------------------------------------------------------------#

    if ax_provided:
        return ax
    else:
        if plotdir is None:
            if len(dirnames) > 1:
                print('Convergence plot directory not specified. '+
                      'Adding to results/'+dirnames[0]+'/figures')
            plotdir = 'results/'+dirnames[0]+'/figures'
        plotname = plotdir+'/'+testname+'_'+error+'.'+format

        fig.savefig(plotname, bbox_extra_artists=(lgd,),
                    bbox_inches='tight', dpi=dpi)

        plt.close()
