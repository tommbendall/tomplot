"""
This file provides some plotting tools
"""
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np
from .plot_decorations import *

def individual_time_series_plot(dirnames, fields, run_ids, diagnostic,
                                testname=None, plotdir=None, override_dirname=False,
                                field_labels=None, field_labels_from=None,
                                figsize=(8,8), colours=None, linestyles=None,
                                linewidth=2, fontsize=24, title=None, ax=None,
                                grid=True, ylabel=None, xlim=None, ylim=None,
                                time_units='seconds', normalise=False, format='png',
                                dpi=None):
    """
    Makes an individual time series plot for fields from a global netCDF
    diagnostics file.

    :arg dirnames:          list of names of the directories with the data
    :arg fields:            a list of all the fields to be plotted.
    :arg run_ids:           a list of lists IDs for the runs to be plotted. This
                            should be a list of run IDs for each field. A list
                            can be replaced by a single integer.
    :arg diagnostic:        a string giving the global diagnostic to be plotted
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
    :arg linestyles:        (Optional) the line styles to use for the best fit
                            lines on the convergence plot. Default is '-'
    :arg linewidth:         (Optional) the widths to use for the best fit lines.
                            Default is 2.
    :arg fontsize:          (Optional) the fontsize to use for axis labels, tick
                            labels, plot titles and the legend labels.
    :arg title:             (Optional) a title of the plot. Default is None.
    :arg ax:                (Optional) a matplotlib.pyplot ax object. If this is
                            provided, then the plotting will be performed on
                            this ax and this will then be returned. This could
                            be useful for making subplots. Otherwise a new fig
                            and ax will be set up and no plot will be returned.
    :arg grid:              (Optional) Boolean determining whether to show the
                            background gridlines or not. Default is True.
    :arg ylabel:            (Optional) the label for the y-axis. If not provided
                            then this will be determined from 'error'.
    :arg xlim:              (Optional) the limits for the x-axis. If not given
                            then just uses the pyplot defaults.
    :arg ylim:              (Optional) the limits for the y-axis. If not given
                            then just uses the pyplot defaults.
    """

    # Each field will represent a single line on the plot

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
        if len(colours) != len(fields):
            raise ValueError('The list of colours should have the same '+
                             'length as the list of fields. Found %d but expected %d' %
                             (len(colours), len(fields)))


    if linestyles is not None:
        if len(linestyles) != len(fields):
            raise ValueError('The list of linestyles should have the same '+
                             'length as the list of fields. Found %d but expected %d' %
                             (len(linestyles), len(fields)))

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

    for k, dirname in enumerate(dirnames):

        if override_dirname:
            filename = dirname
        else:
            filename = 'results/'+dirname+'/global_output.nc'

        data_file = Dataset(filename, 'r')

        counter = 0

        for i, (field, run_id_list) in enumerate(zip(fields, run_ids)):

            for j, run_id in enumerate(run_id_list):

                #------------------------------------------------------------------#
                # Extract data
                #------------------------------------------------------------------#

                # Factor for units of time
                if time_units == 'seconds':
                    time_factor = 1.0
                    time_label = 's'
                elif time_units == 'days':
                    time_factor = 24.*60.*60
                    time_label = 'days'
                else:
                    raise ValueError('time_units %s not recognised' % time_units)

                # Should be an array indexed by run_id
                time_data = data_file.variables['time'][run_id, :] / time_factor

                # Use only the errors from the final time step
                if diagnostic in data_file.groups[field]['global_quantities'].variables.keys():
                    diagnostic_data = data_file.groups[field]['global_quantities'].variables[diagnostic][run_id,:]
                else:
                    diagnostic_data = data_file.groups[field]['errors'].variables[diagnostic][run_id,:]

                if normalise:
                    if abs(diagnostic_data[0]) > 1e-12:
                        diagnostic_data = (diagnostic_data - diagnostic_data[0]) / diagnostic_data[0]
                    else:
                        diagnostic_data = diagnostic_data - diagnostic_data[0]

                #------------------------------------------------------------------#
                # Determine marker colours, shapes and labels
                #------------------------------------------------------------------#

                colour = colours[i] if colours is not None else get_colour(testname, field, counter)
                linestyle = linestyle[i] if linestyles is not None else '-'

                if field_labels is not None:
                    # Label is just read in
                    label = field_labels[counter]
                else:
                    label = get_label(field)

                #------------------------------------------------------------------#
                # Plot diagnostics
                #------------------------------------------------------------------#

                # Plot diagnostics
                ax.plot(time_data, diagnostic_data, color=colour, label=label, linestyle=linestyle)

                counter += 1

        data_file.close()

    #--------------------------------------------------------------------------#
    # Decorations
    #--------------------------------------------------------------------------#

    xlabel = r'$t \ / $ '+time_label
    if ylabel is None:
        ylabel = get_ylabel(diagnostic, 'time series')

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)
    if title is not None:
        ax.set_title(title)

    if grid:
        ax.grid('on')

    #--------------------------------------------------------------------------#
    # Legend
    #--------------------------------------------------------------------------#

    handles, labels = ax.get_legend_handles_labels()
    lgd = ax.legend(handles, labels, loc='upper left',
                    bbox_to_anchor=(1.0,1.0), edgecolor='black')

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
        plotname = plotdir+'/'+testname+'_'+diagnostic+'.'+format

        fig.savefig(plotname, bbox_extra_artists=(lgd,),
                    bbox_inches='tight', dpi=dpi)

        plt.close()
