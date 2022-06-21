"""
This file provides some plotting tools
"""
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np
from .plot_decorations import *

def individual_field_1d_plot(coords, field_data, testname=None, plotname=None,
                             figsize=(8,8), field_name=None,
                             slice_name=None, slice_idx=None,
                             colour='black', linestyle='-', linewidth=2,
                             fontsize=24, title=None, title_method='full', ax=None,
                             grid=True, slice_label=None, time=None, time_method='seconds',
                             xlabel=None, ylabel=None, xlims=None, ylims=None,
                             xticklabels=None, dpi=None):
    """
    Makes an individual 1D plot of a field from a netCDF field output file.

    :arg coords:     a 1D array of coordinates at which to plot the field.
    :arg field_data: a 1D array giving the values of the field to be plotted.
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
    # Plot
    #--------------------------------------------------------------------------#

    ax.plot(coords, field_data, color=colour, linestyle=linestyle,
            linewidth=linewidth, marker='')

    # Add axes labels, set limits and add a title
    field_min = np.amin(field_data)
    field_max = np.amax(field_data)

    axes_limits_labels_and_titles(ax, xlabel=xlabel, xlims=xlims,
                                  xticklabels=xticklabels, ylabel=ylabel, ylims=ylims,
                                  title=title, title_method=title_method, slice_label=slice_label,
                                  time=time, time_method=time_method,
                                  field_min=field_min, field_max=field_max)

    if ax_provided:
        return ax
    else:

        fig.savefig(plotname, bbox_inches='tight', dpi=dpi)

        plt.close()
