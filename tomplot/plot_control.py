"""
This file provides routines for generating multiple plots. There are four types:
- 1D field slices
- 2D field contour plots
- time series plots of global quantities
- convergence plots of global errors
"""
import matplotlib.pyplot as plt
from pyop2.mpi import COMM_WORLD
from netCDF4 import Dataset
import numpy as np
from .extract_field_data import extract_1D_data, extract_2D_data
from .field_1D_plot import individual_field_1d_plot
from .field_contour_plot import individual_field_contour_plot
from .time_series_plot import individual_time_series_plot
from .convergence_plot import individual_convergence_plot

def make_field_plots(dirname, run_id, testname, fields,
                     time_idxs, slices=None, plotdir=None,
                     override_dirname=False, **kwargs):
    """
    A routine for controlling the auto-generation of plots of fields from
    netCDF field data.

    :arg dirname:          the name of the directory with the data
    :arg run_id:           the ID of the run to be plotted
    :arg testname:         a string to name the plots, usually based on the
                           specific transport test case used for the run to be
                           plotted. This is also used for identifying specific
                           decorations for the plots.
    :arg fields:           a list of names of fields to be plotted.
    :arg time_idxs:        a list of the integers that index the points in time
                           at which the data was dumped. If the string 'all' is
                           given it will be plotted at all times.
    :arg slices:           the slices on which to extract the fields. Can be a
                           list, e.g. ['xy', 'xz'], or a list of lists which
                           specifies the indices for the slices,
                           e.g. [['xy', 0, 1], ['yz', -1]]. Slices are not
                           relevant to 1D domains, and must be specified for
                           3D domains.
    :arg plotdir:          (Optional) the name of the directory to output the
                           plots to. If None (the default value) then plotdir
                           will be the same as dirname.
    :arg override_dirname: (Optional) Boolean for determining whether to
                           override the default path to the dirname
    **kwargs               Other keyword arguments to be passed through for
                           producing field plots.
    """

    # Determine whether we are 1D, 2D or 3D
    if override_dirname:
        filename = dirname
    else:
        filename = 'results/'+dirname+'/nc_fields/field_output_'+str(run_id)+'.nc'

    # The same data file is used for all field plots, so open at start
    comm = COMM_WORLD
    if comm.Get_size() > 1:
        raise ValueError('Plotting not implemented for parallel runs')

    data_file = Dataset(filename, 'r')

    topological_dim = data_file.variables['topological_dimension'][0]

    # Check that slices are appropriate
    if ((topological_dim == 1) and (slices is not None)):
        raise ValueError('slices must be None for plotting 1D data')

    elif ((topological_dim == 3) and (slices is None)):
        raise ValueError('to plot data we require slices to be specified')

    # allow 'all' as a time option
    if time_idxs == 'all':
        time_idxs = range(len(data_file['time'][:]))
    elif isinstance(time_idxs, int):
        time_idxs = [time_idxs]

    if isinstance(fields, str):
        fields = [fields]

    if isinstance(slices, str):
        slices = [slices]

    # Flags for adding options to kwargs
    xlabel_added = False
    xlims_added = False
    xticklabels_added = False
    ylabel_added = False
    ylims_added = False
    yticklabels_added = False
    extra_field_added = False

    #--------------------------------------------------------------------------#

    # Now loop through details and pass them to an individual plotter
    for time_idx in time_idxs:
        for field in fields:

            # Consider case of looping through slices
            if slices is not None:
                for slice_type in slices:

                    # Check if slice_type is a list and extract indices if so
                    if isinstance(slice_type, list):
                        slice_name = slice_type[0]
                        slice_idxs = slice_type[1:] # first element should be name
                    else:
                        # If no slice_idxs specified, use 0
                        slice_idxs = [0]
                        slice_name = slice_type

                    for slice_idx in slice_idxs:

                        if plotdir is None:
                            plotdir = 'results/'+dirname+'/figures'
                        plotname = '%s/%s_%s_slice_%s_run_%s_time_%s.png' % (plotdir, testname, field, slice_name,
                                                                             str(run_id), str(time_idx))

                        if slice_name in ['x','y','z']:

                            coords, field_data, time, \
                            coord_label, coord_lims, \
                            coord_ticks, slice_label = extract_1D_data(data_file, field, time_idx,
                                                                       slice_name=slice_name, slice_idx=slice_idx)

                            if 'xlabel' not in kwargs.keys():
                                kwargs['xlabel'] = coord_label
                                xlabel_added = True
                            if 'xlims' not in kwargs.keys():
                                kwargs['xlims'] = coord_lims
                                xlims_added = True
                            if 'xticklabels' not in kwargs.keys():
                                kwargs['xticklabels'] = coord_ticks
                                xticklabels_added = True

                            individual_field_1d_plot(coords, field_data, testname=testname, time=time,
                                                     plotname=plotname, field_name=field,
                                                     slice_name=slice_name, slice_idx=slice_idx,
                                                     slice_label=slice_label, **kwargs)

                            if xlabel_added:
                                kwargs.pop('xlabel')
                                xlabel_added = False
                            if xlims_added:
                                kwargs.pop('xlims')
                                xlims_added = False
                            if xticklabels_added:
                                kwargs.pop('xticklabels')
                                xticklabels_added = False

                        else:

                            coords_X, coords_Y, field_data, time, \
                            coord_labels, coord_lims, coord_ticks, \
                            slice_label = extract_2D_data(data_file, field, time_idx,
                                                          slice_name=slice_name, slice_idx=slice_idx)

                            # Extract second field data if we need it
                            if 'extra_field_name' in kwargs.keys():
                                coords_X, coords_Y, extra_field_data, time, \
                                coord_labels, coord_lims, coord_ticks, \
                                slice_label = extract_2D_data(data_file, kwargs['extra_field_name'], time_idx,
                                                              slice_name=slice_name, slice_idx=slice_idx)
                                kwargs['extra_field_data'] = extra_field_data
                                extra_field_added = True


                            # Add plotting details to kwargs if they were not already there
                            if 'xlabel' not in kwargs.keys():
                                kwargs['xlabel'] = coord_labels[0]
                                xlabel_added = True
                            if 'xlims' not in kwargs.keys():
                                kwargs['xlims'] = coord_lims[0]
                                xlims_added = True
                            if 'xticklabels' not in kwargs.keys():
                                kwargs['xticklabels'] = coord_ticks[0]
                                xticklabels_added = True
                            if 'ylabel' not in kwargs.keys():
                                kwargs['ylabel'] = coord_labels[1]
                                ylabel_added = True
                            if 'ylims' not in kwargs.keys():
                                kwargs['ylims'] = coord_lims[1]
                                ylims_added = True
                            if 'yticklabels' not in kwargs.keys():
                                kwargs['yticklabels'] = coord_ticks[1]
                                yticklabels_added = True

                            individual_field_contour_plot(coords_X, coords_Y, field_data,
                                                          testname=testname, plotname=plotname, time=time,
                                                          field_name=field, slice_name=slice_name,
                                                          slice_idx=slice_idx, slice_label=slice_label,
                                                          **kwargs)

                            # Remove any added options
                            if xlabel_added:
                                kwargs.pop('xlabel')
                                xlabel_added = False
                            if xlims_added:
                                kwargs.pop('xlims')
                                xlims_added = False
                            if xticklabels_added:
                                kwargs.pop('xticklabels')
                                xticklabels_added = False
                            if ylabel_added:
                                kwargs.pop('ylabel')
                                ylabel_added = False
                            if ylims_added:
                                kwargs.pop('ylims')
                                ylims_added = False
                            if yticklabels_added:
                                kwargs.pop('yticklabels')
                                yticklabels_added = False
                            if extra_field_added:
                                kwargs.pop('extra_field_data')
                                extra_field_added = False

            # Slices is None so don't loop through slices
            else:

                if plotdir is None:
                    plotdir = 'results/'+dirname+'/figures'
                plotname = '%s/%s_%s_run_%s_time_%s.png' % (plotdir, testname, field,
                                                            str(run_id), str(time_idx))

                if topological_dim == 1:

                    coords, field_data, time, \
                    coord_label, coord_lims, \
                    coord_ticks, slice_label = extract_1D_data(data_file, field, time_idx,
                                                               slice_name='x', slice_idx=None)

                    # Add plotting details to kwargs if they are not already there
                    if 'xlabel' not in kwargs.keys():
                        kwargs['xlabel'] = coord_label
                        xlabel_added = True
                    if 'xlims' not in kwargs.keys():
                        kwargs['xlims'] = coord_lims
                        xlims_added = True
                    if 'xticklabels' not in kwargs.keys():
                        kwargs['xticklabels'] = coord_ticks
                        xticklabels_added = True

                    individual_field_1d_plot(coords, field_data, testname=testname, time=time,
                                             plotname=plotname, field_name=field,
                                             slice_name=None, slice_idx=None, **kwargs)

                    # Remove any added options to kwargs
                    if xlabel_added:
                        kwargs.pop('xlabel')
                        xlabel_added = False
                    if xlims_added:
                        kwargs.pop('xlims')
                        xlims_added = False
                    if xticklabels_added:
                        kwargs.pop('xticklabels')
                        xticklabels_added = False

                # Then we must be 2D and no slicing is required
                else:
                    coords_X, coords_Y, field_data, time, \
                    coord_labels, coord_lims, coord_ticks, \
                    slice_label = extract_2D_data(data_file, field, time_idx,
                                                  slice_name=None, slice_idx=None)

                    # Extract second field data if we need it
                    if 'extra_field_name' in kwargs.keys():
                        coords_X, coords_Y, extra_field_data, time, \
                        coord_labels, coord_lims, coord_ticks, \
                        slice_label = extract_2D_data(data_file, kwargs['extra_field_name'], time_idx,
                                                      slice_name=None, slice_idx=None)
                        kwargs['extra_field_data'] = extra_field_data
                        extra_field_added = True


                    # Add plotting details to kwargs if they are not already there
                    if 'xlabel' not in kwargs.keys():
                        kwargs['xlabel'] = coord_labels[0]
                        xlabel_added = True
                    if 'xlims' not in kwargs.keys():
                        kwargs['xlims'] = coord_lims[0]
                        xlims_added = True
                    if 'xticklabels' not in kwargs.keys():
                        kwargs['xticklabels'] = coord_ticks[0]
                        xticklabels_added = True
                    if 'ylabel' not in kwargs.keys():
                        kwargs['ylabel'] = coord_labels[1]
                        ylabel_added = True
                    if 'ylims' not in kwargs.keys():
                        kwargs['ylims'] = coord_lims[1]
                        ylims_added = True
                    if 'yticklabels' not in kwargs.keys():
                        kwargs['yticklabels'] = coord_ticks[1]
                        yticklabels_added = True

                    individual_field_contour_plot(coords_X, coords_Y, field_data,
                                                  testname=testname, time=time, plotname=plotname,
                                                  field_name=field, slice_name=None,
                                                  slice_idx=None, **kwargs)

                    # Remove any added options to kwargs
                    if xlabel_added:
                        kwargs.pop('xlabel')
                        xlabel_added = False
                    if xlims_added:
                        kwargs.pop('xlims')
                        xlims_added = False
                    if xticklabels_added:
                        kwargs.pop('xticklabels')
                        xticklabels_added = False
                    if ylabel_added:
                        kwargs.pop('ylabel')
                        ylabel_added = False
                    if ylims_added:
                        kwargs.pop('ylims')
                        ylims_added = False
                    if yticklabels_added:
                        kwargs.pop('yticklabels')
                        yticklabels_added = False
                    if extra_field_added:
                        kwargs.pop('extra_field_data')
                        extra_field_added = False

    data_file.close()


def make_convergence_plots(dirname, variables, fields, run_ids, errors,
                           titles=None, xlabels=None, ylabels=None,
                           xlims=None, ylims=None, **kwargs):
    """
    A routine for controlling the auto-generation of convergence plots
    from errors stored in global netCDF files.

    :arg dirname:           the name of the directory with the data
    :arg variables:         a string or list of strings giving the names of the
                            variable that we measure convergence with respect to.
                            e.g. 'dx' or ['dz', 'dt']
    :arg fields:            a list of the fields to be plotted
    :arg run_ids:           a list of lists IDs for the runs to be plotted. This
                            should be a list of run IDs for each field. A list
                            can be replaced by a single integer.
    :arg errors:            a list of the error diagnostics to be plotted
    :arg best_fit:          (Optional) Boolean to plot best fit lines through
                            the error points. If True then the gradients of
                            these lines will be added to the legend labels.
                            Default is True.
    :arg titles:            (Optional) a list of titles for the plots. Default is None.
                            Must be a list of lists, of length variables by errors.
    :arg xlabels:           (Optional) a list of labels for the x-axis. If not
                            provided then this will be determined from 'variables'.
                            Should be of same length as the list of variables.
    :arg ylabels:           (Optional) a list of labels for the y-axis. If not
                            provided then this will be determined from 'error'.
                            Should be of the same length as the list of errors.
    :arg xlims:             (Optional) a list of limits for the x-axis.
                            If not given then just uses the pyplot defaults.
                            Should be of same length as the list of variables.
    :arg ylims:             (Optional) a list of limits for the y-axis.
                            If not given then just uses the pyplot defaults.
                            Should be of same length as the list of errors.
    **kwargs:               Other optional arguments to pass through to making
                            the individual convergence plot.
    """

    # -------------------------------------------------------------------------#
    # Checks and convert variables into appropriate forms
    #--------------------------------------------------------------------------#

    # Turn errors and variables into list if only one string is provided
    if isinstance(errors, str):
        errors = [errors]
    if isinstance(variables, str):
        variables = [variables]


    # Optional arguments: axis labels and limits
    if xlabels is not None:
        if isinstance(xlabels, str):
            xlabels = [xlabels]
        if len(xlabels) != len(variables):
            raise ValueError('Length %d of xlabels not equal to length %d of variables' %
                             (len(xlabels), len(variables)))

    if xlims is not None:
        if isinstance(xlims, list) and not isinstance(xlims[0], list):
            # We have a list but not a list of lists
            xlims = [xlims]
        if len(xlims) != len(variables):
            raise ValueError('Length %d of xlims not equal to length %d of variables' %
                             (len(xlims), len(variables)))

    if ylabels is not None:
        if isinstance(ylabels, str):
            ylabels = [ylabels]
        if len(ylabels) != len(errors):
            raise ValueError('Length %d of ylabels not equal to length %d of errors' %
                             (len(ylabels), len(errors)))

    if ylims is not None:
        if isinstance(ylims, list) and not isinstance(ylims[0], list):
            # We have a list but not a list of lists
            ylims = [ylims]
        if len(ylims) != len(errors):
            raise ValueError('Length %d of ylims not equal to length %d of errors' %
                             (len(ylims), len(errors)))

    if titles is not None:
        if isinstance(titles, str):
            titles = [[titles]]

        if len(titles) != len(variables):
            raise ValueError('Length %d of titles not equal to length %d of variables' %
                             (len(titles), len(variables)))
        for i, title in enumerate(titles):
            if not isinstance(title, list):
                raise ValueError('Title must be a list of lists or single string')
            if len(title) != len(errors):
                raise ValueError('Length %d of title not equal to length %d of errors' %
                                 (len(title), len(errors)))

    #--------------------------------------------------------------------------#
    # Loop to produce plots
    #--------------------------------------------------------------------------#

    for i, variable in enumerate(variables):
        for j, error in enumerate(errors):

            # Extract appropriate axis labels and limits -- which could be None
            xlabel = xlabels[i] if xlabels is not None else None
            xlim = xlims[i] if xlims is not None else None
            ylabel = ylabels[j] if ylabels is not None else None
            ylim = ylims[j] if ylims is not None else None
            title = titles[i][j] if titles is not None else None

            individual_convergence_plot(dirname, variable, fields, run_ids, error,
                                        title=title, xlabel=xlabel, ylabel=ylabel,
                                        xlim=xlim, ylim=ylim, **kwargs)


def make_time_series_plots(dirname, fields, run_ids, diagnostics,
                           titles=None, ylabels=None, ylims=None, **kwargs):
    """
    A routine for controlling the auto-generation of time series plots from
    diagnostics stored in global netCDF files.

    :arg dirname:           the name of the directory with the data
    :arg fields:            a list of the fields to be plotted
    :arg run_ids:           a list of IDs for the runs to be plotted. This
                            should be of same length as the list of fields.
                            Can be replaced by a single integer.
    :arg diagnostics:       a list of the global diagnostics to be plotted
    :arg titles:            (Optional) a list of titles for the plots. Default is None.
                            Should be of same length as diagnostics.
    :arg ylabels:           (Optional) a list of labels for the y-axis. If not
                            provided then this will be determined from 'diagnostic'.
                            Should be of the same length as the list of diagnostics.
    :arg ylims:             (Optional) a list of limits for the y-axis.
                            If not given then just uses the pyplot defaults.
                            Should be of same length as the list of diagnostics.
    **kwargs                Other keyword arguments that can be passed through
                            to make time series plot.
    """

    # -------------------------------------------------------------------------#
    # Checks and convert variables into appropriate forms
    #--------------------------------------------------------------------------#

    # Turn diagnostics into list if only one string is provided
    if isinstance(diagnostics, str):
        diagnostics = [diagnostics]

    if ylabels is not None:
        if isinstance(ylabels, str):
            ylabels = [ylabels]
        if len(ylabels) != len(diagnostics):
            raise ValueError('Length %d of ylabels not equal to length %d of diagnostics' %
                             (len(ylabels), len(diagnostics)))

    if ylims is not None:
        if isinstance(ylims, list) and not isinstance(ylims[0], list):
            # We have a list but not a list of lists
            ylims = [ylims]
        if len(ylims) != len(diagnostics):
            raise ValueError('Length %d of ylims not equal to length %d of diagnostics' %
                             (len(ylims), len(diagnostics)))

    if titles is not None:
        if isinstance(titles, str):
            titles = [titles]

        if len(titles) != len(diagnostics):
            raise ValueError('Length %d of titles not equal to length %d of diagnostics' %
                             (len(titles), len(diagnostics)))

    #--------------------------------------------------------------------------#
    # Loop to produce plots
    #--------------------------------------------------------------------------#

    for i, diagnostic in enumerate(diagnostics):

        ylabel = ylabels[i] if ylabels is not None else None
        ylim = ylims[i] if ylims is not None else None
        title = titles[i] if titles is not None else None

        individual_time_series_plot(dirname, fields, run_ids, diagnostic,
                                    title=title, ylabel=ylabel, ylim=ylim, **kwargs)
