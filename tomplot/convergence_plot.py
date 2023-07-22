"""
This file provides tools for quickly making a neat convergence plot, to plot
error measurements against some independent variable so as to determine the
order of accuracy of a numerical method.
"""
import numpy as np
import warnings

__all__ = ["plot_convergence", "add_convergence_comparison_line"]


def plot_convergence(ax, input_data, error_data,
                     label=None, marker='+', color='black', markersize=None,
                     # Options relating to logarithms
                     log_by='data', log_base='e',
                     # Options relating to best fit lines
                     best_fit=True, best_fit_deg=1,
                     linestyle='-', linewidth=None,
                     gradient_in_label=True, gradient_format='.2f'):
    u"""
    Plots the values of some error measurements against their independent
    variable to determine a convergence rate.

    An example of how to use this routine is to try to determine how the order
    of accuracy of a numerical method. For instance, if the error of the scheme
    converges as 𝓞(Δx^2), the order of accuracy is the exponent 2. This can be
    determined from error measurements, by plotting the measurements as a
    function of Δx on a log-log plot. The gradient of a line of best fit then
    approximates the order of accuracy.

    Args:
        ax (:class:`AxesSubplot`): the pyplot axes object to plot the data on.
        input_data (`numpy.ndarray`): an array containing the x-axis data. This
            is likely to correspond to a parameter such as Δx or Δt.
        error_data (`numpy.ndarray`): an array containing the y-axis data, which
            corresponds to the error measurements to be plotted.
        label (str, optional): the label for this set of data, to appear in the
            legend. If the `gradient_in_label` argument is True, the gradient of
            a line of best fit will be appended to this label. Defaults to None.
        marker (str, optional): which type of marker to use for this dataset.
            Defaults to '+'.
        color (str, optional): which colour to use for the plotted points and
            any line of best fit. Defaults to 'black'.
        markersize (float, optional): the size of markers to use for plotting
            the data points. Defaults to None.
        log_by (str, optional): how the plot should handle logarithms. Options
            are: "data", in which log(error) is plotted against log(input), or
            "axes" in which the error is plotted against the input but the axes
            take a logarithmic scale. Defaults to "data".
        log_base (float, optional): the base to use for the logarithm. Should be
            a positive float, not equal to 1. Can also be the string "e", which
            is the default value.
        best_fit (bool, optional): whether to plot a line of best fit through
            the data points. Defaults to True.
        best_fit_deg (int, optional): the order of polynomial to use for the
            line of best fit. Defaults to 1.
        linestyle (str, optional): the linestyle to use for the line of best
            fit. Defaults to '-'.
        linewidth (float, optional): the linewidth to use for the line of best
            fit. Defaults to None.
        gradient_in_label (bool, optional): whether to include the gradient of
            the line of best fit in the label shown in the legend. If True, this
            is appended to the `label` argument. Defaults to True.
        gradient_format (str, optional): the format to use for the gradient
            that is shown in the label. Defaults to '.2f'.

    Raises:
        ValueError: if `gradient_in_label` is True, but `best_fit_deg` is not 1.
        NotImplementedError: plotting a line of best fit for a polynomial that
            is not 1 when using the log_by="axes" argument.
    """

    # ------------------------------------------------------------------------ #
    # Checks on arguments
    # ------------------------------------------------------------------------ #

    assert log_by in ['data', 'axes'], 'convergence plot: the "log_by" ' \
        + f'argument can take values of "data" or "axes", and not {log_by}'

    assert (log_base == 'e' or (type(log_base) in [float, int]
                                and log_base > 0 and log_base != 1)), \
        'plot_convergence: the "log_base" argument can take the value of "e" ' \
        + f'or that of a positive real number (not 1). Your value is {log_base}'

    assert type(best_fit_deg) is int, 'plot_convergence: the "best_fit_deg" ' \
        f'must be of integer type, not {type(best_fit_deg)}'

    if best_fit_deg != 1 and gradient_in_label:
        raise ValueError('plot_convergence: cannot use a "best_fit_deg" that is'
                         + 'not 1 and have gradient_in_label == True')

    if best_fit and best_fit_deg != 1 and log_by == 'axes':
        raise NotImplementedError('plot_convergence: cannot plot a line of best'
                                  + ' fit with a polynomial other than 1 when '
                                  + 'using the log_by="axes" argument')

    if log_by == 'axes' and log_base == 'e':
        warnings.warn('plot_convergence: "log_by" is set to "axes" but using a'
                      + 'log with base e. You may prefer to use log_base = 10')

    # ------------------------------------------------------------------------ #
    # Deal with logarithms
    # ------------------------------------------------------------------------ #
    # Turn log_base of "e" into a real number
    base = np.exp(1) if log_base == "e" else log_base
    log_input_data = np.emath.logn(base, input_data)
    log_error_data = np.emath.logn(base, error_data)

    # ------------------------------------------------------------------------ #
    # Line of best fit
    # ------------------------------------------------------------------------ #
    # This is done first to generate the gradient if it's needed in the label
    if best_fit:
        # Line of best fit using standard polynomial fitting
        best_fit_line = np.poly1d(np.polyfit(log_input_data, log_error_data, deg=best_fit_deg))

        # Read off the gradient of if best fit is linear
        if best_fit_deg == 1:
            intercept = best_fit_line[0]
            gradient = best_fit_line[1]

        if log_by == 'data':
            ax.plot(log_input_data, best_fit_line(log_input_data), label=None,
                    linestyle=linestyle, color=color, linewidth=linewidth)
        elif log_by == 'axes':
            # Plotting with actual error values but on log-log scale
            # Assume linear best fit here (otherwise already raised an error)
            # log(error) = gradient*log(input) + intercept
            # => error = base**intercept * input**gradient
            best_fit_data = base**intercept * input_data**gradient
            ax.loglog(input_data, best_fit_data, label=None, base=base,
                      linestyle=linestyle, color=color, linewidth=linewidth)

        # -------------------------------------------------------------------- #
        # Append gradient to label, if specified
        # -------------------------------------------------------------------- #
        if gradient_in_label:
            # Work out string for formatting gradient
            format_str = '{:'+gradient_format+'}'

            # If label is already None, gradient becomes label
            if label is None:
                label = f'{format_str.format(gradient)}'

            # Otherwise, we need to append the gradient to the label
            else:
                label = label + f', {format_str.format(gradient)}'

    # ------------------------------------------------------------------------ #
    # Plot data points, depending on log method
    # ------------------------------------------------------------------------ #

    if log_by == 'data':
        ax.plot(log_input_data, log_error_data, label=label, color=color,
                marker=marker, markersize=markersize, linestyle='')

    elif log_by == 'axes':
        ax.loglog(input_data, error_data, label=label, base=base,
                  color=color, marker=marker, markersize=markersize, linestyle='')

    return


def add_convergence_comparison_line(ax, exponent, x_points=None, y_shift=0.0,
                                    label=None, color='black', linestyle='--',
                                    linewidth=None, log_by='data', log_base='e'):
    """
    Adds a comparison line to a convergence plot.

    Args:
        ax (:class:`AxesSubplot`): the pyplot axes object to plot the data on.
        exponent (float): the exponent of the comparison line to be plotted,
            which corresponds to the gradient when plotting on a log-log plot.
        x_points (`numpy.ndarray`, optional): points through which to plot the
            comparison line. If not provided, the comparison line is chosen to
            go through the centre of the figure. Defaults to None.
        y_shift (float, optional): shift in the y-direction to apply to the
            line. Defaults to 0.0.
        label (str, optional): the label for this line, to appear in the legend
            if specified. Defaults to None.
        color (str, optional): which colour to use for the line. Defaults to
            'black'.
        linestyle (str, optional): the style of line to use. Defaults to '--',
            corresponding to a dashed line.
        log_by (str, optional): how the plot should handle logarithms. Options
            are: "data", in which log(error) is plotted against log(input), or
            "axes" in which the error is plotted against the input but the axes
            take a logarithmic scale. Defaults to "data".
        log_base (float, optional): the base to use for the logarithm. Should be
            a positive float, not equal to 1. Can also be the string "e", which
            is the default value.
    """

    # ------------------------------------------------------------------------ #
    # Checks on arguments
    # ------------------------------------------------------------------------ #

    assert log_by in ['data', 'axes'], 'convergence plot: the "log_by" ' \
        + f'argument can take values of "data" or "axes", and not {log_by}'

    assert (log_base == 'e' or (type(log_base) in [float, int]
                                and log_base > 0 and log_base != 1)), \
        'add_convergence_comparison: "log_base" argument be "e" ' \
        + f'or a positive real number (not 1). Your value is {log_base}'

    if log_by == 'axes' and log_base == 'e':
        warnings.warn('add_convergence_comparison: "log_by" is set to "axes" but'
                      + 'using log with base e. You may prefer log_base = 10')

    # ------------------------------------------------------------------------ #
    # Deal with logarithms
    # ------------------------------------------------------------------------ #
    # Turn log_base of "e" into a real number
    base = np.exp(1) if log_base == "e" else log_base

    # ------------------------------------------------------------------------ #
    # Work out where to plot the line
    # ------------------------------------------------------------------------ #

    if log_by == 'data':
        log_y_lower, log_y_upper = ax.get_ylim()[0], ax.get_ylim()[1]
        log_y_center = 0.5 * (log_y_lower + log_y_upper) + y_shift

        if x_points is None:
            # Get centre of figure
            log_x_lower, log_x_upper = ax.get_xlim()[0], ax.get_xlim()[1]
            log_x_center = 0.5 * (log_x_lower + log_x_upper)
            log_x_points = np.linspace(log_x_lower, log_x_upper, 20)
        else:
            log_x_points = np.emath.logn(base, x_points)
            log_x_center = 0.5 * (log_x_points[0] + log_x_points[-1])

    elif log_by == 'axes':
        y_lower, y_upper = ax.get_ylim()[0], ax.get_ylim()[1]
        log_y_lower = np.emath.logn(base, y_lower)
        log_y_upper = np.emath.logn(base, y_upper)
        log_y_center = 0.5*(log_y_lower + log_y_upper) + y_shift

        if x_points is None:
            # Get centre of figure
            x_lower, x_upper = ax.get_xlim()[0], ax.get_xlim()[1]
            log_x_lower = np.emath.logn(base, x_lower)
            log_x_upper = np.emath.logn(base, x_upper)
            x_points = np.linspace(x_lower, x_upper, 20)
            log_x_points = np.emath.logn(base, x_points)
        else:
            log_x_lower = np.emath.logn(base, x_points[0])
            log_x_upper = np.emath.logn(base, x_points[-1])
            log_x_points = np.emath.logn(base, x_points)
        log_x_center = 0.5*(log_x_lower + log_x_upper)

    # ------------------------------------------------------------------------ #
    # Work out line to plot in log space and plot it
    # ------------------------------------------------------------------------ #

    # Line going through centre of plot in log space
    intercept = - exponent * log_x_center + log_y_center
    log_y_points = exponent * log_x_points + intercept

    if log_by == 'data':
        # Simply plot the logs of the points
        ax.plot(log_x_points, log_y_points, label=label, color=color,
                linestyle=linestyle, linewidth=linewidth)

    elif log_by == 'axes':
        # Need to convert logs of points to curve
        # log(y) = gradient*log(x) + intercept
        # => y = base**intercept * x**gradient
        y_points = base**intercept * x_points**exponent

        ax.loglog(x_points, y_points, label=label, color=color,
                  linestyle=linestyle, linewidth=linewidth, base=base)

    return
