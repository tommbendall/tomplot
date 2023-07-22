"""
This tests the routine to add a comparison line to a convergence plot.
"""

import matplotlib.pyplot as plt
from tomplot import (plot_convergence, add_convergence_comparison_line,
                     tomplot_field_title, set_tomplot_style)
import numpy as np
import pytest

# Set up some data points to use: three different lines with different exponents
np.random.seed(51)
x_points = np.array([2.0, 3.5, 6.0])
y_points_1 = 5.0*x_points + np.random.randn(len(x_points)) + 10
y_points_2 = 5.0*x_points**2 + np.random.randn(len(x_points)) + 10
y_points_3 = 5.0*x_points**3 + np.random.randn(len(x_points)) + 10


@pytest.mark.parametrize("log_by", ['data', 'axes'])
def test_add_convergence_comparison(log_by, plot_setup):

    setup = plot_setup("none", 0.0, 0.0)

    plt.close()

    set_tomplot_style(16)

    _, ax = plt.subplots(1, 1, figsize=(5, 5))

    log_base = 10

    # First add standard convergence lines
    plot_convergence(ax, x_points, y_points_1, label='linear',
                     color='blue', marker='s', log_by=log_by, log_base=log_base)
    plot_convergence(ax, x_points, y_points_2, label='quadratic',
                     color='red', marker='o', log_by=log_by, log_base=log_base)
    plot_convergence(ax, x_points, y_points_3, label='cubic',
                     color='purple', marker='^', log_by=log_by, log_base=log_base)

    add_convergence_comparison_line(ax, 1, label=r'$\Delta x$', color='blue',
                                    log_by=log_by, log_base=log_base)
    add_convergence_comparison_line(ax, 2, label=r'$\Delta x^2$', color='red',
                                    log_by=log_by, log_base=log_base)
    add_convergence_comparison_line(ax, 3, label=r'$\Delta x^3$', color='purple',
                                    log_by=log_by, log_base=log_base)

    plt.legend()

    title = f'convergence comparison log {log_by}'
    plot_name = f'convergence_comparison_log_by_{log_by}.png'

    tomplot_field_title(ax, title)
    setup.make_plots(plot_name)
