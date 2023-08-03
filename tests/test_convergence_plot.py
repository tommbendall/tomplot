"""
This tests the routine to make a convergence plot.
"""

import matplotlib.pyplot as plt
from tomplot import plot_convergence, tomplot_field_title
import numpy as np
import pytest

# Set up some data points to use: three different lines with different exponents
np.random.seed(19)
x_points = np.array([2.0, 3.5, 6.0])
y_points_1 = 10.0*x_points + np.random.randn(len(x_points)) + 5
y_points_2 = 10.0*x_points**2 + np.random.randn(len(x_points)) + 5
y_points_3 = 10.0*x_points**3 + np.random.randn(len(x_points)) + 5


@pytest.mark.parametrize("log_by", ['data', 'axes'])
@pytest.mark.parametrize("log_base", ['e', 10])
def test_convergence_plot(log_by, log_base, plot_setup):

    setup = plot_setup("none", 0.0, 0.0)

    plt.close()

    _, ax = plt.subplots(1, 1, figsize=(5, 5))

    plot_convergence(ax, x_points, y_points_1, label='linear',
                     color='blue', marker='s', log_by=log_by, log_base=log_base)
    plot_convergence(ax, x_points, y_points_2, label='quadratic',
                     color='red', marker='o', log_by=log_by, log_base=log_base)
    plot_convergence(ax, x_points, y_points_3, label='cubic',
                     color='purple', marker='^', log_by=log_by, log_base=log_base)

    plt.legend()

    title = f'convergence log {log_by} with log base {log_base}'
    plot_name = f'convergence_log_{log_base}_by_{log_by}.png'

    tomplot_field_title(ax, title)
    setup.make_plots(plot_name)
