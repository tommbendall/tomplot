"""
This tests the routine to place legends outside of the figure/axes.
"""

import matplotlib.pyplot as plt
from tomplot import tomplot_field_title, tomplot_legend_fig, tomplot_legend_ax
import numpy as np
import pytest

# Set up some data points to use: three different lines with different exponents
x_points = np.linspace(0.0, 1.0, 51)
y_sin = np.sin(2*np.pi*x_points)
y_cos = np.cos(2*np.pi*x_points)
y_exp = np.exp(x_points) - 2
y_sqrt = np.sqrt(x_points)


@pytest.mark.parametrize("add_legend_to", ['figure', 'axes'])
@pytest.mark.parametrize("legend_loc", ['top', 'bottom'])
def test_tomplot_legend(add_legend_to, legend_loc, plot_setup):

    plt.close()

    # Don't care about the setup
    setup = plot_setup("none", 0.0, 0.0)

    if add_legend_to == 'figure':
        fig, axarray = plt.subplots(1, 2, figsize=(7, 3), sharey='row')
    else:
        _, ax = plt.subplots(1, 1, figsize=(5, 5))

    title = f'Legend {add_legend_to} {legend_loc}'

    if add_legend_to == 'figure':
        axarray[0].plot(x_points, y_sin, color='purple', marker='', label='sin')
        axarray[0].plot(x_points, y_cos, color='red', marker='', label='cos')
        axarray[1].plot(x_points, y_exp, color='gold', marker='', label='exp')
        axarray[1].plot(x_points, y_sqrt, color='blue', marker='', label='sqrt')
        # Add legend
        tomplot_legend_fig(fig, legend_loc, ncols=4, title='Functions')
        fig.suptitle(title)

    else:
        ax.plot(x_points, y_sin, color='purple', marker='', label='sin')
        ax.plot(x_points, y_cos, color='red', marker='', label='cos')
        ax.plot(x_points, y_exp, color='gold', marker='', label='exp')
        ax.plot(x_points, y_sqrt, color='blue', marker='', label='sqrt')
        # Add legend
        tomplot_legend_ax(ax, legend_loc, ncols=4, title='Functions')
        tomplot_field_title(ax, title)

    plot_name = f'legend_{add_legend_to}_{legend_loc}.png'
    setup.make_plots(plot_name)
