"""
Routine for plotting a 2D vector-valued field using "quivers" (arrows).
"""
import matplotlib.pyplot as plt
import numpy as np

__all__ = ["plot_field_quivers"]

def plot_field_quivers(ax, coords_X, coords_Y, field_data_X, field_data_Y,
                       projection=None, quiver_npts=1,
                       x_offset=None, y_offset=None,
                       units='xy',
                       scale=None, angles='xy', scale_units='xy',
                       restrict_quivers=False):

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
    # Slice data
    #--------------------------------------------------------------------------#

    if len(np.shape(field_data_X)) > 1:
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

        coords_X_to_plot = coords_X[x_slice, y_slice]
        coords_Y_to_plot = coords_Y[x_slice, y_slice]
        field_data_X_to_plot = field_data_X[x_slice, y_slice]
        field_data_Y_to_plot = field_data_Y[x_slice, y_slice]

    else:
        # TODO: this shouldn't be x_offset
        data_slice = slice(x_offset, None, quiver_npts)
        coords_X_to_plot = coords_X[data_slice]
        coords_Y_to_plot = coords_Y[data_slice]
        field_data_X_to_plot = field_data_X[data_slice]
        field_data_Y_to_plot = field_data_Y[data_slice]

    #--------------------------------------------------------------------------#
    # Plot quivers
    #--------------------------------------------------------------------------#

    if projection is None:
        # separately handle this case as None transform results in no arrows
        qv = ax.quiver(coords_X_to_plot, coords_Y_to_plot,
                       field_data_X_to_plot, field_data_Y_to_plot,
                       units=units, scale=scale, scale_units=scale_units,
                       angles=angles, zorder=3)
    else:
        qv = ax.quiver(coords_X_to_plot, coords_Y_to_plot,
                       field_data_X_to_plot, field_data_Y_to_plot,
                       units=units, scale=scale, scale_units=scale_units,
                       angles=angles, zorder=3, transform=projection)

    #--------------------------------------------------------------------------#
    # Finish
    #--------------------------------------------------------------------------#

    return qv
