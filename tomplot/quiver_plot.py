"""
Routine for plotting a 2D vector-valued field using "quivers" (arrows).
"""
import numpy as np
import pandas as pd

__all__ = ["plot_field_quivers"]

def plot_field_quivers(ax, coords_X, coords_Y, field_data_X, field_data_Y,
                       projection=None, magnitude_filter=None,
                       spatial_filter_step=None, spatial_filter_offset=None,
                       **quiver_kwargs):
    """
    Plots a 2D vector-valued field using quivers. The data can be filtered to
    make plots clearer.

    Keyword arguments can be passed to the underlying quiver plot. There are a
    handful of special keyword arguments whose values will be set by default.
    These are: `units`='xy', `scale_units='xy', `angles`='xy` and `zorder`=3.
    Common kwargs are:
    * scale: Scales the length of the arrow inversely. Number of data units per
      arrow length unit, e.g., m/s per plot width; a smaller scale parameter
      makes the arrow longer. Default is None.
    * minlength: Minimum length as a multiple of shaft width; if an arrow length
      is less than this, plot a dot (hexagon) of this diameter instead.

    Args:
        ax (:class:`AxesSubplot`): the pyplot axes object to plot the field on.
        coords_X (`numpy.ndarray`): an array containing the coordinates of the
            field data, corresponding to the plot's X axis. The shape of this
            array must correspond to that of the field data.
        coords_Y (`numpy.ndarray`): an array containing the coordinates of the
            field data, corresponding to the plot's Y axis.  The shape of this
            array must correspond to that of the field data.
        field_data_X (`numpy.ndarray`): the X-component of the vector field to
            be plotted. The shape of this array must match that of the coords.
        field_data_Y (`numpy.ndarray`): the Y-component of the vector field to
            be plotted. The shape of this array must match that of the coords.
        projection (:class:`Projection`, optional): a cartopy projection object.
            Defaults to None.
        magnitude_filter (float, optional): Filters out vectors whose magnitude
            is below this value. Defaults to None, in which case the filter is
            not applied.
        spatial_filter_step (int or tuple, optional): specifies what step to
            take when applying a spatial filter to the quivers. For instance,
            a step of 3 means that only one in 3 data points will be plotted.
            This can also be a tuple of length 2 for 2D data, so that values
            correspond to the steps in the X and Y direction. Defaults to None.
        spatial_filter_offset (int or tuple, optional): specifies an offset when
            applying a spatial filter to the quivers. This can also be a tuple
            of length 2 for 2D data, so that values correspond to the steps in
            the X and Y direction. Defaults to None.

    Raises:
        TypeError: _description_

    Returns:
        `matplotlib.pyplot.Quiver`: the resulting `Quiver` object.
    """

    #--------------------------------------------------------------------------#
    # Checks
    #--------------------------------------------------------------------------#

    # All arrays must have the same shape. Check this here to make assumptions
    # easier in what follows
    assert np.shape(coords_X) == np.shape(coords_Y), 'Quiver plot: shape of ' \
        + 'coordinate data must be equal'
    assert np.shape(field_data_X) == np.shape(field_data_Y), 'Quiver plot: ' \
        + 'shape of components for field data must be equal'
    assert np.shape(coords_X) == np.shape(field_data_X), 'Quiver plot: shape ' \
        + 'of field data must match shape of coordinate data'

    data_is_1D = (len(np.shape(coords_X)) == 1)

    #--------------------------------------------------------------------------#
    # Handle keyword arguments
    #--------------------------------------------------------------------------#

    # Make copy of keyword dictionary
    local_kwargs = dict(quiver_kwargs)

    # Special keyword arguments that we want to set
    if 'units' in local_kwargs.keys():
        units = local_kwargs['units']
        del local_kwargs['units']
    else:
        units = 'xy'

    if 'angles' in local_kwargs.keys():
        angles = local_kwargs['angles']
        del local_kwargs['angles']
    else:
        angles = 'xy'

    if 'scale_units' in local_kwargs.keys():
        scale_units = local_kwargs['scale_units']
        del local_kwargs['scale_units']
    else:
        scale_units = 'xy'

    if 'zorder' in local_kwargs.keys():
        zorder = local_kwargs['zorder']
        del local_kwargs['zorder']
    else:
        zorder = 3

    #--------------------------------------------------------------------------#
    # Apply spatial filter, if specified
    #--------------------------------------------------------------------------#

    if spatial_filter_step is not None:
        # Filtering 1D data
        if data_is_1D:
            assert type(spatial_filter_step) is int, 'Quiver plot: To apply ' \
                + 'spatial filter to 1D data, "spatial_filter_step" must be ' \
                + f'an integer, not {type(spatial_filter_offset)}'

            if spatial_filter_offset is None:
                spatial_filter_offset = 0

            # Make slice objects for array
            data_slice = slice(spatial_filter_offset, None, spatial_filter_step)

            # Slice data
            coords_X = coords_X[data_slice]
            coords_Y = coords_Y[data_slice]
            field_data_X = field_data_X[data_slice]
            field_data_Y = field_data_Y[data_slice]

        # Filtering 2D data
        else:
            if type(spatial_filter_step) is int:
                spatial_filter_step = (spatial_filter_step, spatial_filter_step)
            assert type(spatial_filter_step) is tuple, 'Quiver plot: To ' \
                + 'apply spatial filter to 2D data, "spatial_filter_step" ' \
                + f'must be a tuple or integer, not {type(spatial_filter_offset)}'

            if spatial_filter_offset is None:
                spatial_filter_offset = (0, 0)
            elif type(spatial_filter_offset) is int:
                spatial_filter_offset = (spatial_filter_offset, spatial_filter_offset)

            # Make slice objects for array
            data_slice_x = slice(spatial_filter_offset[0], None, spatial_filter_step[0])
            data_slice_y = slice(spatial_filter_offset[1], None, spatial_filter_step[1])

            # Slice data
            coords_X = coords_X[data_slice_x, data_slice_y]
            coords_Y = coords_Y[data_slice_x, data_slice_y]
            field_data_X = field_data_X[data_slice_x, data_slice_y]
            field_data_Y = field_data_Y[data_slice_x, data_slice_y]

    #--------------------------------------------------------------------------#
    # Apply filter based on magnitude, if specified
    #--------------------------------------------------------------------------#

    if magnitude_filter is not None:
        # Make a Dataframe object to quickly filter data
        field_data_mag = np.sqrt(field_data_X**2 + field_data_Y**2)

        df = pd.DataFrame({'coords_X': coords_X.flatten(),
                           'coords_Y': coords_Y.flatten(),
                           'field_data_X': field_data_X.flatten(),
                           'field_data_Y': field_data_Y.flatten(),
                           'field_data_mag': field_data_mag.flatten()})

        # Filter data based on vector magnitude
        filtered_df = df[df['field_data_mag'] > magnitude_filter]

        coords_X = filtered_df['coords_X'].values
        coords_Y = filtered_df['coords_Y'].values
        field_data_X = filtered_df['field_data_X'].values
        field_data_Y = filtered_df['field_data_Y'].values

    #--------------------------------------------------------------------------#
    # Plot quivers
    #--------------------------------------------------------------------------#

    if projection is None:
        # separately handle this case as None transform results in no arrows
        qv = ax.quiver(coords_X, coords_Y, field_data_X, field_data_Y,
                       units=units, scale_units=scale_units, angles=angles,
                       zorder=zorder, **local_kwargs)
    else:
        qv = ax.quiver(coords_X, coords_Y, field_data_X, field_data_Y,
                       units=units, scale_units=scale_units, angles=angles,
                       zorder=zorder, transform=projection, **local_kwargs)

    #--------------------------------------------------------------------------#
    # Finish
    #--------------------------------------------------------------------------#

    return qv
