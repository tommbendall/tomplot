import numpy as np
from .plot_lfric_coords import get_lfric_coords_2d

def extract_lfric_2D_data(data_file, field_name, time_idx,
                          extrusion_details=None, slice_name=None,
                          slice_idx=None, num_points=None, central_lon=0.0,
                          plot_coords_1d=None, testname=None):

    #--------------------------------------------------------------------------#
    # Checks
    #--------------------------------------------------------------------------#

    # TODO: fill in checks here

    if slice_name not in ['xy','yz','xz', None]:
        raise ValueError('For 2D plots slice must be xy, yz or xz')

    # Work out what time variable is called
    try:
        time = data_file.variables['time'][time_idx]
        lfric_initial = False
    except KeyError:
        try:
            time = data_file.variables['time_instant'][time_idx]
            lfric_initial = False
        except KeyError:
            # We have no time variable so this is the lfric_initial file
            time = 0
            lfric_initial = True

    #--------------------------------------------------------------------------#
    # Construct coordinate points
    #--------------------------------------------------------------------------#

    # Get dimension names
    if lfric_initial:
        vert_placement = data_file.variables[field_name].dimensions[0]
        # these will begin with "n" so slice to get rid of this
        hori_placement = data_file.variables[field_name].dimensions[1][1:]
    else:
        vert_placement = data_file.variables[field_name].dimensions[1]
        # these will begin with "n" so slice to get rid of this
        hori_placement = data_file.variables[field_name].dimensions[2][1:]

    plot_coords, data_coords, interp_coords_2d, \
    coord_labels, coord_lims, coord_ticks, \
    slice_label = get_lfric_coords_2d(data_file, hori_placement, vert_placement,
                                      slice_name, slice_idx, num_points=num_points,
                                      extrusion_details=extrusion_details,
                                      central_lon=central_lon,
                                      plot_coords_1d=plot_coords_1d)

    #--------------------------------------------------------------------------#
    # Interpolate field data
    #--------------------------------------------------------------------------#

    from scipy.interpolate import griddata

    field_data = np.zeros_like(plot_coords[0])
    if lfric_initial:
        field_full = data_file.variables[field_name][:,:]
    else:
        field_full = data_file.variables[field_name][time_idx,:,:]

    if slice_name == 'xy':
        # Just do interpolation at that level
        field_data = griddata(data_coords[:,slice_idx], field_full[slice_idx], interp_coords_2d, method='linear')
        field_near = griddata(data_coords[:,slice_idx], field_full[slice_idx], interp_coords_2d, method='nearest')
        field_data[np.isnan(field_data)] = field_near[np.isnan(field_data)]

    else:
        # Need to loop through levels
        for level in range(np.shape(field_full)[0]):
            slice_data = griddata(data_coords, field_full[level], interp_coords_2d, method='linear')
            slice_near = griddata(data_coords, field_full[level], interp_coords_2d, method='nearest')
            slice_data[np.isnan(slice_data)] = slice_near[np.isnan(slice_data)]

            field_data[level] = slice_data

    return plot_coords[0], plot_coords[1], field_data, time, coord_labels, \
           coord_lims, coord_ticks, slice_label
