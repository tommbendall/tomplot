"""
Compute L2-error norms from LFRic field data. This is then stored in the
transportdrake global output format.
"""

import numpy as np
from netCDF4 import Dataset
from .diagnostic_info import DiagnosticInfo
from .data_coords import get_lfric_data_coords
from .extrusions import generate_extrusion

all_error_names = ['L2_error']

def compute_field_measures(target_dir, true_source_dir, source_dirs,
                           variables, measures, extrusion_details,
                           time_idx=-1, run_params=None):

    # ------------------------------------------------------------------------ #
    # Checks
    # ------------------------------------------------------------------------ #

    if run_params is None:
        run_params = {}
    elif type(run_params) is not dict:
        raise TypeError(f'run_params should be a dict not type {type(run_params)}')

    if type(source_dirs) is not list:
        source_dirs = [source_dirs]
        # Convert all params into lists
        for param_key, param_value in run_params.items():
            run_params[param_key] = [param_value]
    else:
        for param_key, param_value in run_params.items():
            if type(param_value) is not list:
                raise TypeError(f'If source_dirs argument is a list, each '+
                                f'param in run_params should also be a list, '+
                                f'but {param_key} is type {type(param_value)}')
            elif len(param_value) != len(source_dirs):
                raise ValueError(f'List for {param_key} should be of length '
                                 f'{len(source_dirs)} to match length of source_dirs')
    
    # Want to turn variables into a list if it isn't already
    if type(variables) is str or isinstance(variables, DiagnosticInfo):
        variables = [variables]

    if type(measures) is str:
        measures = [measures]

    # ------------------------------------------------------------------------ #
    # Create global_output.nc
    # ------------------------------------------------------------------------ #

    output_file_name = f'{target_dir}/global_output.nc'
    output_data = Dataset(output_file_name, 'w')

    # Create dimensions
    output_data.createDimension('run_id', None)
    output_data.createDimension('time', None)
    output_data.createVariable('run_id', int, ('run_id',))
    output_data.createVariable('time', float, ('run_id', 'time'))

    # Add run_id variable
    output_data.variables['run_id'][:] = range(len(source_dirs))
    output_data.variables['time'][0] = 1.0

    for param_key, param_value in run_params.items():
        # Assume that all members of param_value are the same type!
        output_data.createVariable(param_key, type(param_value[0]), ('run_id',))
        output_data.variables[param_key][:] = param_value

    # Add field variable groups
    for variable in variables:
        if type(variable) == str:
            variable_name = variable
        elif isinstance(variable, DiagnosticInfo):
            variable_name = variable.name
        else:
            raise ValueError(f'variable should be string or DiagnosticInfo '+
                             f'object, not {type(variable)}')
        output_data.createGroup(variable_name)
        output_data[variable_name].createGroup('errors')
        output_data[variable_name].createGroup('global_quantities')

        for measure in measures:
            if measure in all_error_names:
                output_data[variable_name]['errors'].createVariable(measure, float, ('run_id', 'time'))
            else:
                raise NotImplementedError(f'Measure {measure} not yet implemented, '+
                                          f'and currently no global quantities are implemented')

    output_data.close()

    # ------------------------------------------------------------------------ #
    # Compute quantities
    # ------------------------------------------------------------------------ #

    output_data = Dataset(output_file_name, 'a')

    for variable in variables:

        # -------------------------------------------------------------------- #
        # Extract true data
        # -------------------------------------------------------------------- #

        true_data_file = Dataset(f'{true_source_dir}/lfric_diag.nc', 'r')

        if type(variable) == str:
            var_for_placement = variable
            variable_name = variable
        elif isinstance(variable, DiagnosticInfo):
            # Use first prognostic variable to get staggering information
            var_for_placement = variable.prognostic_variables[0]
            variable_name = variable.name
        else:
            raise ValueError(f'variable should be string or DiagnosticInfo '+
                             f'object, not {type(variable)}')
        
        # Get name of variable placement
        vert_placement = true_data_file.variables[var_for_placement].dimensions[1]
        # these will begin with "n" so slice to get rid of this
        hori_placement = true_data_file.variables[var_for_placement].dimensions[2][1:]

        # Get field data
        if type(variable) == str:
            true_data = true_data_file.variables[variable][:,:,time_idx]
        elif isinstance(variable, DiagnosticInfo):
            variable.extract_lfric_data(true_data_file, time_idx)
            true_data = variable.evaluate()

        else:
            raise ValueError(f'variable should be string or DiagnosticInfo '+
                             f'object, not {type(variable)}')

        # Extract 2D coordinates
        data_coords_2d = get_lfric_data_coords(true_data_file, hori_placement,
                                               extrusion_details)
        # Turn these into 3D field
        # Shape of true_data is (num_levels, num_columns)
        vertical_coords_1d = generate_extrusion(extrusion_details, vert_placement,
                                                np.shape(true_data)[0])
        true_data_coords = np.zeros((3, np.shape(true_data)[0], np.shape(true_data)[1]))
        true_data_coords[0,:,:] = data_coords_2d[0]
        true_data_coords[1,:,:] = data_coords_2d[1]
        true_data_coords[2,:,:] = np.tile(vertical_coords_1d,
                                          (np.shape(data_coords_2d)[1], 1)).T

        true_data_file.close()

        # -------------------------------------------------------------------- #
        # Extract experiment data
        # -------------------------------------------------------------------- #

        for run_id, source_dir in enumerate(source_dirs):

            print(f'Computing error for {source_dir}')
            data_file = Dataset(f'{source_dir}/lfric_diag.nc', 'r')

            # Get field data
            if type(variable) == str:
                field_data = data_file.variables[variable][time_idx,:,:]
            elif isinstance(variable, DiagnosticInfo):
                variable.extract_lfric_data(data_file, time_idx)
                field_data = variable.evaluate()

            else:
                raise ValueError(f'variable should be string or DiagnosticInfo '+
                                f'object, not {type(variable)}')

            # Extract 2D coordinates
            data_coords_2d = get_lfric_data_coords(data_file, hori_placement,
                                                   extrusion_details)
            # Turn these into 3D field
            # Shape of true_data is (num_levels, num_columns)
            vertical_coords_1d = generate_extrusion(extrusion_details, vert_placement,
                                                    np.shape(field_data)[0])
            data_coords = np.zeros((3, np.shape(field_data)[0], np.shape(field_data)[1]))
            data_coords[0,:,:] = data_coords_2d[0]
            data_coords[1,:,:] = data_coords_2d[1]
            data_coords[2,:,:] = np.tile(vertical_coords_1d,
                                         (np.shape(data_coords_2d)[1], 1)).T

            # ---------------------------------------------------------------- #
            # Obtain true field on other mesh
            # ---------------------------------------------------------------- #

            from scipy.interpolate import griddata
            # reshape data to be agnostic to levels
            field_data = np.reshape(field_data, (np.shape(field_data)[0]*np.shape(field_data)[1],))
            data_coords = np.reshape(data_coords, (3, np.shape(data_coords)[1]*np.shape(data_coords)[2])).T
            # true data only needs doing first time
            if run_id == 0:
                true_data = np.reshape(true_data, (np.shape(true_data)[0]*np.shape(true_data)[1],))
                true_data_coords = np.reshape(true_data_coords, (3, np.shape(true_data_coords)[1]*np.shape(true_data_coords)[2])).T


            # interpolate to get coarse data
            # Nearest method much much quicker for 3D data
            true_coarse_data = griddata(true_data_coords, true_data, data_coords, method='nearest')

            # ---------------------------------------------------------------- #
            # Compute errors / measures
            # ---------------------------------------------------------------- #

            for measure in measures:

                if measure == 'L2_error':
                    # TODO: improve this calculation!
                    diff_data = field_data - true_coarse_data
                    value = np.linalg.norm(diff_data) / diff_data.size
                else:
                    raise NotImplementedError(f'Measure {measure} not implemented')

                output_data[variable_name]['errors'][measure][run_id,time_idx] = value

    output_data.close()