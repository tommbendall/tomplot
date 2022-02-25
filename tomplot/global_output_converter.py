"""
Convert LFRic logged data into the transportdrake global output format.
"""

import numpy as np
import pandas as pd
from netCDF4 import Dataset

all_errors = ['Min-initial', 'Max-initial', 'L2-initial', 'L2-final',
              'Rel-L2-error', 'Dissipation', 'Dispersion']

def convert_global_output(target_dir, source_dirs, dxs=None,
                          mode='transport_stats'):

    # ------------------------------------------------------------------------ #
    # Checks
    # ------------------------------------------------------------------------ #

    if type(source_dirs) is not list:
        source_dirs = [source_dirs]
        dxs = [dxs]
    else:
        if type(dxs) is not list:
            raise ValueError(f'If source_dirs argument is a list, dxs should '+
                             f'also be a list and not type {type(dxs)}')
        elif len(dxs) != len(source_dirs):
            raise ValueError(f'List of dxs should be of length {len(dxs)} '+
                             'to match length of source_dirs')

    if mode not in ['transport_stats','gungho_mass']:
        raise ValueError('Converter mode should be "transport_stats" or "gungho"')

    # ------------------------------------------------------------------------ #
    # Read data into pandas dataframe
    # ------------------------------------------------------------------------ #

    if mode == 'transport_stats':
        col_names = ['day','timestamp','log_level','log_file','str1','str2',
                     'str3','str4','time','measure_name','variable','str5',
                     'measure_value']
    elif mode == 'gungho_mass':
        col_names = ['day','timestamp','log_level','log_file','str1','str2',
                     'str3','str4','str5','species_str','timestage1','timestage2',
                     'str6','str7','str8','timestep', 'mass']
    else:
        raise NotImplementedError(f'mode {mode} not implemented')

    data_frame_list = []

    for source_dir in source_dirs:
        file_name = f'{source_dir}/raw_data/output.log'
        print(f'Reading in {file_name}')
        # read in raw data from log file
        # skipinitialspace=True means that bunches of spaces are delimiters
        raw_data = pd.read_csv(file_name, header=None, sep=' ',
                               skipinitialspace=True, names=col_names)

        # -------------------------------------------------------------------- #
        # Data manipulations
        # -------------------------------------------------------------------- #

        if mode == 'gungho_mass':
            # Species are currently a string 'i,'. Lose the comma
            species = [str(species_str[0]) for species_str in raw_data['species_str'].values]
            raw_data['species'] = species
            # Convert mass data to floats
            raw_data['mass'] = raw_data['mass'].astype(float)
            # Combine point in timestep values into single value
            point_in_timestep = [f'{point1}_{point2}' for (point1, point2)
                                 in zip(raw_data['timestage1'].values,
                                        raw_data['timestage2'].values)]
            raw_data['timestage'] = point_in_timestep
            raw_data['time'] = raw_data['timestep'].astype(float)
            # Change values before initial time step to be after time step zero
            raw_data.loc[raw_data.timestage == 'Before_timestep', 'time'] = 0.0
            raw_data.loc[raw_data.timestage == 'Before_timestep', 'timestage'] = 'After_timestep'

        elif mode == 'transport_stats':
            # No manipulations needed for now
            pass
        else:
            raise NotImplementedError(f'mode {mode} not implemented')

        data_frame_list.append(raw_data)

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
    output_data.variables['run_id'][:] = range(len(data_frame_list))

    # Create groups for variables
    if mode == 'gungho_mass':

        # Add time variable
        for i, df in enumerate(data_frame_list):
            for j, t in enumerate(df.time.unique()):
                output_data.variables['time'][i,j] = float(t)

        # Extract unique species and measures
        unique_species = data_frame_list[0].species.unique()
        unique_measures = data_frame_list[0].timestage.unique()

        # Loop through species, adding groups
        for species in unique_species:
            output_data.createGroup(species)
            output_data[species].createGroup('errors')
            output_data[species].createGroup('global_quantities')

            # Loop through measures and add values
            for measure in unique_measures:
                output_data[species]['global_quantities'].createVariable(measure, float, ('run_id', 'time'))

                for i, df in enumerate(data_frame_list):
                    data_table = df[(df.species == species) & (df.timestage == measure)]
                    # Sort by time step and extract mass
                    data = data_table.sort_values('time').mass.values
                    # Fill to the end
                    idx0 = len(df.time.unique()) - len(data)
                    output_data[species]['global_quantities'][measure][i][idx0:] = data

    elif mode == 'transport_stats':

        unique_variables = data_frame_list[0].variable.unique()

        # Loop through variables, adding groups
        for variable in unique_variables:
            output_data.createGroup(variable)

        output_data.createVariable('dx', float, ('run_id',))
        output_data.variables['dx'][:] = dxs

        # Loop through data and separate into error and global data frames
        for i, df in enumerate(data_frame_list):
            error_data = df[df.measure_name.isin(all_errors)]
            quant_data = df[~df.measure_name.isin(all_errors)]

            # ---------------------------------------------------------------- #
            # Error data
            # ---------------------------------------------------------------- #

            unique_variables = error_data.variable.unique()

            # Create groups for the first time
            if i == 0:
                # Errors
                times = error_data.sort_values('time').time.unique()

                output_data.createDimension('error_time', None)
                output_data.createVariable('error_time', float, ('run_id','error_time'))
                output_data['error_time'][i][:] = times

                for variable in unique_variables:
                    output_data[variable].createGroup('errors')

            # Loop through variables
            for variable in unique_variables:
                unique_measures = error_data[error_data.variable == variable].measure_name.unique()

                for measure in unique_measures:

                    # Create variables for the first time
                    if i == 0:
                        output_data[variable]['errors'].createVariable(measure, float, ('run_id', 'error_time'))

                    data_table = error_data[(error_data.variable == variable) &
                                            (error_data.measure_name == measure)]
                    # Sort by time step and extract mass
                    data = data_table.sort_values('time').measure_value.values
                    # Fill to the end
                    idx0 = len(error_data.time.unique()) - len(data)
                    output_data[variable]['errors'][measure][i][idx0:] = data

            # ---------------------------------------------------------------- #
            # Global quantity data
            # ---------------------------------------------------------------- #

            unique_variables = quant_data.variable.unique()
            times = quant_data.sort_values('time').time.unique()
            # TODO: find out how to do this properly
            # I couldn't work out how to get this to work with slices
            for j, t in enumerate(times):
                output_data['time'][i,j] = t

            # Loop through variables
            for variable in unique_variables:
                if i == 0:
                    output_data[variable].createGroup('global_quantities')

                unique_measures = quant_data[quant_data.variable == variable].measure_name.unique()

                for measure in unique_measures:

                    # Create variables for the first time
                    if i == 0:
                        output_data[variable]['global_quantities'].createVariable(measure, float, ('run_id', 'time'))

                    data_table = quant_data[(quant_data.variable == variable) &
                                            (quant_data.measure_name == measure)]
                    # Sort by time step and extract mass
                    data = data_table.sort_values('time').measure_value.values
                    # Fill to the end
                    idx0 = len(quant_data.time.unique()) - len(data)
                    output_data[variable]['global_quantities'][measure][i][idx0:] = data

    else:
        raise NotImplementedError(f'mode {mode} not implemented')

    output_data.close()