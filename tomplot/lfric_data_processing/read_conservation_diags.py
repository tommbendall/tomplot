import numpy as np
import pandas as pd

moisture_dict = {0: 'dry_mass', 1:'vapour_mass', 2:'cloud', 3:'rain', 4:'ice', 5:'snow', 6:'graupel'}
nummr = 6

def make_conservation_dataframe(file_name, data_format, timestep_method='safe'):

    # ------------------------------------------------------------------------ #
    # Declare details of how to extract data for different formats
    # ------------------------------------------------------------------------ #

    # Store columns to extract and their names in a dictionary

    if data_format == 'dry_mass':
        # Data structure is:
        # log_time ":" "Conservation:" "dry" "mass" before/after time_point value
        column_indices = {'before/after':5, 'basic_time_point':6, 'value':7}
        to_compress = None
        value_name = 'dry_mass'

    elif data_format == 'dry_diags':
        # Data structure is:
        # log_time ":" "Conservation:" "dry" variable value
        column_indices = {'variable':4, 'value':5}
        to_compress = 'variable'
        value_name = None

    elif data_format == 'moist_mass':
        # Data structure is:
        # log_time ":" "Conservation:" "mass" "of" "water" "species" species before/after time_point value
        column_indices = {'species_raw':7, 'before/after':8, 'basic_time_point':9, 'value':10}
        to_compress = 'species'
        value_name = moisture_dict

    elif data_format == 'moist_flux':
        # Data structure is:
        # log_time ":" "Conservation:" scheme "moisture" "flux" in/out value
        column_indices = {'scheme':3, 'value':7}
        to_compress = 'scheme'
        value_name = None
    else:
        raise ValueError(f'data_format arg {data_format} not recognised')


    # ------------------------------------------------------------------------ #
    # Read in raw data as a "CSV"
    # ------------------------------------------------------------------------ #
    # key point is skipinitialspace=True, which means that not every space is a delimiter (only bunches of spaces)
    raw_data = pd.read_csv(file_name, header=None, sep=' ', skipinitialspace=True, usecols=column_indices.values())

    # Rename columns to human-readable values
    # First make inverse dictionary of columns
    column_titles = {}
    for key, value in column_indices.items():
        column_titles[value] = key
    raw_data = raw_data.rename(columns=column_titles)

    # ------------------------------------------------------------------------ #
    # Clean up data
    # ------------------------------------------------------------------------ #

    # Prevent floats randomly being read as strings ----------------------------
    print('Converting data values to float')
    raw_data['value'] = raw_data['value'].astype(float)

    # Species are currently a string 'i,'. Convert to integer i ----------------
    if data_format == 'moist_mass':
        print('Cleaning species values')
        species = [int(species_str[0]) for species_str in raw_data['species_raw'].values]
        raw_data.insert(1, 'species', species)
        raw_data = raw_data.drop(columns=['species_raw'])

    # Combine the "point in timestep" to give one value rather than two
    if data_format in ['dry_mass', 'moist_mass']:
        print('Combining point_in_timestep data')
        point_in_timestep = [string1+'_'+string2 for (string1, string2) in zip(raw_data['before/after'].values, raw_data['basic_time_point'].values)]
        raw_data.insert(0, 'point_in_timestep', point_in_timestep)
        raw_data = raw_data.drop(columns=['before/after', 'basic_time_point'])


    # ------------------------------------------------------------------------ #
    # Compress data, either by species or variable
    # ------------------------------------------------------------------------ #

    print('Compressing data')
    # NB: We require the compression to be such that each unique value has the
    # same number of rows
    if to_compress is not None:
        new_columns = np.unique(raw_data[to_compress].values)

        # Generate base data from the first of the unique values ---------------
        tmp_data = raw_data[raw_data[to_compress] == new_columns[0]]
        # Drop the compression column and rename the values
        tmp_data = tmp_data.drop(columns=[to_compress])
        if type(value_name) is dict:
            tmp_data = tmp_data.rename(columns={'value':value_name[new_columns[0]]})
        elif value_name is not None:
            tmp_data = tmp_data.rename(columns={'value':value_name})
        else:
            tmp_data = tmp_data.rename(columns={'value':new_columns[0]})

        # Start building up a dict of new data: ideally want a clean new df ----
        new_data_dict = {}
        for col in tmp_data.columns:
            new_data_dict[col] = tmp_data[col].values

        # Add the data for the other unique values -----------------------------
        for col in new_columns[1:]:
            if type(value_name) is dict:
                new_col_name = value_name[col]
            elif value_name is not None:
                new_col_name = value_name
            else:
                new_col_name = col
            new_data_dict[new_col_name] = raw_data[raw_data[to_compress] == col]['value'].values

        data = pd.DataFrame(new_data_dict)

    else:
        data = raw_data

    # ------------------------------------------------------------------------ #
    # Add time step data
    # ------------------------------------------------------------------------ #
    columns_to_check = list(filter(lambda col: (col != 'value' and col != to_compress), raw_data.columns))

    # This is done last to increase speed
    if timestep_method == 'quick' and 'point_in_timestep' in data.columns:
        print('Warning: generating time steps quickly, which may not be reliable')
        assert len(columns_to_check) > 0, 'Not enough columns to check'
        # Remove all 'Before_timestep' except the first
        # Find difference between two outputs that are the same from different time steps
        # 1) find first row of time step. Discard ones with "Before_timestep"
        print('Finding first row of first time step')
        start_index = -1
        for index, row in data.iterrows():
            if row['point_in_timestep'] != 'Before_timestep':
                start_index = index
                break
        assert start_index > -1, 'Cleaning data: unable to find a valid first row'
        start_row = data.iloc[[start_index]]
        print(f'First row of first time step is row {start_index}')

        print('Finding first row of second time step')
        second_index = -1
        for index, row in data.iterrows():
            if np.all([row[col] == start_row[col] for col in columns_to_check]) and index > start_index:
                second_index = index
                break
        assert second_index > start_index, 'Cleaning data: unable to find a valid second row'
        print(f'First row of second time step is row {second_index}')

        timestep_data = np.zeros(len(data), dtype=int)
        timestep_data[start_index:] = 1
        data.insert(0, 'timestep', timestep_data)

        # Remove all "Before_timestep" calls except the first
        data = data[~((data['point_in_timestep'] == 'Before_timestep') & (data['timestep'] > 0))]

        # Generate proper time step data
        timestep_data = np.zeros(len(data), dtype=int)
        delta_index = second_index - start_index
        assert int(len(data) - start_index) % delta_index == 0
        num_timesteps = int((len(data) - start_index) / delta_index)
        for i in range(num_timesteps):
            timestep_data[i*delta_index+start_index:(i+1)*delta_index+start_index] = i+1
        data['timestep'] = timestep_data

    else:
        columns_to_check = list(filter(lambda col: (col != 'value' and col != to_compress), raw_data.columns))

        # Find difference between two outputs that are the same from different time steps
        # 1) find first row of time step. Discard ones with "Before_timestep"
        print('Finding first row of time step')
        start_index = -1
        if 'point_in_timestep' in data.columns:
            for index, row in data.iterrows():
                if row['point_in_timestep'] != 'Before_timestep':
                    start_index = index
                    break
        else:
            start_index = 0
        assert start_index > -1, 'Cleaning data: unable to find a valid first row'
        start_row = data.iloc[[start_index]]
        print(f'First row of time step is row {start_index}')

        # 2) loop through rows and increment time step when we find a row matching first row
        print('Generating time step data -- may take a while')
        step_ctr = 0
        timestep_data = np.zeros(len(data), dtype=int)
        
        print(f'Comparing against columns {columns_to_check}')
        if len(columns_to_check) > 0:
            for index, row in data.iterrows():
                if np.all([row[col] == start_row[col] for col in columns_to_check]):
                    step_ctr += 1
                timestep_data[index] = step_ctr
        else:
            timestep_data = range(1, len(data)+1)

        data.insert(0, 'timestep', timestep_data)
    print('Data cleaned')

    return data
