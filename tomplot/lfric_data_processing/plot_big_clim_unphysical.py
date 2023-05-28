"""
This script provides a routine for plotting the location of the unphysical min
or max values of certain fields, and how this location evolves with time.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import pandas as pd
import sys
import cartopy.crs as ccrs

# ---------------------------------------------------------------------------- #
# Options
# ---------------------------------------------------------------------------- #
include_u = False

# ---------------------------------------------------------------------------- #
# Routine for manipulating data
# ---------------------------------------------------------------------------- #
def manipulate_data(raw_data):
    """
    Routine to turn DataFrame of raw data into something manageable for
    plotting. This goes through the raw data line-by-line and combines all lines
    for a single time step into a single row of data.

    Mins and maxes of several parameters are logged every time step, but when
    they breach some threshold, the location of the mins and maxes are also
    logged. When the threshold is not reached, we fill the new DataFrame with
    NaNs.
    """
    # Need to reorganise all data to make a neat data frame
    cols = ['Timestep']
    for variable in ['dtheta_fast', 'u_in_w3', 'v_in_w3', 'w_in_wth','dtheta_slow']:
        for stat in ['min', 'max', 'maxabs']:
            cols.append(f'{variable}_{stat}')
        for stat in ['level', 'latitude', 'longitude']:
            cols.append(f'{variable}_min_{stat}')
            cols.append(f'{variable}_max_{stat}')
            cols.append(f'{variable}_maxabs_{stat}')

    new_data = pd.DataFrame(columns=cols)

    last_variable = ''
    timestep = 0
    tmp_data = {col:np.nan for col in cols}
    tmp_data['Timestep'] = 0
    maxabs_is_min = False

    # Loop through all rows of old data
    for _, row in raw_data.iterrows():
        # Sort rows by what is in the second column
        if row[1] == 'of':
            # This is a time step row -- add old data to dataframe
            if row[0] == 'End' or timestep == 0:
                tmp_data = pd.DataFrame(tmp_data, index=[timestep])
                new_data = pd.concat([new_data, tmp_data], ignore_index=True)
                timestep += 1
            else:
                # Make new data dictionary
                tmp_data = {col:np.nan for col in cols}
                tmp_data['Timestep'] = row[3]
        elif row[1] in ['dtheta_fast','u_in_w3','v_in_w3','w_in_wth','dtheta_slow']:
            last_variable = row[1]
            tmp_data[last_variable+'_min'] = row[3]
            tmp_data[last_variable+'_max'] = row[4]
            # Work out whether to save
            if abs(row[3]) > abs(row[4]):
                tmp_data[last_variable+'_maxabs'] = abs(row[3])
                maxabs_is_min = True
            else:
                tmp_data[last_variable+'_maxabs'] = abs(row[4])
                maxabs_is_min = False
            value = tmp_data[last_variable+'_maxabs']
            if type(value) is not float:
                value = float(value)
        elif row[1] in ['level','longitude','latitude']:
            if ((last_variable == 'u_in_w3' and value > 50.0) or
                 last_variable != 'u_in_w3'):
                tmp_data[f'{last_variable}_min_{row[1]}'] = row[3]
                tmp_data[f'{last_variable}_max_{row[1]}'] = row[4]
                # Get min/max location from latest min/max
                if maxabs_is_min:
                    tmp_data[f'{last_variable}_maxabs_{row[1]}'] = row[3]
                else:
                    tmp_data[f'{last_variable}_maxabs_{row[1]}'] = row[4]
            

    return new_data

# ---------------------------------------------------------------------------- #
# Routine for making figure
# ---------------------------------------------------------------------------- #
def make_figure(datapath, plotpath, branch):
    """
    Routine for making a single plot of the locations of unphysical max(abs)
    values of certain fields.

    The plot that is returned shows the evolution of these locations in both
    2D longitude-height slices (which form the top row of the plot), and in
    longitude-latitude slices (the bottom row of the plot).

    Each column of the plot corresponds to a different field.

    A dot is plotted for the initial location that the unphysical value is
    recorded.
    """
    fig_width = 21 if include_u else 14
    fig = plt.figure(figsize=(fig_width,8))
    plt.suptitle("Configuration = "+branch)
    variables = ['u_in_w3', 'w_in_wth', 'dtheta_fast']

    if not include_u:
        variables.pop(0)

    # Load up all data
    fname = datapath + "/" + branch + "/unphysical.txt"
    raw_data = pd.read_csv(fname, header=None, sep=' ',skipinitialspace=True, names=range(5))
    # Manipulate data
    data = manipulate_data(raw_data)

    # Loop through variables -- each will have its own column
    for k, variable in enumerate(variables):

        ax1 = fig.add_subplot(2,len(variables),1+k)
        ax2 = fig.add_subplot(2,len(variables),1+len(variables)+k,projection=ccrs.PlateCarree())

        # Set details of axis
        ax1.set_title(variable)
        ax2.set_title(variable)
        ax1.set_xlim([-180,180])
        ax1.set_ylim([0,85])
        if k == 0:
            ax1.set_ylabel('Level')
        ax2.set_xlim([-180,180])
        ax2.set_ylim([-90,90])
        ax2.set_ylabel('Latitude')
        ax2.coastlines()
        ax1.set_xlabel('Longitude')
        ax2.set_xlabel('Longitude')

        plot_data = data[data[f'{variable}_maxabs_longitude'].notnull()]

        # Only make plots if there is more than one point
        if len(plot_data) > 0:
            initial_time = plot_data['Timestep'].values[0]

            # Initial point in longitude-height slice
            fig.axes[k*2].plot(plot_data[f'{variable}_maxabs_longitude'].values[0],
                               plot_data[f'{variable}_maxabs_level'].values[0],
                               color='black', linestyle='', marker='o')

            # Trajectory in longitude-height slice
            fig.axes[k*2].plot(plot_data[f'{variable}_maxabs_longitude'].values,
                               plot_data[f'{variable}_maxabs_level'].values,
                               color='black', linestyle='-')

            # Initial point in longitude-latitude slice
            fig.axes[k*2+1].plot(plot_data[f'{variable}_maxabs_longitude'].values[0],
                                 plot_data[f'{variable}_maxabs_latitude'].values[0],
                                 color='black', linestyle='', marker='o')

            # Trajectory in longitude-latitude slice
            fig.axes[k*2+1].plot(plot_data[f'{variable}_maxabs_longitude'].values,
                                 plot_data[f'{variable}_maxabs_latitude'].values,
                                 color='black', linestyle='-')

    fig.subplots_adjust(hspace=0.25)
    figname = f'{plotpath}/unphysical_{branch}.png'
    print(f'Plotting to {figname}')
    plt.savefig(figname, bbox_inches='tight')


if __name__ == "__main__":

    try:
        datapath, plotpath, branch = sys.argv[1:4]
    except ValueError:
        print("Usage: {0} <datapath> <plotpath> <branch>".format(sys.argv[0]))
        exit(1)

    make_figure(datapath, plotpath, branch)
