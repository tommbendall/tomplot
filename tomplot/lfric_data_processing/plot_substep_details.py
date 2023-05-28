import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
font = {'size':28}
plt.rc('font', **font)

legend_ncol = 2
legend_bbox = (0.5, 3.85)
results_opts = ['trunk_timing_37721', 'trunk_timing_37424', 'split_timing_37740', 'split_timing_37740_1pt5']
results_labels = ['r37721', 'r37424', 'split 1pt0', 'split 1pt5']
splittings = [False, False, 'vhv', 'vhv']
old_branches = [False, False, False, False]
variables = ['substeps', 'hori_cfl', 'vert_cfl']
colours = ['red', 'black', 'blue', 'purple']
ylabels = {'substeps': 'Number of substeps', 'hori_cfl': 'Horizontal CFL', 'vert_cfl': 'Vertical CFL'}
indices = {'substeps': 11, 'hori_cfl': 6, 'vert_cfl': 6}

for run in ['nwp', 'clim']:
    fig, axarray = plt.subplots(3, 1, sharex='col', sharey='row', figsize=(16, 15))

    for i, (ax, variable) in enumerate(zip(axarray, variables)):
        results_dirs = [f'/data/users/tbendall/results/trunk_timers/{run}_{results_opt}_{variable}.log' for results_opt in results_opts]
        for results_name, results_opt, colour, label, splitting, old_branch in \
                zip(results_dirs, results_opts, colours, results_labels, splittings, old_branches):
            raw_data = pd.read_csv(results_name, header=None, sep=' ', skipinitialspace=True)

            if variable == 'substeps':
                if splitting == 'vhv':
                    vertical_steps = raw_data[indices[variable]].values[::3] + raw_data[indices[variable]].values[2::3]
                    horizontal_steps = raw_data[indices[variable]].values[1::3]
                    total_steps = vertical_steps + horizontal_steps
                    timesteps = np.array(range(len(total_steps))) / 2.0
                    ax.plot(timesteps[:], total_steps, label=label, color=colour)
                    ax.plot(timesteps[:], vertical_steps, linestyle='--', color=colour)
                    ax.plot(timesteps[:], horizontal_steps, linestyle=':', color=colour)
                elif splitting == False:
                    horizontal_steps = raw_data[indices[variable]].values[:]
                    vertical_steps = raw_data[indices[variable]].values[:]
                    total_steps = vertical_steps + horizontal_steps
                    timesteps = np.array(range(len(total_steps))) / 2.0
                    ax.plot(timesteps[:], raw_data[indices[variable]].values[:], label=label, color=colour)
                else:
                    raise NotImplementedError
                print(f'{results_opt} {run} job: {np.sum(horizontal_steps)} horizontal substeps, {np.sum(vertical_steps)} vertical substeps')
            else:
                timesteps = np.array(range(len(raw_data[indices[variable]].values))) / 2.0
                if old_branch:
                    timesteps = timesteps / 2.0
                    step = 2
                else:
                    step = 1
                ax.plot(timesteps[::step], raw_data[indices[variable]].values[::step], label=label, color=colour)


        ax.set_ylabel(ylabels[variable])
        if i == len(axarray) - 1:
            ax.set_xlabel('Time step')

    handles, lgd_labels = axarray[-1].get_legend_handles_labels()
    lgd = axarray[-1].legend(handles, lgd_labels, loc='upper center', ncol=legend_ncol,
                        bbox_to_anchor=legend_bbox, edgecolor='black')

    figname = f'/data/users/tbendall/results/trunk_timers/figures/{run}_transport_steps.png'
    print(f'Saving plot to {figname}')
    fig.savefig(figname)