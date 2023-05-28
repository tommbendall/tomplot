import numpy as np
import matplotlib.pyplot as plt
import sys

def make_figure(datapath, plotpath, branch, months):

    fig, ((ax, bx), (cx, dx), (ex, fx)) = plt.subplots(3,2,figsize=(15,6))
    plt.suptitle("Configuration = "+branch)

    # Load up CFL-h data
    fname = f'{datapath}/CFL_H.txt'
    cfl_h = np.loadtxt(fname)
    nper_time = 4
    t = np.arange(len(cfl_h))/nper_time

    ax.plot(t, cfl_h, color='black', linestyle='-')
    ax.set(ylabel='CFL H')
    ax.set(ylim=[0,15])

    # Load up CFL-v data
    fname = f'{datapath}/CFL_V.txt'
    cfl_v = np.loadtxt(fname)
    nper_time = 4
    t = np.arange(len(cfl_v))/nper_time

    bx.plot(t, cfl_v, color='black', linestyle='-')
    bx.set(ylabel='CFL V')
    bx.set(ylim=[0,50])

    # Load up dtheta-slow data
    fname = f'{datapath}/dtheta_slow.txt'
    dt_s = np.loadtxt(fname)
    nper_time = 1
    t = np.arange(len(dt_s))/nper_time

    cx.plot(t, dt_s[:,0], color='black', linestyle='-')
    cx.set(ylabel='dtheta slow (Min)')
    cx.set(ylim=[-8,0])
    dx.plot(t, dt_s[:,1], color='black', linestyle='-')
    dx.set(ylabel='dtheta slow (Max)')
    dx.set(ylim=[0,30])

    # Load up dtheta-fast data
    fname = f'{datapath}/dtheta_fast.txt'
    try:
        dt_f = np.loadtxt(fname)
        nper_time = 2
        t = np.arange(len(dt_f))/nper_time

        ex.plot(t, dt_f[:,0], color='black', linestyle='-')
        ex.set(ylabel='dtheta fast (Min)')
        ex.set(ylim=[-25,0])
        fx.plot(t, dt_f[:,1], color='black', linestyle='-')
        fx.set(ylabel='dtheta fast (Max)')
        fx.set(ylim=[0,12])

    except ValueError:
        print(f'Unable to read branch {branch}')

    figname = f'{plotpath}/sensitivity_{branch}.png'
    print(f'Plotting to {figname}')
    plt.savefig(figname)

if __name__ == "__main__":

    try:
        datapath, plotpath, branch, months = sys.argv[1:5]
    except ValueError:
        print("Usage: {0} <datapath> <plotpath> <branch> <months>".format(sys.argv[0]))
        exit(1)

    month_list = months.split(':')
    make_figure(datapath, plotpath, branch, month_list)
