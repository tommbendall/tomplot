"""
Makes plots for the cylindrical convergence test case
"""

from netCDF4 import Dataset
import numpy as np
from tomplot import make_quiver_plots, make_convergence_plots

# ---------------------------------------------------------------------------- #
# Things that can be altered and parameters for the test case
# ---------------------------------------------------------------------------- #

results_dirname = 'demo_5_cylindrical_curly'
plot_times = 'all'
field_labels = ['RTCF2']
plot_testname = 'cylinder_curly'
cbar_label = r'$|F| \ / $ m s$^{-1}$'

# ---------------------------------------------------------------------------- #
# Derived things from options
# ---------------------------------------------------------------------------- #

data = Dataset('results/'+results_dirname+'/global_output.nc','r')
run_ids = data['run_id'][:]
num_setups = len(field_labels)
data.close()

# ---------------------------------------------------------------------------- #
# Field plots
# ---------------------------------------------------------------------------- #

print('Making quiver plots')
for run_id in run_ids:
    for i in range(num_setups):
        print('Field plot %d %d' % (run_id, i))
        field_info = ('F_'+str(i), 'F_'+str(i)+'_zonal', 'F_'+str(i)+'_meridional')

        make_quiver_plots(results_dirname, run_id, plot_testname,
                          field_info, plot_times, 'xy', colour_scheme='OrRd',
                          contours=np.arange(0.0, 1.4, 0.1), contour_method='magnitude',
                          extend_cmap=False, quiver_npts=2, scale=0.5, cbar_label=cbar_label,
                          ylabelpad=-30, xlabel=r'$\varrho \phi \ / $ m')
