"""
Makes quiver plots for the planar demo case
"""

import numpy as np
from tomplot import make_quiver_plots

# ---------------------------------------------------------------------------- #
# Things that can be altered and parameters for the test case
# ---------------------------------------------------------------------------- #

plot_times = 'all'
run_id = 2
field_info  = ('F_0', 'F_0_x', 'F_0_y')

# ---------------------------------------------------------------------------- #
# Field plots
# ---------------------------------------------------------------------------- #

print('Making quiver plots')

results_dirname = 'demo_1_RTCF1_planar'

make_quiver_plots(results_dirname, run_id, 'planar_demo',
                  field_info, plot_times, 'xy',
                  contours=np.arange(0.0, 2.4, 0.2),
                  cbar_label=r'$|F| \ /$ m s$^{-1}$',
                  contour_method='magnitude', quiver_npts=2, scale=0.5)

results_dirname = 'demo_2_RTCF2_planar'
make_quiver_plots(results_dirname, run_id, 'planar_demo',
                  field_info, plot_times, 'xy',
                  contours=np.arange(0.0, 2.4, 0.2),
                  cbar_label=r'$|F| \ /$ m s$^{-1}$',
                  contour_method='magnitude', quiver_npts=1, scale=0.5)

