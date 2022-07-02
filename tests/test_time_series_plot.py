"""
In this file we test the routines for generating time series plots.
"""

from tomplot import *
import os
import numpy as np

def test_time_series_plot():

    have_i_passed = False

    # Things the same for all plots
    variables = 'dx'

    #--------------------------------------------------------------------------#
    # Get NetCDF file                                                          #
    #--------------------------------------------------------------------------#

    file_name = 'tests/data/tomdata_sphere/global_output.nc'
    plotdir = 'tests/figures'

    if not os.path.exists(file_name):
        raise IOError('Unable to find the model file for this run')

    #--------------------------------------------------------------------------#
    # Plot mins and maxes
    #--------------------------------------------------------------------------#

    diagnostics = ['min', 'max']
    fields = ['v_Vec_DG1', 'v_Vec_DG1', 'w', 'w']
    run_ids = [0, 2, 0, 2]

    make_time_series_plots(file_name, fields, run_ids, diagnostics,
                           testname='test', plotdir=plotdir, override_dirname=True)

    #--------------------------------------------------------------------------#
    # Plot L2
    #--------------------------------------------------------------------------#

    diagnostics = 'L2'
    fields = ['v', 'w', 'v', 'w']
    run_ids = [0, 0, 1, 1]

    make_time_series_plots(file_name, fields, run_ids, diagnostics,
                           testname='test', plotdir=plotdir, override_dirname=True)

    have_i_passed = True

    assert have_i_passed
