"""
In this file we test the routines for generating convergence plots.
"""

from tomplot import *
import os
import numpy as np
import pytest

def test_convergence_plot():

    have_i_passed = False

    # Things the same for all plots
    variables = 'dx'

    #--------------------------------------------------------------------------#
    # Get NetCDF file                                                          #
    #--------------------------------------------------------------------------#

    file_name = 'tests/data/global_output.nc'
    plotdir = 'tests/figures'

    if not os.path.exists(file_name):
        raise IOError('Unable to find the model file for this run')

    #--------------------------------------------------------------------------#
    # Plot dissipation errors                                                  #
    #--------------------------------------------------------------------------#

    errors = 'dissipation_error'
    fields = ['v', 'w', 'v', 'w']
    field_labels = [r'RTCF1\_vector\_invariant', r'RTCF1\_dg\_advection',
                    r'Vec\_DG1\_vector\_invariant', r'Vec\_DG1\_dg\_advection']
    run_ids = [[0, 1], [0, 1], [2, 3], [2, 3]]

    make_convergence_plots(file_name, variables, fields, run_ids, errors,
                           testname='test', plotdir=plotdir, override_dirname=True)

    #--------------------------------------------------------------------------#
    # Plot two types of error for v                                            #
    #--------------------------------------------------------------------------#

    errors = ['dispersion_error', 'L2_error']
    fields = ['v', 'v']
    run_ids = [[0, 1], [2, 3]]
    comparison_lines = [1, 1.5]

    make_convergence_plots(file_name, variables, fields, run_ids, errors,
                           testname='test', plotdir=plotdir, override_dirname=True,
                           comparison_lines=comparison_lines)

    have_i_passed = True

    assert have_i_passed
