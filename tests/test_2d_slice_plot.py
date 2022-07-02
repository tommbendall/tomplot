"""
In this file we test the routines for generating 2D slice plots.
"""

from tomplot import *
import os
import numpy as np

def test_2d_slice_plot():

    have_i_passed = False

    #--------------------------------------------------------------------------#
    # Get NetCDF file                                                          #
    #--------------------------------------------------------------------------#

    file_name = 'tests/data/tomdata_sphere/field_output_0.nc'
    plotdir = 'tests/figures'

    if not os.path.exists(file_name):
        raise IOError('Unable to find the model file for this run')

    #--------------------------------------------------------------------------#
    # Plot fields
    #--------------------------------------------------------------------------#

    fields = ['v_zonal']
    time_idxs = 'all'
    slices = ['xy']
    run_id = 0

    make_field_plots(file_name, run_id, 'test',
                     fields, time_idxs, slices=slices,
                     plotdir=plotdir, override_dirname=True)


    have_i_passed = True

    assert have_i_passed
