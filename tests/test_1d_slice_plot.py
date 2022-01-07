"""
In this file we test the routines for generating 1D slice plots.
"""

from tomplot import *
import os
import numpy as np
import pytest

def test_1d_slice_plot():

    have_i_passed = False

    #--------------------------------------------------------------------------#
    # Get NetCDF file                                                          #
    #--------------------------------------------------------------------------#

    file_name = 'tests/data/nc_fields/field_output_0.nc'
    plotdir = 'tests/figures'

    if not os.path.exists(file_name):
        raise IOError('Unable to find the model file for this run')

    #--------------------------------------------------------------------------#
    # Plot fields
    #--------------------------------------------------------------------------#

    fields = ['v_zonal']
    time_idxs = 'all'
    slices = [['x','midpoint'], ['y',25]]
    run_id = 0
    field_labels_from = ['space', 'scheme']

    make_field_plots(file_name, run_id, 'test',
                     fields, time_idxs, slices=slices,
                     plotdir=plotdir, override_dirname=True)


    have_i_passed = True

    assert have_i_passed
