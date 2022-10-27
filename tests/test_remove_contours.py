"""
In this file we test the function for making a 2D field plot with a removed contour
"""

from firedrake import Function, SpatialCoordinate, sin, pi
from tomplot import *
from transportdrake import Outputting, build_mesh, build_spaces
import os
import numpy as np

def test_remove_contour():

    have_i_passed = False

    domain = 'plane'
    n = 20
    L = 10
    file_name = 'tests/data/test_remove_contour/field_output_0.nc'

    if not os.path.exists(file_name):
        raise IOError('Unable to find the model file for this run')

    #--------------------------------------------------------------------------#
    # Plot field                                                               #
    #--------------------------------------------------------------------------#

    run_id = 0
    make_field_plots(file_name, run_id, 'test_unremoved_contour',
                     'f', 'all', remove_contour=False,
                     plotdir='tests/figures', override_dirname=True)


    have_i_passed = True

    assert have_i_passed
