"""
In this file we test the function for making a 2D field plot with a removed contour
"""

from firedrake import Function, SpatialCoordinate, sin, pi
from tomplot import *
from transportdrake import Outputting, build_mesh, build_spaces
import os
import numpy as np
import pytest

def test_remove_contour():

    have_i_passed = False

    domain = 'plane'
    n = 20
    L = 10
    data_name_string = 'test_remove_contour'

    #------------------------------------------------------------------------------#
    # Make NetCDF file                                                             #
    #------------------------------------------------------------------------------#
    # First remove the results directory if it already exists
    if os.path.exists('results/'+data_name_string):
        os.system('rm -r '+'results/'+data_name_string)

    outputting = Outputting(data_name_string, field_method='nc')

    #--------------------------------------------------------------------------#
    # Set-up of meshes, function spaces and fields                             #
    #--------------------------------------------------------------------------#
    mesh_parameters = {'Lx':L, 'Ly':L, 'nx': n, 'ny':n, 'cell_shape':'triangle'}
    mesh, mesh_metadata = build_mesh(domain, mesh_parameters)
    x, y = SpatialCoordinate(mesh)

    # make function spaces
    spaces_dict = build_spaces(mesh, space_names=[])

    V = spaces_dict['DG0']
    f = Function(V, name='f')

    #--------------------------------------------------------------------------#
    # Set-up initial conditions                                                #
    #--------------------------------------------------------------------------#

    # Something equally positive and negative so we have a clear contour to remove
    f.interpolate(sin(2*pi*x/L)*sin(4*pi*y/L))

    #--------------------------------------------------------------------------#
    # Set-up things for dumping output and dump first time step                #
    #--------------------------------------------------------------------------#

    run_metadata = {'mesh': mesh_metadata}
    run_id = outputting.register_run(run_metadata, spaces_dict, diagnostic_spaces=['DG0'])

    outputting.register_prognostic_field(name='f', field=f)

    outputting.dump_output_fields(0)

    #--------------------------------------------------------------------------#
    # Plot field                                                               #
    #--------------------------------------------------------------------------#

    make_field_plots(data_name_string, run_id, 'test_unremoved_contour',
                     'f', 'all', remove_contour=False,
                     plotdir='tests/figures')


    have_i_passed = True

    assert have_i_passed
