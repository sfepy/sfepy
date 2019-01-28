from __future__ import absolute_import
input_name = '../examples/navier_stokes/navier_stokes2d_iga.py'
output_name = 'test_navier_stokes2d_iga.vtk'

try:
    from igakit import igalib
except ImportError:
    from tests_basic import TestDummy
    class Test(TestDummy):
        pass
else:
    from tests_basic import TestInput
    class Test(TestInput):
        pass
