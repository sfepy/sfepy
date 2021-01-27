from __future__ import absolute_import
input_name = '../examples/diffusion/poisson_iga.py'
output_name = 'test_poisson_iga.vtk'

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
