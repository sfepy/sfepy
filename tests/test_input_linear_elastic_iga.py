from __future__ import absolute_import
input_name = '../examples/linear_elasticity/linear_elastic_iga.py'
output_name = 'test_linear_elastic_iga.vtk'

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
