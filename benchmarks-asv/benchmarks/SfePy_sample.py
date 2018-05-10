from __future__ import absolute_import

# Used for 'asv dev' only (benchmark development).
import sys
sys.path.append('/Users/kejzlar/Work/SfePy/SfePy')

from sfepy.mesh.mesh_generators import gen_block_mesh
from sfepy.discrete.fem.domain import FEDomain


class SfePySample:
    """
    Sample SfePy time/memory benchmarking.
    """
    version = None
    params = [10, 50, 100]
    param_names = ['Shape dim']
    timeout = 120

    def setup(self, dim):
        self.shape = [dim] * 3

    def teardown(self, dim):
        pass

    def domain_creation(self, dim):
        dims = [1, 1, 1][:dim]
        centre = [0, 0, 0][:dim]

        mesh = gen_block_mesh(dims, self.shape, centre, name='block')
        domain = FEDomain('domain', mesh)

        return domain

    def mem_domain_creation(self, dim):
        return SfePySample.domain_creation(self, dim)

    def peakmem_domain_creation(self, dim):
        return SfePySample.domain_creation(self, dim)

    def time_domain_creation(self, dim):
        SfePySample.domain_creation(self, dim)
