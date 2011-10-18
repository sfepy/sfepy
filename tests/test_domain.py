import os.path as op

from sfepy.base.testing import TestCommon
from sfepy import data_dir
from sfepy.fem import Mesh, Domain

def refine(domain, out_dir, level=3):
    for ii in range(3):
        domain = domain.refine()
        filename = op.join(out_dir, 'refine_' + domain.mesh.name + '.mesh')
        domain.mesh.write(filename, io='auto')

    return domain

class Test(TestCommon):

    @staticmethod
    def from_conf(conf, options):
        mesh = Mesh('mesh_tetra',
                    data_dir + '/meshes/various_formats/small3d.mesh')
        domain = Domain('domain', mesh)

        return Test(conf=conf, options=options, domain=domain)

    def test_facets(self):
        ok = True

        ed = self.domain.ed
        _ok = ed.n_unique == 26
        self.report('unique edges: %s' % _ok)
        ok = ok and _ok

        fa = self.domain.fa
        _ok = fa.n_unique == 30
        self.report('unique faces: %s' % _ok)
        ok = ok and _ok

        return ok

    def test_refine_tetra(self):
        refine(self.domain, self.options.out_dir)

        return True

    def test_refine_hexa(self):
        mesh = Mesh('mesh_hexa',
                    data_dir + '/meshes/various_formats/abaqus_hex.inp')
        domain = Domain('domain', mesh)

        refine(domain, self.options.out_dir)

        return True
