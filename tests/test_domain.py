import os.path as op

from sfepy.base.testing import TestCommon

class Test(TestCommon):

    @staticmethod
    def from_conf(conf, options):
        from sfepy import data_dir
        from sfepy.fem import Mesh, Domain

        mesh = Mesh('mesh', data_dir + '/meshes/various_formats/small3d.mesh')
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

    def test_refine(self):
        domain = self.domain

        for ii in range(3):
            domain = domain.refine()
            filename = op.join(self.options.out_dir,
                               'refine_' + domain.mesh.name + '.mesh')
            domain.mesh.write(filename, io='auto')

        return True
