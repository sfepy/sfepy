import numpy as nm

from sfepy.base.conf import transform_functions
from sfepy.base.testing import TestCommon

def get_nodes(coors, domain=None):
    x, z = coors[:,0], coors[:,2]

    return nm.where((z < 0.1) & (x < 0.1))[0]

def get_elements(coors, domain=None):
    return {0 : [1, 4, 5]}

class Test(TestCommon):

    @staticmethod
    def from_conf( conf, options ):
        from sfepy import data_dir
        from sfepy.fem import Mesh, Domain, Functions

        mesh = Mesh('test mesh',
                    data_dir + '/meshes/various_formats/abaqus_tet.inp')
        domain = Domain('test domain', mesh)

        conf_functions = {
            'get_nodes' : (get_nodes,),
            'get_elements' : (get_elements,),
        }
        functions = Functions.from_conf(transform_functions(conf_functions))

        test = Test(conf=conf, options=options,
                    domain=domain, functions=functions)
        return test

    def test_selectors(self):
        """
        Test basic region selectors.
        """
        selectors = [
            'all',
            'nodes of surface',
            'nodes of group 0',
            'nodes in (z < 0.1) & (x < 0.1)',
            'nodes by get_nodes',
            'node 0, 1, 2',
            'elements of group 0',
            'elements by get_elements',
            'element 1, 4, 5',
            'element (0, 1), (0, 4), (0, 5)'
        ]

        all_vertices = [
            [0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12],
            [0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12],
            [0,  1,  3,  7],
            [1,  2,  3,  4,  5,  9, 11],
            [1,  2,  3,  4,  5,  9, 11],
            [0,  1,  2],
            [0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12],
            [0,  1,  2,  3,  4,  5,  6,  8],
            [0,  1,  2,  3,  4,  5,  6,  8],
            [0,  1,  2,  3,  4,  5,  6,  8],
        ]

        ok = True
        for ii, sel in enumerate(selectors):
            self.report('select:', sel)
            reg = self.domain.create_region('r', sel, functions=self.functions)

            _ok = (reg.all_vertices == all_vertices[ii]).all()
            self.report('  all_vertices:', _ok)

            ok = ok and _ok

        return ok

    def test_operators(self):
        """
        Test operators in region selectors.
        """
        ok = True

        r1 = self.domain.create_region('r1', 'nodes of surface')

        sel = 'r.r1 -n nodes of group 0'
        self.report('select:', sel)
        reg = self.domain.create_region('reg', sel)
        av = [2,  4,  5,  6,  8,  9, 10, 11, 12]
        _ok = (reg.all_vertices == nm.array(av)).all()
        self.report('  all_vertices:', _ok)
        ok = ok and _ok

        sel = 'node 0, 1, 2 +n nodes of group 0'
        self.report('select:', sel)
        reg = self.domain.create_region('reg', sel)
        av = [0,  1,  2,  3,  7]
        _ok = (reg.all_vertices == nm.array(av)).all()
        self.report('  all_vertices:', _ok)
        ok = ok and _ok

        sel = 'node 0, 1, 2 *n nodes of group 0'
        self.report('select:', sel)
        reg = self.domain.create_region('reg', sel)
        av = [0,  1]
        _ok = (reg.all_vertices == nm.array(av)).all()
        self.report('  all_vertices:', _ok)
        ok = ok and _ok

        sel = 'r.r1 -e element 1, 4, 5'
        self.report('select:', sel)
        reg = self.domain.create_region('reg', sel)
        _ok = (nm.setdiff1d(r1.cells[0], [1, 4, 5]) == reg.cells[0]).all()
        self.report('  cells:', _ok)
        ok = ok and _ok

        sel = 'element 8, 3 +e element 1, 4, 5'
        self.report('select:', sel)
        reg = self.domain.create_region('reg', sel)
        cells = [1,  3,  4,  5,  8]
        _ok = (reg.cells[0] == nm.array(cells)).all()
        self.report('  cells:', _ok)
        ok = ok and _ok

        sel = 'element 8, 3, 2 *e element 8, 4, 2, 7'
        self.report('select:', sel)
        reg = self.domain.create_region('reg', sel)
        cells = [2,  8]
        _ok = (reg.cells[0] == nm.array(cells)).all()
        self.report('  cells:', _ok)
        ok = ok and _ok

        return ok
