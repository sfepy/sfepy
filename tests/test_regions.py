from __future__ import absolute_import
import numpy as nm

from sfepy.base.conf import transform_functions
from sfepy.base.testing import TestCommon

def get_vertices(coors, domain=None):
    x, z = coors[:,0], coors[:,2]

    return nm.where((z < 0.1) & (x < 0.1))[0]

def get_cells(coors, domain=None):
    return nm.where(coors[:, 0] < 0)[0]

class Test(TestCommon):

    @staticmethod
    def from_conf( conf, options ):
        from sfepy import data_dir
        from sfepy.discrete.fem import Mesh, FEDomain
        from sfepy.discrete import Functions

        mesh = Mesh.from_file(data_dir
                              + '/meshes/various_formats/abaqus_tet.inp')
        mesh.nodal_bcs['set0'] = [0, 7]
        domain = FEDomain('test domain', mesh)

        conf_functions = {
            'get_vertices' : (get_vertices,),
            'get_cells' : (get_cells,),
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
            ['all', 'cell'],
            ['vertices of surface', 'facet'],
            ['vertices of group %d' % self.domain.cmesh.vertex_groups[0],
             'facet'],
            ['vertices of set set0', 'vertex'],
            ['vertices in (z < 0.1) & (x < 0.1)', 'facet'],
            ['vertices by get_vertices', 'cell'],
            ['vertex 0, 1, 2', 'vertex'],
            ['vertex in r.r6', 'vertex'],
            ['cells of group %d' % self.domain.cmesh.cell_groups[0], 'cell'],
            # ['cells of set 0', 'cell'], not implemented...
            ['cells by get_cells', 'cell'],
            ['cell 1, 4, 5', 'cell'],
            ['copy r.r5', 'cell'],
            ['r.r5', 'cell'],
        ]

        vertices = [
            [0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12],
            [0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12],
            [0,  1,  3,  7],
            [0,  7],
            [1,  2,  3,  4,  5,  9, 11],
            [1,  2,  3,  4,  5,  9, 11],
            [0,  1,  2],
            [0],
            [0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10],
            [0,  1,  2,  3,  4,  5,  6,  9, 10, 11],
            [0,  1,  2,  3,  4,  5,  6,  8],
            [1,  2,  3,  4,  5,  9, 11],
            [1,  2,  3,  4,  5,  9, 11],
        ]

        ok = True
        for ii, sel in enumerate(selectors):
            self.report('select:', sel)
            reg = self.domain.create_region('r%d' % ii, sel[0], kind=sel[1],
                                            functions=self.functions)
            _ok = ((len(reg.vertices) == len(vertices[ii]))
                   and (reg.vertices == vertices[ii]).all())
            self.report('  vertices:', _ok)

            ok = ok and _ok

        return ok

    def test_operators(self):
        """
        Test operators in region selectors.
        """
        ok = True

        r1 = self.domain.create_region('r1', 'all')

        sel = 'r.r1 -v vertices of group 0'
        self.report('select:', sel)
        reg = self.domain.create_region('reg', sel, kind='vertex')
        av = [2,  4,  5,  6,  8,  9, 10, 11, 12]
        _ok = (reg.vertices == nm.array(av)).all()
        self.report('  vertices:', _ok)
        ok = ok and _ok

        sel = 'vertex 0, 1, 2 +v vertices of group 0'
        self.report('select:', sel)
        reg = self.domain.create_region('reg', sel, kind='vertex')
        av = [0,  1,  2,  3,  7]
        _ok = (reg.vertices == nm.array(av)).all()
        self.report('  vertices:', _ok)
        ok = ok and _ok

        sel = 'vertex 0, 1, 2 *v vertices of group 0'
        self.report('select:', sel)
        reg = self.domain.create_region('reg', sel, kind='vertex')
        av = [0,  1]
        _ok = (reg.vertices == nm.array(av)).all()
        self.report('  vertices:', _ok)
        ok = ok and _ok

        sel = 'vertex 20 *v vertices of group 0'
        self.report('select:', sel, ', allow_empty == False')
        _ok = False
        try:
            reg = self.domain.create_region('reg', sel, kind='vertex',
                                            allow_empty=False)
        except ValueError:
            _ok = True
        self.report('  exception raised:', _ok)
        ok = ok and _ok

        sel = 'vertex 20 *v vertices of group 0'
        self.report('select:', sel, ', allow_empty == True')
        reg = self.domain.create_region('reg', sel, kind='vertex',
                                        allow_empty=True)
        av = []
        _ok = (reg.vertices == nm.array(av)).all()
        self.report('  vertices:', _ok)
        ok = ok and _ok

        sel = 'r.r1 -c cell 1, 4, 5'
        self.report('select:', sel)
        reg = self.domain.create_region('reg', sel)
        _ok = (nm.setdiff1d(r1.cells[0], [1, 4, 5]) == reg.cells[0]).all()
        self.report('  cells:', _ok)
        ok = ok and _ok

        sel = 'cell 8, 3 +c cell 1, 4, 5'
        self.report('select:', sel)
        reg = self.domain.create_region('reg', sel)
        cells = [1,  3,  4,  5,  8]
        _ok = (reg.cells == nm.array(cells)).all()
        self.report('  cells:', _ok)
        ok = ok and _ok

        sel = 'cell 8, 3, 2 *c cell 8, 4, 2, 7'
        self.report('select:', sel)
        reg = self.domain.create_region('reg', sel)
        cells = [2,  8]
        _ok = (reg.cells == nm.array(cells)).all()
        self.report('  cells:', _ok)
        ok = ok and _ok

        return ok
