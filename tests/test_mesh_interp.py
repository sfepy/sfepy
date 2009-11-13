import os.path as op

from sfepy.base.base import *
from sfepy.base.conf import transform_variables, transform_fields
from sfepy.base.testing import TestCommon

fields = {
#    'f' : ((1,1), 'real', 'Omega', {'Omega' : '3_4_P1'}),
    'f' : ((1,1), 'real', 'Omega', {'Omega' : '3_8_Q1'}),
}

variables = {
    'u'       : ('unknown field', 'f', 0),
    'v'       : ('test field',    'f', 'u'),
}

def in_dir(adir):
    return lambda x: op.join(adir, x)

class Test(TestCommon):

    @staticmethod
    def from_conf(conf, options):
        test = Test(conf=conf, options=options)
        return test

    def test_interpolation(self):
        from sfepy.fem import Mesh, Domain, Fields, Variables
        from sfepy.base.la import make_axis_rotation_matrix

#        m1 = Mesh('original mesh', 'database/simple.mesh')
#        m1 = Mesh('original mesh', 'database/kostka_medium_tetra.mesh')
        m1 = Mesh('original mesh', 'database/tests/block.mesh')
#        m1 = Mesh('original mesh', '../sfepy-projects/screw_press/meshes/screw_press_3.mesh')

        for ii, angle in enumerate(nm.linspace(0.0, nm.pi, 11)):
            print '>>>>>>>>>>>>>>>>', ii, angle
            shift = nm.array([2.0, 0.5, 0.5]) * (ii / 10.0)
#            shift = [0.0, 0.0, 0.0]
            mtx = make_axis_rotation_matrix([0, 1, 0], angle)

            m2 = m1.copy('rotated mesh')
            m2.transform_coors(mtx)
            m2.coors += shift

            d1 = Domain('d1', m1)

            omega1 = d1.create_region('Omega', 'all')

            ff = Fields.from_conf(transform_fields(fields), d1.regions)
            field1 = ff[0]

            vv = Variables.from_conf(transform_variables(variables), ff)
            vv.setup_dof_info()
            u1 = vv['u']

            bbox = m1.get_bounding_box()
            nx = bbox[1,0] - bbox[0,0]
            data1 = nm.sin(4.0 * nm.pi * m1.coors[:,0:1] / nx)
            u1.set_from_mesh_vertices(data1)

            d2 = Domain('d2', m2)
            omega2 = d2.create_region('Omega', 'all')

            ff2 = Fields.from_conf(transform_fields(fields), d2.regions)
            field2 = ff2[0]

            vv2 = Variables.from_conf(transform_variables(variables), ff2)
            vv2.setup_dof_info()
            u2 = vv2['u']

            # Performs interpolation, if other field differs from self.field
            # or, in particular, is defined on a different mesh.
            u2.set_from_other(u1, strategy='interpolation')
##             u2.set_from_other(u1, strategy='crawl')

            fname = in_dir(self.options.out_dir)

            if ii == 0:
                u1.save_as_mesh(fname('test_mesh_interp_u1.vtk'))

            u2.save_as_mesh(fname('test_mesh_interp_u2.%03d.vtk' % ii))
        
        return True
