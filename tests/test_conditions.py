import os.path as op

import numpy as nm

from sfepy import data_dir

from sfepy.base.testing import TestCommon

def init_vec(variables):
    return nm.random.rand(variables.di.ptr[-1])

class Test(TestCommon):

    @staticmethod
    def from_conf(conf, options):
        from sfepy.discrete import FieldVariable, Variables, Problem
        from sfepy.discrete.fem import Mesh, FEDomain, Field

        mesh = Mesh.from_file(data_dir + '/meshes/2d/square_unit_tri.mesh')
        domain = FEDomain('domain', mesh)

        omega = domain.create_region('Omega', 'all')
        domain.create_region('Left',
                             'vertices in (x < -0.499)',
                             'facet')
        domain.create_region('LeftStrip',
                             'vertices in (x < -0.499)'
                             ' & (y > -0.199) & (y < 0.199)',
                             'facet')
        domain.create_region('LeftFix',
                             'r.Left -v r.LeftStrip',
                             'facet')
        domain.create_region('Right',
                             'vertices in (x > 0.499)',
                             'facet')
        domain.create_region('RightStrip',
                             'vertices in (x > 0.499)'
                             ' & (y > -0.199) & (y < 0.199)',
                             'facet')
        domain.create_region('RightFix',
                             'r.Right -v r.RightStrip',
                             'facet')

        fu = Field.from_args('fu', nm.float64, 'vector', omega, approx_order=2)
        u = FieldVariable('u', 'unknown', fu)

        fp = Field.from_args('fp', nm.float64, 'scalar', omega, approx_order=2)
        p = FieldVariable('p', 'unknown', fp)

        pb = Problem('test', domain=domain, fields=[fu, fp],
                     auto_conf=False, auto_solvers=False)

        test = Test(problem=pb, variables=Variables([u, p]),
                    conf=conf, options=options)
        return test

    def test_ebcs(self):
        from sfepy.discrete.conditions import Conditions, EssentialBC
        from sfepy.discrete.common.dof_info import expand_nodes_to_equations

        variables = self.variables
        regions = self.problem.domain.regions

        all_ebcs = []
        all_ebcs.append(EssentialBC('fix_u1', regions['LeftFix'],
                                    {'u.all' : nm.array([0.0, 1.0])}))
        all_ebcs.append(EssentialBC('fix_u2', regions['LeftStrip'],
                                    {'u.0' : 0.0, 'u.1' : 1.0}))
        all_ebcs.append(EssentialBC('fix_p1', regions['RightFix'],
                                    {'p.all' : 0.0}))
        all_ebcs.append(EssentialBC('fix_p2', regions['RightStrip'],
                                    {'p.0' : 0.0}))
        all_ebcs.append([EssentialBC('fix_p3', regions['Right'],
                                     {'p.0' : 0.0}),
                         EssentialBC('fix_u3', regions['Left'],
                                     {'u.0' : 0.0, 'u.1' : 1.0})])

        ok = True
        for ii, bcs in enumerate(all_ebcs):
            if not isinstance(bcs, list): bcs = [bcs]

            ebcs = Conditions(bcs)
            variables.equation_mapping(ebcs=ebcs, epbcs=None,
                                       ts=None, functions=None)
            vec = init_vec(variables)
            variables.apply_ebc(vec)

            for var_name, var_bcs in ebcs.group_by_variables().iteritems():
                var = variables[var_name]
                for bc in var_bcs:
                    bc.canonize_dof_names(var.dofs)
                    self.report('%d: %s %s: %s %s'
                                % (ii, var.name,
                                   bc.name, bc.region.name, bc.dofs[0]))
                    nods = var.field.get_dofs_in_region(bc.region)
                    eq = expand_nodes_to_equations(nods, bc.dofs[0], var.dofs)

                    off = variables.di.indx[var_name].start
                    n_nod = len(nods)
                    for bdof, dof_name in enumerate(bc.dofs[0]):
                        idof = var.dofs.index(dof_name)
                        eqs = eq[n_nod * bdof : n_nod * (bdof + 1)]

                        _ok = nm.allclose(vec[off + eqs], idof,
                                          atol=1e-14, rtol=0.0)
                        if not _ok:
                            self.report(' %s: failed! (all of %s == %f)'
                                        % (dof_name, vec[off + eqs], idof))
                        ok = ok and _ok

        return ok

    def test_epbcs(self):
        from sfepy.discrete import Function, Functions
        from sfepy.discrete.conditions import Conditions, PeriodicBC
        from sfepy.discrete.common.dof_info import expand_nodes_to_equations
        from sfepy.discrete.fem.periodic import match_y_line

        variables = self.variables
        regions = self.problem.domain.regions

        match_y_line = Function('match_y_line', match_y_line)
        pbc = PeriodicBC('pbc', [regions['LeftStrip'], regions['RightStrip']],
                         {'u.[1,0]' : 'u.[0,1]'}, match='match_y_line')

        functions = Functions([match_y_line])

        epbcs = Conditions([pbc])
        variables.equation_mapping(ebcs=None, epbcs=epbcs,
                                   ts=None, functions=functions)

        vec = init_vec(variables)
        variables.apply_ebc(vec)

        var = variables['u']
        var_bcs = epbcs.group_by_variables()['u']
        bc = var_bcs['pbc']
        bc.canonize_dof_names(var.dofs)

        nods0 = var.field.get_dofs_in_region(bc.regions[0])
        nods1 = var.field.get_dofs_in_region(bc.regions[1])

        coors0 = var.field.get_coor(nods0)
        coors1 = var.field.get_coor(nods1)
        i0, i1 = match_y_line(coors0, coors1)

        eq0 = expand_nodes_to_equations(nods0[i0], bc.dofs[0], var.dofs)
        eq1 = expand_nodes_to_equations(nods1[i1], bc.dofs[1], var.dofs)

        ok = True

        _ok = len(nm.setdiff1d(eq0, var.eq_map.master)) == 0
        if not _ok:
            self.report('master equations mismatch! (set(%s) == set(%s))'
                        % (eq0, var.eq_map.master))
        ok = ok and _ok

        _ok = len(nm.setdiff1d(eq1, var.eq_map.slave)) == 0
        if not _ok:
            self.report('slave equations mismatch! (set(%s) == set(%s))'
                        % (eq1, var.eq_map.slave))
        ok = ok and _ok

        off = variables.di.indx['u'].start
        _ok = nm.allclose(vec[off + eq0], vec[off + eq1], atol=1e-14, rtol=0.0)
        if not _ok:
            self.report('periodicity test failed! (%s == %s)'
                        % (vec[off + eq0], vec[off + eq0]))
        ok = ok and _ok

        return ok

    def test_save_ebc(self):
        name = op.join(self.options.out_dir,
                       op.splitext(op.basename(__file__))[0])
        self.problem.save_ebc(name + '_ebcs_f.vtk', force=True)
        self.problem.save_ebc(name + '_ebcs.vtk', default=-1, force=False)

        return True
