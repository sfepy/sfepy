from __future__ import print_function
from __future__ import absolute_import
from sfepy import data_dir

filename_mesh = data_dir + '/meshes/2d/special/circle_in_square.mesh'

dim = 2

fields = {
    'scalar_field' : ('real', 'scalar', 'Omega', 1),
    'vector_field' : ('real', 'vector', 'Omega', 1),
}

variables = {
    'us'  : ('unknown field',   'scalar_field', 0),
    'ts'  : ('test field',      'scalar_field', 'us'),
    'ps1' : ('parameter field', 'scalar_field', 'us'),
    'ps2' : ('parameter field', 'scalar_field', 'us'),
    'uv'  : ('unknown field',   'vector_field', 1),
    'tv'  : ('test field',      'vector_field', 'uv'),
    'pv1' : ('parameter field', 'vector_field', 'uv'),
    'pv2' : ('parameter field', 'vector_field', 'uv'),
}

regions = {
    'Omega' : 'all',
    'Left' : ('vertices in (x < -0.499)', 'facet'),
}

integrals = {
    'i' : 2,
}

materials = {
    'm' : 'get_pars',
    'm2' : ({'K' : [[3.0, 0.1], [0.3, 1.0]]},),
}

equations = {
    'eq' : """dw_diffusion.i.Omega( m2.K, ts, us ) = 0"""
}

def get_pars(ts, coor, mode=None, term=None, **kwargs):
    if mode == 'qp':
        n_nod, dim = coor.shape
        sym = (dim + 1) * dim // 2

        if 'biot' in term.name:
            val = nm.zeros((sym, 1), dtype=nm.float64)
            val[:dim] = 0.132
            val[dim:sym] = 0.092
        elif 'volume_dot' in term.name:
            val = 1.0 / nm.array([3.8], dtype=nm.float64)
        elif 'diffusion' in term.name:
            val = nm.eye(dim, dtype=nm.float64)
        else:
            raise ValueError

        return {'val' : nm.tile(val, (coor.shape[0], 1, 1))}

functions = {
    'get_pars' : (get_pars,),
}

# (eval term prefix, parameter corresponding to test variable, 'd' variables,
# 'dw' variables (test must be paired with unknown, which should be at
# index 2!), mat mode)
test_terms = [
    ('%s_biot.i.Omega( m.val, %s, %s )',
     ('dw', 'ps1', ('pv1', 'ps1'), ('pv1', 'ts', 'us', 'uv', 'tv'))),
    ('%s_biot.i.Omega( m.val, %s, %s )',
     ('dw', 'pv1', ('pv1', 'ps1'), ('tv', 'ps1', 'uv', 'us', 'ts'))),
    ('%s_diffusion.i.Omega( m.val, %s, %s )',
     ('dw', 'ps1', ('ps1', 'ps2'), ('ts', 'ps1', 'us'))),
    ('%s_volume_dot.i.Omega( m.val, %s, %s )',
     ('dw', 'ps1', ('ps1', 'ps2'), ('ts', 'ps1', 'us'))),
]

import numpy as nm
from sfepy.base.testing import TestCommon

def _integrate(var, val_qp):
    from sfepy.discrete import Integral
    from sfepy.discrete.common.mappings import get_jacobian

    integral = Integral('i', 2)
    det = get_jacobian(var.field, integral)
    val = (val_qp * det).sum(axis=1) / det.sum(axis=1)

    return val

class Test(TestCommon):

    @staticmethod
    def from_conf(conf, options):
        from sfepy.discrete import Problem

        problem = Problem.from_conf(conf, init_equations=False)
        test = Test(problem=problem,
                    conf=conf, options=options)
        return test

    def test_consistency_d_dw(self):
        from sfepy.discrete import Variables

        ok = True
        pb = self.problem
        for aux in test_terms:
            term_template, (prefix, par_name, d_vars, dw_vars) = aux
            print(term_template, prefix, par_name, d_vars, dw_vars)

            term1 = term_template % ((prefix,) + d_vars)

            variables = Variables.from_conf(self.conf.variables, pb.fields)

            for var_name in d_vars:
                var = variables[var_name]
                n_dof = var.field.n_nod * var.field.shape[0]
                aux = nm.arange(n_dof, dtype=nm.float64)
                var.set_data(aux)

            if prefix == 'd':
                val1 = pb.evaluate(term1, var_dict=variables.as_dict())

            else:
                val1 = pb.evaluate(term1, call_mode='d_eval',
                                   var_dict=variables.as_dict())

            self.report('%s: %s' % (term1, val1))

            term2 = term_template % (('dw',) + dw_vars[:2])

            vec, vv = pb.evaluate(term2, mode='weak',
                                  var_dict=variables.as_dict(),
                                  ret_variables=True)

            pvec = vv.get_state_part_view(vec, dw_vars[2])
            val2 = nm.dot(variables[par_name](), pvec)
            self.report('%s: %s' % (term2, val2))

            err = nm.abs(val1 - val2) / nm.abs(val1)
            _ok = err < 1e-12
            self.report('relative difference: %e -> %s' % (err, _ok))

            ok = ok and _ok

        return ok

    def test_eval_matrix(self):
        problem = self.problem

        problem.set_equations()
        problem.time_update(ebcs={}, epbcs={})

        var = problem.get_variables()['us']

        vec = nm.arange(var.n_dof, dtype=var.dtype)

        var.set_data(vec)

        val1 = problem.evaluate('dw_diffusion.i.Omega( m2.K, us, us )',
                                mode='eval', us=var)

        mtx = problem.evaluate('dw_diffusion.i.Omega( m2.K, ts, us )',
                               mode='weak', dw_mode='matrix')

        val2 = nm.dot(vec, mtx * vec)

        ok = (nm.abs(val1 - val2) / nm.abs(val1)) < 1e-14
        self.report('eval: %s, weak: %s, ok: %s' % (val1, val2, ok))

        return ok

    def test_vector_matrix(self):
        problem = self.problem

        problem.set_equations()
        problem.time_update()

        state = problem.create_state()
        state.apply_ebc()

        aux1 = problem.evaluate("dw_diffusion.i.Omega( m2.K, ts, us )",
                                mode='weak', dw_mode='vector')

        problem.time_update(ebcs={}, epbcs={})

        mtx = problem.evaluate("dw_diffusion.i.Omega( m2.K, ts, us )",
                               mode='weak', dw_mode='matrix')
        aux2g = mtx * state()
        problem.time_update(ebcs=self.conf.ebcs,
                            epbcs=self.conf.epbcs)
        aux2 = problem.equations.strip_state_vector(aux2g, follow_epbc=True)

        ret = self.compare_vectors(aux1, aux2,
                                   label1='vector mode',
                                   label2='matrix mode')
        if not ret:
            self.report('failed')

        return ret

    def test_surface_evaluate(self):
        from sfepy.discrete import FieldVariable
        problem = self.problem

        us = problem.get_variables()['us']
        vec = nm.empty(us.n_dof, dtype=us.dtype)
        vec[:] = 1.0
        us.set_data(vec)

        expr = 'ev_surface_integrate.i.Left( us )'
        val = problem.evaluate(expr, us=us)
        ok1 = nm.abs(val - 1.0) < 1e-15
        self.report('with unknown: %s, value: %s, ok: %s'
                    % (expr, val, ok1))

        ps1 = FieldVariable('ps1', 'parameter', us.get_field(),
                            primary_var_name='(set-to-None)')
        ps1.set_data(vec)

        expr = 'ev_surface_integrate.i.Left( ps1 )'
        val = problem.evaluate(expr, ps1=ps1)
        ok2 = nm.abs(val - 1.0) < 1e-15
        self.report('with parameter: %s, value: %s, ok: %s'
                    % (expr, val, ok2))
        ok2 = True

        return ok1 and ok2

    def test_ev_grad(self):
        problem = self.problem

        var = problem.create_variables(['us'])['us']
        val = nm.arange(var.n_dof, dtype=var.dtype)
        var.set_data(val)

        val1 = problem.evaluate('ev_grad.i.Omega( us )', us=var, mode='el_avg')
        self.report('ev_grad(el_avg): min, max:', val1.min(), val1.max())

        aux = problem.evaluate('ev_grad.i.Omega( us )', us=var, mode='qp')
        val2 = _integrate(var, aux)
        val2.shape = val1.shape
        self.report('ev_grad(qp): min, max:', val2.min(), val2.max())

        ok = self.compare_vectors(val1, val2,
                                  label1='de mode',
                                  label2='dq mode')
        if not ok:
            self.report('failed')

        return ok

    def test_ev_div(self):
        problem = self.problem

        var = problem.create_variables(['uv'])['uv']
        val = nm.arange(var.n_dof, dtype=var.dtype)
        var.set_data(val)

        val1 = problem.evaluate('ev_div.i.Omega( uv )', uv=var, mode='el_avg')
        self.report('ev_div(el_avg): min, max:', val1.min(), val1.max())

        aux = problem.evaluate('ev_div.i.Omega( uv )', uv=var, mode='qp')
        val2 = _integrate(var, aux)
        val2.shape = val1.shape
        self.report('ev_div(qp): min, max:', val2.min(), val2.max())

        ok = self.compare_vectors(val1, val2,
                                  label1='de mode',
                                  label2='dq mode')
        if not ok:
            self.report('failed')

        return ok
