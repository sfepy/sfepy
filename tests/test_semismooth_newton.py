from __future__ import absolute_import
import sympy as sm
import numpy as nm
import scipy.sparse as sps

from sfepy.base.base import dict_to_struct
from sfepy.base.testing import TestCommon
import six
from six.moves import range

conf = {
    'name' : 'semismooth_newton',
    'kind' : 'nls.semismooth_newton',

    'semismooth' : True,

    'i_max'      : 10,
    'eps_a'      : 1e-8,
    'eps_r'      : 1e-2,
    'macheps'   : 1e-16,
    'lin_red'    : 1e-2, # Linear system error < (eps_a * lin_red).
    'ls_red_reg' : 0.1,
    'ls_red_alt' : 0.01,
    'ls_red_warp' : 0.001,
    'ls_on'      : 2.0,
    'ls_min'     : 1e-10,
    ## 'log'        : {'plot' : 'aux.png'},
}

ls_conf = {
    'name' : 'ls',
    'kind' : 'ls.scipy_direct',
}

alpha = 30.0 * nm.pi / 180.0
ca = nm.cos(alpha)
sa = nm.sin(alpha)

# Coefficient of friction.
fc = 0.01

# Spring (bar) stiffnesses.
ks = nm.ones((7,), dtype=nm.float64)

# Load vector.
fs = nm.array([0, 0, 1, 1, -1, -1], dtype=nm.float64)


def eval_matrix(mtx, **kwargs):

    num_mtx = nm.zeros(mtx.shape, dtype=nm.float64)

    mtx = mtx.subs(kwargs)
    for ir in range(mtx.shape[0]):
        for ic in range(mtx.shape[1]):
            num_mtx[ir,ic] = float(mtx[ir,ic])

    return num_mtx

def convert_to_csr(m_in):
    m_out = {}
    for key, mtx in six.iteritems(m_in):
        m_out[key] = sps.csr_matrix(mtx)

    return m_out

def define_matrices():
    e = sm.zeros(7, 1)
    e[0] = sm.sympify('u1')
    e[1] = sm.sympify('c * u3 + s * v3')
    e[2] = sm.sympify('u2 - u1')
    e[3] = sm.sympify('v4 + u2n')
    e[4] = sm.sympify('u4 - u3')
    e[5] = sm.sympify('v3 + u1n')
    e[6] = sm.sympify('c * u4 + s * v4 - c * u1 + s * u1n')
    k = sm.Matrix(7, 1, lambda i, j: sm.Symbol('k%d' % i))

    psi = (k.T * e.applyfunc(lambda x: x**2))[0]
    work = sm.sympify('u1 * g1 + u2 * g2')

    w = ['u1', 'u2', 'u3', 'u4', 'v3', 'v4']
    wbar = w[:2]
    g = ['g1', 'g2']

    mtx_a = sm.zeros(6, 6)
    for ir in range(mtx_a.shape[0]):
        for ic in range(mtx_a.shape[1]):
            mtx_a[ir,ic] = sm.diff(sm.diff(psi, w[ir]), w[ic])

    mtx_b_bar = sm.zeros(2, 6)
    for ir in range(mtx_b_bar.shape[0]):
        for ic in range(mtx_b_bar.shape[1]):
            mtx_b_bar[ir,ic] = sm.diff(sm.diff(work, g[ir]), w[ic])

    mtx_b = sm.zeros(2, 2)
    for ir in range(mtx_b.shape[0]):
        for ic in range(mtx_b.shape[1]):
            mtx_b[ir,ic] = sm.diff(sm.diff(work, wbar[ir]), g[ic])

    mtx_c = sm.eye(2)

    sn1 = sm.diff(psi, 'u1n').subs('u1n', 0)
    sn2 = sm.diff(psi, 'u2n').subs('u2n', 0)
    sn = sm.Matrix([sn1, sn2])

    mtx_d = sm.zeros(2, 6)
    for ir in range(mtx_d.shape[0]):
        for ic in range(mtx_d.shape[1]):
            mtx_d[ir,ic] = sm.diff(sn[ir], w[ic])

    subsd = {'c' : ca, 's' : sa}
    for ii, _k in enumerate(k):
        subsd['k%d' % ii] = ks[ii]

    m = {
        'A' : eval_matrix(mtx_a, **subsd),
        'Bb' : eval_matrix(mtx_b_bar),
        'B' : eval_matrix(mtx_b),
        'C' : eval_matrix(mtx_c),
    }

    ms = convert_to_csr(m)

    m['sn'] = sn.subs(subsd)
    m['D'] = mtx_d.subs(subsd)
    m['fs'] = fs

    return m, ms, w

class Test(TestCommon):

    @staticmethod
    def from_conf(conf, options):

        test = Test(conf=conf, options=options)

        test.m, test.ms, test.w_names = define_matrices()

        return test

    def test_semismooth_newton(self):
        import numpy as nm
        from sfepy.solvers import Solver

        ns = [0, 6, 2, 2]

        offsets = nm.cumsum(ns)
        nx = offsets[-1]

        iw = slice(offsets[0], offsets[1])
        ig = slice(offsets[1], offsets[2])
        il = slice(offsets[2], offsets[3])

        def fun_smooth(vec_x):
            xw = vec_x[iw]
            xg = vec_x[ig]
            xl = vec_x[il]

            m = self.ms
            rw = m['A'] * xw - m['Bb'].T * xg - self.m['fs'] 
            rg = m['Bb'] * xw + xl * xg

            rwg = nm.r_[rw, rg]

            return rwg

        def fun_smooth_grad(vec_x):
            xw = vec_x[iw]
            xl = vec_x[il]
            xg = vec_x[ig]

            m = self.m

            mzl = nm.zeros((6, 2), dtype=nm.float64)

            mw = nm.c_[m['A'], -m['Bb'].T, mzl]
            mg = nm.c_[m['Bb'], nm.diag(xl), nm.diag(xg)]

            mx = nm.r_[mw, mg]

            mx = sps.csr_matrix(mx)
            return mx

        def fun_a(vec_x):
            xw = vec_x[iw]
            xg = vec_x[ig]

            subsd = {}
            for ii, key in enumerate(self.w_names):
                subsd[key] = xw[ii]

            sn = eval_matrix(self.m['sn'], **subsd).squeeze()

            ra = nm.abs(xg) - fc * nm.abs(sn)

            return -ra

        def fun_a_grad(vec_x):
            xw = vec_x[iw]
            xg = vec_x[ig]
            xl = vec_x[il]

            subsd = {}
            for ii, key in enumerate(self.w_names):
                subsd[key] = xw[ii]

            md = eval_matrix(self.m['D'], **subsd)
            sn = eval_matrix(self.m['sn'], **subsd).squeeze()

            ma = nm.zeros((xl.shape[0], nx), dtype=nm.float64)
            ma[:,iw] = - fc * nm.sign(sn)[:,None] * md
            ma[:,ig] = nm.sign(xg)[:,None] * self.m['C']

            return -sps.csr_matrix(ma)

        def fun_b(vec_x):
            xl = vec_x[il]

            return xl

        def fun_b_grad(vec_x):
            xl = vec_x[il]

            mb = nm.zeros((xl.shape[0], nx), dtype=nm.float64)
            mb[:,il] = self.m['C']

            return sps.csr_matrix(mb)

        vec_x0 = 0.1 * nm.ones((nx,), dtype=nm.float64)

        lin_solver = Solver.any_from_conf(dict_to_struct(ls_conf))
        status = {}
        solver = Solver.any_from_conf(dict_to_struct(conf),
                                      fun_smooth=fun_smooth,
                                      fun_smooth_grad=fun_smooth_grad,
                                      fun_a=fun_a,
                                      fun_a_grad=fun_a_grad,
                                      fun_b=fun_b,
                                      fun_b_grad=fun_b_grad,
                                      lin_solver=lin_solver,
                                      status=status)

        vec_x = solver(vec_x0)

        xw = vec_x[iw]
        xg = vec_x[ig]
        xl = vec_x[il]

        self.report('x:', xw)
        self.report('g:', xg)
        self.report('l:', xl)


        subsd = {}
        for ii, key in enumerate(self.w_names):
            subsd[key] = xw[ii]

        sn = eval_matrix(self.m['sn'], **subsd).squeeze()
        self.report('sn:', sn)

        ok = status['condition'] == 0

        return ok
