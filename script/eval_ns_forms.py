#!/usr/bin/env python
"""
Operators present in the FE discretization of (adjoint) Navier-Stokes terms.
"""
from __future__ import print_function
from __future__ import absolute_import
import sympy as s
from six.moves import range

def create_scalar(name, n_ep):
    vec = s.zeros(n_ep, 1)
    for ip in range(n_ep):
        vec[ip,0] = '%s%d' % (name, ip)
    return vec

def create_vector(name, n_ep, dim):
    """ordering is DOF-by-DOF"""
    vec = s.zeros(dim * n_ep, 1)
    for ii in range(dim):
        for ip in range(n_ep):
            vec[n_ep*ii+ip,0] = '%s%d%d' % (name, ii, ip)
    return vec

def create_scalar_base(name, n_ep):
    phi = s.zeros(1, n_ep)
    for ip in range(n_ep):
        phi[0,ip] = '%s%d' % (name, ip)
    return phi

def create_vector_base(name, phic, dim):
    n_ep = phic.shape[1]
    phi = s.zeros(dim, dim * n_ep)
    indx = []
    for ii in range(dim):
        phi[ii,n_ep*ii:n_ep*(ii+1)] = phic
        indx.append(ii)
    return phi, indx

def create_scalar_base_grad(name, phic, dim):
    n_ep = phic.shape[1]
    gc = s.zeros(dim, n_ep)
    for ii in range(dim):
        for ip in range(n_ep):
            gc[ii,ip] = '%s%d%d' % (name, ii, ip)
    return gc

def create_vector_base_grad(name, gc, transpose=False):
    dim, n_ep = gc.shape
    g = s.zeros(dim * dim, dim * n_ep)
    indx = []
    if transpose:
        for ir in range(dim):
            for ic in range(dim):
                g[dim*ir+ic,n_ep*ic:n_ep*(ic+1)] = gc[ir,:]
                indx.append((ic, ir))
    else:
        for ir in range(dim):
            for ic in range(dim):
                g[dim*ir+ic,n_ep*ir:n_ep*(ir+1)] = gc[ic,:]
                indx.append((ir, ic))
    return g, indx

def create_u_operator(u, transpose=False):
    dim = u.shape[0]
    op_u = s.zeros(dim * dim, dim)
    if transpose:
        for ir in range(dim):
            for ic in range(dim):
                op_u[dim*ir+ic,ic] = u[ir]
    else:
        for ii in range(dim):
            op_u[dim*ii:dim*(ii+1),ii] = u
    return op_u

def grad_vector_to_matrix(name, gv):
    dim2 = gv.shape[0]
    dim = int(s.sqrt(dim2))
    gm = s.zeros(dim, dim)
    for ir in range(dim):
        for ic in range(dim):
            gm[ir,ic] = gv[dim*ir+ic,0]
    return gm

def substitute_continuous(expr, names, u, phi):
    pu = phi * u
    for ii in range(phi.rows):
        expr = expr.subs(pu[ii,0], names[ii])
    return expr

def create_vector_var_data(name, phi, vindx, g, gt, vgindx, u):
    gu = g * u
    gum = grad_vector_to_matrix('gum', gu)
    print('g %s:\n' % name, gum)

    gut = gt * u
    gutm = grad_vector_to_matrix('gutm', gut)
    print('gt %s:\n' % name, gutm)

    pu = phi * u
    names = ['c%s%d' % (name, indx) for indx in vindx]
    cu = substitute_continuous(pu, names, u, phi)
    print('continuous %s:\n' % name, cu)

    gnames = ['cg%s%d_%d' % (name, indx[0], indx[1]) for indx in vgindx]
    cgu = substitute_continuous(gu, gnames, u, g)
    cgum = grad_vector_to_matrix('gum', cgu)
    print('continuous g %s:\n' % name, cgum)

    cgut = substitute_continuous(gut, gnames, u, g)
    cgutm = grad_vector_to_matrix('gutm', cgut)
    print('continuous gt %s:\n' % name, cgutm)

    op_u = create_u_operator(cu)
    print(op_u)

    op_ut = create_u_operator(cu, transpose=True)
    print(op_ut)

    out = {
        'g' : gu,
        'g_m' : gum,
        'q' : pu,
        'c' : cu,
        'cg' : cgu,
        'cg_m' : cgum,
        'cgt' : cgut,
        'cgt_m' : cgutm,
        'op' : op_u,
        'opt' : op_ut,
        'names' : names,
        'gnames' : gnames,
    }

    return out

def create_scalar_var_data(name, phi, g, u):
    gu = g * u

    pu = phi * u
    names = ['c%s' % name]
    cu = substitute_continuous(pu, names, u, phi)
    print('continuous %s:\n' % name, cu)

    gnames = ['cg%s_%d' % (name, ii) for ii in range(g.shape[0])]
    cgu = substitute_continuous(gu, gnames, u, g)
    print('continuous g %s:\n' % name, cgu)

    op_gu = create_u_operator(cgu)
    print(op_gu)

    out = {
        'g' : gu,
        'q' : pu,
        'c' : cu,
        'cg' : cgu,
        'gop' : op_gu,
        'names' : names,
        'gnames' : gnames,
    }

    return out

def main():
    n_ep = 3
    dim = 2

    u = create_vector('u', n_ep, dim)
    v = create_vector('v', n_ep, dim)
    b = create_vector('b', n_ep, dim)
    p = create_scalar('p', n_ep)
    q = create_scalar('q', n_ep)
    r = create_scalar('r', n_ep)

    print(u)
    print(v)

    phic = create_scalar_base('phic', n_ep)
    phi, vindx = create_vector_base('phi', phic, dim)
    gc = create_scalar_base_grad('gc', phic, dim)
    g, vgindx = create_vector_base_grad('g', gc)
    gt, aux = create_vector_base_grad('gt', gc, transpose=True)

    print(phi)
    print(phic)
    print(gc)
    print(g)
    print(gt)

    ud = create_vector_var_data('u', phi, vindx, g, gt, vgindx, u)
    vd = create_vector_var_data('v', phi, vindx, g, gt, vgindx, v)
    bd = create_vector_var_data('b', phi, vindx, g, gt, vgindx, b)
    pd = create_scalar_var_data('p', phic, gc, p)
    qd = create_scalar_var_data('q', phic, gc, q)
    rd = create_scalar_var_data('r', phic, gc, r)
    print(list(ud.keys()))

    assert bool(bd['op'].T * g == bd['opt'].T * gt)
    assert bool(bd['opt'].T * g == bd['op'].T * gt)
    assert bool(bd['cgt_m'] == bd['cg_m'].T)

    print('((b * grad) u), v)')
    form1 = vd['c'].T * bd['op'].T * ud['cg']
    form2 = vd['c'].T * bd['opt'].T * ud['cgt']
    print(form1)
    print(form2)
    print(bool(form1 == form2))

    print('((v * grad) u), b)')
    form1 = vd['c'].T * bd['op'].T * ud['cgt']
    form2 = vd['c'].T * bd['opt'].T * ud['cg']
    print(form1)
    print(form2)
    print(bool(form1 == form2))

    print('((u * grad) v), b)')
    form1 = vd['cgt'].T * bd['op'] * ud['c']
    form2 = vd['cg'].T * bd['opt'] * ud['c']
    print(form1)
    print(form2)
    print(bool(form1 == form2))

    print('((b * grad) v), u)')
    form1 = vd['cg'].T * bd['op'] * ud['c']
    form2 = vd['cgt'].T * bd['opt'] * ud['c']
    print(form1)
    print(form2)
    print(bool(form1 == form2))

    print('((v * grad) b), u)')
    form1 = vd['c'].T * bd['cgt_m'] * ud['c']
    form2 = vd['c'].T * bd['cg_m'].T * ud['c']
    print(form1)
    print(form2)
    print(bool(form1 == form2))

    print('((b * grad) u), (b * grad) v)')
    form1 = vd['cg'].T * bd['op'] * bd['op'].T * ud['cg']
    print(form1)

    print('((u * grad) b), (b * grad) v)')
    form1 = vd['cg'].T * bd['op'] * bd['cg_m'] * ud['c']
    print(form1)

    print('(grad p, (b * grad) v)')
    form1 = vd['cg'].T * bd['op'] * pd['cg']
    print(form1)

    print('(grad q, (b * grad) u)')
    form1 = qd['cg'].T * bd['op'].T * ud['cg']
    print(form1)

    print('(grad q, (u * grad) b)')
    form1 = qd['cg'].T * bd['cg_m'] * ud['c']
    print(form1)

    print('(grad r, (u * grad) v)')
    form1 = vd['cgt'].T * rd['gop'] * ud['c']
    print(form1)

    return ud, vd, bd, pd, qd, rd

if __name__ == '__main__':
    ud, vd, bd, pd, qd, rd = main()
