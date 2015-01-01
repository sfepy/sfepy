#!/usr/bin/env python
"""
Modal analysis of a linear elastic block in 2D or 3D.

The dimension of the problem is determined by the length of the vector
in ``--dims`` option.

Examples
--------

- Run with the default arguments, show results (color = strain)::

    python examples/standalone/interactive/modal_analysis.py --show

- Clamp bottom surface of the domain, show 9 eigen-shapes::

    python examples/standalone/interactive/modal_analysis.py -b clamped -n 9 --show

- Increase mesh resolution::

    python examples/standalone/interactive/modal_analysis.py -s 31,31 -n 9 --show

- Use 3D domain::

    python examples/standalone/interactive/modal_analysis.py -d 1,1,1 -c 0,0,0 -s 8,8,8 --show
"""
import sys
sys.path.append('.')
from optparse import OptionParser

import numpy as nm
import scipy.sparse.linalg as sla

from sfepy.base.base import assert_, output, Struct
from sfepy.discrete import (FieldVariable, Material, Integral, Integrals,
                            Equation, Equations, Problem)
from sfepy.discrete.fem import FEDomain, Field
from sfepy.terms import Term
from sfepy.discrete.conditions import Conditions, EssentialBC
from sfepy.mechanics.matcoefs import stiffness_from_youngpoisson
from sfepy.mesh.mesh_generators import gen_block_mesh

usage = '%prog [options]\n' + __doc__.rstrip()

helps = {
    'dims' :
    'dimensions of the block [default: %default]',
    'centre' :
    'centre of the block [default: %default]',
    'shape' :
    'numbers of vertices along each axis [default: %default]',
    'bc_kind' :
    'kind of Dirichlet boundary conditions on the bottom surface, one of:'
    ' free, clamped [default: %default]',
    'young' : "the Young's modulus [default: %default]",
    'poisson' : "the Poisson's ratio [default: %default]",
    'order' : 'displacement field approximation order [default: %default]',
    'n_eigs' : 'the number of eigenvalues to compute [default: %default]',
    'show' : 'show the results figure',
}

def main():
    parser = OptionParser(usage=usage, version='%prog')
    parser.add_option('-d', '--dims', metavar='dims',
                      action='store', dest='dims',
                      default='[1.0, 1.0]', help=helps['dims'])
    parser.add_option('-c', '--centre', metavar='centre',
                      action='store', dest='centre',
                      default='[0.0, 0.0]', help=helps['centre'])
    parser.add_option('-s', '--shape', metavar='shape',
                      action='store', dest='shape',
                      default='[11, 11]', help=helps['shape'])
    parser.add_option('-b', '--bc-kind', metavar='kind',
                      action='store', dest='bc_kind',
                      choices=['free', 'clamped'],
                      default='free', help=helps['bc_kind'])
    parser.add_option('--young', metavar='float', type=float,
                      action='store', dest='young',
                      default=2000.0, help=helps['young'])
    parser.add_option('--poisson', metavar='float', type=float,
                      action='store', dest='poisson',
                      default=0.4, help=helps['poisson'])
    parser.add_option('--order', metavar='int', type=int,
                      action='store', dest='order',
                      default=1, help=helps['order'])
    parser.add_option('-n', '--n-eigs', metavar='int', type=int,
                      action='store', dest='n_eigs',
                      default=6, help=helps['order'])
    parser.add_option('', '--show',
                      action="store_true", dest='show',
                      default=False, help=helps['show'])
    options, args = parser.parse_args()

    assert_((0.0 < options.poisson < 0.5),
            "Poisson's ratio must be in ]0, 0.5[!")
    assert_((0 < options.order),
            'displacement approximation order must be at least 1!')

    dims = nm.array(eval(options.dims), dtype=nm.float64)
    dim = len(dims)
    centre = nm.array(eval(options.centre), dtype=nm.float64)[:dim]
    shape = nm.array(eval(options.shape), dtype=nm.int32)[:dim]

    output('dimensions:', dims)
    output('centre:    ', centre)
    output('shape:     ', shape)
    output('using values:')
    output("  Young's modulus:", options.young)
    output("  Poisson's ratio:", options.poisson)

    # Build the problem definition.
    mesh = gen_block_mesh(dims, shape, centre, name='mesh')
    domain = FEDomain('domain', mesh)

    bbox = domain.get_mesh_bounding_box()
    min_y, max_y = bbox[:, 1]
    eps = 1e-8 * (max_y - min_y)
    omega = domain.create_region('Omega', 'all')
    bottom = domain.create_region('Bottom',
                                  'vertices in (y < %.10f)' % (min_y + eps),
                                  'facet')

    field = Field.from_args('fu', nm.float64, 'vector', omega,
                            approx_order=options.order)

    u = FieldVariable('u', 'unknown', field)
    v = FieldVariable('v', 'test', field, primary_var_name='u')

    mtx_d = stiffness_from_youngpoisson(dim, options.young, options.poisson)

    m = Material('m', D=mtx_d)

    integral = Integral('i', order=2*options.order)

    t1 = Term.new('dw_lin_elastic(m.D, v, u)',
                  integral, omega, m=m, v=v, u=u)
    t2 = Term.new('dw_volume_dot(v, u)',
                  integral, omega, v=v, u=u)
    eq1 = Equation('balance', t1)
    eq2 = Equation('mass', t2)
    lhs_eqs = Equations([eq1, eq2])

    pb = Problem('modal', equations=lhs_eqs)

    if options.bc_kind == 'free':
        pb.time_update()
        n_rbm = dim * (dim + 1) / 2

    else:
        fixed_b = EssentialBC('FixedB', bottom, {'u.all' : 0.0})
        pb.time_update(ebcs=Conditions([fixed_b]))
        n_rbm = 0

    pb.update_materials()

    # Assemble stiffness and mass matrices.
    mtx_k = eq1.evaluate(mode='weak', dw_mode='matrix', asm_obj=pb.mtx_a)
    mtx_m = mtx_k.copy()
    mtx_m.data[:] = 0.0
    mtx_m = eq2.evaluate(mode='weak', dw_mode='matrix', asm_obj=mtx_m)

    try:
        eigs, svecs = sla.eigsh(mtx_k, k=options.n_eigs + n_rbm, M=mtx_m,
                                which='SM',
                                maxiter=10000)
    except sla.ArpackNoConvergence as ee:
        eigs = ee.eigenvalues
        svecs = ee.eigenvectors
        output('only %d eigenvalues converged!' % len(eigs))

    eigs = eigs[n_rbm:]
    svecs = svecs[:, n_rbm:]

    output('eigenvalues:', eigs)
    output('eigen-frequencies:', nm.sqrt(eigs))

    # Make full eigenvectors (add DOFs fixed by boundary conditions).
    variables = pb.get_variables()

    vecs = nm.empty((variables.di.ptr[-1], svecs.shape[1]),
                    dtype=nm.float64)
    for ii in xrange(svecs.shape[1]):
        vecs[:, ii] = variables.make_full_vec(svecs[:, ii])

    # Save the eigenvectors.
    out = {}
    state = pb.create_state()
    for ii in xrange(eigs.shape[0]):
        state.set_full(vecs[:, ii])
        aux = state.create_output_dict()
        strain = pb.evaluate('ev_cauchy_strain.i.Omega(u)',
                             integrals=Integrals([integral]),
                             mode='el_avg', verbose=False)
        out['u%03d' % ii] = aux.popitem()[1]
        out['strain%03d' % ii] = Struct(mode='cell', data=strain)

    pb.save_state('eigenshapes.vtk', out=out)
    pb.save_regions_as_groups('regions')

    if options.show:
        # Show the solution. If the approximation order is greater than 1, the
        # extra DOFs are simply thrown away.
        from sfepy.postprocess.viewer import Viewer
        from sfepy.postprocess.domain_specific import DomainSpecificPlot

        scaling = 0.05 * dims.max() / nm.abs(vecs).max()

        ds = {}
        for ii in xrange(eigs.shape[0]):
            pd = DomainSpecificPlot('plot_displacements',
                                    ['rel_scaling=%s' % scaling,
                                     'color_kind="tensors"',
                                     'color_name="strain%03d"' % ii])
            ds['u%03d' % ii] = pd

        view = Viewer('eigenshapes.vtk')
        view(domain_specific=ds, only_names=sorted(ds.keys()),
             is_scalar_bar=False, is_wireframe=True)

if __name__ == '__main__':
    main()
