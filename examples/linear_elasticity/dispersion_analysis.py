#!/usr/bin/env python
"""
Dispersion analysis of a heterogeneous finite scale periodic cell.

The periodic cell mesh has to contain two subdomains Y1, Y2, so that different
material properties can be defined in each of the subdomains (see `--pars`
option).
"""
from __future__ import absolute_import
import os
import sys
sys.path.append('.')
import functools
from copy import copy
from argparse import ArgumentParser, RawDescriptionHelpFormatter

import numpy as nm

from sfepy.base.base import import_file, output, Struct
from sfepy.base.conf import dict_from_string, ProblemConf
from sfepy.base.ioutils import ensure_path, remove_files_patterns, save_options
from sfepy.base.log import Log
from sfepy.discrete.fem import MeshIO
from sfepy.mechanics.matcoefs import stiffness_from_youngpoisson as stiffness
import sfepy.mechanics.matcoefs as mc
from sfepy.mechanics.units import apply_unit_multipliers
import sfepy.discrete.fem.periodic as per
from sfepy.homogenization.utils import define_box_regions
from sfepy.discrete import Problem
from sfepy.solvers import Solver
from sfepy.solvers.ts import TimeStepper

def apply_units_le(pars, unit_multipliers):
    new_pars = apply_unit_multipliers(pars,
                                      ['stress', 'one', 'density',
                                       'stress', 'one' ,'density'],
                                      unit_multipliers)
    return new_pars

def define_le(filename_mesh, pars, approx_order, refinement_level, solver_conf,
              plane='strain'):
    io = MeshIO.any_from_filename(filename_mesh)
    bbox = io.read_bounding_box()
    dim = bbox.shape[1]
    size = (bbox[1] - bbox[0]).max()

    options = {
        'absolute_mesh_path' : True,
        'refinement_level' : refinement_level,
    }

    fields = {
        'displacement': ('complex', dim, 'Omega', approx_order),
    }

    young1, poisson1, density1, young2, poisson2, density2 = pars
    materials = {
        'm' : ({
            'D' : {'Y1' : stiffness(dim, young=young1, poisson=poisson1,
                                    plane=plane),
                   'Y2' : stiffness(dim, young=young2, poisson=poisson2,
                                    plane=plane)},
            'density' : {'Y1' : density1, 'Y2' : density2},
        },),
        'wave' : ({
            '.vec' : [1] * dim,
        },),
    }

    variables = {
        'u' : ('unknown field', 'displacement', 0),
        'v' : ('test field', 'displacement', 'u'),
    }

    regions = {
        'Omega' : 'all',
        'Y1': 'cells of group 1',
        'Y2': 'cells of group 2',
    }
    regions.update(define_box_regions(dim,
                                      bbox[0], bbox[1], 1e-8))

    ebcs = {
    }

    if dim == 3:
        epbcs = {
            'periodic_x' : (['Left', 'Right'], {'u.all' : 'u.all'},
                            'match_x_plane'),
            'periodic_y' : (['Near', 'Far'], {'u.all' : 'u.all'},
                            'match_y_plane'),
            'periodic_z' : (['Top', 'Bottom'], {'u.all' : 'u.all'},
                            'match_z_plane'),
        }
    else:
        epbcs = {
            'periodic_x' : (['Left', 'Right'], {'u.all' : 'u.all'},
                            'match_y_line'),
            'periodic_y' : (['Bottom', 'Top'], {'u.all' : 'u.all'},
                            'match_x_line'),
        }

    per.set_accuracy(1e-8 * size)
    functions = {
        'match_x_plane' : (per.match_x_plane,),
        'match_y_plane' : (per.match_y_plane,),
        'match_z_plane' : (per.match_z_plane,),
        'match_x_line' : (per.match_x_line,),
        'match_y_line' : (per.match_y_line,),
    }

    integrals = {
        'i' : 2 * approx_order,
    }

    equations = {
        'K' : 'dw_lin_elastic.i.Omega(m.D, v, u)',
        'S' : 'dw_elastic_wave.i.Omega(m.D, wave.vec, v, u)',
        'R' : """dw_elastic_wave_cauchy.i.Omega(m.D, wave.vec, u, v)
               - dw_elastic_wave_cauchy.i.Omega(m.D, wave.vec, v, u)""",
        'M' : 'dw_volume_dot.i.Omega(m.density, v, u)',
    }

    solver_0 = solver_conf.copy()
    solver_0['name'] = 'eig'

    return locals()

def set_wave_dir_le(materials, wdir):
    wave_mat = materials['wave']
    wave_mat.datas['special']['vec'] = wdir

def _max_diff_csr(mtx1, mtx2):
    aux = nm.abs((mtx1 - mtx2).data)
    return aux.max() if len(aux) else 0.0

def save_eigenvectors(filename, svecs, pb):
    if svecs is None: return

    variables = pb.get_variables()
    # Make full eigenvectors (add DOFs fixed by boundary conditions).
    vecs = nm.empty((variables.di.ptr[-1], svecs.shape[1]),
                    dtype=svecs.dtype)
    for ii in range(svecs.shape[1]):
        vecs[:, ii] = variables.make_full_vec(svecs[:, ii])

    # Save the eigenvectors.
    out = {}
    state = pb.create_state()
    for ii in range(svecs.shape[1]):
        state.set_full(vecs[:, ii])
        aux = state.create_output_dict()
        out.update({key + '%03d' % ii : aux[key] for key in aux})

    pb.save_state(filename, out=out)

helps = {
    'pars' :
    'material parameters in Y1, Y2 subdomains in basic units'
    ' [default: %(default)s]',
    'conf' :
    'if given, an alternative problem description file with apply_units() and'
    ' define() functions [default: %(default)s]',
    'mesh_size' :
    'desired mesh size (max. of bounding box dimensions) in basic units'
    ' - the input periodic cell mesh is rescaled to this size'
    ' [default: %(default)s]',
    'unit_multipliers' :
    'basic unit multipliers (time, length, mass) [default: %(default)s]',
    'plane' :
    'for 2D problems, plane strain or stress hypothesis selection'
    ' [default: %(default)s]',
    'wave_dir' : 'the wave vector direction (will be normalized)'
    ' [default: %(default)s]',
    'mode' : 'solution mode: omega = solve a generalized EVP for omega,'
    ' kappa = solve a quadratic generalized EVP for kappa'
    ' [default: %(default)s]',
    'range' : 'the wave vector magnitude / frequency range'
    ' (like numpy.linspace) depending on the mode option'
    ' [default: %(default)s]',
    'order' : 'displacement field approximation order [default: %(default)s]',
    'refine' : 'number of uniform mesh refinements [default: %(default)s]',
    'n_eigs' : 'the number of eigenvalues to compute [default: %(default)s]',
    'eigs_only' : 'compute only eigenvalues, not eigenvectors',
    'solver_conf' : 'eigenvalue problem solver configuration options'
    ' [default: %(default)s]',
    'save_materials' : 'save material parameters into'
    ' <output_directory>/materials.vtk',
    'log_std_waves' : 'log also standard pressure dilatation and shear waves',
    'silent' : 'do not print messages to screen',
    'clear' :
    'clear old solution files from output directory',
    'output_dir' :
    'output directory [default: %(default)s]',
    'mesh_filename' :
    'input periodic cell mesh file name [default: %(default)s]',
}

def main():
    # Aluminium and epoxy.
    default_pars = '70e9,0.35,2.799e3, 3.8e9,0.27,1.142e3'
    default_solver_conf = ("kind='eig.scipy',method='eigh',tol=1.0e-5,"
                           "maxiter=1000,which='LM',sigma=0.0")

    parser = ArgumentParser(description=__doc__,
                            formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument('--pars', metavar='young1,poisson1,density1'
                        ',young2,poisson2,density2',
                        action='store', dest='pars',
                        default=default_pars, help=helps['pars'])
    parser.add_argument('--conf', metavar='filename',
                        action='store', dest='conf',
                        default=None, help=helps['conf'])
    parser.add_argument('--mesh-size', type=float, metavar='float',
                        action='store', dest='mesh_size',
                        default=None, help=helps['mesh_size'])
    parser.add_argument('--unit-multipliers',
                        metavar='c_time,c_length,c_mass',
                        action='store', dest='unit_multipliers',
                        default='1.0,1.0,1.0', help=helps['unit_multipliers'])
    parser.add_argument('--plane', action='store', dest='plane',
                        choices=['strain', 'stress'],
                        default='strain', help=helps['plane'])
    parser.add_argument('--wave-dir', metavar='float,float[,float]',
                        action='store', dest='wave_dir',
                        default='1.0,0.0,0.0', help=helps['wave_dir'])
    parser.add_argument('--mode', action='store', dest='mode',
                        choices=['omega', 'kappa'],
                        default='omega', help=helps['mode'])
    parser.add_argument('--range', metavar='start,stop,count',
                        action='store', dest='range',
                        default='10,100,10', help=helps['range'])
    parser.add_argument('--order', metavar='int', type=int,
                        action='store', dest='order',
                        default=1, help=helps['order'])
    parser.add_argument('--refine', metavar='int', type=int,
                        action='store', dest='refine',
                        default=0, help=helps['refine'])
    parser.add_argument('-n', '--n-eigs', metavar='int', type=int,
                        action='store', dest='n_eigs',
                        default=6, help=helps['n_eigs'])
    parser.add_argument('--eigs-only',
                        action='store_true', dest='eigs_only',
                        default=False, help=helps['eigs_only'])
    parser.add_argument('--solver-conf', metavar='dict-like',
                        action='store', dest='solver_conf',
                        default=default_solver_conf, help=helps['solver_conf'])
    parser.add_argument('--save-materials',
                        action='store_true', dest='save_materials',
                        default=False, help=helps['save_materials'])
    parser.add_argument('--log-std-waves',
                        action='store_true', dest='log_std_waves',
                        default=False, help=helps['log_std_waves'])
    parser.add_argument('--silent',
                        action='store_true', dest='silent',
                        default=False, help=helps['silent'])
    parser.add_argument('-c', '--clear',
                        action='store_true', dest='clear',
                        default=False, help=helps['clear'])
    parser.add_argument('-o', '--output-dir', metavar='path',
                        action='store', dest='output_dir',
                        default='output', help=helps['output_dir'])
    parser.add_argument('mesh_filename', default='',
                        help=helps['mesh_filename'])
    options = parser.parse_args()

    output_dir = options.output_dir

    output.set_output(filename=os.path.join(output_dir,'output_log.txt'),
                      combined=options.silent == False)

    if options.conf is not None:
        mod = import_file(options.conf)
        apply_units = mod.apply_units
        define = mod.define
        set_wave_dir = mod.set_wave_dir

    else:
        apply_units = apply_units_le
        define = define_le
        set_wave_dir = set_wave_dir_le

    options.pars = [float(ii) for ii in options.pars.split(',')]
    options.unit_multipliers = [float(ii)
                                for ii in options.unit_multipliers.split(',')]
    options.wave_dir = [float(ii)
                        for ii in options.wave_dir.split(',')]
    aux = options.range.split(',')
    options.range = [float(aux[0]), float(aux[1]), int(aux[2])]
    options.solver_conf = dict_from_string(options.solver_conf)

    if options.clear:
        remove_files_patterns(output_dir,
                              ['*.h5', '*.vtk', '*.txt'],
                              ignores=['output_log.txt'],
                              verbose=True)

    filename = os.path.join(output_dir, 'options.txt')
    ensure_path(filename)
    save_options(filename, [('options', vars(options))])

    pars = apply_units(options.pars, options.unit_multipliers)
    output('material parameters with applied unit multipliers:')
    output(pars)

    if options.mode == 'omega':
        rng = copy(options.range)
        rng[:2] = apply_unit_multipliers(options.range[:2],
                                         ['wave_number', 'wave_number'],
                                         options.unit_multipliers)
        output('wave number range with applied unit multipliers:', rng)

    else:
        rng = copy(options.range)
        rng[:2] = apply_unit_multipliers(options.range[:2],
                                         ['frequency', 'frequency'],
                                         options.unit_multipliers)
        output('frequency range with applied unit multipliers:', rng)

    define_problem = functools.partial(define,
                                       filename_mesh=options.mesh_filename,
                                       pars=pars,
                                       approx_order=options.order,
                                       refinement_level=options.refine,
                                       solver_conf=options.solver_conf,
                                       plane=options.plane)

    conf = ProblemConf.from_dict(define_problem(), sys.modules[__name__])

    pb = Problem.from_conf(conf)
    dim = pb.domain.shape.dim

    if dim != 2:
        options.plane = 'strain'

    wdir = nm.asarray(options.wave_dir[:dim], dtype=nm.float64)
    wdir = wdir / nm.linalg.norm(wdir)

    stepper = TimeStepper(rng[0], rng[1], dt=None, n_step=rng[2])

    bbox = pb.domain.mesh.get_bounding_box()
    size = (bbox[1] - bbox[0]).max()
    scaling0 = apply_unit_multipliers([1.0], ['length'],
                                      options.unit_multipliers)[0]
    scaling = scaling0
    if options.mesh_size is not None:
        scaling *= options.mesh_size / size
    output('scaling factor of periodic cell mesh coordinates:', scaling)
    output('new mesh size with applied unit multipliers:', scaling * size)
    pb.domain.mesh.coors[:] *= scaling
    pb.set_mesh_coors(pb.domain.mesh.coors, update_fields=True)

    bzone = 2.0 * nm.pi / (scaling * size)
    output('1. Brillouin zone size:', bzone * scaling0)
    output('1. Brillouin zone size with applied unit multipliers:', bzone)

    pb.time_update()
    pb.update_materials()

    if options.save_materials or options.log_std_waves:
        stiffness = pb.evaluate('ev_integrate_mat.2.Omega(m.D, u)',
                            mode='el_avg', copy_materials=False, verbose=False)
        young, poisson = mc.youngpoisson_from_stiffness(stiffness,
                                                        plane=options.plane)
        density = pb.evaluate('ev_integrate_mat.2.Omega(m.density, u)',
                            mode='el_avg', copy_materials=False, verbose=False)

    if options.save_materials:
        out = {}
        out['young'] = Struct(name='young', mode='cell',
                              data=young[..., None, None])
        out['poisson'] = Struct(name='poisson', mode='cell',
                                data=poisson[..., None, None])
        out['density'] = Struct(name='density', mode='cell', data=density)
        materials_filename = os.path.join(output_dir, 'materials.vtk')
        pb.save_state(materials_filename, out=out)

    # Set the normalized wave vector direction to the material(s).
    set_wave_dir(pb.get_materials(), wdir)

    conf = pb.solver_confs['eig']
    eig_solver = Solver.any_from_conf(conf)

    # Assemble the matrices.
    mtx_m = pb.mtx_a.copy()
    eq_m = pb.equations['M']
    mtx_m = eq_m.evaluate(mode='weak', dw_mode='matrix', asm_obj=mtx_m)
    mtx_m.eliminate_zeros()

    mtx_k = pb.mtx_a.copy()
    eq_k = pb.equations['K']
    mtx_k = eq_k.evaluate(mode='weak', dw_mode='matrix', asm_obj=mtx_k)
    mtx_k.eliminate_zeros()

    mtx_s = pb.mtx_a.copy()
    eq_s = pb.equations['S']
    mtx_s = eq_s.evaluate(mode='weak', dw_mode='matrix', asm_obj=mtx_s)
    mtx_s.eliminate_zeros()

    mtx_r = pb.mtx_a.copy()
    eq_r = pb.equations['R']
    mtx_r = eq_r.evaluate(mode='weak', dw_mode='matrix', asm_obj=mtx_r)
    mtx_r.eliminate_zeros()

    output('symmetry checks of real blocks:')
    output('M - M^T:', _max_diff_csr(mtx_m, mtx_m.T))
    output('K - K^T:', _max_diff_csr(mtx_k, mtx_k.T))
    output('S - S^T:', _max_diff_csr(mtx_s, mtx_s.T))
    output('R + R^T:', _max_diff_csr(mtx_r, -mtx_r.T))

    n_eigs = options.n_eigs
    if options.n_eigs > mtx_k.shape[0]:
        options.n_eigs = mtx_k.shape[0]
        n_eigs = None

    if options.mode == 'omega':
        eigenshapes_filename = os.path.join(output_dir,
                                            'frequency-eigenshapes-%s.vtk'
                                            % stepper.suffix)

        extra = []
        extra_plot_kwargs = []
        if options.log_std_waves:
            lam, mu = mc.lame_from_youngpoisson(young, poisson,
                                                plane=options.plane)
            alam = nm.average(lam)
            amu = nm.average(mu)
            adensity = nm.average(density)

            cp = nm.sqrt((alam + 2.0 * amu) / adensity)
            cs = nm.sqrt(amu / adensity)
            output('average p-wave speed:', cp)
            output('average shear wave speed:', cs)

            extra = [r'$\omega_p$', r'$\omega_s$']
            extra_plot_kwargs = [{'ls' : '--', 'color' : 'k'},
                                 {'ls' : '--', 'color' : 'gray'}]

        log = Log([[r'$\lambda_{%d}$' % ii for ii in range(options.n_eigs)],
                   [r'$\omega_{%d}$'
                    % ii for ii in range(options.n_eigs)] + extra],
                  plot_kwargs=[{}, [{}] * options.n_eigs + extra_plot_kwargs],
                  yscales=['linear', 'linear'],
                  xlabels=[r'$\kappa$', r'$\kappa$'],
                  ylabels=[r'eigenvalues $\lambda_i$',
                           r'frequencies $\omega_i$'],
                  log_filename=os.path.join(output_dir, 'frequencies.txt'),
                  aggregate=1000, sleep=0.1)

        for iv, wmag in stepper:
            output('step %d: wave vector %s' % (iv, wmag * wdir))

            mtx_a = mtx_k + wmag**2 * mtx_s + (1j * wmag) * mtx_r
            mtx_b = mtx_m

            output('A - A^H:', _max_diff_csr(mtx_a, mtx_a.H))

            if options.eigs_only:
                eigs = eig_solver(mtx_a, mtx_b, n_eigs=n_eigs,
                                  eigenvectors=False)
                svecs = None

            else:
                eigs, svecs = eig_solver(mtx_a, mtx_b, n_eigs=options.n_eigs,
                                         eigenvectors=True)
            omegas = nm.sqrt(eigs)

            output('eigs, omegas:\n', nm.c_[eigs, omegas])

            out = tuple(eigs) + tuple(omegas)
            if options.log_std_waves:
                out = out + (cp * wmag, cs * wmag)
            log(*out, x=[wmag, wmag])

            save_eigenvectors(eigenshapes_filename % iv, svecs, pb)

        log(save_figure=os.path.join(output_dir, 'frequencies.png'))
        log(finished=True)

    else:
        import scipy.sparse as sps
        from sksparse.cholmod import cholesky

        eigenshapes_filename = os.path.join(output_dir,
                                            'wave-number-eigenshapes-%s.vtk'
                                            % stepper.suffix)

        factor = cholesky(mtx_s)
        perm = factor.P()
        ir = nm.arange(len(perm))
        mtx_p = sps.coo_matrix((nm.ones_like(perm), (ir, perm)))
        mtx_l = mtx_p.T * factor.L()
        mtx_eye = sps.eye(mtx_l.shape[0], dtype=nm.float64)

        output('S - LL^T:', _max_diff_csr(mtx_s, mtx_l * mtx_l.T))

        log = Log([[r'$\kappa_{%d}$' % ii for ii in range(options.n_eigs)]],
                  plot_kwargs=[{'ls' : 'None', 'marker' : 'o'}],
                  yscales=['linear'],
                  xlabels=[r'$\omega$'],
                  ylabels=[r'wave numbers $\kappa_i$'],
                  log_filename=os.path.join(output_dir, 'wave-numbers.txt'),
                  aggregate=1000, sleep=0.1)
        for io, omega in stepper:
            output('step %d: frequency %s' % (io, omega))

            mtx_a = sps.bmat([[mtx_k - omega**2 * mtx_m, None],
                              [None, mtx_eye]])
            mtx_b = sps.bmat([[1j * mtx_r, mtx_l],
                              [mtx_l.T, None]])

            output('A - A^T:', _max_diff_csr(mtx_a, mtx_a.T))
            output('A - A^H:', _max_diff_csr(mtx_a, mtx_a.T))
            output('B - B^H:', _max_diff_csr(mtx_b, mtx_b.H))

            if options.eigs_only:
                eigs = eig_solver(mtx_a, mtx_b, n_eigs=n_eigs,
                                  eigenvectors=False)
                svecs = None

            else:
                eigs, svecs = eig_solver(mtx_a, mtx_b, n_eigs=options.n_eigs,
                                         eigenvectors=True)
            kappas = eigs

            output('kappas:\n', kappas[:, None])

            out = tuple(kappas)
            log(*out, x=[omega])

            save_eigenvectors(eigenshapes_filename % io, svecs, pb)

        log(save_figure=os.path.join(output_dir, 'wave-numbers.png'))
        log(finished=True)

if __name__ == '__main__':
    main()
